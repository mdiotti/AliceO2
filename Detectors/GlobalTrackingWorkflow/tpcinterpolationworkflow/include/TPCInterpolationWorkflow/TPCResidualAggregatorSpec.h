// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef O2_TPC_RESIDUALAGGREGATORSPEC_H
#define O2_TPC_RESIDUALAGGREGATORSPEC_H

/// \file   TPCResidualAggregatorSpec.h
/// \brief DPL device for collecting and binning TPC cluster residuals
/// \author Ole Schmidt

#include "DetectorsCalibration/Utils.h"
#include "SpacePoints/TrackResiduals.h"
#include "SpacePoints/ResidualAggregator.h"
#include "CommonUtils/MemFileHelper.h"
#include "Framework/Task.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/ControlService.h"
#include "Framework/WorkflowSpec.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "Framework/RawDeviceService.h"
#include <fairmq/Device.h>
#include <chrono>

using namespace o2::framework;
using GID = o2::dataformats::GlobalTrackID;

namespace o2
{
namespace calibration
{

class ResidualAggregatorDevice : public o2::framework::Task
{
 public:
  ResidualAggregatorDevice(std::shared_ptr<o2::base::GRPGeomRequest> req, bool trackInput, bool ctpInput, bool writeOutput, bool writeUnbinnedResiduals, bool writeBinnedResiduals, bool writeTrackData, std::shared_ptr<o2::globaltracking::DataRequest> dataRequest) : mCCDBRequest(req), mTrackInput(trackInput), mCTPInput(ctpInput), mWriteOutput(writeOutput), mWriteUnbinnedResiduals(writeUnbinnedResiduals), mWriteBinnedResiduals(writeBinnedResiduals), mWriteTrackData(writeTrackData), mDataRequest(dataRequest) {}

  void init(o2::framework::InitContext& ic) final
  {
    o2::base::GRPGeomHelper::instance().setRequest(mCCDBRequest);
    int minEnt = ic.options().get<int>("min-entries");
    auto slotLength = ic.options().get<uint32_t>("tf-per-slot");
    if (slotLength == 0) {
      slotLength = o2::calibration::INFINITE_TF;
    }
    auto updateInterval = ic.options().get<uint32_t>("updateInterval");
    auto delay = ic.options().get<uint32_t>("max-delay");

    std::string outputDirConf = ic.options().get<std::string>("output-dir");
    std::string outputDir = o2::utils::Str::rectifyDirectory(outputDirConf);
    std::string metaFileDir = ic.options().get<std::string>("meta-output-dir");
    bool storeMetaFile = false;
    if (metaFileDir != "/dev/null") {
      metaFileDir = o2::utils::Str::rectifyDirectory(metaFileDir);
      storeMetaFile = true;
    }
    if (!mWriteOutput && (outputDirConf != "none" || storeMetaFile)) {
      LOGF(info, "File output is disabled, but output directory %s was specified and meta file storage is set to %i. No output file will be written", outputDirConf, storeMetaFile);
      storeMetaFile = false;
    }
    mAggregator = std::make_unique<o2::tpc::ResidualAggregator>(minEnt);
    mAggregator->setOutputDir(outputDir);
    if (storeMetaFile) {
      mAggregator->setMetaFileOutputDir(metaFileDir);
    }
    if (!mWriteOutput) {
      mAggregator->disableFileWriting();
    }
    int autosave = ic.options().get<int>("autosave-interval");
    mAggregator->setAutosaveInterval(autosave);
    // TODO mAggregator should get an option to set the binning externally (expose TrackResiduals::setBinning methods to user? as command line option?)
    mAggregator->setSlotLength(slotLength);
    mAggregator->setMaxSlotsDelay(delay);
    mAggregator->setCheckIntervalInfiniteSlot(updateInterval);
    mAggregator->setWriteBinnedResiduals(mWriteBinnedResiduals);
    mAggregator->setWriteUnbinnedResiduals(mWriteUnbinnedResiduals);
    mAggregator->setWriteTrackData(mWriteTrackData);
    mAggregator->setCompression(ic.options().get<int>("compression"));
  }

  void finaliseCCDB(o2::framework::ConcreteDataMatcher& matcher, void* obj) final
  {
    o2::base::GRPGeomHelper::instance().finaliseCCDB(matcher, obj);
  }

  void run(o2::framework::ProcessingContext& pc) final
  {
    auto runStartTime = std::chrono::high_resolution_clock::now();
    o2::globaltracking::RecoContainer recoCont;
    recoCont.collectData(pc, *mDataRequest);
    updateTimeDependentParams(pc);
    std::chrono::duration<double, std::milli> ccdbUpdateTime = std::chrono::high_resolution_clock::now() - runStartTime;

    auto residualsData = pc.inputs().get<gsl::span<o2::tpc::UnbinnedResid>>("unbinnedRes");
    // track data input is optional
    const gsl::span<const o2::tpc::TrackData>* trkDataPtr = nullptr;
    using trkDataType = std::decay_t<decltype(pc.inputs().get<gsl::span<o2::tpc::TrackData>>(""))>;
    std::optional<trkDataType> trkData;
    if (mTrackInput) {
      trkData.emplace(pc.inputs().get<gsl::span<o2::tpc::TrackData>>("trkData"));
      trkDataPtr = &trkData.value();
    }
    // CTP lumi input (optional)
    const o2::ctp::LumiInfo* lumi = nullptr;
    using lumiDataType = std::decay_t<decltype(pc.inputs().get<o2::ctp::LumiInfo>(""))>;
    std::optional<lumiDataType> lumiInput;
    if (mCTPInput) {
      recoCont.getCTPLumi();
      lumiInput = recoCont.getCTPLumi();
      lumi = &lumiInput.value();
    }

    auto data = std::make_pair<gsl::span<const o2::tpc::TrackData>, gsl::span<const o2::tpc::UnbinnedResid>>(std::move(*trkData), std::move(residualsData));
    o2::base::TFIDInfoHelper::fillTFIDInfo(pc, mAggregator->getCurrentTFInfo());
    LOG(info) << "Processing TF " << mAggregator->getCurrentTFInfo().tfCounter << " with " << trkData->size() << " tracks and " << residualsData.size() << " unbinned residuals associated to them";
    mAggregator->process(data, lumi);
    std::chrono::duration<double, std::milli> runDuration = std::chrono::high_resolution_clock::now() - runStartTime;
    LOGP(info, "Duration for run method: {} ms. From this taken for time dependent param update: {} ms",
         std::chrono::duration_cast<std::chrono::milliseconds>(runDuration).count(),
         std::chrono::duration_cast<std::chrono::milliseconds>(ccdbUpdateTime).count());
  }

  void endOfStream(o2::framework::EndOfStreamContext& ec) final
  {
    LOG(info) << "Finalizing calibration for end of stream";
    mAggregator->checkSlotsToFinalize(o2::calibration::INFINITE_TF);
    mAggregator.reset(); // must invoke destructor manually here, otherwise we get a segfault
  }

 private:
  void updateTimeDependentParams(ProcessingContext& pc)
  {
    o2::base::GRPGeomHelper::instance().checkUpdates(pc);
    static bool initOnceDone = false;
    if (!initOnceDone) {
      initOnceDone = true;
      mAggregator->setDataTakingContext(pc.services().get<DataTakingContext>());
    }
  }
  std::unique_ptr<o2::tpc::ResidualAggregator> mAggregator; ///< the TimeSlotCalibration device
  std::shared_ptr<o2::base::GRPGeomRequest> mCCDBRequest;
  std::shared_ptr<o2::globaltracking::DataRequest> mDataRequest; ///< optional CTP input
  bool mTrackInput{false};             ///< flag whether to expect track data as input
  bool mCTPInput{false};               ///< flag whether to expect luminosity input from CTP
  bool mWriteOutput{true};             ///< if false, no output file will be written
  bool mWriteBinnedResiduals{false};   ///< flag, whether to write binned residuals to output file
  bool mWriteUnbinnedResiduals{false}; ///< flag, whether to write unbinned residuals to output file
  bool mWriteTrackData{false};         ///< flag, whether to write track data to output file
};

} // namespace calibration

namespace framework
{

DataProcessorSpec getTPCResidualAggregatorSpec(bool trackInput, bool ctpInput, bool writeOutput, bool writeUnbinnedResiduals, bool writeBinnedResiduals, bool writeTrackData)
{
  std::shared_ptr<o2::globaltracking::DataRequest> dataRequest = std::make_shared<o2::globaltracking::DataRequest>();
  if (ctpInput) {
    dataRequest->requestClusters(GID::getSourcesMask("CTP"), false);
  }
  auto& inputs = dataRequest->inputs;
  inputs.emplace_back("unbinnedRes", "GLO", "UNBINNEDRES");
  if (trackInput) {
    inputs.emplace_back("trkData", "GLO", "TRKDATA");
  }
  auto ccdbRequest = std::make_shared<o2::base::GRPGeomRequest>(true,                           // orbitResetTime
                                                                true,                           // GRPECS=true
                                                                false,                          // GRPLHCIF
                                                                false,                          // GRPMagField
                                                                false,                          // askMatLUT
                                                                o2::base::GRPGeomRequest::None, // geometry
                                                                inputs);
  return DataProcessorSpec{
    "residual-aggregator",
    inputs,
    Outputs{},
    AlgorithmSpec{adaptFromTask<o2::calibration::ResidualAggregatorDevice>(ccdbRequest, trackInput, ctpInput, writeOutput, writeUnbinnedResiduals, writeBinnedResiduals, writeTrackData, dataRequest)},
    Options{
      {"tf-per-slot", VariantType::UInt32, 6'000u, {"number of TFs per calibration time slot (put 0 for infinite slot length)"}},
      {"updateInterval", VariantType::UInt32, 6'000u, {"update interval in number of TFs in case slot length is infinite"}},
      {"max-delay", VariantType::UInt32, 1u, {"number of slots in past to consider"}},
      {"min-entries", VariantType::Int, 0, {"minimum number of entries on average per voxel"}},
      {"compression", VariantType::Int, 101, {"ROOT compression setting for output file (see TFile documentation for meaning of this number)"}},
      {"output-dir", VariantType::String, "none", {"Output directory for residuals, must exist"}},
      {"meta-output-dir", VariantType::String, "/dev/null", {"Residuals metadata output directory, must exist (if not /dev/null)"}},
      {"autosave-interval", VariantType::Int, 0, {"Write output to file for every n-th TF. 0 means this feature is OFF"}}}};
}

} // namespace framework
} // namespace o2

#endif // O2_TPC_RESIDUALAGGREGATORSPEC_H
