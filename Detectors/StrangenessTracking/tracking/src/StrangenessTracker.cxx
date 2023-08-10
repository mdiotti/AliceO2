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
/// \file StrangenessTracker.cxx
/// \brief

#include <numeric>
#include "StrangenessTracking/StrangenessTracker.h"
#include "ITStracking/IOUtils.h"

namespace o2
{
namespace strangeness_tracking
{

bool StrangenessTracker::loadData(const o2::globaltracking::RecoContainer& recoData)
{
  clear();

  mInputV0tracks = recoData.getV0s();
  mInputCascadeTracks = recoData.getCascades();
  mInputITStracks = recoData.getITSTracks();
  mInputITSidxs = recoData.getITSTracksClusterRefs();

  auto compClus = recoData.getITSClusters();
  auto clusPatt = recoData.getITSClustersPatterns();
  auto pattIt = clusPatt.begin();
  mInputITSclusters.reserve(compClus.size());
  mInputClusterSizes.resize(compClus.size());
  o2::its::ioutils::convertCompactClusters(compClus, pattIt, mInputITSclusters, mDict);
  auto pattIt2 = clusPatt.begin();
  getClusterSizes(mInputClusterSizes, compClus, pattIt2, mDict);

  mITSvtxBrackets.resize(mInputITStracks.size());
  for (int i = 0; i < mInputITStracks.size(); i++) {
    mITSvtxBrackets[i] = {-1, -1};
  }

<<<<<<< HEAD
  // build time bracket for each ITS and kink tracks

=======
  // build time bracket for each ITS and kink track
>>>>>>> ba066281d2e3d3164616381d8262de1f3971633a
  auto trackIndex = recoData.getPrimaryVertexMatchedTracks(); // Global ID's for associated tracks
  auto vtxRefs = recoData.getPrimaryVertexMatchedTrackRefs(); // references from vertex to these track IDs

  std::unordered_map<GIndex, int> kinkMap;

  int nv = vtxRefs.size();
  for (int iv = 0; iv < nv; iv++) {
    const auto& vtref = vtxRefs[iv];
    int it = vtref.getFirstEntry(), itLim = it + vtref.getEntries();
    for (; it < itLim; it++) {
      auto tvid = trackIndex[it];
      if (!recoData.isTrackSourceLoaded(tvid.getSource())) {
        continue;
      }

      if (tvid.getSource() == GIndex::ITS) {
<<<<<<< HEAD

=======
>>>>>>> ba066281d2e3d3164616381d8262de1f3971633a
        if (mITSvtxBrackets[tvid.getIndex()].getMin() == -1) {
          mITSvtxBrackets[tvid.getIndex()].setMin(iv);
          mITSvtxBrackets[tvid.getIndex()].setMax(iv);
        } else {
          mITSvtxBrackets[tvid.getIndex()].setMax(iv);
        }
        continue;
      }

      // if (!mStrParams->mKinkFinder || tvid.getSource() != GIndex::TPC) { // exclude TPC only tracks for the time being (?)
      if (!mStrParams->mKinkFinder) {
        continue;
      }

      if (tvid.isAmbiguous()) {

        if (kinkMap.find(tvid) != kinkMap.end()) { // track already in the map, update the time bracket
          mKinkTracks[kinkMap[tvid]].vtxBracket.setMax(iv);
          continue;
        }
      }

      kinkTrackHelper kinkHelper;
      kinkHelper.track = recoData.getTrackParam(tvid);
      kinkHelper.index = tvid;
      kinkHelper.vtxBracket = {iv, iv};
      kinkHelper.itsRef = recoData.getITSContributorGID(tvid);
      if (kinkHelper.itsRef.getSource() == GIndex::ITSAB) {
        mKinkTracks.push_back(kinkHelper);
        kinkMap[tvid] = mKinkTracks.size() - 1;
      }
    }
  }

  if (mMCTruthON) {
    mITSClsLabels = recoData.mcITSClusters.get();
    mITSTrkLabels = recoData.getITSTracksMCLabels();
  }

  LOG(debug) << "V0 tracks size: " << mInputV0tracks.size();
  LOG(debug) << "Cascade tracks size: " << mInputCascadeTracks.size();
  LOG(debug) << "ITS tracks size: " << mInputITStracks.size();
  LOG(debug) << "ITS idxs size: " << mInputITSidxs.size();
  LOG(debug) << "ITS clusters size: " << mInputITSclusters.size();
  LOG(debug) << "VtxRefs size: " << vtxRefs.size();
  LOG(debug) << "Kink tracks size: " << mKinkTracks.size();

  return true;
}

void StrangenessTracker::prepareITStracks() // sort tracks by eta and phi and select only tracks with vertex matching
{

  for (int iTrack{0}; iTrack < mInputITStracks.size(); iTrack++) {
    if (mITSvtxBrackets[iTrack].getMin() == -1) {
      continue;
    }
    mSortedITStracks.push_back(mInputITStracks[iTrack]);
    mSortedITSindexes.push_back(iTrack);
  }

  mTracksIdxTable.resize(mUtils.mPhiBins * mUtils.mEtaBins + 1);
  std::sort(mSortedITStracks.begin(), mSortedITStracks.end(), [&](o2::its::TrackITS& a, o2::its::TrackITS& b) { return mUtils.getBinIndex(a.getEta(), a.getPhi()) < mUtils.getBinIndex(b.getEta(), b.getPhi()); });
  std::sort(mSortedITSindexes.begin(), mSortedITSindexes.end(), [&](int i, int j) { return mUtils.getBinIndex(mInputITStracks[i].getEta(), mInputITStracks[i].getPhi()) < mUtils.getBinIndex(mInputITStracks[j].getEta(), mInputITStracks[j].getPhi()); });

  for (auto& track : mSortedITStracks) {
    mTracksIdxTable[mUtils.getBinIndex(track.getEta(), track.getPhi())]++;
  }
  std::exclusive_scan(mTracksIdxTable.begin(), mTracksIdxTable.begin() + mUtils.mPhiBins * mUtils.mEtaBins, mTracksIdxTable.begin(), 0);
  mTracksIdxTable[mUtils.mPhiBins * mUtils.mEtaBins] = mSortedITStracks.size();
}

void StrangenessTracker::process()
{
<<<<<<< HEAD

  std::vector<int> ITSidx = {19167, 61925, 98281, 98329, 113062, 74660, 85521, 103831, 24592, 76456, 131501, 36667, 37472, 97235, 16289, 24907, 64634, 21184, 13105, 14897, 108210, 100322, 141805, 143252, 59511, 16964, 52173, 53271, 112165, 5870, 10629, 101034, 22340, 42871, 111376, 133087, 136063, 48765, 26155, 121965, 30387, 31582, 33364, 77546, 129909, 110429, 128618, 23110, 109749, 76782, 95490, 91254, 92699, 33867, 10804, 67972};
  std::vector<int> ITSTPCidx = {108229, 97511, 97731, 97731, 103540, 92531, 95391, 99081, 106058, 102899, 94105, 91842, 103727, 95834, 98241, 96521, 96349, 101056, 90870, 92382, 90464, 106512, 107629, 105819, 96777, 102940, 91378, 95267, 93804, 94456, 102674, 95229, 100501, 93771, 94140, 97601, 100731, 102347, 93713, 103804, 98648, 99297, 98271, 93526, 92622, 96470, 96078, 97388, 99124, 103006, 97821, 90098, 92486, 97994, 95952, 97787};
  std::vector<int> firstIdx = {767, 2383, 3804, 3807, 4428, 3034, 3486, 4218, 959, 3020, 5240, 1502, 1513, 3906, 681, 1014, 2621, 877, 557, 638, 4304, 3924, 5560, 5624, 2305, 677, 2083, 2132, 4412, 216, 433, 4150, 852, 1644, 4398, 5292, 5416, 1944, 1058, 4872, 1187, 1254, 1343, 3040, 5161, 4412, 5174, 931, 4502, 3070, 3843, 3652, 3676, 1386, 414, 2707};
  std::vector<int> lastIdx = {767, 2383, 3804, 3807, 4435, 3034, 3489, 4218, 961, 3020, 5240, 1502, 1513, 3906, 681, 1019, 2621, 881, 557, 638, 4304, 3929, 5560, 5624, 2305, 677, 2083, 2132, 4412, 216, 433, 4150, 852, 1644, 4403, 5297, 5416, 1950, 1058, 4872, 1193, 1254, 1343, 3054, 5167, 4412, 5183, 931, 4502, 3075, 3843, 3659, 3676, 1386, 414, 2707};
  std::vector<int> firstIdxITSTPC = {763, 2377, 3801, 3801, 4432, 3030, 3486, 4216, 959, 3015, 5236, 1497, 1512, 3904, 680, 1014, 2618, 878, 554, 638, 4301, 3924, 5559, 5624, 2301, 675, 2080, 2132, 4411, 205, 431, 4146, 850, 1643, 4398, 5294, 5415, 1947, 1056, 4869, 1188, 1251, 1342, 3043, 5165, 4409, 5174, 929, 4499, 3071, 3841, 3653, 3672, 1383, 411, 2698};
  std::vector<int> lastIdxITSTPC = {766, 2384, 3806, 3806, 4433, 3032, 3489, 4222, 960, 3024, 5244, 1502, 1517, 3914, 681, 1019, 2619, 880, 559, 639, 4305, 3929, 5562, 5626, 2303, 679, 2083, 2132, 4417, 210, 434, 4151, 852, 1649, 4404, 5297, 5419, 1951, 1060, 4874, 1193, 1256, 1345, 3048, 5167, 4415, 5179, 929, 4506, 3074, 3843, 3659, 3675, 1386, 418, 2702};

=======
  /*
>>>>>>> ba066281d2e3d3164616381d8262de1f3971633a
  // Loop over V0s
  mDaughterTracks.resize(2); // resize to 2 prongs: first positive second negative

  for (int iV0{0}; iV0 < mInputV0tracks.size(); iV0++) {
    LOG(debug) << "Analysing V0: " << iV0 + 1 << "/" << mInputV0tracks.size();
    auto& v0 = mInputV0tracks[iV0];
    mV0dauIDs[kV0DauPos] = v0.getProngID(kV0DauPos), mV0dauIDs[kV0DauNeg] = v0.getProngID(kV0DauNeg);
    auto posTrack = v0.getProng(kV0DauPos);
    auto negTrack = v0.getProng(kV0DauNeg);
    auto alphaV0 = calcV0alpha(v0);
    alphaV0 > 0 ? posTrack.setAbsCharge(2) : negTrack.setAbsCharge(2);
    V0 correctedV0; // recompute V0 for Hypertriton

    if (!recreateV0(posTrack, negTrack, correctedV0)) {
      continue;
    }

    mStrangeTrack.mPartType = dataformats::kStrkV0;

    auto v0R2 = v0.calcR2();
    auto iBinsV0 = mUtils.getBinRect(correctedV0.getEta(), correctedV0.getPhi(), mStrParams->mEtaBinSize, mStrParams->mPhiBinSize);
    for (int& iBinV0 : iBinsV0) {
      for (int iTrack{mTracksIdxTable[iBinV0]}; iTrack < TMath::Min(mTracksIdxTable[iBinV0 + 1], int(mSortedITStracks.size())); iTrack++) {
        mStrangeTrack.mMother = (o2::track::TrackParCovF)correctedV0;
        mDaughterTracks[kV0DauPos] = correctedV0.getProng(kV0DauPos);
        mDaughterTracks[kV0DauNeg] = correctedV0.getProng(kV0DauNeg);
        mITStrack = mSortedITStracks[iTrack];
        auto& ITSindexRef = mSortedITSindexes[iTrack];
<<<<<<< HEAD
        if (mITSvtxBrackets[ITSindexRef].getMin() > v0.getVertexID() || mITSvtxBrackets[ITSindexRef].getMax() < v0.getVertexID()) {
=======

        if (mStrParams->mVertexMatching && (mITSvtxBrackets[ITSindexRef].getMin() > v0.getVertexID() ||
                                            mITSvtxBrackets[ITSindexRef].getMax() < v0.getVertexID())) {
>>>>>>> ba066281d2e3d3164616381d8262de1f3971633a
          continue;
        }

        if (matchDecayToITStrack(sqrt(v0R2))) {

          auto propInstance = o2::base::Propagator::Instance();
          o2::track::TrackParCov decayVtxTrackClone = mStrangeTrack.mMother; // clone track and propagate to decay vertex
          if (!propInstance->propagateToX(decayVtxTrackClone, mStrangeTrack.mDecayVtx[0], getBz(), o2::base::PropagatorImpl<float>::MAX_SIN_PHI, o2::base::PropagatorImpl<float>::MAX_STEP, mCorrType)) {
            LOG(debug) << "Mother propagation to decay vertex failed";
            continue;
          }

          decayVtxTrackClone.getPxPyPzGlo(mStrangeTrack.mDecayMom);
          auto p2mom = decayVtxTrackClone.getP2();
          auto p2pos = mFitter3Body.getTrack(kV0DauPos).getP2(); // positive V0 daughter
          auto p2neg = mFitter3Body.getTrack(kV0DauNeg).getP2(); // negative V0 daughter
          if (alphaV0 > 0) {
            mStrangeTrack.mMasses[0] = calcMotherMass(p2mom, p2pos, p2neg, PID::Helium3, PID::Pion); // Hypertriton invariant mass at decay vertex
            mStrangeTrack.mMasses[1] = calcMotherMass(p2mom, p2pos, p2neg, PID::Alpha, PID::Pion);   // Hyperhydrogen4Lam invariant mass at decay vertex
          } else {
            mStrangeTrack.mMasses[0] = calcMotherMass(p2mom, p2neg, p2pos, PID::Helium3, PID::Pion); // Anti-Hypertriton invariant mass at decay vertex
            mStrangeTrack.mMasses[1] = calcMotherMass(p2mom, p2neg, p2pos, PID::Alpha, PID::Pion);   // Anti-Hyperhydrogen4Lam invariant mass at decay vertex
          }

          LOG(debug) << "ITS Track matched with a V0 decay topology ....";
          LOG(debug) << "Number of ITS track clusters attached: " << mITStrack.getNumberOfClusters();
          mStrangeTrack.mDecayRef = iV0;
          mStrangeTrack.mITSRef = mSortedITSindexes[iTrack];
          mStrangeTrackVec.push_back(mStrangeTrack);
          mClusAttachments.push_back(mStructClus);
          if (mMCTruthON) {
            auto lab = getStrangeTrackLabel();
            mStrangeTrackLabels.push_back(lab);
          }
        }
      }
    }
  }

  // Loop over Cascades
  mDaughterTracks.resize(3); // resize to 3 prongs: first bachelor, second V0 pos, third V0 neg

  for (int iCasc{0}; iCasc < mInputCascadeTracks.size(); iCasc++) {
    LOG(debug) << "Analysing Cascade: " << iCasc + 1 << "/" << mInputCascadeTracks.size();
<<<<<<< HEAD
=======

>>>>>>> ba066281d2e3d3164616381d8262de1f3971633a
    auto& casc = mInputCascadeTracks[iCasc];
    auto& cascV0 = mInputV0tracks[casc.getV0ID()];
    mV0dauIDs[kV0DauPos] = cascV0.getProngID(kV0DauPos);
    mV0dauIDs[kV0DauNeg] = cascV0.getProngID(kV0DauNeg);

    mStrangeTrack.mPartType = dataformats::kStrkCascade;
    // first: bachelor, second: V0 pos, third: V0 neg
    auto cascR2 = casc.calcR2();
    auto iBinsCasc = mUtils.getBinRect(casc.getEta(), casc.getPhi(), mStrParams->mEtaBinSize, mStrParams->mPhiBinSize);
    for (int& iBinCasc : iBinsCasc) {
      for (int iTrack{mTracksIdxTable[iBinCasc]}; iTrack < TMath::Min(mTracksIdxTable[iBinCasc + 1], int(mSortedITStracks.size())); iTrack++) {
        mStrangeTrack.mMother = (o2::track::TrackParCovF)casc;
        mDaughterTracks[kV0DauPos] = cascV0.getProng(kV0DauPos);
        mDaughterTracks[kV0DauNeg] = cascV0.getProng(kV0DauNeg);
        mDaughterTracks[kBach] = casc.getBachelorTrack();
        mITStrack = mSortedITStracks[iTrack];
        auto& ITSindexRef = mSortedITSindexes[iTrack];
        LOG(debug) << "----------------------";
        LOG(debug) << "CascV0: " << casc.getV0ID() << ", Bach ID: " << casc.getBachelorID() << ", ITS track ref: " << mSortedITSindexes[iTrack];

        if (mITSvtxBrackets[ITSindexRef].getMin() > casc.getVertexID() || mITSvtxBrackets[ITSindexRef].getMax() < casc.getVertexID()) {
          LOG(debug) << "Vertex ID mismatch: " << mITSvtxBrackets[ITSindexRef].getMin() << " < " << casc.getVertexID() << " < " << mITSvtxBrackets[ITSindexRef].getMax();
          continue;
        }

        if (matchDecayToITStrack(sqrt(cascR2))) {

          auto propInstance = o2::base::Propagator::Instance();
          o2::track::TrackParCov decayVtxTrackClone = mStrangeTrack.mMother; // clone track and propagate to decay vertex
          if (!propInstance->propagateToX(decayVtxTrackClone, mStrangeTrack.mDecayVtx[0], getBz(), o2::base::PropagatorImpl<float>::MAX_SIN_PHI, o2::base::PropagatorImpl<float>::MAX_STEP, mCorrType)) {
            LOG(debug) << "Mother propagation to decay vertex failed";
            continue;
          }
          decayVtxTrackClone.getPxPyPzGlo(mStrangeTrack.mDecayMom);
          auto p2mom = decayVtxTrackClone.getP2();
          auto p2V0 = mFitter3Body.getTrack(0).getP2();                                           // V0 momentum at decay vertex
          auto p2bach = mFitter3Body.getTrack(1).getP2();                                         // bachelor momentum at decay vertex
          mStrangeTrack.mMasses[0] = calcMotherMass(p2mom, p2V0, p2bach, PID::Lambda, PID::Pion); // Xi invariant mass at decay vertex
          mStrangeTrack.mMasses[1] = calcMotherMass(p2mom, p2V0, p2bach, PID::Lambda, PID::Kaon); // Omega invariant mass at decay vertex

          LOG(debug) << "ITS Track matched with a Cascade decay topology ....";
          LOG(debug) << "Number of ITS track clusters attached: " << mITStrack.getNumberOfClusters();

          mStrangeTrack.mDecayRef = iCasc;
          mStrangeTrack.mITSRef = mSortedITSindexes[iTrack];
          mStrangeTrackVec.push_back(mStrangeTrack);
          mClusAttachments.push_back(mStructClus);
          if (mMCTruthON) {
            auto lab = getStrangeTrackLabel();
            mStrangeTrackLabels.push_back(lab);
          }
        }
      }
    }
  }
<<<<<<< HEAD

=======
  */
>>>>>>> ba066281d2e3d3164616381d8262de1f3971633a
  if (!mStrParams->mKinkFinder) {
    return;
  }

  // loop over Kinks

  LOG(debug) << "Analysing " << mKinkTracks.size() << " Kinks ...";

  for (int iKink{0}; iKink < mKinkTracks.size(); iKink++) {
    LOG(debug) << "--------------------------------------------";
    LOG(debug) << "Analysing Kink: " << iKink + 1 << "/" << mKinkTracks.size();
    auto kink = mKinkTracks[iKink].track;
    auto vrtxTrackIdx = mKinkTracks[iKink].index;
    auto kinkVrtxIDmin = mKinkTracks[iKink].vtxBracket.getMin();
    auto kinkVrtxIDmax = mKinkTracks[iKink].vtxBracket.getMax();
    auto kinkITSid = mKinkTracks[iKink].itsRef;
    LOG(debug) << "Kink source: " << vrtxTrackIdx.getSourceName();
    LOG(debug) << "Kink Vtx ID: " << kinkVrtxIDmin << " - " << kinkVrtxIDmax;
    if (kinkITSid.isSourceSet() == false) {
      LOG(debug) << "ITS source not set";
    } else {
      LOG(debug) << "Kink ITS ref: " << kinkITSid.getIndex();
      LOG(debug) << "get source: " << int(kinkITSid.getSource());
    }

    auto iBinskink = mUtils.getBinRect(kink.getEta(), kink.getPhi(), mStrParams->mEtaBinSize, mStrParams->mPhiBinSize);
    LOG(debug) << "Kink phi: " << kink.getPhi() << ", eta: " << kink.getEta();

    for (int& iBinKink : iBinskink) {
      for (int iTrack{mTracksIdxTable[iBinKink]}; iTrack < TMath::Min(mTracksIdxTable[iBinKink + 1], int(mSortedITStracks.size())); iTrack++) {
        mITStrack = mSortedITStracks[iTrack];

        auto& ITSindexRef = mSortedITSindexes[iTrack];

<<<<<<< HEAD
        if (mITStrack.getCharge() != kink.getCharge())
          continue;

        for (int idx = 0; idx < ITSidx.size(); idx++) {
          if (ITSidx[idx] == ITSindexRef && ITSTPCidx[idx] == vrtxTrackIdx.getIndex()) {
            // match
            LOG(info) << "Match found!";
            LOG(info) << "Kink Vertexes:";
            LOG(info) << "ITS: " << mITSvtxBrackets[ITSindexRef].getMin() << " - " << mITSvtxBrackets[ITSindexRef].getMax();
            LOG(info) << "ITSTPC: " << kinkVrtxIDmin << " - " << kinkVrtxIDmax;
            LOG(info) << "Macro Vertexes:";
            LOG(info) << "ITS: " << firstIdx[idx] << " - " << lastIdx[idx];
            LOG(info) << "ITSTPC: " << firstIdxITSTPC[idx] << " - " << lastIdxITSTPC[idx];
          }
        }
=======
        mKinkTrack.likeSign = (mITStrack.getCharge() == kink.getCharge());

        if(mITStrack.getCharge() != mITStrack.getParamOut().getCharge())
          continue;
>>>>>>> ba066281d2e3d3164616381d8262de1f3971633a

        if (mITSvtxBrackets[ITSindexRef].getMax() < kinkVrtxIDmin || mITSvtxBrackets[ITSindexRef].getMin() > kinkVrtxIDmax) {
          continue;
        }

        if (matchKinkToITSTrack(kink)) {
          LOG(debug) << "----------------------";
          LOG(debug) << "ITS phi: " << mITStrack.getPhi() << ", eta: " << mITStrack.getEta();
          LOG(debug) << "ITS track ref: " << mSortedITSindexes[iTrack];
          LOG(debug) << "ITS Track matched with a Kink topology ....";
          LOG(debug) << "Matching chi2: " << mKinkTrack.mChi2Match;

          mKinkTrack.mNLayers = mITStrack.getNumberOfClusters();
          auto trackClusters = getTrackClusters();
          auto& lastClus = trackClusters[0];
          mKinkTrack.mChi2Match = getMatchingChi2(kink, mITStrack);

          // compute angle between kink and ITS track
          std::array<float, 3> kinkMom;
          kink.getPxPyPzGlo(kinkMom);
          auto kinkMomNorm = TMath::Sqrt(kinkMom[0] * kinkMom[0] + kinkMom[1] * kinkMom[1] + kinkMom[2] * kinkMom[2]);
          std::array<float, 3> itsMom;
          mITStrack.getPxPyPzGlo(itsMom);
          auto itsMomNorm = TMath::Sqrt(itsMom[0] * itsMom[0] + itsMom[1] * itsMom[1] + itsMom[2] * itsMom[2]);
          auto dotProd = kinkMom[0] * itsMom[0] + kinkMom[1] * itsMom[1] + kinkMom[2] * itsMom[2];
          auto angle = TMath::ACos(dotProd / (kinkMomNorm * itsMomNorm)) * TMath::RadToDeg();
          LOG(debug) << "Angle between kink and ITS track: " << angle;

          auto motherMass = calcKinkMotherMass(mKinkTrack.mMotherP, mKinkTrack.mDaughterP, PID::Triton, PID::PI0); // just the Hypertriton case for now
          mKinkTrack.mMasses[0] = motherMass;
<<<<<<< HEAD
=======
          mKinkTrack.mMasses[1] = calcKinkMotherMass(mKinkTrack.mMotherP, mKinkTrack.mDaughterP, PID::Pion, PID::Proton);  //Sigma-, Proton because there is not PID::Neutron
>>>>>>> ba066281d2e3d3164616381d8262de1f3971633a

          mKinkTrack.mTrackIdx = vrtxTrackIdx;
          mKinkTrack.mITSRef = ITSindexRef;
          mKinkTrackVec.push_back(mKinkTrack);
        }
      }
    }
  }
}

bool StrangenessTracker::matchKinkToITSTrack(o2::track::TrackParCovF daughterTrack)
{
  int nCand;

  try {
    nCand = mFitterKink.process(mITStrack.getParamOut(), daughterTrack);
  } catch (std::runtime_error& e) {
    LOG(debug) << "Fitterkink failed: " << e.what();
    return false;
  }

  if (!mFitterKink.propagateTracksToVertex() || !nCand)
    return false;

  auto chi2 = mFitterKink.getChi2AtPCACandidate();
  std::array<float, 3> R = mFitterKink.getPCACandidatePos();

  LOG(debug) << "Chi 2: " << chi2;
  if (chi2 < 0 || chi2 > mStrParams->mMaxChi2)
    return false;

  double recR = sqrt(R[0] * R[0] + R[1] * R[1]);
  LOG(debug) << "R: " << recR;
  if (recR < 18)
    return false;

  mKinkTrack.mChi2Vertex = chi2;
  mKinkTrack.mDecayVtx = R;

  mFitterKink.getTrack(0).getPxPyPzGlo(mKinkTrack.mMotherP);
  mFitterKink.getTrack(1).getPxPyPzGlo(mKinkTrack.mDaughterP);

  // computer cluster sizes
  std::vector<int> motherClusSizes;
  auto trackClusSizes = getTrackClusterSizes();
  for (int iClus{0}; iClus < trackClusSizes.size(); iClus++) {
    auto& compClus = trackClusSizes[iClus];
    motherClusSizes.push_back(compClus);
  }
  mKinkTrack.mITSClusSize = float(std::accumulate(motherClusSizes.begin(), motherClusSizes.end(), 0)) / motherClusSizes.size();

  return true;
}

bool StrangenessTracker::matchDecayToITStrack(float decayR)
{
  auto geom = o2::its::GeometryTGeo::Instance();
  auto trackClusters = getTrackClusters();
  auto trackClusSizes = getTrackClusterSizes();
  auto& lastClus = trackClusters[0];
  mStrangeTrack.mMatchChi2 = getMatchingChi2(mStrangeTrack.mMother, mITStrack);

  auto radTol = decayR < 4 ? mStrParams->mRadiusTolIB : mStrParams->mRadiusTolOB;
  auto nMinClusMother = trackClusters.size() < 4 ? 2 : mStrParams->mMinMotherClus;

  std::vector<ITSCluster> motherClusters;
  std::vector<int> motherClusSizes;
  std::array<unsigned int, 7> nAttachments;
  nAttachments.fill(-1); // fill arr with -1

  int nUpdates = 0;
  bool isMotherUpdated = false;

  for (int iClus{0}; iClus < trackClusters.size(); iClus++) {
    auto& clus = trackClusters[iClus];
    auto& compClus = trackClusSizes[iClus];
    int nUpdOld = nUpdates;
    double clusRad = sqrt(clus.getX() * clus.getX() - clus.getY() * clus.getY());
    auto diffR = decayR - clusRad;
    auto relDiffR = diffR / decayR;
    // Look for the Mother if the Decay radius allows for it, within a tolerance
    LOG(debug) << "decayR: " << decayR << ", diffR: " << diffR << ", clus rad: " << clusRad << ", radTol: " << radTol;
    if (relDiffR > -radTol) {
      LOG(debug) << "Try to attach cluster to Mother, layer: " << geom->getLayer(clus.getSensorID());
      if (updateTrack(clus, mStrangeTrack.mMother)) {
        motherClusters.push_back(clus);
        motherClusSizes.push_back(compClus);
        nAttachments[geom->getLayer(clus.getSensorID())] = 0;
        isMotherUpdated = true;
        nUpdates++;
        LOG(debug) << "Cluster attached to Mother";
        continue; // if the cluster is attached to the mother, skip the rest of the loop
      }
    }

    // if Mother is not found, check for V0 daughters compatibility
    if (relDiffR < radTol && !isMotherUpdated) {
      bool isDauUpdated = false;
      LOG(debug) << "Try to attach cluster to Daughters, layer: " << geom->getLayer(clus.getSensorID());
      for (int iDau{0}; iDau < mDaughterTracks.size(); iDau++) {
        auto& dauTrack = mDaughterTracks[iDau];
        if (updateTrack(clus, dauTrack)) {
          nAttachments[geom->getLayer(clus.getSensorID())] = iDau + 1;
          isDauUpdated = true;
          break;
        }
      }
      if (!isDauUpdated) {
        break; // no daughter track updated, stop the loop
      }
      nUpdates++;
    }
    if (nUpdates == nUpdOld) {
      break; // no track updated, stop the loop
    }
  }

  if (nUpdates < trackClusters.size() || motherClusters.size() < nMinClusMother) {
    return false;
  }

  o2::track::TrackParCov motherTrackClone = mStrangeTrack.mMother; // clone and reset covariance for final topology refit
  motherTrackClone.resetCovariance();

  LOG(debug) << "Clusters attached, starting inward-outward refit";

  std::reverse(motherClusters.begin(), motherClusters.end());

  for (auto& clus : motherClusters) {
    if (!updateTrack(clus, motherTrackClone)) {
      break;
    }
  }

  // compute mother average cluster size
  mStrangeTrack.mITSClusSize = float(std::accumulate(motherClusSizes.begin(), motherClusSizes.end(), 0)) / motherClusSizes.size();

  LOG(debug) << "Inward-outward refit finished, starting final topology refit";
  // final Topology refit

  int cand = 0; // best V0 candidate
  int nCand;

  // refit cascade
  if (mStrangeTrack.mPartType == dataformats::kStrkCascade) {
    V0 cascV0Upd;
    if (!recreateV0(mDaughterTracks[kV0DauPos], mDaughterTracks[kV0DauNeg], cascV0Upd)) {
      LOG(debug) << "Cascade V0 refit failed";
      return false;
    }
    try {
      nCand = mFitter3Body.process(cascV0Upd, mDaughterTracks[kBach], motherTrackClone);
    } catch (std::runtime_error& e) {
      LOG(debug) << "Fitter3Body failed: " << e.what();
      return false;
    }
    if (!nCand || !mFitter3Body.propagateTracksToVertex()) {
      LOG(debug) << "Fitter3Body failed: propagation to vertex failed";
      return false;
    }
  }

  // refit V0
  else if (mStrangeTrack.mPartType == dataformats::kStrkV0) {
    try {
      nCand = mFitter3Body.process(mDaughterTracks[kV0DauPos], mDaughterTracks[kV0DauNeg], motherTrackClone);
    } catch (std::runtime_error& e) {
      LOG(debug) << "Fitter3Body failed: " << e.what();
      return false;
    }
    if (!nCand || !mFitter3Body.propagateTracksToVertex()) {
      LOG(debug) << "Fitter3Body failed: propagation to vertex failed";
      return false;
    }
  }

  mStrangeTrack.mDecayVtx = mFitter3Body.getPCACandidatePos();
  mStrangeTrack.mTopoChi2 = mFitter3Body.getChi2AtPCACandidate();
  mStructClus.arr = nAttachments;

  return true;
}

bool StrangenessTracker::updateTrack(const ITSCluster& clus, o2::track::TrackParCov& track)
{
  auto geom = o2::its::GeometryTGeo::Instance();
  auto propInstance = o2::base::Propagator::Instance();
  float alpha = geom->getSensorRefAlpha(clus.getSensorID()), x = clus.getX();
  int layer{geom->getLayer(clus.getSensorID())};

  if (!track.rotate(alpha)) {
    return false;
  }

  if (!propInstance->propagateToX(track, x, getBz(), o2::base::PropagatorImpl<float>::MAX_SIN_PHI, o2::base::PropagatorImpl<float>::MAX_STEP, mCorrType)) {
    return false;
  }

  if (mCorrType == o2::base::PropagatorF::MatCorrType::USEMatCorrNONE) {
    float thick = layer < 3 ? 0.005 : 0.01;
    constexpr float radl = 9.36f; // Radiation length of Si [cm]
    constexpr float rho = 2.33f;  // Density of Si [g/cm^3]
    if (!track.correctForMaterial(thick, thick * rho * radl)) {
      return false;
    }
  }
  auto chi2 = std::abs(track.getPredictedChi2(clus)); // abs to be understood
  LOG(debug) << "Chi2: " << chi2;
  if (chi2 > mStrParams->mMaxChi2 || chi2 < 0) {
    return false;
  }

  if (!track.update(clus)) {
    return false;
  }

  return true;
}

} // namespace strangeness_tracking
} // namespace o2