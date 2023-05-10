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

std::vector<int> ITSidx = {7748, 11802, 1084, 1304, 1304, 6719, 848, 1842, 7963, 11535, 11630, 10321, 10664, 5063, 7923, 4931, 4466, 11529, 1846, 8370, 5151, 5439, 5433, 12036, 418, 12747, 722, 7691, 11942, 4049, 2244, 8058, 9391, 8904, 7671, 7499, 3865, 9330, 676, 2188, 2188, 2387, 270, 8910, 8434, 2687, 8022, 9139, 7329, 11362, 12406, 8596, 1558, 2241, 9677, 3121, 9197, 9460, 2838, 984, 1355, 11578, 5283, 2663, 11584, 22, 3911, 2448, 5279, 2887, 796, 2289, 2619, 891, 9154, 1709, 7513, 783, 1152, 4613, 7743, 2058, 9560, 11705, 10649, 9852, 936, 11027, 8333, 10205, 10854, 6262, 6671, 1926, 1136, 10564, 8771, 2228, 11192, 5777, 5203, 7401, 9137, 5265, 7482, 9075, 3890, 11869, 5062, 2764, 3916, 8507, 7254, 5086, 6921, 8672, 10338, 1415};
std::vector<int> ITSTPCidx = {100109, 108229, 97511, 97731, 97731, 103146, 91152, 93847, 99968, 103540, 103635, 97789, 98132, 92531, 95391, 92399, 99081, 106058, 96375, 102899, 98238, 98526, 98520, 105123, 94105, 103867, 91842, 98811, 103727, 95834, 97734, 103562, 98241, 97754, 96521, 96349, 92715, 101056, 90870, 92382, 92382, 92581, 90464, 107400, 106924, 101177, 106512, 107629, 105819, 105898, 106942, 103132, 96094, 96777, 97884, 94710, 100786, 101049, 92617, 94182, 94553, 106352, 100057, 97437, 102940, 91378, 95267, 93804, 97677, 95285, 94456, 95949, 96279, 94551, 102674, 95229, 100501, 93771, 94140, 97601, 100731, 94845, 102347, 104492, 103436, 101776, 93713, 103804, 96776, 98648, 99297, 99165, 98271, 93526, 94554, 103982, 102189, 92622, 101885, 96470, 96078, 97388, 99124, 99196, 101413, 103006, 97821, 105800, 92396, 90098, 92486, 97077, 97994, 95952, 97787, 102645, 104311, 95388};
std::vector<int> firstIdxITS = {1718, 767, 2383, 3804, 3807, 4948, 4004, 428, 1122, 4428, 4813, 1301, 2415, 3034, 3486, 4788, 4218, 959, 2251, 3020, 624, 3292, 3599, 4708, 5240, 840, 1502, 3156, 1513, 3906, 4828, 5484, 681, 969, 1014, 2621, 3151, 877, 557, 638, 638, 3669, 4304, 317, 1407, 3778, 3924, 5560, 5624, 451, 451, 1158, 1961, 2305, 1928, 903, 2252, 3544, 3370, 2735, 2799, 153, 881, 4048, 677, 2083, 2132, 4412, 1230, 2842, 216, 2421, 3840, 5356, 433, 4150, 852, 1644, 4398, 5292, 5416, 73, 1944, 3000, 4818, 2369, 1058, 4872, 633, 1187, 1254, 3108, 1343, 3040, 665, 2018, 5199, 5161, 3181, 4412, 5174, 931, 4502, 847, 1612, 3070, 3843, 5177, 2344, 3652, 3676, 4855, 1386, 414, 2707, 391, 5231, 5392};
std::vector<int> lastIdxITS = {1722, 767, 2383, 3804, 3807, 4952, 4004, 433, 1122, 4435, 4813, 1301, 2415, 3034, 3489, 4788, 4218, 961, 2251, 3020, 624, 3300, 3599, 4708, 5240, 840, 1502, 3156, 1513, 3906, 4828, 5484, 681, 969, 1019, 2621, 3151, 881, 557, 643, 638, 3669, 4304, 317, 1407, 3778, 3929, 5560, 5624, 451, 451, 1158, 1961, 2305, 1928, 903, 2252, 3544, 3370, 2741, 2799, 162, 881, 4048, 677, 2083, 2132, 4412, 1230, 2842, 216, 2421, 3840, 5356, 433, 4150, 852, 1644, 4403, 5297, 5416, 73, 1950, 3000, 4821, 2377, 1058, 4872, 633, 1193, 1254, 3108, 1343, 3054, 679, 2018, 5199, 5167, 3181, 4412, 5183, 931, 4502, 847, 1612, 3075, 3843, 5177, 2344, 3659, 3676, 4855, 1386, 414, 2707, 391, 5231, 5392};
std::vector<int> firstIdxITSTPC = {1719, 763, 2377, 3801, 3801, 4942, 4002, 430, 1122, 4432, 4813, 1296, 2413, 3030, 3486, 4786, 4216, 959, 2250, 3015, 623, 3290, 3597, 4706, 5236, 840, 1497, 3153, 1512, 3904, 4827, 5483, 680, 967, 1014, 2618, 3148, 878, 554, 638, 638, 3669, 4301, 313, 1405, 3778, 3924, 5559, 5624, 447, 447, 1156, 1959, 2301, 1928, 896, 2247, 3544, 3370, 2729, 2799, 149, 878, 4044, 675, 2080, 2132, 4411, 1229, 2839, 205, 2418, 3838, 5352, 431, 4146, 850, 1643, 4398, 5294, 5415, 73, 1947, 2999, 4818, 2372, 1056, 4869, 630, 1188, 1251, 3107, 1342, 3043, 671, 2017, 5198, 5165, 3181, 4409, 5174, 929, 4499, 844, 1611, 3071, 3841, 5175, 2340, 3653, 3672, 4855, 1383, 411, 2698, 389, 5225, 5392};
std::vector<int> lastIdxITSTPC = {1722, 766, 2384, 3806, 3806, 4947, 4007, 433, 1123, 4433, 4816, 1301, 2415, 3032, 3489, 4791, 4222, 960, 2252, 3024, 626, 3294, 3602, 4707, 5244, 841, 1502, 3158, 1517, 3914, 4831, 5485, 681, 969, 1019, 2619, 3154, 880, 559, 639, 639, 3675, 4305, 317, 1407, 3786, 3929, 5562, 5626, 452, 452, 1159, 1966, 2303, 1929, 902, 2249, 3544, 3370, 2733, 2805, 153, 881, 4048, 679, 2083, 2132, 4417, 1233, 2845, 210, 2424, 3842, 5356, 434, 4151, 852, 1649, 4404, 5297, 5419, 75, 1951, 3001, 4821, 2377, 1060, 4874, 634, 1193, 1256, 3112, 1345, 3048, 676, 2021, 5200, 5167, 3181, 4415, 5179, 929, 4506, 849, 1617, 3074, 3843, 5181, 2346, 3659, 3675, 4856, 1386, 418, 2702, 395, 5231, 5396};
std::vector<int> timeFrame = {0, 2, 2, 2, 2, 2, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 7, 7, 7, 8, 8, 8, 8, 9, 11, 11, 11, 12, 12, 13, 14, 15, 15, 15, 15, 15, 16, 17, 17, 17, 17, 17, 18, 18, 18, 18, 18, 18, 19, 19, 19, 19, 19, 20, 21, 21, 21, 22, 23, 23, 24, 24, 24, 25, 25, 25, 25, 26, 26, 27, 27, 27, 27, 28, 28, 29, 29, 29, 29, 29, 30, 30, 30, 30, 31, 32, 32, 33, 33, 33, 35, 36, 36, 37, 37, 37, 38, 39, 39, 40, 41, 41, 42, 42, 42, 42, 42, 43, 43, 45, 45, 46, 48, 48, 49, 49, 49};

//int TFcount = 0;

bool StrangenessTracker::loadData(const o2::globaltracking::RecoContainer& recoData)
{
  //TFcount++;
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

  // build time bracket for each ITS and kink track
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
        auto itstpcidx = kinkHelper.index;
        for(int i = 0; i< ITSTPCidx.size(); i++){
          //if(ITSTPCidx[i] == itstpcidx.getIndex() && TFcount == timeFrame[i]){
            if(ITSTPCidx[i] == itstpcidx.getIndex()){
            LOG(info) <<"ITSTPC track : " << itstpcidx.getIndex() << " with source " << itstpcidx.getSourceName();
            LOG(info) <<"Momentum : " << kinkHelper.track.getP();
            LOG(info) <<"TF : " << timeFrame[i];
          }
        }
        
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

        if (mStrParams->mVertexMatching && (mITSvtxBrackets[ITSindexRef].getMin() > v0.getVertexID() ||
                                            mITSvtxBrackets[ITSindexRef].getMax() < v0.getVertexID())) {
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
          auto p2pos = mFitter3Body.getTrack(kV0DauPos).getP2();                                     // positive V0 daughter
          auto p2neg = mFitter3Body.getTrack(kV0DauNeg).getP2();                                     // negative V0 daughter
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

        if (mStrParams->mVertexMatching && (mITSvtxBrackets[ITSindexRef].getMin() > casc.getVertexID() ||
                                            mITSvtxBrackets[ITSindexRef].getMax() < casc.getVertexID())) {
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
  if (!mStrParams->mKinkFinder) {
    return;
  }
  // loop over Kinks

  LOG(debug) << "Analysing " << mKinkTracks.size() << " Kinks ...";
  LOG(info) << "There are " << mKinkTracks.size() << " kinks in the event";

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

        for (int idx = 0; idx < ITSidx.size(); idx++) {
          if (ITSidx[idx] == ITSindexRef && ITSTPCidx[idx] == vrtxTrackIdx.getIndex()) {
            // match
            LOG(info) << "Match found!";
            LOG(info) << "Track info: ITS " << ITSidx[idx] << " ITSTPC " << ITSTPCidx[idx] << " Source " << vrtxTrackIdx.getSourceName();
            LOG(info) << "Kink Vertexes:";
            LOG(info) << "ITS: " << mITSvtxBrackets[ITSindexRef].getMin() << " - " << mITSvtxBrackets[ITSindexRef].getMax();
            LOG(info) << "ITSTPC: " << kinkVrtxIDmin << " - " << kinkVrtxIDmax;
            LOG(info) << "Macro Vertexes:";
            LOG(info) << "ITS: " << firstIdxITS[idx] << " - " << lastIdxITS[idx];
            LOG(info) << "ITSTPC: " << firstIdxITSTPC[idx] << " - " << lastIdxITSTPC[idx];
          }
        }

        if (mITStrack.getCharge() != kink.getCharge())
          continue;

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

          mKinkTrack.mTrackIdx = vrtxTrackIdx;
          mKinkTrack.mITSRef = mSortedITSindexes[iTrack];
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