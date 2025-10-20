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

#include "TRDQC/RawDataManager.h"

#include <RtypesCore.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>

#include <iterator>
#include <set>
#include <utility>
#include <boost/range/distance.hpp>
#include <boost/range/iterator_range_core.hpp>

#include "TRDQC/CoordinateTransformer.h"
#include "TRDBase/Geometry.h"
#include "TRDBase/GeometryFlat.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "Framework/Logger.h"

using namespace o2::trd;
constexpr bool debugprint = false;
/// comparison function to order digits by det / row / MCM / -channel
bool comp_digit(const o2::trd::Digit& a, const o2::trd::Digit& b)
{
  if (a.getDetector() != b.getDetector()) {
    return a.getDetector() < b.getDetector();
  }

  if (a.getPadRow() != b.getPadRow()) {
    return a.getPadRow() < b.getPadRow();
  }

  if (a.getROB() != b.getROB()) {
    return a.getROB() < b.getROB();
  }

  if (a.getMCM() != b.getMCM()) {
    return a.getMCM() < b.getMCM();
  }

  // sort channels in descending order, to ensure ordering of pad columns
  if (a.getChannel() != b.getChannel()) {
    return a.getChannel() > b.getChannel();
  }

  return true;
}

/// comparison function to order tracklets by det / row / MCM / channel
bool comp_tracklet(const o2::trd::Tracklet64& a, const o2::trd::Tracklet64& b)
{
  // upper bits of hcid and padrow from Tracklet64 word
  const uint64_t det_row_mask = 0x0ffde00000000000;

  // lowest bit of hcid (side), MCM col and pos from Tracklet64 word
  const uint64_t col_pos_mask = 0x00011fff00000000;

  auto a_det_row = a.getTrackletWord() & det_row_mask;
  auto b_det_row = b.getTrackletWord() & det_row_mask;

  if (a_det_row != b_det_row) {
    return a_det_row < b_det_row;
  }

  auto a_col_pos = a.getTrackletWord() & col_pos_mask;
  auto b_col_pos = b.getTrackletWord() & col_pos_mask;

  return a_col_pos < b_col_pos;
};

bool comp_spacepoint(const ChamberSpacePoint& a, const ChamberSpacePoint& b)
{
  if (a.getDetector() != b.getDetector()) {
    return a.getDetector() < b.getDetector();
  }

  if (a.getPadRow() != b.getPadRow()) {
    return a.getPadRow() < b.getPadRow();
  }

  if (a.getPadCol() != b.getPadCol()) {
    return a.getPadCol() < b.getPadCol();
  }
  return true;
}

bool comp_tracksegments(const TrackSegment& a, const TrackSegment& b)
{
  if(a.getTriggerTime() != b.getTriggerTime()){
    return a.getTriggerTime() < b.getTriggerTime();
  }
  if (a.getDetector() != b.getDetector()) {
    return a.getDetector() < b.getDetector();
  }

  if (a.getPadRow() != b.getPadRow()) {
    return a.getPadRow() < b.getPadRow();
  }

  if (a.getPadCol() != b.getPadCol()) {
    return a.getPadCol() < b.getPadCol();
  }

  return true;
}

void RawDataSpan::sort()
{
  std::stable_sort(std::begin(digits), std::end(digits), comp_digit);
  std::stable_sort(std::begin(tracklets), std::end(tracklets), comp_tracklet);
  std::stable_sort(std::begin(hits), std::end(hits), comp_spacepoint);
  // how to handle multiple mcm per track ?
  // do i do 6 tracks and 1 per layer ?
  std::stable_sort(std::begin(tracks_itstpc_seg), std::end(tracks_itstpc_seg), comp_tracksegments);
}

template <typename keyfunc>
std::vector<RawDataSpan> RawDataSpan::iterateBy()
{
  // an map for keeping track which ranges correspond to which key
  std::map<uint32_t, RawDataSpan> spanmap;

  // sort digits and tracklets
  sort();

  // add all the digits to a map
  for (auto cur = digits.begin(); cur != digits.end(); /* noop */) {
    // calculate the key of the current (first unprocessed) digit
    auto key = keyfunc::key(*cur);
    // find the first digit with a different key
    auto nxt = std::find_if(cur, digits.end(), [key](auto x) { return keyfunc::key(x) != key; });
    // store the range cur:nxt in the map
    spanmap[key].digits = boost::make_iterator_range(cur, nxt);
    // continue after this range
    cur = nxt;
  }

  // add tracklets to the map
  for (auto cur = tracklets.begin(); cur != tracklets.end(); /* noop */) {
    auto key = keyfunc::key(*cur);
    auto nxt = std::find_if(cur, tracklets.end(), [key](auto x) { return keyfunc::key(x) != key; });
    spanmap[key].tracklets = boost::make_iterator_range(cur, nxt);
    cur = nxt;
  }

  // spanmap contains all TRD data - either digits or tracklets. Now we insert hit information into these spans. The
  // tricky part is that space points or hits can belong to more than one MCM, i.e. they could appear in two spans.
  // We keep the begin iterator for each key in a map
  std::map<uint32_t, std::vector<HitPoint>::iterator> firsthit;
  for (auto cur = hits.begin(); cur != hits.end(); ++cur) {
    // calculate the keys for this hit
    auto keys = keyfunc::keys(*cur);
    // if we are not yet aware of this key, register the current hit as the first hit
    for (auto key : keys) {
      firsthit.insert({key, cur});
    }
    // remote the keys from the firsthit map that are no longer found in the hits
    for (auto it = firsthit.cbegin(); it != firsthit.cend(); /* no increment */) {
      if (keys.find(it->first) == keys.end()) {
        spanmap[it->first].hits = boost::make_iterator_range(it->second, cur);
        it = firsthit.erase(it);
      } else {
        ++it;
      }
    }
  }
  // now for tracks .....
  //
  //  spanmap contains all TRD data - either digits or tracklets. Now we insert tracking information into these spans. The
  //  tricky part is that space points or hits can belong to more than one MCM, i.e. they could appear in two spans.
  //  We keep the begin iterator for each key in a map
  //  ITSTPC, ITSTPCTRD first.
  //  1. track all the tracks through the detector, finding the detector, and the entry and exit points of the drift and ionisation regions, 3 points, and the curvature.
  //  2. store said points as a structure TrackSegment.
  // std::vector<TrackSegment> mITSTPCTracks_segments{0};
  // std::map<uint32_t, std::vector<o2::dataformats::TrackTPCITS>::iterator> firsttrack_itstpc;
  std::map<uint32_t, std::vector<TrackSegment>::iterator> firsttrack_itstpc_segment;
  std::map<uint32_t, std::vector<TrackSegment>::iterator> firsttrack_itstpctrd_segment;
  // std::map<uint32_t, std::vector<o2::trd::TrackTRD>::iterator> firsttrack_itstpctrd;
  for (auto cur = tracks_itstpc_seg.begin(); cur != tracks_itstpc_seg.end(); ++cur) {
    // calculate the keys for this track
    auto keys = keyfunc::keys(*cur);
    // if we are not yet aware of this key, register the current track as the first track
    for (auto key : keys) {
      firsttrack_itstpc_segment.insert({key, cur});
    }
    // remove the keys from the firsthit map that are no longer found in the hits
    for (auto it = firsttrack_itstpc_segment.begin(); it != firsttrack_itstpc_segment.cend(); /* no increment */) {
      if (keys.find(it->first) == keys.end()) {
        spanmap[it->first].tracks_itstpc_seg = boost::make_iterator_range(it->second, cur);
        it = firsttrack_itstpc_segment.erase(it);
      } else {
        ++it;
      }
    }
  }

  // convert the map of spans into a vector, as we do not need the access by key
  // any longer, and having a vector makes the looping by the user easier.
  std::vector<RawDataSpan> spans;
  transform(spanmap.begin(), spanmap.end(), back_inserter(spans), [](auto const& pair) { return pair.second; });

  return spans;
}

/// PadRowID is a struct to calculate unique identifiers per pad row.
/// The struct can be passed as a template parameter to the RawDataSpan::IterateBy
/// method to split the data span by pad row and iterate over the pad rows.
struct PADROW_ID {
  /// The static `key` method calculates a padrow ID for digits and tracklets
  template <typename T>
  static uint32_t key(const T& x)
  {
    return 100 * x.getDetector() + x.getPadRow();
  }

  static std::set<uint32_t> keys(const o2::trd::ChamberSpacePoint& x)
  {
    uint32_t key = 100 * x.getDetector() + x.getPadRow();
    return {key};
  }

  static std::set<uint32_t> keys(const o2::trd::TrackSegment& x)
  {
    uint32_t key = 100 * x.getDetector() + x.getStartPoint().getPadRow();
    return {key};
  }

  static bool match(const uint32_t key, const o2::trd::ChamberSpacePoint& x)
  {
    return key == 100 * x.getDetector() + x.getPadRow();
  }
  static bool match(const uint32_t key, const o2::trd::TrackSegment& x)
  {
    return key == 100 * x.getDetector() + x.getStartPoint().getPadRow();
  }
  static int getDetector(uint32_t k) { return k / 1000; }
  // static int getPadRow(key) {return (key%1000) / 8;}
  static int getMcmRowCol(uint32_t k) { return k % 1000; }
};

// instantiate the template to iterate by padrow
template std::vector<RawDataSpan> RawDataSpan::iterateBy<PADROW_ID>();
// non-template wrapper function to keep PadRowID within the .cxx file
std::vector<RawDataSpan> RawDataSpan::iterateByPadRow() { return iterateBy<PADROW_ID>(); }

/// A struct that can be used to calculate unique identifiers for MCMs, to be
/// used to split ranges by MCM.
struct MCM_ID {
  template <typename T>
  static uint32_t key(const T& x)
  {
    return 1000 * x.getDetector() + 8 * x.getPadRow() + 4 * (x.getROB() % 2) + x.getMCM() % 4;
  }

  static std::set<uint32_t> keys(const o2::trd::ChamberSpacePoint& x)
  {
    uint32_t detrow = 1000 * x.getDetector() + 8 * x.getPadRow();
    uint32_t mcmcol = uint32_t(x.getPadCol() / float(o2::trd::constants::NCOLMCM));

    // float c = x.getPadCol() - float(mcmcol * o2::trd::constants::NCOLMCM);
    float c = x.getMCMChannel(mcmcol);

    if (c >= 19.0 && mcmcol >= 1) {
      return {detrow + mcmcol - 1, detrow + mcmcol};
    } else if (c <= 1.0 && mcmcol <= 6) {
      return {detrow + mcmcol, detrow + mcmcol + 1};
    } else {
      return {detrow + mcmcol};
    }
  }
  static std::set<uint32_t> keys(const o2::trd::TrackSegment& x)
  {
    return keys(x.getStartPoint());
  }

  static int getDetector(uint32_t k) { return k / 1000; }
  // static int getPadRow(key) {return (key%1000) / 8;}
  static int getMcmRowCol(uint32_t k) { return k % 1000; }
};

// template instantion and non-template wrapper function
template std::vector<RawDataSpan> RawDataSpan::iterateBy<MCM_ID>();
std::vector<RawDataSpan> RawDataSpan::iterateByMCM() { return iterateBy<MCM_ID>(); }

struct DET_ID {
  template <typename T>
  static uint32_t key(const T& x)
  {
    return x.getDetector();
  }

  static std::set<uint32_t> keys(const o2::trd::ChamberSpacePoint& x)
  {
    uint32_t det = x.getDetector();
    return {det};
  }
  static std::set<uint32_t> keys(const o2::trd::TrackSegment& x)
  {
    return keys(x.getStartPoint());
  }
  /*
     static bool match(const uint32_t key, const o2::trd::ChamberSpacePoint& x)
     {
       return key == x.getDetector();
     }*/
};

template std::vector<RawDataSpan> RawDataSpan::iterateBy<DET_ID>();
std::vector<RawDataSpan> RawDataSpan::iterateByDetector() { return iterateBy<DET_ID>(); }

std::vector<TrackSegment> RawDataSpan::makeMCTrackSegments()
{
  // define a struct to keep track of the first and last MC hit of one track in one chamber
  struct SegmentInfo {
    // The first hit is the hit closest to the anode region, i.e. with the largest x coordinate.
    size_t firsthit{0}; //
    // The last hit is the hit closest to the radiator, i.e. with the smallest x coordinate.
    size_t lasthit{0};
    float start{-999.9}; // local x cordinate of the first hit, init value ensures any hit updates
    float end{999.9};    // local x cordinate of the last hit, init value ensures any hit updates
  };
  // Keep information about found track segments in a map indexed by track ID and detector number.
  // If the span only covers (part of) a detector, the detector information is redundant, but in
  // the case of processing a whole event, the distinction by detector will be needed.
  std::map<std::pair<int, int>, SegmentInfo> trackSegmentInfo;

  for (int iHit = 0; iHit < hits.size(); ++iHit) {
    auto hit = hits[iHit];

    // in the following, we will look for track segments using hits in the drift region
    if (hit.isFromDriftRegion()) {
      // The first hit is the hit closest to the anode region, i.e. with the largest x coordinate.
      auto id = std::make_pair(hit.getID(), hit.getDetector());
      if (hit.getX() > trackSegmentInfo[id].start) {
        trackSegmentInfo[id].firsthit = iHit;
        trackSegmentInfo[id].start = hit.getX();
      }
      // The last hit is the hit closest to the radiator, i.e. with the smallest x coordinate.
      if (hit.getX() < trackSegmentInfo[id].end) {
        trackSegmentInfo[id].lasthit = iHit;
        trackSegmentInfo[id].end = hit.getX();
      }
    }
  } // hit loop

  std::vector<TrackSegment> trackSegments;
  for (auto x : trackSegmentInfo) {
    auto trackid = x.first.first;
    auto detector = x.first.second;
    auto firsthit = hits[x.second.firsthit];
    auto lasthit = hits[x.second.lasthit];
    trackSegments.emplace_back(firsthit, lasthit, trackid);
  }
  return trackSegments;
}

/// The RawDataManager constructor: connects all data files and sets up trees, readers etc.
RawDataManager::RawDataManager(std::filesystem::path dir)
{

  if (!std::filesystem::exists(dir) || !std::filesystem::is_directory(dir)) {
    O2ERROR("'%s' is not a directory", dir.c_str());
    return;
  }

  // We allways need the trigger records, which are stored in trdtracklets.root.
  // While at it, let's also set up reading the tracklets.
  if (!std::filesystem::exists(dir / "trdtracklets.root")) {
    O2ERROR("'tracklets.root' not found in directory '%s'", dir.c_str());
    return;
  }

  mMainFile = new TFile((dir / "trdtracklets.root").c_str());
  mMainFile->GetObject("o2sim", mDataTree);

  // set up the branches we want to read
  mDataTree->SetBranchAddress("Tracklet", &mTracklets);
  mDataTree->SetBranchAddress("TrackTrg", &mTrgRecords);

  if (std::filesystem::exists(dir / "trddigits.root")) {
    mDataTree->AddFriend("o2sim", (dir / "trddigits.root").c_str());
    mDataTree->SetBranchAddress("TRDDigit", &mDigits);
  }

  if (std::filesystem::exists(dir / "o2match_itstpc.root")) {
    mDataTree->AddFriend("matchTPCITS", (dir / "o2match_itstpc.root").c_str());
    mDataTree->SetBranchAddress("TPCITS", &mITSTPCTracks);
  } else {
    LOGP(info, " o2match_itstpc.root not found!");
  }
  if (std::filesystem::exists(dir / "o2match_tof_itstpc.root")) {
    mDataTree->AddFriend("matchTOF", (dir / "o2match_tof_itstpc.root").c_str());
    mDataTree->SetBranchAddress("TOFMatchInfo", &mITSTPCTracks_TOF);
  }
  if (std::filesystem::exists(dir / "o2match_tof_itstpctrd.root")) {
    mDataTree->AddFriend("matchTOF", (dir / "o2match_tof_itstpctrd.root").c_str());
    mDataTree->SetBranchAddress("TOFMatchInfo", &mITSTPCTRDTracks_TOF);
  }
  if (std::filesystem::exists(dir / "trdmatches_itstpc.root")) {
    mDataTree->AddFriend("tracksTRD", (dir / "trdmatches_itstpc.root").c_str());
    mDataTree->SetBranchAddress("tracks", &mITSTPCTracks_TRD);
    mDataTree->SetBranchAddress("trgrec", &mITSTPCTracksTrigRec_TRD);
  }
  if (std::filesystem::exists(dir / "trdmatches_itstpctof.root")) {
    mDataTree->AddFriend("tracksTRD", (dir / "trdmatches_itstpctof.root").c_str());
    mDataTree->SetBranchAddress("tracks", &mITSTPCTOFTracks_TRD);
    mDataTree->SetBranchAddress("trgrec", &mITSTPCTOFTracksTrigRec_TRD);
  }
  if (std::filesystem::exists(dir / "trdmatches_tpc.root")) {
    mDataTree->AddFriend("tracksTRD", (dir / "trdmatches_tpc.root").c_str());
    mDataTree->SetBranchAddress("tracks", &mTPCTracks_TRD);
    mDataTree->SetBranchAddress("trgrec", &mTPCTracksTrigRec_TRD);
  }

  // For data, we need info about time frames to match ITS and TPC tracks to trigger records.
  if (std::filesystem::exists(dir / "o2_tfidinfo.root")) {
    TFile fInTFID((dir / "o2_tfidinfo.root").c_str());
    mTFIDs = (std::vector<o2::dataformats::TFIDInfo>*)fInTFID.Get("tfidinfo");
  }

  // For MC, we first read the collision context
  if (std::filesystem::exists(dir / "collisioncontext.root")) {
    TFile fInCollCtx((dir / "collisioncontext.root").c_str());
    mCollisionContext = (o2::steer::DigitizationContext*)fInCollCtx.Get("DigitizationContext");
    // mCollisionContext->printCollisionSummary();
  }

  // We create the MC TTree using event header and tracks from the kinematics file
  if (std::filesystem::exists(dir / "o2sim_Kine.root")) {
    mMCFile = new TFile((dir / "o2sim_Kine.root").c_str());
    mMCFile->GetObject("o2sim", mMCTree);
    mMCTree->SetBranchAddress("MCEventHeader.", &mMCEventHeader);
    mMCTree->SetBranchAddress("MCTrack", &mMCTracks);
  }

  // We then add the TRD hits to the MC tree
  if (mMCFile && std::filesystem::exists(dir / "o2sim_HitsTRD.root")) {
    mMCTree->AddFriend("o2sim", (dir / "o2sim_HitsTRD.root").c_str());
    mMCTree->SetBranchAddress("TRDHit", &mHits);
  }
  // instantiate the class that handles all the data access
  //--------------------------------
  // load up the geometry
  // //TODO move this to only run once
  mGeo = o2::trd::Geometry::instance();
  //  auto mGeogpu = GPUTRDGeometry();

  //-------- init geometry and field --------//
  o2::base::GeometryManager::loadGeometry();
  o2::base::Propagator::initFieldFromGRP();
  mProp = o2::base::Propagator::Instance();
  // auto geo = o2::trd::Geometry::instance();
  mGeo->createPadPlaneArray();
  mGeo->createClusterMatrixArray();
  const o2::trd::GeometryFlat geoFlat(*mGeo);

  /*
    auto& elParam = o2::tpc::ParameterElectronics::Instance();
    auto& gasParam = o2::tpc::ParameterGas::Instance();
    auto& detParam = o2::tpc::ParameterDetector::Instance();
    auto tpcTBinMUS = elParam.ZbinWidth;
    auto tpcVdrift = gasParam.DriftV;
    auto tpcTOffset = detParam.DriftTimeOffset;
   */

  // GPUTRDGeometry mGeogpu;     // TRD geometry
  if (mGeo == nullptr) {
    LOGP(error, " mGeo is null");
  } else {
    LOGP(info, " mGeo is not null");
  }
  //
  // obtain average radius of TRD chambers
  o2::math_utils::Transform3D matrix = mGeo->getMatrixT2L(0);
  // o2::math_utilgpu::Transform3D matrix = geo->getMatrixT2L(0);
  // double loc[3] = {geo->anodePos(), 0.f, 0.f};
  // double glb[3] = {0.f, 0.f, 0.f};
  std::array<double, 3> loc = {mGeo->anodePos(), 0.f, 0.f};
  std::array<double, 3> driftstart = {mGeo->anodePos() - mGeo->camHght() / 2 - mGeo->craHght(), 0.f, 0.f}; // position of the start of the drift region
  std::array<double, 3> driftend = {mGeo->anodePos() - mGeo->camHght() / 2, 0.f, 0.f};                     // position of the end of the drift region and start of the amplification region
  std::array<double, 3> amplificationend = {mGeo->anodePos() + mGeo->camHght() / 2, 0.f, 0.f};             // position of the end of amplification region
  std::array<double, 3> glb = {0.f, 0.f, 0.f};
  std::array<double, 3> glbdriftstart = {0.f, 0.f, 0.f};
  std::array<double, 3> glbdriftend = {0.f, 0.f, 0.f};
  std::array<double, 3> glbamplificationend = {0.f, 0.f, 0.f};
  for (int32_t iDet = 0; iDet < o2::trd::constants::NCHAMBER; ++iDet) {
    // matrix = mGeo->getMatrixT2L(iDet);
    // matrix.print();
    // ApplyInverse(loc,glb,matrix);
    std::array<double, 12> m{0};
    mGeo->getMatrixT2L(iDet).GetComponents(m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8], m[9], m[10], m[11]);
    enum Transform3DMatrixIndex { kXX = 0,
                                  kXY = 1,
                                  kXZ = 2,
                                  kDX = 3,
                                  kYX = 4,
                                  kYY = 5,
                                  kYZ = 6,
                                  kDY = 7,
                                  kZX = 8,
                                  kZY = 9,
                                  kZZ = 10,
                                  kDZ = 11 };
    const double tmp[3] = {loc[0] - m[kDX], loc[1] - m[kDY], loc[2] - m[kDZ]};
    glb[0] = m[kXX] * tmp[0] + m[kYX] * tmp[1] + m[kZX] * tmp[2];
    glb[1] = m[kXY] * tmp[0] + m[kYY] * tmp[1] + m[kZY] * tmp[2];
    glb[2] = m[kXZ] * tmp[0] + m[kYZ] * tmp[1] + m[kZZ] * tmp[2];
    const double tmpds[3] = {driftstart[0] - m[kDX], driftstart[1] - m[kDY], driftstart[2] - m[kDZ]};
    glbdriftstart[0] = m[kXX] * tmpds[0] + m[kYX] * tmpds[1] + m[kZX] * tmpds[2];
    glbdriftstart[1] = m[kXY] * tmpds[0] + m[kYY] * tmpds[1] + m[kZY] * tmpds[2];
    glbdriftstart[2] = m[kXZ] * tmpds[0] + m[kYZ] * tmpds[1] + m[kZZ] * tmpds[2];
    const double tmpde[3] = {loc[0] - m[kDX], loc[1] - m[kDY], loc[2] - m[kDZ]};
    glbdriftend[0] = m[kXX] * tmpde[0] + m[kYX] * tmpde[1] + m[kZX] * tmpde[2];
    glbdriftend[1] = m[kXY] * tmpde[0] + m[kYY] * tmpde[1] + m[kZY] * tmpde[2];
    glbdriftend[2] = m[kXZ] * tmpde[0] + m[kYZ] * tmpde[1] + m[kZZ] * tmpde[2];
    const double tmpae[3] = {amplificationend[0] - m[kDX], amplificationend[1] - m[kDY], amplificationend[2] - m[kDZ]};
    glbamplificationend[0] = m[kXX] * tmpae[0] + m[kYX] * tmpae[1] + m[kZX] * tmpae[2];
    glbamplificationend[1] = m[kXY] * tmpae[0] + m[kYY] * tmpae[1] + m[kZY] * tmpae[2];
    glbamplificationend[2] = m[kXZ] * tmpae[0] + m[kYZ] * tmpae[1] + m[kZZ] * tmpae[2];
    // matrix.LocalToMaster(loc, glb);
    // LOGP(info,"{} ---> {} default {} {} {} {} {} {}",loc[0],glb[0],300.2f, 312.8f, 325.4f, 338.0f, 350.6f, 363.2f);
    mR[iDet] = glb[0];
    mRdriftstart[iDet] = glbdriftstart[0];
    mRdriftend[iDet] = glbdriftend[0];
    mRamplificationend[iDet] = glbamplificationend[0];
    //LOGP(info, "iDet:{} R:{} driftstart:{} driftend:{} ampend:{}", iDet, glb[0], glbdriftstart[0], glbdriftend[0], glbamplificationend[0]);
    glb.fill(0.f);
    glbdriftstart.fill(0.f);
    glbdriftend.fill(0.f);
    glbamplificationend.fill(0.f);
    // static constexpr float ANODEPOS = CRAH + CDRH + CAMH / 2.0 - CHSV / 2.0;
    // mRadiusOfDrift[iDet] = glb[0];
  }

  TFile* vdriftexbf = TFile::Open("o2-trd-CalVdriftExB.root");
  if (vdriftexbf == nullptr) {
  }
  o2::trd::CalVdriftExB* calvdriftexb{0};
  vdriftexbf->GetObject("ccdb_object", calvdriftexb);
  if (calvdriftexb == nullptr) {
  }
  mTransformer.init();
  mTransformer.setCalVdriftExB(calvdriftexb);
LOGP(info,"Setup root tree for output");
   mfile = new TFile("track2tracklet.root","RECREATE");
   moutputtree = new TTree("t","tracklets and tracksegments");

}
/***********************************************************************************************/

float RawDataManager::getTriggerTime(const o2::trd::TriggerRecord& trig, const std::vector<o2::dataformats::TFIDInfo>& tfids, const int timeframeno)
{
  o2::dataformats::TFIDInfo tfid;
  if (tfids.size() != 0 && timeframeno != 0) {
    tfid = tfids.at(timeframeno);
  } else {
    tfid = o2::dataformats::TFIDInfo();
  }

  if (tfid.isDummy()) {
    // LOGP(info,"Triggertime : dummy {}",trig.getBCData().bc2ns() * 1e-3);
    return trig.getBCData().bc2ns() * 1e-3;
  } else {
    o2::InteractionRecord intrec = {0, tfid.firstTForbit};
    //    std::cout << "returning time diff of "  << trig.getBCData().differenceInBCMUS(intrec) << " trig time : " << trig.getBCData().toLong() << "  from trig : " << tfid.tfCounter << "  and firstorbit : " <<  tfid.firstTForbit << "\n";
    // LOGP(info,"Triggertime : {}", trig.getBCData().differenceInBCMUS(intrec));
    return trig.getBCData().differenceInBCMUS(intrec);
  }
}

bool RawDataManager::trackMatchesCollision(const float& triggertime, const float& tracktime)
{

  float timestamperror = 1.0;
  float timestampedeviations;
  float timewindow = 2.0;
  float max = tracktime + timewindow; // timestampeerror*timestampdeviation+2.5; // in us
  float min = tracktime - timewindow; // timestampeerror*timestampdeviation-2.5; // in us
  if (triggertime > min && triggertime < max) {
    //   LOGP(info,"{} < {} < {} ----- ",min,triggertime,max);
    return true;
  }
  //  LOGP(info,"{} < {} < {} !!!! ",min,triggertime,max);
  return false;
}

void RawDataManager::prepareTracking() // std::vector<o2::trd::TriggerRecord>& trdTriggerRecords, std::vector<o2::trd::Tracklet64>& trdTracklets, std::array<int32_t,540*mMaxTriggers>& mTrackletIndexArray)//,  trdTrigRecMask)
{
  //--------------------------------------------------------------------
  // Prepare tracklet index array and if requested calculate space points
  //--------------------------------------------------------------------
  for (uint32_t iColl = 0; iColl < mTrgRecords->size(); /*trdTriggerRecords.size();*/ ++iColl) {
    int32_t nTrklts = 0;
    int32_t idxOffset = 0;
    idxOffset = (*mTrgRecords)[iColl].getFirstTracklet();
    // auto idxOffset1 = trdTriggerRecords[iColl].getFirstTracklet();
    // nTrklts = (iColl < trdTriggerRecords.size() - 1) ? trdTriggerRecords[iColl + 1].getFirstTracklet() - trdTriggerRecords[iColl].getFirstTracklet() : trdTracklets.size() - trdTriggerRecords[iColl].getFirstTracklet();
    nTrklts = (iColl < mTrgRecords->size() - 1) ? (*mTrgRecords)[iColl + 1].getFirstTracklet() - (*mTrgRecords)[iColl].getFirstTracklet() : mTrgRecords->size() - (*mTrgRecords)[iColl].getFirstTracklet();
    const o2::trd::Tracklet64* tracklets = &((*mTracklets)[idxOffset]);
    int32_t* trkltIndexArray = &mTrackletIndexArray[iColl * (540 + 1) + 1];
    trkltIndexArray[-1] = 0;
    int32_t currDet = 0;
    int32_t nextDet = 0;
    int32_t trkltCounter = 0;
    for (int32_t iTrklt = 0; iTrklt < nTrklts; ++iTrklt) {
      if (tracklets[iTrklt].getDetector() > currDet) {
        nextDet = tracklets[iTrklt].getDetector();
        for (int32_t iDet = currDet; iDet < nextDet; ++iDet) {
          trkltIndexArray[iDet] = trkltCounter;
        }
        currDet = nextDet;
      }
      ++trkltCounter;
    }
    for (int32_t iDet = currDet; iDet <= o2::trd::constants::NCHAMBER; ++iDet) {
      trkltIndexArray[iDet] = trkltCounter;
    }
    /*    if (!CalculateSpacePoints(iColl)) {
        }*/
  }
  //  mNEvents++;
}

bool RawDataManager::propagateToLayerX(o2::dataformats::TrackTPCITS& track, float xToGo, float e, float maxStep)
{
  if(debugprint) LOGP(info,"{} at line {} track.X:{} xToGo:{}",__func__,__LINE__,track.getX(),xToGo);
#define UseMaterialCorr 1
  //----------------------------------------------------------------
  //
  // Propagates the track to the plane X=xk (cm)
  // taking into account all the three components of the magnetic field
  // and correcting for the crossed material.
  //
  // maxStep  - maximal step for propagation
  // tofInfo  - optional container for track length and PID-dependent TOF integration
  //
  // matCorr  - material correction type, it is up to the user to make sure the pointer is attached (if LUT is requested)
  //----------------------------------------------------------------
  // LOGP(info,"propagate track to xToGo:{} track at : {}",xToGo,track.getX());
  auto dx = xToGo - track.getX();
  //  LOGP(info,"a dx : {}",dx);
  int dir = dx > 0.f ? 1 : -1;

  std::array<float, 3> b{};
  while (std::abs(dx) > 0.00001) {
    auto step = std::min(std::abs(dx), maxStep);
    if (dir < 0) {
      step = -step;
    }
    auto x = track.getX() + step;
    auto xyz0 = track.getXYZGlo();
    std::array<float, 3> b{};
    mProp->getFieldXYZ(xyz0, &b[0]);

    // commented out to ignore material budget for now.
    /*    auto correct = [&track, &xyz0, tofInfo, matCorr, signCorr, this]() {
          bool res = true;
          if (matCorr != MatCorrType::USEMatCorrNONE) {
            auto xyz1 = track.getXYZGlo();
            auto mb = this->getMatBudget(matCorr, xyz0, xyz1);
            if (!track.correctForMaterial(mb.meanX2X0, mb.getXRho(signCorr))) {
              res = false;
            }
            if (tofInfo) {
              tofInfo->addStep(mb.length, track.getQ2P2()); // fill L,ToF info using already calculated step length
              tofInfo->addX2X0(mb.meanX2X0);
              tofInfo->addXRho(mb.getXRho(signCorr));
            }
          } else if (tofInfo) { // if tofInfo filling was requested w/o material correction, we need to calculate the step lenght
            auto xyz1 = track.getXYZGlo();
            math_utils::Vector3D<value_type> stepV(xyz1.X() - xyz0.X(), xyz1.Y() - xyz0.Y(), xyz1.Z() - xyz0.Z());
            tofInfo->addStep(stepV.R(), track.getQ2P2());
          }
          return res;
        };*/

    if (!track.propagateTo(x, b)) {
     // if (debugprint)
      //  LOGP(info, " returning false trying to propagate track to x:{} ", x);
      return false;
    }
     //   LOGP(info, " returning true trying to propagate track to x:{} ", x);

    /*    if (maxSnp > 0 && math_utils::detail::abs<value_type>(track.getSnp()) >= maxSnp) {
          correct();
          return false;
        }*/
    /*    if (!correct()) {
          return false;
        }*/
    dx = xToGo - track.getX();
//       LOGP(info,"b dx : {} xToGo:{} trackX:{}",dx,xToGo, track.getX());
  }
  track.setX(xToGo);
  // LOGP(info,"XXXXXXXXXXXXXXXXXXX x moved by {} to ",xToGo,track.getX());
 //   LOGP(info,"{} at line {}",__func__,__LINE__);
  return true;
}

int32_t RawDataManager::getSector(float alpha)
{
  //--------------------------------------------------------------------
  // TRD sector number for reference system alpha
  //--------------------------------------------------------------------
  if (alpha < 0) {
    alpha += 2.f * std::numbers::pi;
  } else if (alpha >= 2.f * std::numbers::pi) {
    alpha -= 2.f * std::numbers::pi;
  }
  return (int32_t)(alpha * (float)18 / (2.f * std::numbers::pi));
}

float RawDataManager::getAlphaOfSector(const int32_t sec)
{
  //--------------------------------------------------------------------
  // rotation angle for TRD sector sec
  //--------------------------------------------------------------------
  float alpha = 2.0f * std::numbers::pi / (float)6 * ((float)sec + 0.5f);
  if (alpha > std::numbers::pi) {
    alpha -= 2 * std::numbers::pi;
  }
  // LOGP(info,"alpha for sector {} is {} ",sec,alpha);
  return alpha;
}

int32_t RawDataManager::getDetectorNumber(const float zPos, const float alpha, const int32_t layer)
{
  //--------------------------------------------------------------------
  // if track position is within chamber return the chamber number
  // otherwise return -1
  //--------------------------------------------------------------------
  int32_t stack = mGeo->getStack(zPos, layer);
  if (stack < 0) {
    // LOGP(info, " get detector number :: stack : {}, zpos: {} layer: {}",stack,zPos, layer);
    return -1;
  }
  int32_t sector = getSector(alpha);
  // LOGP(info, " sector : {}",sector);
  // LOGP(info, " det number : : {}",sector,mGeo->getDetector(layer, stack, sector));

  return mGeo->getDetector(layer, stack, sector);
}

bool RawDataManager::isGeoFindable(o2::dataformats::TrackTPCITS& track, const int32_t layer, const float alpha, const float zShiftTrk)
{
  //--------------------------------------------------------------------
  // returns true if track position inside active area of the TRD
  // and not too close to the boundaries
  //--------------------------------------------------------------------

  float zTrk = track.getZ() + zShiftTrk;

  int32_t det = getDetectorNumber(zTrk, alpha, layer);
  // LOGP(info, " Det : {} zShiftTrk {} getZ() {}, zTrk: {}",det,zShiftTrk,track.getZ(),zTrk);
  //  reject tracks between stacks
  if (det < 0) {
    return false;
  }

  // reject tracks in PHOS hole and for non existent chamber 17_4_4
  if (!mGeo->chamberInGeometry(det)) {
    // LOGP(info, " Chamber not in gemoetry for Det : {}",det);
    return false;
  }

  const o2::trd::PadPlane* pp = mGeo->getPadPlane(det);
  float yMax = pp->getColEnd();
  float zMax = pp->getRow0();
  float zMin = pp->getRowEnd();

  float epsY = 5.f;
  float epsZ = 5.f;

  // reject tracks closer than epsY cm to pad plane boundary
  if (yMax - std::abs(track.getY()) < epsY) {
    // LOGP(info, " Track too close to plane boundary in Y : {} - {} < {}",yMax,std::abs(track.getY()),epsY);
    return false;
  }
  // reject tracks closer than epsZ cm to stack boundary
  if (!((zTrk > zMin + epsZ) && (zTrk < zMax - epsZ))) {
    // LOGP(info, " Track too close to plane boundary in Z : {} > {} +{} && {} < {}= {}",zTrk,zMin,epsZ,zTrk,zMax,epsZ);
    return false;
  }
  // LOGP(info, "Is findable");

  return true;
}

void RawDataManager::findChambersInRoad(o2::dataformats::TrackTPCITS& track, const float roadY, const float roadZ, const int32_t iLayer, std::array<int, 18>& det, const float zMax, const float alpha, const float zShiftTrk)
{
  //--------------------------------------------------------------------
  // determine initial chamber where the track ends up
  // add more chambers of the same sector or (and) neighbouring
  // stack if track is close the edge(s) of the chamber
  //--------------------------------------------------------------------

  const float yMax = std::abs(mGeo->getCol0(iLayer));
  float zTrk = track.getZ() + zShiftTrk;

  int32_t currStack = mGeo->getStack(zTrk, iLayer);
  int32_t currSec = getSector(alpha);
  int32_t currDet;

  int32_t nDets = 0;

  if (currStack > -1) {
    // chamber unambiguous
    currDet = mGeo->getDetector(iLayer, currStack, currSec);
    // LOGP(info,"Adding for det {} with currstack {} for det {}",nDets,currStack, mGeo->getDetector(iLayer, currStack, currSec));
    det[nDets++] = currDet;
    const o2::trd::PadPlane* pp = mGeo->getPadPlane(iLayer, currStack);
    int32_t lastPadRow = mGeo->getRowMax(iLayer, currStack, 0);
    float zCenter = pp->getRowPos(lastPadRow / 2);
    if ((zTrk + roadZ) > pp->getRow0() || (zTrk - roadZ) < pp->getRowEnd()) {
      int32_t addStack = zTrk > zCenter ? currStack - 1 : currStack + 1;
      if (addStack < 5 && addStack > -1) {
        //   LOGP(info,"Adding for det {} with stack {} for det {}",nDets,addStack, mGeo->getDetector(iLayer, addStack, currSec));
        det[nDets++] = mGeo->getDetector(iLayer, addStack, currSec);
      }
    }
  } else {
    if (std::abs(zTrk) > zMax) {
      // shift track in z so it is in the TRD acceptance
      if (zTrk > 0) {
        currDet = mGeo->getDetector(iLayer, 0, currSec);
      } else {
        currDet = mGeo->getDetector(iLayer, o2::trd::constants::NSTACK - 1, currSec);
      }
      det[nDets++] = currDet;
      // LOGP(info,"Adding for det {} with currDet {} for currStack {}",nDets,currDet, mGeo->getStack(currDet));
      currStack = mGeo->getStack(currDet);
    } else {
      // track in between two stacks, add both surrounding chambers
      // gap between two stacks is 4 cm wide
      currDet = getDetectorNumber(zTrk + 4.0f, alpha, iLayer);
      if (currDet != -1) {
        det[nDets++] = currDet;
      }
      currDet = getDetectorNumber(zTrk - 4.0f, alpha, iLayer);
      if (currDet != -1) {
        det[nDets++] = currDet;
      }
    }
  }
  // add chamber(s) from neighbouring sector in case the track is close to the boundary
  if ((std::abs(track.getY()) + roadY) > yMax) {
    const int32_t nStacksToSearch = nDets;
    int32_t newSec;
    if (track.getY() > 0) {
      newSec = (currSec + 1) % 18;
    } else {
      newSec = (currSec > 0) ? currSec - 1 : 18 - 1;
    }
    for (int32_t idx = 0; idx < nStacksToSearch; ++idx) {
      currStack = mGeo->getStack(det[idx]);
      det[nDets++] = mGeo->getDetector(iLayer, currStack, newSec);
    }
  }
  // skip PHOS hole and non-existing chamber 17_4_4
  /*for (int32_t iDet = 0; iDet < nDets; iDet++) {
    if (!mGeo->chamberInGeometry(det[iDet])) {
      det[iDet] = -1;
    }
  }
  */
  // LOGP(info," {} {} ndets:{}",__func__,__LINE__,nDets);
}

bool RawDataManager::getYZAt(float xk, float b, float& y, float& z, o2::dataformats::TrackTPCITS& track)
{
  //----------------------------------------------------------------
  // estimate Y,Z in tracking frame at given X
  //----------------------------------------------------------------
  float almost0 = 0x1.0p-126f;
  float almost1 = 1.f - 1.0e-6f;
  float dx = xk - track.getX();
  y = track.getY(); // mP[kY];
  z = track.getZ(); // mP[kZ];
  if (std::abs(dx) < almost0) {
    return true;
  }
  float crv = track.getCurvature(b);
  float x2r = crv * dx;
  float f1 = track.getSnp(), f2 = f1 + x2r;
  if ((std::abs(f1) > almost1) || (std::abs(f2) > almost1)) {
    return false;
  }
  float r1 = std::sqrt((1.f - f1) * (1.f + f1));
  if (std::abs(r1) < almost0) {
    return false;
  }
  float r2 = std::sqrt((1.f - f2) * (1.f + f2));
  if (std::abs(r2) < almost0) {
    return false;
  }
  double dy2dx = (f1 + f2) / (r1 + r2);
  y += dx * dy2dx;
  if (std::abs(x2r) < 0.05f) {
    z += dx * (r2 + f2 * dy2dx) * track.getTgl();
  } else {
    // for small dx/R the linear apporximation of the arc by the segment is OK,
    // but at large dx/R the error is very large and leads to incorrect Z propagation
    // angle traversed delta = 2*asin(dist_start_end / R / 2), hence the arc is: R*deltaPhi
    // The dist_start_end is obtained from sqrt(dx^2+dy^2) = x/(r1+r2)*sqrt(2+f1*f2+r1*r2)
    //    double chord = dx*TMath::Sqrt(1+dy2dx*dy2dx);   // distance from old position to new one
    //    double rot = 2*TMath::ASin(0.5*chord*crv); // angular difference seen from the circle center
    //    track1 += rot/crv*track3;
    //
    float rot = std::sin(r1 * f2 - r2 * f1);        // more economic version from Yura.
    if (f1 * f1 + f2 * f2 > 1.f && f1 * f2 < 0.f) { // special cases of large rotations or large abs angles
      if (f2 > 0.f) {
        rot = std::numbers::pi - rot; //
      } else {
        rot = -std::numbers::pi - rot;
      }
    }
    z += track.getTgl() / crv * rot;
  }
  return true;
}

int RawDataManager::findNearestTracklet(o2::trd::TrackSegment& tracksegment)
{
  //loop through tracklets and find closest.
  //
  //
  //
  ////// WE ARE HERE !!!!!!!
  float mindistance=2.0;
  int closesttracklet=-1;
  auto& trgrec = mTriggerRecord;
  auto triggertime = getTriggerTime(trgrec, *mTFIDs, mTimeFrameNo - 1);
  //if (!trackMatchesCollision(triggertime, tracksegment.getTrackTime())) {
  //  continue;
 // }
  for(int trklt=trgrec.getFirstTracklet(); trklt<trgrec.getFirstTracklet()+trgrec.getNumberOfTracklets();++trklt){
    auto tracklet=(*mTracklets)[trklt]; 
    if(!(tracklet.getDetector()==tracksegment.getDetector())){
      continue;
    }
    LOGP(info,"Comparing tracklet padrow {} to segment padrow {}",tracklet.getPadRow(),tracksegment.getPadRow());
   // LOGP(info,"Comparing tracklet det {} to segment det {}",tracklet.getDetector(),tracksegment.getDetector());
    if(!(tracklet.getPadRow()==tracksegment.getPadRow())){
      continue;
    }
    float distance = std::abs(tracklet.getPadCol() - tracksegment.getPadColAtTimeBin());
    if(distance<mindistance){
      mindistance=distance;
      closesttracklet=trklt;
    }
  }
  return closesttracklet;
}


int RawDataManager::propagateTrack(o2::dataformats::TrackTPCITS& track, float e, float maxStep, float triggertime, int& glbTrkltIdxOffset, int collisionId)
{
  if (debugprint)
    LOGP(info, "{} {} Propagating track with track pos {:.2f} {:.2f} {:.2f} pt:{:.4f} its:{} tpc:{}", __func__, __LINE__, track.getX(), track.getY(), track.getZ(), track.getPt(), (int)track.getRefITS(), (int)track.getRefTPC());
  // Propagates the track to the plane X=xk (cm) for each respective layer
  //  float zShiftTrk = (mTrackAttribs[iTrk].mTime - GetConstantMem()->ioPtrs.trdTriggerTimes[collisionId]) * mTPCVdrift * mTrackAttribs[iTrk].mSide;
  int layerCount = 0;
  const int32_t nMaxChambersToSearch = 18;
  int32_t trkltIdxOffset = collisionId * (o2::trd::constants::NCHAMBER + 1); // offset for accessing mTrackletIndexArray for given collision
  std::vector<o2::trd::TrackSegment> btracksegments;
  std::vector<o2::trd::Tracklet64> btracklets;
  int timeframe;
  int eventno;
  int matched;
  float trackalpha,tracksnp,pT,trackxstart;
  float tracktime,trdtriggertime;
  int itsindex,tpcindex;
  int layers;
  moutputtree->Branch("segment", &btracksegments);
  moutputtree->Branch("tracklet", &btracklets);
  moutputtree->Branch("timeframe",&timeframe);
  //moutputtree->Branch("trackxstart",&trackxstart);
  moutputtree->Branch("alpha",&trackalpha);
  moutputtree->Branch("snp",&tracksnp);
  moutputtree->Branch("pt",&pT);
  moutputtree->Branch("itsindex",&itsindex);
  moutputtree->Branch("tpcindex",&tpcindex);
  moutputtree->Branch("event",&eventno);
  moutputtree->Branch("matched",&matched);
  moutputtree->Branch("triggertime",&trdtriggertime);
  moutputtree->Branch("tracktime",&tracktime);
  moutputtree->Branch("layers",&layers);
                                                                             //
  for (int32_t iLayer = 0; iLayer < 6; ++iLayer) {
    //   nCurrHypothesis = 0;
    if (debugprint)
      LOGP(info, " $$$$ Layer : {} layercount: {}", iLayer,layerCount);
    const o2::trd::PadPlane* pad = mGeo->getPadPlane(iLayer, 0);
    float tilt = std::tan(std::numbers::pi / 180.f * pad->getTiltingAngle());
    const float zMaxTRD = pad->getRow0();

    std::array<int32_t, 18> det{-1}; // TRD chambers to be searched for tracklets
    det.fill(-1);
    // propagate track to average radius of TRD layer iLayer (sector 0, stack 2 is chosen as a reference)
    if (!propagateToLayerX(track, mRdriftstart[2 * 6 + iLayer] - 1.0f, e, maxStep)) { // unwind by 1cm so inside the radiator.

      // LOGP(info,"Track propagation failed for in layer {} (pt={}, x={}, mR[layer]={})",iLayer, track.getPt(), track.getX(), mR[2 * 6 + iLayer]);
      continue;
    }
    // LOGP(info," {} {} track.x={}",__func__,__LINE__,track.getX());
    /*
    // rotate track in new sector in case of sector crossing
    if (!AdjustSector(prop, trkWork)) {
      if (ENABLE_INFO) {
        GPUInfo("Adjusting sector failed for track %i candidate %i in layer %i", iTrk, iCandidate, iLayer);
      }
      continue;
    }
*/
    // check if track is findable
    float mTPCVdrift = 2.58f;
    float side = 1.0f;
    float zShiftTrk = (track.getTimeMUS().getTimeStamp() - triggertime) * mTPCVdrift * side;

    if (!isGeoFindable(track, iLayer, track.getAlpha(), zShiftTrk)) {
   //        LOGP(info,"Track not geofindable");
      continue;
    }
    layerCount++;

    // define search window
    float roadY = 7.f * std::sqrt(track.getSigmaY2() + 0.1f * 0.1f) + 4; // Param().rec.trd.extraRoadY; // add constant to the road to account for uncertainty due to radial deviations (few mm)
    // roadZ = 7.f * CAMath::Sqrt(trkWork->getSigmaZ2() + 9.f * 9.f / 12.f); // take longest pad length
    // mRoadZ=18 is the default, extraRoadZ is 10 by default
    float roadZ = 18 + 10; // Param().rec.trd.extraRoadZ; // simply twice the longest pad length -> efficiency 99.996%
    //
    if (std::abs(track.getZ() + zShiftTrk) - roadZ >= zMaxTRD) {
           LOGP(info,"Track out of TRD acceptance with z={} in layer {} (eta={})", track.getZ() + zShiftTrk, iLayer, track.getEta());
      continue;
    }

    findChambersInRoad(track, roadY, roadZ, iLayer, det, zMaxTRD, track.getAlpha(), zShiftTrk);
    int detcounter = 0;
    int totaltrkltcounter = 0;

    std::string s = std::accumulate(
      det.begin(), det.end(), std::string{},
      [](std::string acc, int x) { return acc + std::format("[{}] ", x); }
      //[](std::string acc, int x) { return acc + (x>-1)?std::format("[{}] ", x): std::format("[{}] ", 0); }
    );
    std::string aa;
    for (auto& x : det) {
      aa += std::format("[{}]", x);
    }
    //     LOGP(info,"chambers to search : {}",s);
    //    LOGP(info,"chambers to search : {}",aa);

    for (int32_t iDet = 0; iDet < nMaxChambersToSearch; iDet++) {
      int32_t currDet = det[iDet];
      if (currDet == -1) {
        break; // continue;
      }
      pad = mGeo->getPadPlane(currDet);
      int32_t currSec = mGeo->getSector(currDet);

      if (currSec != getSector(track.getAlpha())) {
        if (!track.rotate(getAlphaOfSector(currSec))) {
          if (debugprint)
            LOGP(warn, "Track could not be rotated in tracklet coordinate system currSec:{} alpha:{}", currSec, getAlphaOfSector(currSec));
          break;
        }
      }
      if (currSec != getSector(track.getAlpha())) {
        if (debugprint)
          LOGP(info, "Track is in sector {} and we are in sector {}", getSector(track.getAlpha()), currSec);
        continue;
      }
      // propagate track to start radius of chamber drift start
      const PadPlane* pp = mGeo->getPadPlane(currDet);
      if (propagateToLayerX(track, mRdriftstart[currDet], 0.8f, 0.2f)) { // prop.propagateToX(mR[currDet], .8f, .2f)) {
        // we are at the start of a layer now propagate to the outter radius and build a tracksegment for the voxel of an mcm. TODO voxel of a padrow
        if (debugprint)
          LOGP(info, "{} {} Propagating to start track in det:{} with track pos {:.2f} {:.2f} {:.2f} pt:{:.4f} its:{} tpc:{}", __func__, __LINE__, currDet, track.getX(), track.getY(), track.getZ(), track.getPt(), (int)track.getRefITS(), (int)track.getRefTPC());
        float projY, projZ;
        auto ctrans = o2::trd::CoordinateTransformer::instance();
        TrackSegment tracksegment;
        std::array<float, 3> b{};
        auto xyz0 = track.getXYZGlo();
        math_utils::Point3D<float> trackxyz = {track.getX(), track.getY(), track.getZ()};
        // mProp->getFieldXYZ(xyz0, &b[0]);
        std::array<float, 3> rct{};
        std::array<float, 3> rcts{};
        //-------------------------------------/
        auto l2gmatrix = mGeo->getMatrixL2G(currDet);
        auto t2gmatrix = mGeo->getMatrixT2G(currDet);
        auto localpointG = t2gmatrix * trackxyz;
        auto localpoint = l2gmatrix * localpointG;
        // LOGP(info,"{} {} start of track local pos: {:.2f} {:.2f} {:.2f} pt:{:.4f} its:{} tpc:{}",__func__,__LINE__, localpoint.X(),localpoint.Y(),localpoint.Z(),track.getPt(),(int)track.getRefITS(),(int)track.getRefTPC());
        //-------------------------------------/
        // convert global to local ROC x,y,z:
        rct = ctrans->RecalculateRCT(currDet, localpoint.X(), localpoint.Y(), localpoint.Z(), ctrans->GetT0(), ctrans->GetVdrift(), ctrans->GetExB());
        //ChamberSpacePoint a(track.getRefTPC(),currDet,  localpoint.X(), localpoint.Y(), localpoint.Z(), rct, false);
        ChamberSpacePoint a(track.getRefTPC(),currDet,  localpoint.X(), localpoint.Y(), localpoint.Z(), rct, false);

        tracksegment.setStartPoint(a);
        // correction for tilted pads (only applied if deltaZ < lPad && track z err << lPad)
        float tiltCorr = tilt * (track.getZ() - projZ);
        float lPad = pad->getRowSize(a.getPadRow());
        if (!((std::abs(track.getZ() - projZ) < lPad) && (track.getSigmaZ2() < (lPad * lPad / 12.f)))) {
          tiltCorr = 0.f; // will be zero also for TPC tracks which are shifted in z
        }
        // correction for mean z position of tracklet (is not the center of the pad if track eta != 0)
        float mZCorrCoefNRC = 1.4f;
        // float zPosCorr -= zShiftTrk; // shift tracklet instead of track in order to avoid having to do a re-fit for each collision
        propagateToLayerX(track, mRamplificationend[currDet], 0.8f, 0.2f);
        if (debugprint)
          LOGP(info, "{} {} Propagating to end track in det:{} with track pos {:.2f} {:.2f} {:.2f} pt:{:.4f} its:{} tpc:{}", __func__, __LINE__, currDet, track.getX(), track.getY(), track.getZ(), track.getPt(), (int)track.getRefITS(), (int)track.getRefTPC());
        math_utils::Point3D<float> trackxyzend = {track.getX(), track.getY(), track.getZ()};
        xyz0 = track.getXYZGlo();
        // mProp->getFieldXYZ(xyz0, &b[0]);
        auto localpointendG = t2gmatrix * trackxyzend;
        auto localpointend = l2gmatrix * localpointendG;
        if (debugprint)
          LOGP(info, "{} {} end of track local pos: {:.2f} {:.2f} {:.2f} pt:{:.4f} its:{} tpc:{} rct:{}:{}:{}", __func__, __LINE__, localpointend.X(), localpointend.Y(), localpointend.Z(), track.getPt(), (int)track.getRefITS(), (int)track.getRefTPC(),rct[0],rct[1],rct[2]);
        rcts = ctrans->RecalculateRCT(currDet, localpointend.X(), localpointend.Y(), localpointend.Z(), ctrans->GetT0(), ctrans->GetVdrift(), ctrans->GetExB());
        ChamberSpacePoint ae(track.getRefTPC(),currDet, localpointend.X(), localpointend.Y(), localpointend.Z(), rcts, false);
        tracksegment.setEndPoint(ae);
        tracksegment.setCollisionId(collisionId);
        tracksegment.setTriggerTime(triggertime);
        tracksegment.setTrackTime(track.getTimeMUS().getTimeStamp());
        tracksegment.setRefTPCId(track.getRefTPC());
        tracksegment.setRefITSId(track.getRefITS());
        tracksegment.setPhi(track.getPhi());
        tracksegment.setSnp(track.getSnp());
        tracksegment.setPt(track.getPt());
        tracksegment.setSnp(track.getSnp());
        tracksegment.setAlpha(track.getAlpha());
        //TODO what to do if the tracksegment spans a padrow or mcm ?
        //postprocess the tracksegment and split it up?

       // if (debugprint)
          LOGP(info, "TrackSegment padrow:padcol:timebin {:.2f}:{:.2f}:{:.2f} --> {:.2f}:{:.2f}:{:.2f} for det:{} and padrow:{}", tracksegment.getStartPoint().getPadRowF(), tracksegment.getStartPoint().getPadCol(), tracksegment.getStartPoint().getTimeBin(), tracksegment.getEndPoint().getPadRowF(), tracksegment.getEndPoint().getPadCol(), tracksegment.getEndPoint().getTimeBin(),currDet,tracksegment.getPadRow());
        if (debugprint)
          LOGP(info, "TrackSegment  x:y:z {:.2f}:{:.2f}:{:.2f} --> {:.2f}:{:.2f}:{:.2f}", localpoint.X(), localpoint.Y(), localpoint.Z(), localpointend.X(), localpointend.Y(), localpointend.Z());
        if (debugprint)
          LOGP(info, "TrackSegmentG x:y:z {:.2f}:{:.2f}:{:.2f} --> {:.2f}:{:.2f}:{:.2f}", localpointG.X(), localpointG.Y(), localpointG.Z(), localpointendG.X(), localpointendG.Y(), localpointendG.Z());
        mITSTPCTracks_segments.push_back(tracksegment);
        btracksegments.push_back(tracksegment);
        //find nearest tracklet to this track segment.
        int trackletindex=findNearestTracklet(tracksegment);
        btracklets.push_back((*mTracklets)[trackletindex]);
      } else {
        LOGP(info, "Track could not be propagated to radius of chamber {} which is layer: {}", currDet, iLayer);
      }
    } // chamber loop

    // add no update to hypothesis list
    //} // end candidate loop
  } // end of layer loop
  // if(layerCount>0) LOGP(info,"Layer count : {}",layerCount);
  /****************************************************************************************/
  //now write the tree event ..
      timeframe=getTimeFrameNumber();
      eventno=getEventNumber();
      trackalpha=track.getAlpha();
      tracksnp=track.getSnp();
      pT=track.getPt();
      itsindex=track.getRefITS();
      tpcindex=track.getRefTPC();
      tracktime=track.getTimeMUS().getTimeStamp();
      trdtriggertime=triggertime;;   
      layers=layerCount;

 // if (layerCount > 2)
 //   return true;
 //LOGP(info,"propagateTrack returning with layercount of {}",layerCount);
  return layerCount;
}
/***********************************************************************************************/

bool RawDataManager::buildTrackSegments(bool onlydigits)
{
  // take the track, step it through the relevant detectors/mcm/padrow etc. and build the tracksegment objects.
  //  just ITSTPC for now.
  //
  //LOGP(info,"{} {}",__func__,__LINE__);
  bool timeframehasdigits = false;
  bool timeframehastracks = false;

  //clear the tracksegments
  mITSTPCTracks_segments.clear();
  if (onlydigits) {
  //LOGP(info,"{} {}",__func__,__LINE__);
    for (auto& trig : *mTrgRecords) {
      if (trig.getNumberOfDigits() > 0) {
        LOGP(info, "time frame has {} digits", trig.getNumberOfDigits());
        timeframehasdigits = true;
      }
    }
  //LOGP(info,"{} {}",__func__,__LINE__);
  }
  else {
  //LOGP(info,"{} {}",__func__,__LINE__);
    timeframehasdigits=true;
  }

  //LOGP(info,"{} {}",__func__,__LINE__);
  if ((*mITSTPCTracks).size() > 0) {
    timeframehastracks = true;
  //LOGP(info,"{} {}",__func__,__LINE__);
  }
  if (!timeframehasdigits) {
    LOGP(info, "time frame does not have digits");
  //LOGP(info,"{} {}",__func__,__LINE__);
    return false;
  }
  if (!timeframehastracks) {
    LOGP(info, "time frame does not have ITSTPC tracks");
  //LOGP(info,"{} {}",__func__,__LINE__);
    return false;
  }

  LOGP(info, "itstpc track count : {} ", (*mITSTPCTracks).size());

  // --------------------------------------------------------------------
  // loop over timeframes
  int trackcounter = 0;
  int tracktimewindowcounter = 0;
  int tracklowptcounter = 0;
  // ntimeframe=2;
  // buildsegments is called per timeframe, so dont loop over timeframe
  //  for(int timeframe=1;timeframe<ntimeframe;++timeframe) {
  mTrackletIndexArray.fill(-1);
  prepareTracking(); //,  trdTrigRecMask)

  int goodtrackcounter = 0;
  int badtrackcounter = 0;
  // std::vector<o2::trd::CalibratedTracklet> trackletscal;
  int numberoftracklets = mTracklets->size();
  // LOGP(info,"Timeframe {}  number of tracklets {}",timeframe, numberoftracklets);
  if (mTrgRecords->size() > mMaxTriggers) {
    LOGP(error, "We have a problem with too many trigger records : {} > {}", mTrgRecords->size(), mMaxTriggers);
    exit(1);
  }
  int counttracksegments=0;
  for (auto& track : *mITSTPCTracks) {
    // auto propagator = o2::base::Propagator::Instance();
    auto tpcid = track.getRefTPC();
    auto itsid = track.getRefITS();
    int trigcount = 0;
    int collisionId = 0;
    int totalNumberOfLayers=0;
    int numberOfTimeTrackMatched=0;
    for (auto& trdtrig : *mTrgRecords) {

      auto triggertime = getTriggerTime(trdtrig, *mTFIDs, mTimeFrameNo - 1);
      auto ntracklets = trdtrig.getNumberOfTracklets();
      auto trackletstart = trdtrig.getFirstTracklet();
      // LOGP(info,"Triggertime : {} collisionId {}",triggertime,collisionId);

      //      std::vector<o2::dataformats::TrackTPCITS> TracksForThisEvent;
      int trackcounter = 0;
      goodtrackcounter = 0;
      badtrackcounter = 0;
      auto ttrack = track;
      trackcounter++;
      if (trackMatchesCollision(triggertime, ttrack.getTimeMUS().getTimeStamp())) {
        continue;
      }
      trackcounter++;
      int numberOfLayers=0;
      auto difftime = ttrack.getTimeMUS().getTimeStamp() - triggertime;
      auto pt = track.getPt();
      float timeWindow = 4.0; // time is within 20us
      if (difftime < timeWindow && pt>2.0) {
        // TracksForThisEvent.push_back(track);
        //  Track_pad_row_timebin.push_back(GeneratePadRowTimeBin(track));
        numberOfTimeTrackMatched++;
        tracktimewindowcounter++;
        // if(pt>1.0){
        // LOGP(info,"$$$$ Propagating track ..... {} {} track.x={} trigtime:{}   tracktime:{}  its:{}  tpc:{} trackcounter:{} collionsId:{} timeframe:{} eventno:{}",__func__,__LINE__,ttrack.getX(),triggertime,ttrack.getTimeMUS().getTimeStamp(),
        //    (int)ttrack.getRefITS(),(int)ttrack.getRefTPC(),trackcounter,collisionId,mTimeFrameNo,mEventNo);
        numberOfLayers=propagateTrack(ttrack, .8f, 2.f, triggertime, trackletstart, collisionId);
        totalNumberOfLayers+=numberOfLayers;
        if (numberOfLayers > -1) {
          // LOGP(info,"  track propagated..... collid:{} timeframe:{} eventno:{}",collisionId,mTimeFrameNo,mEventNo);
          goodtrackcounter++; // we have at least a singular layer 0,1,2,3,4,5;
        } else {
          badtrackcounter++;
          // LOGP(info," track failed to propagate ..... collid:{} timeframe:{} eventno:{}",collisionId,mTimeFrameNo,mEventNo);
        }
        //LOGP(info,"$$$$  Finished Propagating track ..... with track.x={} collid:{} timeframe:{} eventno:{} layerscount:{} layersfromprop:{}, trackpt:{}",ttrack.getX(),collisionId,mTimeFrameNo,mEventNo,totalNumberOfLayers, numberOfLayers, ttrack.getPt());
        //}
      }

      collisionId++;
    }
    //LOGP(info,"Tracks in collisionid {} this event {} out of a total of {} good tracks from {} tracks, with {}% good events, track matched {} times numberoflayers:{}", collisionId, goodtrackcounter, trackcounter, (*mITSTPCTracks).size(), (float)goodtrackcounter/(float)trackcounter*100,numberOfTimeTrackMatched,totalNumberOfLayers);
    /********************************************************************************************/

    /********************************************************************************************/

    // now we have tracks and tracklets for a given "event" based on the time window, pair off and calculate distance.
  }
  LOGP(info, "FOR TIMEFRAME {} totaltracks {} goodtracks {} nonlowptctracks {} timewindowtracks {}::{:.2}% lowpttracks {}::{:.2}%", mTimeFrameNo, trackcounter, goodtrackcounter, trackcounter - tracklowptcounter, tracktimewindowcounter, ((float)tracktimewindowcounter / (float)(trackcounter - tracklowptcounter)) * 100.0, tracklowptcounter, ((float)tracklowptcounter / (float)trackcounter) * 100.0);
  return true;
}

bool RawDataManager::nextTimeFrame(bool onlydigits)
{
  if (!mDataTree->GetEntry(mTimeFrameNo)) {
    // loading time frame will fail at end of file
    return false;
  }

  mEventNo = 0;
  mTimeFrameNo++;

  LOGP(info, "Loaded data for time frame #{} with {} TRD trigger records, {} digits and {} tracklets",
       mTimeFrameNo, mTrgRecords->size(), mDigits->size(), mTracklets->size());
  //  for(auto& track : *mITSTPCTracks){
  //    if((int)track.getRefITS()==73 && (int)track.getRefTPC()>121000) LOGP(info,">>> {} {} track.x={} tracktime:{}  its:{}  tpc:{}",__func__,__LINE__,track.getX(),track.getTimeMUS().getTimeStamp(),(int)track.getRefITS(),(int)track.getRefTPC());
  //  }
  //count tracks about a pt value:
  //auto highpttracks=std::count_if(mITSTPCTracks->begin(),mITSTPCTracks->end(),[](const Track &t){return t.getPt()>1.0;});
  //auto higherpttracks=std::count_if(mITSTPCTracks->begin(),mITSTPCTracks->end(),[](const Track &t){return t.getPt()>1.0;});
  int highpttracks=0;
  int higherpttracks=0;
  for(auto& tmptrack : *mITSTPCTracks){
    if(tmptrack.getPt()>1.0) highpttracks++;
    if(tmptrack.getPt()>2.0) higherpttracks++;
  }
  //LOGP(info, "Building track segments for time frame {} that has {} tracks, with pt>1.0 {} with pt>2.0 {}", mTimeFrameNo, mITSTPCTracks->size(),highpttracks,higherpttracks);
  auto tracksegmentstart = std::chrono::high_resolution_clock::now(); // measure total processing time
  if(mTimeFrameNo>6)buildTrackSegments(onlydigits);
  //buildTrackSegments(onlydigits);
  //LOGP(info,"sorting ITSTPC track segments with size : {} mTimeFrameNo : {}",mITSTPCTracks_segments.size(),mTimeFrameNo);
  std::stable_sort(mITSTPCTracks_segments.begin(),mITSTPCTracks_segments.end(),comp_tracksegments);
  //tracksegements are now trd trigger order.
  auto tracksegmenttime = std::chrono::high_resolution_clock::now() - tracksegmentstart;
  LOGP(info, "Built {} track segments for time frame {} in {} ms ITSTPCTracks {} ITSTPCTracks(>1GeV): {} ITSTPCTracks(>2GeV): {}", mITSTPCTracks_segments.size(), mTimeFrameNo, std::chrono::duration_cast<std::chrono::milliseconds>(tracksegmenttime).count(),mITSTPCTracks->size(),highpttracks,higherpttracks);
  std::array<int,75> collissionsegmentcount; collissionsegmentcount.fill(-1);
  for(const auto& tg : mITSTPCTracks_segments ){
    collissionsegmentcount[tg.getCollisionId()]++;
  }
  int counter=0;
  for(auto& a : collissionsegmentcount){
    if(a==-1) break;
    LOGP(info,"collision [{}] has {} tracksegments",counter++,a);

  }
  // sort by collision_id;
  return true;
}

bool RawDataManager::nextEvent()
{
  // get the next trigger record
  if (mEventNo >= mTrgRecords->size()) {
    return false;
  }
  mTriggerRecord = mTrgRecords->at(mEventNo);
  LOGP(info, "Processing event: orbit {} bc {:04d} with {} digits, {} tracklets, and ITSTPC matched track segments {}",
       mTriggerRecord.getBCData().orbit, mTriggerRecord.getBCData().bc,
       mTriggerRecord.getNumberOfDigits(), mTriggerRecord.getNumberOfTracklets(), mITSTPCTracks_segments.size());

  if (mCollisionContext) {

    // clear MC data
    mHitPoints.clear();

    for (int i = 0; i < mCollisionContext->getNCollisions(); ++i) {
      auto evrec = mCollisionContext->getEventRecords()[i];
      if (std::abs(mTriggerRecord.getBCData().differenceInBCNS(evrec)) <= 3000) {
        // if (mMCReader) {
        mMCTree->GetEntry(i);
        // }

        LOGP(info, "Loaded matching MC event #{} with time offset {:.2f} ns and {} hits",
               i, mTriggerRecord.getBCData().differenceInBCNS(evrec), mHits->size());

        // convert hits to spacepoints
        auto ctrans = o2::trd::CoordinateTransformer::instance();
        for (auto& hit : *mHits) {
          mHitPoints.emplace_back(ctrans->MakeSpacePoint(hit), hit.GetCharge());
        }
      }
    }
  }
  if (mITSTPCTracks_segments.size() > 0) {
    // we have track segments for this timeframe.
    //  they are time sorted so find the first one with in the
    // if(mITSTPCTracks_segments.)
  }

  mEventNo++;
  return true;
}

RawDataSpan RawDataManager::getEvent()
{
  RawDataSpan ev;

  ev.digits = boost::make_iterator_range_n(mDigits->begin() + mTriggerRecord.getFirstDigit(), mTriggerRecord.getNumberOfDigits());
  ev.tracklets = boost::make_iterator_range_n(mTracklets->begin() + mTriggerRecord.getFirstTracklet(), mTriggerRecord.getNumberOfTracklets());

  ev.hits = boost::make_iterator_range(mHitPoints.begin(), mHitPoints.end());

  auto evtime = getTriggerTime();
  // if (tpctracks) {
  //   for (auto &track : *mTpcTracks) {
  //     //   // auto tracktime = track.getTimeMUS().getTimeStamp();
  //     auto dtime = track.getTime0() / 5.0 - evtime;
  //     if (dtime > mMatchTimeMinTPC && dtime < mMatchTimeMaxTPC) {
  //       ev.mTpcTracks.push_back(track);
  //     }
  //   }
  // }
  // find first tracksegment for this event.
  //std::vector<TrackSegment>::iterator first;
  //std::vector<TrackSegment>::iterator last;
  int first=-1,last=9999999;
  int counter=0;
  //LOGP(info,"mITSTPCTracks_segments.size() {} first {}  last {} before finding first and last",mITSTPCTracks_segments.size(),first,last);
  for(auto& ts : mITSTPCTracks_segments ){
    //LOGP(info,"comparing tracksegments for extraction {} {} first {} last {} counter {}",ts.getTriggerTime(),evtime,first,last,counter);
    
    if(first!=-1 && std::abs(ts.getTriggerTime()-evtime)>1.5){
     // LOGP(info,"setting last tracksegments for extraction {} {} first {} last {} counter {}",ts.getTriggerTime(),evtime,first,last,counter);
      last = counter;
      break;
    }
    
    if(std::abs(ts.getTriggerTime()-evtime)<1.5 && first==-1){
    //LOGP(info,"setting first tracksegments for extraction {} {} first {} last {} counter {}",ts.getTriggerTime(),evtime,first,last,counter);
      first = counter;
    }
    ++counter;
  }
  // find last tracksegment for this event.
  if(first !=-1 && last != 9999999) ev.tracks_itstpc_seg= 
    boost::make_iterator_range(mITSTPCTracks_segments.begin()+first, mITSTPCTracks_segments.begin()+last);

  LOGP(info,"mITSTPCTracks_segments {} ",last-first);
  // ev.trackpoints.begin() = ev.evtrackpoints.begin();
  // ev.trackpoints.end() = ev.evtrackpoints.end();

  return ev;
}

o2::dataformats::TFIDInfo RawDataManager::getTimeFrameInfo()
{
  if (mTFIDs) {
    // LOGP(info,"mTFIDs is valid !");
    return mTFIDs->at(mTimeFrameNo - 1);
  } else {
    LOGP(info, "mTFIDs is invalid !");
    return o2::dataformats::TFIDInfo();
  }
}

float RawDataManager::getTriggerTime()
{
  auto tfid = getTimeFrameInfo();

  if (tfid.isDummy()) {
    return mTriggerRecord.getBCData().bc2ns() * 1e-3;
  } else {
    o2::InteractionRecord intrec = {0, tfid.firstTForbit};
    return mTriggerRecord.getBCData().differenceInBCMUS(intrec);
  }
}

std::string RawDataManager::describeFiles()
{
  std::ostringstream out;
  if (!mMainFile) {
    out << "RawDataManager is not connected to any files" << std::flush;
    return out.str();
  }
  if (!mDataTree) {
    out << "ERROR: main datatree not connected" << std::flush;
    return out.str();
  }
  out << "Main file:" << mMainFile->GetPath() << " has " << mDataTree->GetEntries() << " time frames " << std::endl;
  if (mDataTree->GetFriend("TRDDigit")) {
    out << "digits" << std::endl;
  }
  /*  if (mDataTree->GetFriend("TPCITS")) {
      out << "tpc its matches" << std::endl;
    }*/

  if (mTFIDs) {
    out << mTFIDs->size() << " TFIDs were read from o2_tfidinfo.root" << std::flush;
  }
  return out.str();
}

std::string RawDataManager::describeTimeFrame()
{
  std::ostringstream out;
  out << "## Time frame " << mTimeFrameNo << ": ";
  // out << mDatareader->GetEntries() << "";
  return out.str();
}

std::string RawDataManager::describeEvent()
{
  std::ostringstream out;
  out << "## TF:Event " << mTimeFrameNo << ":" << mEventNo << ":  "
      //  << hits->getsize() << " hits   "
      << mTriggerRecord.getNumberOfDigits() << " digits and "
      << mTriggerRecord.getNumberOfTracklets() << " tracklets";
  return out.str();
}
