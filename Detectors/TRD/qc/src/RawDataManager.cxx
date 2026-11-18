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
#include <TBranch.h>
#include <boost/range/distance.hpp>
#include <boost/range/iterator_range_core.hpp>
#include <iterator>
#include "TRDQC/CoordinateTransformer.h"
#include "Framework/Logger.h"

#include <map>
#include <set>
#include <utility>
#include <algorithm>
#include <cmath>
#include <TGeoManager.h>
#include "TRDBase/Geometry.h"
#include "MathUtils/Utils.h"

using namespace o2::trd;

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

  return false;
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

  return false;
}

void RawDataSpan::sort()
{
  std::stable_sort(std::begin(digits), std::end(digits), comp_digit);
  std::stable_sort(std::begin(tracklets), std::end(tracklets), comp_tracklet);
  std::stable_sort(std::begin(hits), std::end(hits), comp_spacepoint);
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

  for (const auto& [key, begin] : firsthit) {
    spanmap[key].hits = boost::make_iterator_range(begin, hits.end());
  }

  std::multimap<uint32_t, const RecoTrackOverlay*> overlaysByRow;
  for (const auto& track : recoTracks) {
    overlaysByRow.emplace(100 * track.detector + track.row, &track);
  }
  for (auto& [key, span] : spanmap) {
    const auto [first, last] = overlaysByRow.equal_range(keyfunc::overlayRowKey(key));
    for (auto it = first; it != last; ++it) {
      if (keyfunc::matchesOverlay(key, *it->second)) {
        span.recoTracks.push_back(*it->second);
      }
    }
  }

  // convert the map of spans into a vector, as we do not need the access by key
  // and longer, and having a vector makes the looping by the user easier.
  std::vector<RawDataSpan> spans;
  transform(spanmap.begin(), spanmap.end(), back_inserter(spans), [](auto const& pair) { return pair.second; });

  return spans;
}

/// PadRowID is a struct to calculate unique identifiers per pad row.
/// The struct can be passed as a template parameter to the RawDataSpan::IterateBy
/// method to split the data span by pad row and iterate over the pad rows.
struct PadRowID {
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

  static uint32_t overlayRowKey(uint32_t key) { return key; }
  static bool matchesOverlay(uint32_t, const RecoTrackOverlay&) { return true; }

  static bool match(const uint32_t key, const o2::trd::ChamberSpacePoint& x)
  {
    return key == 100 * x.getDetector() + x.getPadRow();
  }
};

// instantiate the template to iterate by padrow
template std::vector<RawDataSpan> RawDataSpan::iterateBy<PadRowID>();

// non-template wrapper function to keep PadRowID within the .cxx file
std::vector<RawDataSpan> RawDataSpan::iterateByPadRow() { return iterateBy<PadRowID>(); }

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

  static uint32_t overlayRowKey(uint32_t key) { return 100 * (key / 1000) + (key % 1000) / 8; }

  static bool matchesOverlay(uint32_t key, const RecoTrackOverlay& x)
  {
    if (x.trackletIndex >= 0) {
      return key == MCM_ID::key(x.tracklet);
    }
    const int firstPad = (key % 1000 % 8) * constants::NCOLMCM - 1;
    const int lastPad = firstPad + constants::NCOLMCM + 3;
    return std::max(x.pad0, x.pad1) >= firstPad && std::min(x.pad0, x.pad1) <= lastPad;
  }

  static int getDetector(uint32_t k) { return k / 1000; }
  // static int getPadRow(key) {return (key%1000) / 8;}
  static int getMcmRowCol(uint32_t k) { return k % 1000; }
};

// template instantion and non-template wrapper function
template std::vector<RawDataSpan> RawDataSpan::iterateBy<MCM_ID>();
std::vector<RawDataSpan> RawDataSpan::iterateByMCM() { return iterateBy<MCM_ID>(); }

// I started to implement a struct to iterate by detector, but did not finish this
// struct DetectorID {
//   /// The static `key` method calculates a padrow ID for digits and tracklets
//   template <typename T>
//   static uint32_t key(const T& x)
//   {
//     return x.getDetector();
//   }

//   static std::vector<uint32_t> keys(const o2::trd::ChamberSpacePoint& x)
//   {
//     uint32_t key = x.getDetector();
//     return {key};
//   }

//   static bool match(const uint32_t key, const o2::trd::ChamberSpacePoint& x)
//   {
//     return key == x.getDetector();
//   }
// };

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

  if (!mDataTree) {
    return;
  }

  // set up the branches we want to read
  mDataTree->SetBranchAddress("Tracklet", &mTracklets);
  mDataTree->SetBranchAddress("TrackTrg", &mTrgRecords);

  if (std::filesystem::exists(dir / "trddigits.root")) {
    mDataTree->AddFriend("o2sim", (dir / "trddigits.root").c_str());
    mDataTree->SetBranchAddress("TRDDigit", &mDigits);
  }

  mITSTPCTree = openTrackTree(dir / "o2match_itstpc.root", "matchTPCITS");
  mTRDTrackTree = openTrackTree(dir / "trdmatches_itstpc.root", "tracksTRD");
  mITSTPCTOFTree = openTrackTree(dir / "o2match_tof_itstpc.root", "matchTOF");
  mTRDTOFTree = openTrackTree(dir / "o2match_tof_itstpctrd.root", "matchTOF");
  if ((mITSTPCTOFTree && !mITSTPCTree) || (mTRDTOFTree && !mTRDTrackTree)) {
    LOG(warning) << "A TOF match file is missing its parent track file";
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
}

RawDataManager::~RawDataManager() = default;

TTree* RawDataManager::openTrackTree(const std::filesystem::path& path, const char* name)
{
  if (!std::filesystem::exists(path)) {
    return nullptr;
  }
  std::unique_ptr<TFile> file(TFile::Open(path.c_str()));
  TTree* tree = nullptr;
  if (file && !file->IsZombie()) {
    file->GetObject(name, tree);
  }
  if (!tree) {
    LOG(warning) << "Cannot read tree " << name << " in " << path.string();
    return nullptr;
  }
  if (tree->GetEntries() != mDataTree->GetEntries()) {
    LOG(warning) << "Skipping " << path.string() << ": timeframe entry count differs from trdtracklets.root";
    return nullptr;
  }
  mTrackFiles.push_back(std::move(file));
  return tree;
}

namespace
{
template <typename T>
void readTrackBranch(TTree* tree, const char* name, Long64_t entry, std::vector<T>& data)
{
  data.clear();
  if (!tree || entry >= tree->GetEntries()) {
    return;
  }
  auto* address = &data;
  TBranch* branch = nullptr;
  if (tree->SetBranchAddress(name, &address, &branch) < 0 || !branch) {
    LOG(warning) << "Cannot bind track branch " << name << " in " << tree->GetName();
    return;
  }
  const auto bytes = branch->GetEntry(entry);
  tree->ResetBranchAddress(branch);
  if (bytes <= 0) {
    data.clear();
  }
}
} // namespace

bool RawDataManager::nextTimeFrame()
{
  if (!mDataTree || mDataTree->GetEntry(mTimeFrameNo) <= 0 || !mTracklets || !mTrgRecords) {
    // loading time frame will fail at end of file
    return false;
  }

  readTrackBranch(mITSTPCTree, "TPCITS", mTimeFrameNo, mTracks);
  readTrackBranch(mTRDTrackTree, "tracks", mTimeFrameNo, mTRDTracks);
  readTrackBranch(mTRDTrackTree, "trgrec", mTimeFrameNo, mTRDTrackTriggers);
  readTrackBranch(mITSTPCTOFTree, "TOFMatchInfo", mTimeFrameNo, mITSTPCTOFMatches);
  readTrackBranch(mTRDTOFTree, "TOFMatchInfo", mTimeFrameNo, mTRDTOFMatches);
  mDisplayTracklets = *mTracklets;
  mRecoTracks.clear();
  mEventNo = 0;
  mTimeFrameNo++;

  O2INFO("Loaded data for time frame #%d with %d TRD trigger records, %d digits and %d tracklets",
         mTimeFrameNo, mTrgRecords->size(), mDigits ? mDigits->size() : 0, mTracklets->size());

  return true;
}

bool RawDataManager::nextEvent()
{
  // get the next trigger record
  if (!mTrgRecords || mEventNo >= mTrgRecords->size()) {
    return false;
  }
  mTriggerRecord = mTrgRecords->at(mEventNo);
  O2INFO("Processing event: orbit %d bc %04d with %d digits and %d tracklets",
         mTriggerRecord.getBCData().orbit, mTriggerRecord.getBCData().bc,
         mTriggerRecord.getNumberOfDigits(), mTriggerRecord.getNumberOfTracklets());

  mHitPoints.clear();
  if (mCollisionContext && mMCTree) {

    for (int i = 0; i < mCollisionContext->getNCollisions(); ++i) {
      auto evrec = mCollisionContext->getEventRecords()[i];
      if (abs(mTriggerRecord.getBCData().differenceInBCNS(evrec)) <= 3000) {
        // if (mMCReader) {
        if (mMCTree->GetEntry(i) <= 0 || !mHits) {
          continue;
        }
        // }

        O2INFO("Loaded matching MC event #%d with time offset %f ns and %d hits",
               i, mTriggerRecord.getBCData().differenceInBCNS(evrec), mHits->size());

        // convert hits to spacepoints
        auto ctrans = o2::trd::CoordinateTransformer::instance();
        for (auto& hit : *mHits) {
          mHitPoints.emplace_back(ctrans->MakeSpacePoint(hit), hit.GetCharge());
        }
      }
    }
  }

  buildTrackOverlays();
  mEventNo++;
  return true;
}

namespace
{
using TrackID = o2::dataformats::GlobalTrackID;
using MatchInfo = o2::dataformats::MatchInfoTOF;
using ChamberPoint = ROOT::Math::Impl::Transform3D<double>::Point;

float tofTimeMUS(const MatchInfo& match)
{
  return (match.getSignal() - match.getLTIntegralOut().getTOF(o2::track::PID::Pion)) * 1.e-6;
}

std::vector<const MatchInfo*> indexTOFMatches(const std::vector<MatchInfo>& matches,
                                             TrackID::Source parent, size_t size)
{
  std::vector<const MatchInfo*> byTrack(size, nullptr);
  for (const auto& match : matches) {
    auto ref = match.getTrackRef();
    if (ref.getSource() == parent && ref.getIndex() < size) {
      auto*& selected = byTrack[ref.getIndex()];
      if (!selected || match.getChi2() < selected->getChi2()) {
        selected = &match;
      }
    }
  }
  return byTrack;
}

void projectTrack(const o2::track::TrackParCov& state, float bz, const RecoTrackOverlay& info,
                  std::vector<RecoTrackOverlay>& output)
{
  using namespace o2::trd::constants;
  auto* geo = Geometry::instance();
  auto* coordinates = CoordinateTransformer::instance();
  constexpr int nSamples = 13; // 2.5 mm spacing across the 3 cm drift gap
  const double driftStart = -Geometry::camHght() / 2.;
  for (int layer = 0; layer < NLAYER; ++layer) {
    bool ok = false;
    auto nominal = state.getXYZGloAt(Geometry::getTime0(layer), bz, ok);
    if (!ok) {
      continue;
    }
    if (!std::isfinite(nominal.X()) || !std::isfinite(nominal.Y())) {
      continue;
    }
    const int centralSector = o2::math_utils::angle2Sector(std::atan2(nominal.Y(), nominal.X()));
    // Include neighbouring sectors at chamber edges.
    for (int offset = -1; offset <= 1; ++offset) {
      const int sector = (centralSector + offset + NSECTOR) % NSECTOR;
      const float alpha = o2::math_utils::sector2Angle(sector);
      auto track = state;
      if (!track.rotate(alpha)) {
        continue;
      }
      for (int stack = 0; stack < NSTACK; ++stack) {
        const int det = geo->getDetector(layer, stack, sector);
        if (!geo->chamberInGeometry(det)) {
          continue;
        }
        const auto& matrix = geo->getMatrixT2L(det);
        const auto* pads = geo->getPadPlane(det);
        int previousRow = -1;
        std::array<float, 3> previous{};
        for (int sample = 0; sample < nSamples; ++sample) {
          const double driftX = driftStart - Geometry::cdrHght() * sample / (nSamples - 1);
          const double chamberX = Geometry::anodePos() + driftX;
          double trackX = (matrix ^ ChamberPoint(chamberX, 0., 0.)).X();
          ChamberPoint local;
          bool converged = false;
          // Match the aligned chamber plane rather than the nominal sector X.
          for (int iteration = 0; iteration < 4; ++iteration) {
            float y, z;
            if (!track.getYZAt(trackX, bz, y, z)) {
              break;
            }
            local = matrix(ChamberPoint(trackX, y, z));
            const double dx = chamberX - local.X();
            if (std::abs(dx) < 0.002) {
              converged = true;
              break;
            }
            trackX += dx;
          }
          if (!converged || local.Z() > pads->getRow0ROC() || local.Z() <= pads->getRowEndROC()) {
            previousRow = -1;
            continue;
          }
          const int row = pads->getPadRowNumberROC(local.Z());
          if (row < 0 || pads->getPadRow(local.Z()) >= pads->getNrows()) {
            previousRow = -1;
            continue;
          }
          // Stay on the drift-region side of the transform's piecewise time map.
          const float x = std::min(float(driftX), std::nextafter(float(driftStart), -1.f));
          auto point = coordinates->Local2RCT(det, x, local.Y(), local.Z());
          if (!std::isfinite(point[1]) || point[1] < 0.f || point[1] >= pads->getNcols()) {
            previousRow = -1;
            continue;
          }
          if (row == previousRow) {
            auto piece = info;
            piece.detector = det;
            piece.row = row;
            piece.pad0 = previous[1];
            piece.time0 = previous[2];
            piece.pad1 = point[1];
            piece.time1 = point[2];
            output.push_back(piece);
          }
          previous = point;
          previousRow = row;
        }
      }
    }
  }
}
} // namespace

void RawDataManager::buildTrackOverlays()
{
  mRecoTracks.clear();
  if (mTracks.empty() && mTRDTracks.empty()) {
    return;
  }
  const auto tfid = getTimeFrameInfo();
  const bool hasTriggerTime = !tfid.isDummy();
  const float triggerTime = hasTriggerTime ? getTriggerTime() : 0.f;
  auto* geo = Geometry::instance();
  const bool canProject = std::isfinite(mTrackBz) && gGeoManager &&
                          (geo->isBuilt() || geo->createClusterMatrixArray());
  auto makeInfo = [triggerTime, hasTriggerTime](size_t index, TrackID::Source source,
                                               const std::vector<MatchInfo>& matches,
                                               const MatchInfo* match, float time) {
    RecoTrackOverlay info;
    info.track = TrackID(index, source);
    info.triggerTimeResidual = hasTriggerTime ? time - triggerTime : std::numeric_limits<float>::quiet_NaN();
    if (match) {
      info.hasTOF = true;
      info.tofMatch = TrackID(match - matches.data(),
                              source == TrackID::ITSTPC ? TrackID::ITSTPCTOF : TrackID::ITSTPCTRDTOF);
    }
    return info;
  };

  if (canProject && hasTriggerTime) {
    const auto tofByTrack = indexTOFMatches(mITSTPCTOFMatches, TrackID::ITSTPC, mTracks.size());
    for (size_t i = 0; i < mTracks.size(); ++i) {
      const auto& track = mTracks[i];
      const float time = track.getTimeMUS().getTimeStamp();
      const auto* tof = tofByTrack[i];
      if (tof && !(std::abs(tofTimeMUS(*tof) - triggerTime) <= mTrackTimeWindow)) {
        tof = nullptr;
      }
      const bool parentAssociated = std::abs(time - triggerTime) <= mTrackTimeWindow;
      if (tof || parentAssociated) {
        auto info = makeInfo(i, TrackID::ITSTPC, mITSTPCTOFMatches, tof, tof ? tofTimeMUS(*tof) : time);
        info.parentAssociated = parentAssociated;
        projectTrack(track.getParamOut(), mTrackBz, info, mRecoTracks);
      }
    }
  }

  const auto tofByTrack = indexTOFMatches(mTRDTOFMatches, TrackID::ITSTPCTRD, mTRDTracks.size());
  for (const auto& record : mTRDTrackTriggers) {
    if (record.getBCData() != mTriggerRecord.getBCData()) {
      continue;
    }
    const size_t first = record.getFirstTrack();
    const size_t count = record.getNumberOfTracks();
    if (first > mTRDTracks.size() || count > mTRDTracks.size() - first) {
      continue;
    }
    for (size_t i = first; i < first + count; ++i) {
      const auto& track = mTRDTracks[i];
      const auto* tof = tofByTrack[i];
      const auto info = makeInfo(i, TrackID::ITSTPCTRD, mTRDTOFMatches, tof,
                                 tof ? tofTimeMUS(*tof) : triggerTime);
      if (canProject) {
        projectTrack(track.getOuterParam(), mTrackBz, info, mRecoTracks);
      }
      for (int layer = 0; layer < constants::NLAYER; ++layer) {
        const int index = track.getTrackletIndex(layer);
        const int firstTracklet = mTriggerRecord.getFirstTracklet();
        if (index < 0 || index < firstTracklet || index - firstTracklet >= mTriggerRecord.getNumberOfTracklets() ||
            size_t(index) >= mTracklets->size()) {
          continue;
        }
        auto attached = info;
        attached.trackletIndex = index;
        attached.tracklet = (*mTracklets)[index];
        attached.detector = attached.tracklet.getDetector();
        attached.row = attached.tracklet.getPadRow();
        mRecoTracks.push_back(attached);
      }
    }
  }
}

RawDataSpan RawDataManager::getEvent()
{
  RawDataSpan ev;

  ev.digits = boost::make_iterator_range(mEmptyDigits.begin(), mEmptyDigits.end());
  const size_t firstDigit = mTriggerRecord.getFirstDigit();
  const size_t nDigits = mTriggerRecord.getNumberOfDigits();
  if (mDigits && firstDigit <= mDigits->size() && nDigits <= mDigits->size() - firstDigit) {
    ev.digits = boost::make_iterator_range_n(mDigits->begin() + firstDigit, nDigits);
  }
  ev.tracklets = boost::make_iterator_range(mDisplayTracklets.end(), mDisplayTracklets.end());
  const size_t firstTracklet = mTriggerRecord.getFirstTracklet();
  const size_t nTracklets = mTriggerRecord.getNumberOfTracklets();
  if (firstTracklet <= mDisplayTracklets.size() && nTracklets <= mDisplayTracklets.size() - firstTracklet) {
    ev.tracklets = boost::make_iterator_range_n(mDisplayTracklets.begin() + firstTracklet, nTracklets);
  }

  ev.hits = boost::make_iterator_range(mHitPoints.begin(), mHitPoints.end());

  ev.recoTracks = mRecoTracks;
  return ev;
}

o2::dataformats::TFIDInfo RawDataManager::getTimeFrameInfo()
{
  if (mTFIDs && mTimeFrameNo > 0 && mTimeFrameNo <= mTFIDs->size()) {
    return mTFIDs->at(mTimeFrameNo - 1);
  } else {
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
  if (mDataTree->GetBranch("TRDDigit")) {
    out << "digits" << std::endl;
  }
  for (const auto& file : mTrackFiles) {
    out << file->GetName() << '\n';
  }

  if (mTFIDs) {
    out << mTFIDs->size() << " TFIDs were read from o2_tfidinfo.root\n";
  } else if (mITSTPCTree) {
    out << "ITS-TPC timing overlays require o2_tfidinfo.root\n";
  }
  if ((mITSTPCTree || mTRDTrackTree) && (!std::isfinite(mTrackBz) || !gGeoManager)) {
    out << "Projected tracks require run geometry and nominal Bz\n";
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
