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
#include "DataFormatsTRD/HelperMethods.h"
#include "TRDQC/CoordinateTransformer.h"

#include "Framework/Logger.h"
#include <TGeoManager.h>
#include <TVirtualMC.h>

#include <RtypesCore.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <boost/range/distance.hpp>
#include <boost/range/iterator_range_core.hpp>
#include <iterator>
#include <set>
#include <utility>

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

bool comp_trackrefs(const o2::TrackReference& a, const o2::TrackReference& b)
{
  if (a.getTrackID() != b.getTrackID()) {
    return a.getTrackID() < b.getTrackID();
  }

  if (a.getUserId()>>2 != b.getUserId()>>2) {
    return a.getUserId()>>2 < b.getUserId()>>2; // this gets the detector/stack/layer
  }

  if (a.getTime() != b.getTime()) {
    return a.getTime() < b.getTime(); // time of flight, so order from interaction region
  }

  return true;
}

bool comp_trackrefshits(const MCTrackletSegmentInfo& a, const MCTrackletSegmentInfo& b)
{
  if (a.mEnter.getTime() != b.mEnter.getTime()) {
    return a.mEnter.getTime() < b.mEnter.getTime(); // time of flight, so order from interaction region
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
  /*if (a.mTrackId != b.mTrackId) {
    return a.mTrackId < b.mTrackId;
  }

  if (a.mDet>>2 != b.mDet>>2) {
    return (a.mDet>>2) <  (b.mDet>>2); // this gets the detector/stack/layer
  }

  if (a.mEnter.getTime() != b.mEnter.getTime()) {
    return a.mEnter.getTime() < b.mEnter.getTime(); // time of flight, so order from interaction region
  }
*/
  return true;
}

void RawDataSpan::sort()
{
  std::stable_sort(std::begin(digits), std::end(digits), comp_digit);
  std::stable_sort(std::begin(tracklets), std::end(tracklets), comp_tracklet);
  std::stable_sort(std::begin(hits), std::end(hits), comp_spacepoint);
  std::stable_sort(std::begin(trackrefsegments), std::end(trackrefsegments), comp_trackrefshits);
}

template <typename keyfunc>
std::vector<RawDataSpan> RawDataSpan::iterateBy()
{
  // an map for keeping track which ranges correspond to which key
  std::map<uint32_t, RawDataSpan> spanmap;
  std::vector<uint32_t> foundkeys;

  //TODO come back and try assign the trackrefhits to particular keys.
  //for now we just add all of them to each key :-(

  // sort digits and tracklets
  sort();

  // add all the digits to a map
  for (auto cur = digits.begin(); cur != digits.end(); /* noop */) {
    // calculate the key of the current (first unprocessed) digit
    auto key = keyfunc::key(*cur,true);
    foundkeys.push_back(key);
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
    foundkeys.push_back(key);
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
    auto keys = keyfunc::keys(*cur,true);
    // if we are not yet aware of this key, register the current hit as the first hit
    for (auto key : keys) {
      firsthit.insert({key, cur});
    }
    // remove the keys from the firsthit map that are no longer found in the hits
    for (auto it = firsthit.cbegin(); it != firsthit.cend(); /* no increment */) {
      if (keys.find(it->first) == keys.end()) {
        spanmap[it->first].hits = boost::make_iterator_range(it->second, cur);
        it = firsthit.erase(it);
      } else {
        ++it;
      }
    }
  }
  // trackreferences underpinning the tracksegments are in time order, match the times with the 
  // given a key, find the first and last tracklet, take the first and last time.
  // advance to the first time in the trackrefsegments
  auto uniquekeys=std::unique(foundkeys.begin(),foundkeys.end());
  foundkeys.erase(uniquekeys,foundkeys.end());
  for( auto &key : foundkeys) {
    //spanmap[key].trackrefsegments= boost::make_iterator_range(trackrefsegments.begin(),trackrefsegments.end());
  //  std::cout << "key : " << key << "\n";
  }
  int count=0;
  for (auto cur = trackrefsegments.begin(); cur != trackrefsegments.end(); /* noop */) {
    //LOGP(info,"********************** tracksegment start:  ***********************************");
    //std::cout << "index: " << count++ << "\n";


    auto key = keyfunc::key(*cur,true);
    auto it= std::find(foundkeys.begin(),foundkeys.end(),key);
    if(it == foundkeys.end()){
      LOGP(info, "tracksegment key is not present so error key:{}", key);
      
    }
    else{
      LOGP(info, "tracksegment key is present so error key:{}", key);
    }
    
    //LOGP(info,"!!!!!!!!!!!!!!!!!!!!!! tracksegment find_if:  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    //std::cout << "!!!!!!!!!!!!!!!!!!!!!! cur " << *cur << "\n";
    auto nxt = std::find_if(cur, trackrefsegments.end(), [key](auto x) { return keyfunc::key(x) != key; });
    //std::cout << "!!!!!!!!!!!!!!!!!!!!!! nxt " << *nxt << "\n";
    //LOGP(info,"!!!!!!!!!!!!!!!!!!!!!! tracksegment find_if:  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    spanmap[key].trackrefsegments = boost::make_iterator_range(cur, nxt);
    cur = nxt;
    //LOGP(info,"********************** tracksegment end:  ***********************************");
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
  static uint32_t key(const T& x, bool print=false)
  {
    return 100 * x.getDetector() + x.getPadRow();
  }

  static std::set<uint32_t> keys(const o2::trd::ChamberSpacePoint& x, bool print=false)
  {
    uint32_t key = 100 * x.getDetector() + x.getPadRow();
    return {key};
  }

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
  static uint32_t key(const T& x, bool print=false)
  {
    if(print)std::cout << " calculate key : " << 1000 * x.getDetector() + 8 * x.getPadRow() + 4 * (x.getROB() % 2) + x.getMCM() % 4 << " for det:" << x.getDetector() << " padrow:" << x.getPadRow() << " rob:" << x.getROB() << " mcm:" << x.getMCM() << std::endl;
    return 1000 * x.getDetector() + 8 * x.getPadRow() + 4 * (x.getROB() % 2) + x.getMCM() % 4;
  }

  static std::set<uint32_t> keys(const o2::trd::ChamberSpacePoint& x, bool print=false)
  {
    uint32_t detrow = 1000 * x.getDetector() + 8 * x.getPadRow();
    uint32_t mcmcol = uint32_t(x.getPadCol() / float(o2::trd::constants::NCOLMCM));

    // float c = x.getPadCol() - float(mcmcol * o2::trd::constants::NCOLMCM);
    float c = x.getMCMChannel(mcmcol);
    if(print && x.getDetector()==9)std::cout << " calculate key : " << 1000 * x.getDetector() + 8 * x.getPadRow() + 4 * (x.getROB() % 2) + x.getMCM() % 4 << " for det:" << x.getDetector() << " padrow:" << x.getPadRow() << " rob:" << x.getROB() << " mcm:" << x.getMCM() << std::endl;
    if(print && x.getDetector()==9)std::cout << " calculates key : detrow : " << detrow << " mcmcol : " << mcmcol << std::endl;

    if (c >= 19.0 && mcmcol >= 1) {

      if(print&& x.getDetector()==9)std::cout << "returning detrow+mcmcol-1 : " << detrow+mcmcol-1 << " detrow+mcmcol : " << detrow+mcmcol << std::endl;
      return {detrow + mcmcol - 1, detrow + mcmcol};
    } else if (c <= 1.0 && mcmcol <= 6) {
      if(print&& x.getDetector()==9)std::cout << "returning detrow+mcmcol : " << detrow+mcmcol << " detrow+mcmcol+1 : " << detrow+mcmcol+1 << std::endl;
      return {detrow + mcmcol, detrow + mcmcol + 1};
    } else {
      if(print&& x.getDetector()==9)std::cout << "returning detrow+mcmcol : " << detrow+mcmcol << std::endl;
      return {detrow + mcmcol};
    }
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


std::vector<o2::TrackReference> RawDataSpan::makeMCTrackReferences()
{
  // define a struct to keep track of the first and last MC hit of one track in one chamber
  struct TrackReferencesInfo {
    // The first reference is hopefully the entry point
    size_t firstref{0}; //
    // The last reference should be the exit point of the amplification region
    size_t lastref{0};
    float trackid{-999.9}; // local x cordinate of the first trackef
    float start{-100.0};
    float end{999.9};      // local x cordinate of the last trackref
  };
  // Keep information about found track references in a map indexed by track ID and detector number.
  // If the span only covers (part of) a detector, the detector information is redundant, but in
  // the case of processing a whole event, the distinction by detector will be needed.
  std::map<std::pair<int, int>, TrackReferencesInfo> trackReferenceInfo;

/*  for (int iTrackRef = 0; iTrackRef < trackrefs.size(); ++iTrackRef) {
    auto ref = trackrefs[iTrackRef];

    // we look for track references classified as entering the drift region
    if ((ref.getUserId() & 0x3) == 0x1) {
      // The first hit is the hit closest to the anode region, i.e. with the largest x coordinate.
      auto id = std::make_pair(ref.getTrackID(), ref.getUserId() >> 2);
      if (ref.X() > trackReferenceInfo[id].start) {
        trackReferenceInfo[id].firstref = iTrackRef;
        trackReferenceInfo[id].firstref = ref.X();
      }
      // The last hit is the hit closest to the radiator, i.e. with the smallest x coordinate.
      if (ref.X() < trackReferenceInfo[id].end) {
        trackReferenceInfo[id].trackid = iTrackRef;
        trackReferenceInfo[id].end = ref.X();
      }
    }
  } // trackreference loop
*/
  std::vector<o2::TrackReference> trackReferences;
  for (auto x : trackReferenceInfo) {
    auto trackid = x.first.first;
    auto detector = x.first.second;
    auto firstref = hits[x.second.firstref];
    auto lastref = hits[x.second.lastref];
    //trackReferences.emplace_back(firstref, lastref, trackid);
  }
  return trackReferences;
}

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
//    if (hit.isFromDriftRegion()) {
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
 //   }
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

std::vector<TrackSegment> RawDataSpan::makeMCTrackSegmentsGeant()
{
  std::vector<TrackSegment> trackSegments;
// go through the TrackletSegmentInfo and match on time and return the relevant segemnts to this time period.
  std::cout << " in " << __func__ << " trackrefhits has : " << trackrefsegments.size() << " entries " << std::endl;
  for(auto &segmentinfo : trackrefsegments ){
    //convert these to TrackSegments
     
  }
 return trackSegments;

}
std::vector<TrackSegment> RawDataSpan::makeMCTrackSegmentsEntryExit()
{
  // define a struct to keep track of the first and last MC hit of one track in one chamber
  struct SegmentInfo {
    //enering and exiting trackreferences.
    size_t entertrackref{0}; //
    size_t exittrackref{0};
    float start{-999.9}; // local x cordinate of the first trackref.
    float end{999.9};    // local x cordinate of the last trackref.
  };
  // Keep information about found track segments in a map indexed by track ID and detector number.
  // If the span only covers (part of) a detector, the detector information is redundant, but in
  // the case of processing a whole event, the distinction by detector will be needed.
  std::map<std::pair<int, int>, SegmentInfo> trackRefSegmentInfo;
  int trackrefcounter=0; 
 // std::cout << " trackrefs has : " << trackrefs.size() << " entries " << std::endl;
  std::cout << " trackrefs has : " << trackrefsegments.size() << " entries " << std::endl;
  
  for (auto &trackref : trackrefsegments) {
    //   LOGP(info, "trackrefs  det id : {} ", trackref.getDetectorId());
    // only copy the trd trackrefs.
    // a TrackID is for a specific geant track, so it and detectorId uniquely identifies the trackref entry and exit points, defined below. 
    auto id = std::make_pair(trackref.mEnter.getTrackID(), trackref.mEnter.getUserId()>>2);
      //case 0x1:
      // direction entering
      trackRefSegmentInfo[id].entertrackref=trackrefcounter;
       //   std::cout << " Now to look for corresponding hit of trackref entering" << std::endl;
      for (int iHit = 0; iHit < hits.size(); ++iHit) {
        auto hit = hits[iHit];
        float distance = sqrt(pow((hit.getX() - trackref.mEnter.X()),2) +  pow(( hit.getY()- trackref.mEnter.Y()),2)+ pow( (hit.getZ() - trackref.mEnter.Z()),2));
        //std::cout << " distance of " << distance << " for " <<hit.getX()<<":"<<hit.getY()<<":"<<hit.getZ()<<" == "<< trackref.mEnter.X()<<":"<<trackref.mEnter.Y()<<":"<<trackref.mEnter.Z() << std::endl;

      std::array<float,3> rct;
      Hit a=convertTrackReferenceToHit(trackref.mEnter,trackref.mEnter.getTrackID(),trackref.mEnter.getUserId()>>2);
      //ChamberSpacePoint(;
       // std::cout << " distance of " << distance << " for " <<hit.getX()<<":"<<hit.getY()<<":"<<hit.getZ()<<" == "<< trackref.mEnter.X()<<":"<<trackref.mEnter.Y()<<":"<<trackref.mEnter.Z() << std::endl;
       // std::cout << " distance of " << distance << " for " <<hit.getX()<<":"<<hit.getY()<<":"<<hit.getZ()<<" == "<< a.GetX()<<":"<<a.GetY()<<":"<<a.GetZ() << std::endl;
       // std::cout << " distance of " << distance << " for " <<hit.getX()<<":"<<hit.getY()<<":"<<hit.getZ()<<" == "<< a.getLocalC()<<":"<<a.getLocalR()<<":"<<a.getLocalT() << std::endl;
        if(hit.getX() == trackref.mEnter.X() && hit.getY()== trackref.mEnter.Y() && hit.getZ() == trackref.mEnter.Z()) {
          auto hit = hits[iHit];
          trackRefSegmentInfo[id].entertrackref=iHit;
         // std::cout << " found corresponding hit for the in at" << std::endl;
        }
      }
      trackRefSegmentInfo[id].exittrackref=trackrefcounter;
      //exitunmatched.push_back(o2::math_utils::Point3D<float>(
      //   trackref.X(), trackref.Y(), trackref.Z()));
      // direction exiting
      //std::cout << " Now to look for corresponding hit of trackref exiting" << std::endl;
      for (int iHit = 0; iHit < hits.size(); ++iHit) {
        auto hit = hits[iHit];
        float distance = sqrt(pow((hit.getX() - trackref.mExit.X()),2) +  pow(( hit.getY()- trackref.mExit.Y()),2)+ pow( (hit.getZ() - trackref.mExit.Z()),2));
       // std::cout << " distance of " << distance << " for " <<hit.getX()<<":"<<hit.getY()<<":"<<hit.getZ()<<" == "<< trackref.mExit.X()<<":"<<trackref.mExit.Y()<<":"<<trackref.mExit.Z() << std::endl;
      std::array<float,3> rct;
      Hit a=convertTrackReferenceToHit(trackref.mEnter,trackref.mEnter.getTrackID(),trackref.mEnter.getUserId()>>2);
        //std::cout << " distance of " << distance << " for " <<hit.getX()<<":"<<hit.getY()<<":"<<hit.getZ()<<" == "<< a.getLocalC()<<":"<<a.getLocalR()<<":"<<a.getLocalT() << std::endl;
        if(hit.getX() == trackref.mExit.X() && hit.getY()== trackref.mExit.Y() && hit.getZ() == trackref.mExit.Z()) {
          auto hit = hits[iHit];
          trackRefSegmentInfo[id].exittrackref=iHit;
        //  std::cout << " found corresponding hit for the out at" << std::endl;
        }
      }
    trackrefcounter++;
    }
//  for (const auto &trackref : trackrefs) {
//    //   LOGP(info, "trackrefs  det id : {} ", trackref.getDetectorId());
//    if (trackref.getDetector() == 2) {
//    // only copy the trd trackrefs.
//    // a TrackID is for a specific geant track, so it and detectorId uniquely identifies the trackref entry and exit points, defined below. 
//      auto id = std::make_pair(trackref.getTrackID(), trackref.getDetector());
//      std::cout << " trackref userid : " << (trackref.getUserId() & 0x3)<< std::endl;
//      switch (trackref.getUserId() & 0x3) {
//        case 0x1:
//        // direction entering
//        trackRefSegmentInfo[id].entertrackref=trackrefcounter;
//          //enterunmatched.push_back(o2::math_utils::Point3D<float>(
//           // trackref.X(), trackref.Y(), trackref.Z()));
//            std::cout << " Now to look for corresponding hit of trackref entering" << std::endl;
//        for (int iHit = 0; iHit < hits.size(); ++iHit) {
//          auto hit = hits[iHit];
//          float distance = sqrt(pow((hit.getX() - trackref.getX()),2) +  pow(( hit.getY()- trackref.getY()),2)+ pow( (hit.getZ() - trackref.getZ()),2));
//          std::cout << " distance of " << distance << " for " <<hit.getX()<<":"<<hit.getY()<<":"<<hit.getZ()<<" == "<< trackref.getX()<<":"<<trackref.getY()<<":"<<trackref.getZ() << std::endl;
//
//          if(hit.getX() == trackref.getX() && hit.getY()== trackref.getY() && hit.getZ() == trackref.getZ()) {
//            auto hit = hits[iHit];
//            trackRefSegmentInfo[id].entertrackref=iHit;
//            std::cout << " found corresponding hit for the in at" << std::endl;
//          }
//        }
//        break;
//      case 0x2:
//        trackRefSegmentInfo[id].exittrackref=trackrefcounter;
//        //exitunmatched.push_back(o2::math_utils::Point3D<float>(
//        //   trackref.X(), trackref.Y(), trackref.Z()));
//        // direction exiting
//        std::cout << " Now to look for corresponding hit of trackref exiting" << std::endl;
//        for (int iHit = 0; iHit < hits.size(); ++iHit) {
//          auto hit = hits[iHit];
//          float distance = sqrt(pow((hit.getX() - trackref.getX()),2) +  pow(( hit.getY()- trackref.getY()),2)+ pow( (hit.getZ() - trackref.getZ()),2));
//          std::cout << " distance of " << distance << " for " <<hit.getX()<<":"<<hit.getY()<<":"<<hit.getZ()<<" == "<< trackref.getX()<<":"<<trackref.getY()<<":"<<trackref.getZ() << std::endl;
//          if(hit.getX() == trackref.getX() && hit.getY()== trackref.getY() && hit.getZ() == trackref.getZ()) {
//            auto hit = hits[iHit];
//            trackRefSegmentInfo[id].exittrackref=iHit;
//            std::cout << " found corresponding hit for the out at" << std::endl;
//          }
//        }
//        break;
//      }
//    }
//    trackrefcounter++;
//  }
  //we want hits not trackreferences, to make our lives easier.
  std::vector<TrackSegment> trackSegments;

  auto ctrans = o2::trd::CoordinateTransformer::instance();
  for (auto x : trackRefSegmentInfo) {
    HitPoint enter,exit;
 //   enter=ctrans->MakeSpacePoint(trackrefs[x.second.entertrackref]);
  //  exit=ctrans->MakeSpacePoint(trackrefs[x.second.exittrackref]);
    auto trackid = x.first.first;
    auto detector = x.first.second;
//    o2::trd::Hit entering= convertTrackReferenceToHit(enter,trackid,detector);
//    o2::trd::Hit leaving= convertTrackReferenceToHit(exit,trackid,detector);
    //find the corresponding hit to each.
    //// this is massively expensive but will fix later.
    auto firsthit = hits[x.second.entertrackref];
    auto lasthit = hits[x.second.exittrackref];
    trackSegments.emplace_back(firsthit, lasthit, trackid);
    }
  return trackSegments;
}

Hit RawDataSpan::convertTrackReferenceToHit(o2::TrackReference &ref, 
                                               int trackid, int detector)
{
  //std::cout << __func__ << " " << ref << " trackid: " << trackid << " det: " << detector << std::endl;
  std::cout << std::dec;
  float tof = ref.getTime() * 1e6; // The time of flight in micro-seconds
  double pos[3] = {ref.X(),ref.Y(),ref.Z()};
  double loc[3] = {-99, -99, -99};
  // gGeoManager->Export("geometry.root"); // is there a corresponding import ?
  if(gGeoManager == nullptr){
    gGeoManager = new TGeoManager();
    gGeoManager->Import("o2sim_geometry.root");
  }
  auto node = gGeoManager->FindNode(ref.X(), ref.Y(), ref.Z());
  gGeoManager->MasterToLocal(pos, loc); // Go to the local coordinate system (locR, locC, locT)
 // if (!node) {
 //   std::cout << "node name : unknown\n";
 // }
  auto userid=ref.getUserId();
  auto enterleave=userid&0x3;
 // std::cout << "node name : " << node->GetVolume()->GetName() <<  " userid: " << userid << " ref enter/leave : " << enterleave << "\n";
  //std::cout << "node name : " << node->GetVolume()->GetName() <<  " ref enter/leave : " << (ref.getUserId()&0x3) << "\n";
  // find node
  float locC = loc[0], locR = loc[1], locT = loc[2];
  locT = locT - 0.5 * (Geometry::drThick() + Geometry::amThick()); // Relative to middle of amplification region
  float charge=0.;// nominal error charge as this is not a hit.
  //std::cout << "converTrackReferenceToHit : " << pos[0] <<"," << pos[1] <<"," << pos[2] <<" :: " << locC << "," << locR << "," << locT <<  " local : " << ref.LocalX() << ":" << ref.LocalY() << std::endl;
  /********/
  /*int sector = Geometry::getSector(ref.getUserId()>>2);
  int layer = Geometry::getLayer(detector);
  int stack = Geometry::getStack(detector);
  float phi = 2.0 * TMath::Pi() / (float)o2::trd::constants::NSECTOR * ((float)sector + 0.5);
  double loc1[3] = {-99, -99, -99};
  double pos1[3] = {ref.X(),ref.Y(),ref.Z()};
  loc1[0] = pos1[0] * TMath::Cos(phi) - pos1[1] * TMath::Sin(phi);
  loc1[1] = -1.0*pos1[0] * TMath::Sin(phi) + pos1[1] * TMath::Cos(phi);
  loc1[2] = pos1[2];*/
 // std::cout << fmt::format("converTrackReferenceToHit manual phi:{} : {:.5}:{:.5}:{:.5} :: {:.5}:{:.5}:{:.5}",phi,pos1[0],pos1[1],pos1[2],loc1[0],loc1[1],loc1[2]);
/*  double rmin = mgeo->getTime0(layer);
  double rmax = mgeo->getTime0(layer) - mgeo->drThick() - mgeo->amThick();
  double ymin = -mgeo->getChamberWidth(layer) / 2;
  double ymax = mgeo->getChamberWidth(layer) / 2;
  double zmin = -mgeo->getChamberLength(layer, stack) / 2;
  double zmax = mgeo->getChamberLength(layer, stack) / 2;
 */
  /********/


  o2::trd::Hit thehit(ref.X(),ref.Y(),ref.Z(),locC,locR,locT, tof, charge,  trackid, detector,true);
  return thehit;
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
    mDataTree->SetBranchAddress("TPCITS", &mTracks);
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
    mCollisionContext->printCollisionSummary();
  }

  // We create the MC TTree using event header and tracks from the kinematics file
  if (std::filesystem::exists(dir / "o2sim_Kine.root")) {
    mMCFile = new TFile((dir / "o2sim_Kine.root").c_str());
    mMCFile->GetObject("o2sim", mMCTree);
    mMCTree->SetBranchAddress("MCEventHeader.", &mMCEventHeader);
    mMCTree->SetBranchAddress("MCTrack", &mMCTracks);
    mMCTree->SetBranchAddress("TrackRefs", &mMCTrackReferences);
  }
  //process trackreferences into MCTrackletSegmentInfo
  processTrackReferences();
  // We then add the TRD hits to the MC tree
  if (mMCFile && std::filesystem::exists(dir / "o2sim_HitsTRD.root")) {
    mMCTree->AddFriend("o2sim", (dir / "o2sim_HitsTRD.root").c_str());
    mMCTree->SetBranchAddress("TRDHit", &mHits);
  }
}

void RawDataManager::processTrackReferences()
{
  //Trackreferences come in ordered by :  trackid, then time.
  std::map<std::pair<uint32_t,uint32_t>, MCTrackletSegmentInfo> trackreferences;
  int nev = mMCTree->GetEntries();
  // take the incoming trackreferences, remove those not applicable to TRD, and order them by det, trackid, and time.
  for (int iev = 0; iev < nev; ++iev) {
    mMCTree->GetEvent(iev); // all trackreferences are in a single branch.
    //  sort track ref by time and enter before exit.
    for (const auto &trackref : *mMCTrackReferences) {
      if (trackref.getDetectorId() == 2) {
//      LOGP(info, "trackrefs  det id : {} trackid: {} time : {} ", trackref.getDetectorId(), trackref.getTrackID(),trackref.getTime());
        uint32_t trackid=trackref.getTrackID();
        uint32_t det;
        auto key = std::make_pair(trackref.getUserId()>>2,trackref.getTrackID());
        // only copy the trd trackrefs.
        // std::map<int, o2::TrackReference> inner;
        // inner.insert(std::make_pair(det, trackref));
        // trackreferences.insert(trackref.getTrackID(), inner);
        if(!trackreferences.contains(key)){
        trackreferences[key]=MCTrackletSegmentInfo();  // map[key].   trackrefsvectorwithextra.push_back(trackref);
        if((trackref.getUserId() & 0x3) == 0x2) {
            // we have an exit but no preceding enter ??
          std::cout << "Exit but no entry key : " << key.first << ":"<< key.second << " Setting entry point with entry point of : " <<  trackreferences[key].mEnter << " and corresponding exit of : " << trackreferences[key].mExit << std::endl;
          }
        }
        switch (trackref.getUserId() & 0x3) {
        case 0x1:
          // direction entering
          //trackreferences[key].setEntry(trackref.X(),trackref.Y(),trackref.Z(),trackid, det);
          trackreferences[key].mEnter=trackref;//setEntry(trackref.X(),trackref.Y(),trackref.Z(),trackid, det);
          //std::cout << " key : " << key.first << ":"<< key.second << " Setting entry point with entry point of : " <<  trackreferences[key].mEnter << " and corresponding exit of : " << trackreferences[key].mExit << std::endl;
          break;
        case 0x2:
         // trackreferences[key].setExit(trackref.X(),trackref.Y(),trackref.Z(), trackid,det);
          trackreferences[key].mExit=trackref;
           //std::cout << " key : " << key.first << ":"<< key.second << " Setting exit point with exit point of : " <<  trackreferences[key].mExit << " and corresponding enter of : " << trackreferences[key].mEnter << std::endl;
          // direction exiting
          break;
        }
      }
    }
  }
// now walk through the map and insert the entry/exit ref points into std::vector<MCTrackletSegmentInfo> mMCTrackletSegmentInfo{0};
  // Iterate using C++17 facilities
  for (const auto& [key, value] : trackreferences){
      //std::cout << '[' << key << "] = " << value << "; "; 
    auto det = key.first;
    auto trackid = key.second;
    if(trackreferences[key].mEnter.getTrackID() == 0){
      std::cout << " something wrong, trackid for key.1 "<< key.first << " key.2" << key.second << " is : " << trackreferences[key].mEnter.getTrackID() << " " << trackreferences[key].mExit.getTrackID() << std::endl;
    }
    auto ctrans = o2::trd::CoordinateTransformer::instance();
    trackreferences[key].setMidPoint();
   // trackreferences[key].setMidSpace(ctrans->MakeSpacePoint(mExit), 0.0,mExit.getTrackID());
    //mMCTrackletSegmentInfo.emplace_back(trackreferences[key].mEnter,trackreferences[key].mExit,det,trackid);
    if(trackreferences[key].isGood()){
      o2::TrackReference ent=trackreferences[key].mEnter;
      o2::TrackReference ext=trackreferences[key].mExit;
      trackreferences[key].setMidSpacePoint();
      ChamberSpacePoint enter=ctrans->MakeSpacePoint(ent);
      ChamberSpacePoint exit=ctrans->MakeSpacePoint(ext);

      //ChamberSpacePoint mid=ctrans->MakeSpacePoint(trackreferences[key].mMidSpace);
      
      //std::cout << "ZZZ det: " << exit.getDetector() << "\n";
      trackreferences[key].setEnterSpace(enter, 0.0,trackreferences[key].mEnter.getTrackID());
      trackreferences[key].setExitSpace(exit, 0.0,trackreferences[key].mExit.getTrackID());
      trackreferences[key].setMidSpace(mid, 0.0,trackreferences[key].mExit.getTrackID());

      if(trackreferences[key].getDetector()==9){
        LOGP(info,"ZZ Det : 9 Trackid : {} : {}:{}:{} padrow {} padcol {} mcm {} rob {}",trackreferences[key].mEnter.getTrackID(),
             trackreferences[key].mMidSpace.getX(),trackreferences[key].mMidSpace.getY(),trackreferences[key].mMidSpace.getZ(),trackreferences[key].mMidSpace.getPadRow(), trackreferences[key].mMidSpace.getPadCol(),trackreferences[key].mMidSpace.getMCM(),trackreferences[key].mMidSpace.getROB());
      }
      // remove those that dont have a enter or exit
    mMCTrackletSegmentInfo.emplace_back(value);
     // std::cout << "GOOD: " << trackreferences[key] << "\n";
    }
    //else {
    //  std::cout << "BAD: " << trackreferences[key] << "\n";
   // }
    //mMCTrackletSegmentInfo.emplace_back(value.mEnter,value.mExit,det,trackid);
//    std::cout << "trackid for key.1 "<< key.first << " key.2" << key.second << "\n Enter: " << trackreferences[key].mEnter  << "\n Exit: " << trackreferences[key].mExit << "\n";
//    std::cout << "    det:" << trackreferences[key].getDetector() << " padrow:"<< trackreferences[key].getPadRow() <<" padcol:"<< trackreferences[key].getPadRow()<<" rob:"<< trackreferences[key].getROB()<<" mcm:"<< trackreferences[key].getMCM() << "\n";
  }
}

bool RawDataManager::nextTimeFrame()
{
  if (!mDataTree->GetEntry(mTimeFrameNo)) {
    // loading time frame will fail at end of file
    return false;
  }

  mEventNo = 0;
  mTimeFrameNo++;

  O2INFO("Loaded data for time frame #%d with %d TRD trigger records, %d digits and %d tracklets",
         mTimeFrameNo, mTrgRecords->size(), mDigits->size(), mTracklets->size());

  return true;
}

bool RawDataManager::nextEvent()
{
  // get the next trigger record
  if (mEventNo >= mTrgRecords->size()) {
    return false;
  }
  mTriggerRecord = mTrgRecords->at(mEventNo);
  O2INFO("Processing event: orbit %d bc %04d with %d digits and %d tracklets",
         mTriggerRecord.getBCData().orbit, mTriggerRecord.getBCData().bc,
         mTriggerRecord.getNumberOfDigits(), mTriggerRecord.getNumberOfTracklets());

  if (mCollisionContext) {

    // clear MC data
    mHitPoints.clear();

    for (int i = 0; i < mCollisionContext->getNCollisions(); ++i) {
      auto evrec = mCollisionContext->getEventRecords()[i];
      if (abs(mTriggerRecord.getBCData().differenceInBCNS(evrec)) <= 3000) {
        // if (mMCReader) {
        mMCTree->GetEntry(i);
        // }

        O2INFO("Loaded matching MC event #%d with time offset %f ns and %d hits and %d trackreferences ",
               i, mTriggerRecord.getBCData().differenceInBCNS(evrec), mHits->size(), mMCTrackReferences->size());
        O2INFO("Loaded matching MC event #%d with time offset %f ns and %d hits and %d trackreferences ",
               i, mTriggerRecord.getBCData().differenceInBCNS(evrec), mHits->size(), mMCTrackReferences->size());

        // convert hits to spacepoints
        auto ctrans = o2::trd::CoordinateTransformer::instance();
        for (auto& hit : *mHits) {
         // if(hit.getDetector()==9){
            if(hit.GetTrackID() == 19179) LOGP(info,"ZZ det 9 trackid : {} hit : {} {} {} : {} {} {} ",hit.GetTrackID(),hit.GetX(),hit.GetY(),hit.GetZ(),hit.getLocalC(),hit.getLocalR(),hit.getLocalT());
          //}
          mHitPoints.emplace_back(ctrans->MakeSpacePoint(hit), hit.GetCharge(),hit.GetTrackID());
         // std::cout << " building mHitPoints trackid : " << hit.GetTrackID() << std::endl;
        }
       /* for (auto& trackref : *mMCTrackReferences) {
          std::cout << "added track ref at : " << trackref.X()<<":"<<trackref.Y()<<":"<<trackref.Z() << std::endl;
          mTrackReferences.emplace_back(ctrans->MakeSpacePoint(trackref));
        }*/
      }
    }
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
  
  //ev.trackrefs = boost::make_iterator_range(mTrackReferences.begin(), mTrackReferences.end());
  //ev.trackrefhits = boost::make_iterator_range(mMCTrackReferences.begin(), mMCTrackReferences.end());
  ev.trackrefsegments = boost::make_iterator_range(mMCTrackletSegmentInfo.begin(), mMCTrackletSegmentInfo.end());
  
  auto evtime = getTriggerTime();
  std::cout << "event time : " << evtime << std::endl;
  std::vector<uint32_t> trackids;
  for(auto hit : ev.hits){
    trackids.push_back(hit.getTrackID());

  }
  int idcountpre=trackids.size();
  sort(trackids.begin(),trackids.end());
  auto uniqueids=std::unique(trackids.begin(),trackids.end());
  trackids.erase(uniqueids,trackids.end());
  int idcountpost=trackids.size();
  std::cout << " pre trackids:" << idcountpre << " post trackids:" << idcountpost << std::endl;

  // if (tpctracks) {
  //   for (auto &track : *mTpcTracks) {
  //     //   // auto tracktime = track.getTimeMUS().getTimeStamp();
  //     auto dtime = track.getTime0() / 5.0 - evtime;
  //     if (dtime > mMatchTimeMinTPC && dtime < mMatchTimeMaxTPC) {
  //       ev.mTpcTracks.push_back(track);
  //     }
  //   }
  // }

  if (mTracks) {
    for (auto& track : *mTracks) {
      //   // auto tracktime = track.getTimeMUS().getTimeStamp();
      // auto dtime = track.getTimeMUS().getTimeStamp() - evtime;
      // if (dtime > mMatchTimeMinTPC && dtime < mMatchTimeMaxTPC) {
      //   ev.tracks.push_back(track);

      // for(int ly=0; ly<6; ++ly) {
      //   auto point = extra.extrapolate(track.getParamOut(), ly);
      //   if (point.isValid()) {
      //     ev.evtrackpoints.push_back(point);
      //   }
      // }
      // }
    }
  }

  // ev.trackpoints.begin() = ev.evtrackpoints.begin();
  // ev.trackpoints.end() = ev.evtrackpoints.end();
  std::cout << " returning event" << std::endl;
  return ev;
}

o2::dataformats::TFIDInfo RawDataManager::getTimeFrameInfo()
{
  if (mTFIDs) {
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
  if (mDataTree->GetFriend("TRDDigit")) {
    out << "digits" << std::endl;
  }
  if (mDataTree->GetFriend("TPCITS")) {
    out << "tpc its matches" << std::endl;
  }

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




int MCTrackletSegmentInfo::getPadRow() const
{
  auto enterlocalx=mEnter.LocalX();
  auto exitlocalx=mExit.LocalX();
  //: mID(id), mDetector(detector), mX(x), mY(y), mZ(z), mPadrow(rct[0]), mPadcol(rct[1]), mTimebin(rct[2]), mInDrift(inDrift){};
  double pos[3] = {mEnter.X(),mEnter.Y(),mEnter.Z()};
  double posexit[3] = {mExit.X(),mExit.Y(),mExit.Z()};
  double loc1[3] = {-99, -99, -99};
  double loc[3] = {-99, -99, -99};
  // gGeoManager->Export("geometry.root"); // is there a corresponding import ?
  if(gGeoManager == nullptr){
    gGeoManager = new TGeoManager();
    gGeoManager->Import("o2sim_geometry.root");
  }
  auto node = gGeoManager->FindNode(mMidPoint[0],mMidPoint[1],mMidPoint[2]);
  gGeoManager->MasterToLocal(posexit, loc); // Go to the local coordinate system (locR, locC, locT)
  //std::cout <<__func__ <<"aapos : " << pos[0] << ":"<<pos[1]<<":"<<pos[2] << " :: " <<loc[0] << ":" << loc[1] << ":" <<loc[2] << "  midpoint:" << mMidPoint[0] << ":" << mMidPoint[1] <<":" << mMidPoint[2] <<" " << node->GetVolume()->GetName() << "\n";
  
  //return (int)loc[1];
  return mExitSpace.getPadRow();
};

int MCTrackletSegmentInfo::getPadCol() const
{
  double pos[3] = {mEnter.X(),mEnter.Y(),mEnter.Z()};
  double posexit[3] = {mExit.X(),mExit.Y(),mExit.Z()};
  double loc[3] = {-99, -99, -99};
  // gGeoManager->Export("geometry.root"); // is there a corresponding import ?
  if(gGeoManager == nullptr){
    gGeoManager = new TGeoManager();
    gGeoManager->Import("o2sim_geometry.root");
  }
  auto node = gGeoManager->FindNode(mMidPoint[0],mMidPoint[1],mMidPoint[2]);
  gGeoManager->MasterToLocal(posexit, loc); // Go to the local coordinate system (locR, locC, locT)
  float locC = loc[0], locR = loc[1], locT = loc[2];
  //std::cout <<__func__ <<"pos : " << pos[0] << ":"<<pos[1]<<":"<<pos[2] << " :: " <<loc[0] << ":" << loc[1] << ":" <<loc[2] << "  midpoint:" << mMidPoint[0] << ":" << mMidPoint[1] <<":" << mMidPoint[2] <<" " << node->GetVolume()->GetName() << "\n";
  //return (int)locC;
  return mMidSpace.getPadCol();///(int)locC;
}

float MCTrackletSegmentInfo::getMCMf() const
{
  return o2::trd::HelperMethods::getMCMfromPad(getPadRow(),getPadCol());
};

float MCTrackletSegmentInfo::getROBf() const
{
  return o2::trd::HelperMethods::getROBfromPad(getPadRow(),getPadCol());
};

int MCTrackletSegmentInfo::getMCM() const
{
  auto padrow = getPadRow();
  auto padcol = getPadCol();
  //std::cout << __func__ << " padrow:" << padrow << "[0:"<< constants::NROWC1<<"]  padcol:" << padcol << "[0:"<< constants::NCOLUMN << "] \n";
  //std::cout <<__func__ << " midpoint:" << mMidPoint[0] << ":" << mMidPoint[1] <<":" << mMidPoint[2] <<" " << gGeoManager->FindNode(mMidPoint[0],mMidPoint[1],mMidPoint[2])->GetVolume()->GetName() << "\n";
  return mMidSpace.getMCM();//o2::trd::HelperMethods::getMCMfromPad(getPadRow(),getPadCol());
};

int MCTrackletSegmentInfo::getROB() const
{
  return mMidSpace.getROB();//o2::trd::HelperMethods::getROBfromPad(getPadRow(),getPadCol());
};


namespace o2::trd
{

std::ostream& operator<<(std::ostream& os, const MCTrackletSegmentInfo& p)
{
  int sector = p.getDetector() / 30;
  int stack = (p.getDetector() % 30) / 6;
  int layer = p.getDetector() % 6;
  os << fmt::format("mctrackletsegmentinfo : ({:.3}:{:.3}:{:.3}--{:.3}:{:.3}:{:.3} : mid:{:.3}:{:.3}:{:.3})",p.mEnter.X(),p.mEnter.Y(),p.mEnter.Z(),p.mExit.X(),p.mExit.Y(),p.mExit.Z(), p.mMidPoint[0],p.mMidPoint[1],p.mMidPoint[2] );
  return os;
}

}
