// Copyright 2019-2023 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef ALICEO2_TRD_RawDataManager_H_
#define ALICEO2_TRD_RawDataManager_H_

///
/// \file   RawDataManager.h
/// \author Thomas Dietel, tom@dietel.net
///

#include "DataFormatsTRD/Digit.h"
#include "DataFormatsTRD/Tracklet64.h"
#include "DataFormatsTRD/TriggerRecord.h"
#include "DataFormatsTRD/Hit.h"
#include "DataFormatsTRD/TrackTriggerRecord.h"
#include "DataFormatsTRD/TrackTRD.h"
#include "TRDBase/TrackletTransformer.h"

#include "TRDQC/CoordinateTransformer.h"

#include "CommonDataFormat/TFIDInfo.h"
#include "ReconstructionDataFormats/TrackTPCITS.h"
#include "ReconstructionDataFormats/MatchInfoTOF.h"
#include "SimulationDataFormat/DigitizationContext.h"
#include "SimulationDataFormat/MCEventHeader.h"
#include "SimulationDataFormat/MCTrack.h"
#include "DetectorsBase/Propagator.h"

// #include <TTreeReaderArray.h>

#include <boost/range/iterator_range_core.hpp>
#include <vector>
#include <filesystem>
#include <boost/range.hpp>

class TFile;
class TTree;

static constexpr  int32_t mMaxTriggers=100;
namespace o2::trd
{

/// RawDataSpan holds ranges pointing to parts of other containers
/// This class helps with restricting the view of digits/trackets/... to part
/// of a timeframe, e.g. data from one event, detector, padrow or MCM.
/// It also provides methods to help with iterating over the data
/// The main data members are public, and can be accessed without setters/getters.
struct RawDataSpan {
 public:
  boost::iterator_range<std::vector<o2::trd::Digit>::iterator> digits;
  boost::iterator_range<std::vector<o2::trd::Tracklet64>::iterator> tracklets;
  boost::iterator_range<std::vector<HitPoint>::iterator> hits;
 /* boost::iterator_range<std::vector<o2::dataformats::TrackTPCITS>::iterator> tracks_itstpc;
  boost::iterator_range<std::vector<o2::dataformats::MatchInfoTOF>::iterator> tracks_itstpc_tof;
  boost::iterator_range<std::vector<o2::dataformats::MatchInfoTOF>::iterator> tracks_itstpctrd_tof;
  boost::iterator_range<std::vector<o2::trd::TrackTRD>::iterator> tracks_itstpctof_trd;
  boost::iterator_range<std::vector<o2::trd::TrackTRD>::iterator> tracks_itstpc_trd;
  boost::iterator_range<std::vector<o2::trd::TrackTRD>::iterator> tracks_tpc_trd;
*/

  boost::iterator_range<std::vector<TrackSegment>::iterator> tracks_itstpc_seg;
  boost::iterator_range<std::vector<TrackSegment>::iterator> tracks_itstpc_tof_seg;
  boost::iterator_range<std::vector<TrackSegment>::iterator> tracks_itstpctrd_tof_seg;
  boost::iterator_range<std::vector<TrackSegment>::iterator> tracks_itstpctof_trd_seg;
  boost::iterator_range<std::vector<TrackSegment>::iterator> tracks_itstpc_trd_seg;
  boost::iterator_range<std::vector<TrackSegment>::iterator> tracks_tpc_trd_seg;

  /// Sort digits, tracklets and space points by detector, pad row, column
  /// The digits, tracklets, hits and other future data members must be sorted
  /// for the IterateBy method.
  void sort();

  /// Return a vector with data spans, split according to the keyfunc struct
  /// The keyfunc struct must have a method `key` to calculate a key for tracklets
  /// and digits that uniquely identifies the span this digit or tracklet will be
  /// part of. The key must be monotonically increasing for the digits and tracklets
  // stored in the raw data span.
  template <typename keyfunc>
  std::vector<RawDataSpan> iterateBy();

  std::vector<RawDataSpan> iterateByMCM();
  std::vector<RawDataSpan> iterateByPadRow();
  std::vector<RawDataSpan> iterateByDetector();

  std::vector<TrackSegment> makeMCTrackSegments();

  //   pair<int, int> getMaxADCsumAndChannel();
  //   int getMaxADCsum(){ return getMaxADCsumAndChannel().first; }

  //   int getDetector() { if (mDetector == -1) calculateCoordinates(); return mDetector; }
  //   int getPadRow() { if (mPadRow == -1) calculateCoordinates(); return mPadRow; }

  // protected:
  //   // The following variables cache the calculations of raw coordinates:
  //   //   non-negative values indicate the actual position
  //   //   -1 indicates that the values has not been calculated yet
  //   //   -2 indicates that the value is not unique in this span
  //   int mDetector{-1};
  //   int mPadRow{-1};

  //   void calculateCoordinates();
};

/// RawDataManager: read raw, MC and reconstruced data files and loop over them
///
/// The input for the RawDataManager are directories with raw and, optionally,
/// reconstructed and Monte-Carlo files. It scans for available files, reads
/// them and provides access functions to loop over time frames and events.
/// Supported data files:
///   - trdtracklets.root (required): TRD tracklets and trigger records
///   - trddigits.root: TRD digits - might only be available for some events
///   - o2_tfidinfo.root (reconstructed data only): time frame information
///   - o2sim_HitsTRD.root: TRD MC hits with global and TRD chamber coordinates
///
/// The following file types had some support in previous private branches,
/// but need further cleanup before integration into O2:
///   -  o2match_itstpc.root: ITS-TPC tracks
///   -  tpctracks.root: TPC-only tracks
///   -  tpctracks.root: TPC-only tracks
class RawDataManager
{

 public:
  RawDataManager(std::filesystem::path dir = ".");
  // RawDataManager(std::string dir = "./");

  // void SetMatchWindowTPC(float min, float max)
  // { mMatchTimeMinTPC=min; mMatchTimeMaxTPC=max; }

  bool nextTimeFrame(bool onlydigits=false);
  bool nextEvent();

  /// access time frame info
  o2::dataformats::TFIDInfo getTimeFrameInfo();

  // TTreeReaderArray<o2::tpc::TrackTPC> *GetTimeFrameTPCTracks() {return mTpcTracks; }
  std::vector<o2::dataformats::TrackTPCITS>* getTimeFrameTracks() { return mITSTPCTracks; }

  // access event info
  RawDataSpan getEvent();
  float getTriggerTime();

  size_t getTimeFrameNumber() { return mTimeFrameNo; }
  size_t getEventNumber() { return mEventNo; }

  o2::steer::DigitizationContext* getCollisionContext() { return mCollisionContext; }

  std::string describeFiles();
  std::string describeTimeFrame();
  std::string describeEvent();
  //walk through the available tracks and build track segments through the detectors with 10 points per detector.
 // void buildTrackSegments();
 private:
  // access to TRD digits and tracklets
  TFile* mMainFile{0}; // the main trdtracklets.root file
  TTree* mDataTree{0}; // tree and friends from digits, tracklets files
                       //  TTreeReader* mDataReader{0};

  std::vector<o2::trd::Digit>* mDigits{0};
  std::vector<o2::trd::Tracklet64>* mTracklets{0};
  std::vector<o2::trd::TriggerRecord>* mTrgRecords{0};

  // access tracks
  std::array<int32_t,mMaxTriggers*o2::trd::constants::NCHAMBER> mTrackletIndexArray;
  std::vector<o2::dataformats::TrackTPCITS>* mITSTPCTracks{0};
  std::vector<o2::dataformats::MatchInfoTOF>* mITSTPCTracks_TOF{0};
  std::vector<o2::dataformats::MatchInfoTOF>* mITSTPCTRDTracks_TOF{0};
  std::vector<o2::trd::TrackTRD>* mITSTPCTOFTracks_TRD{0};
  std::vector<o2::trd::TrackTriggerRecord>* mITSTPCTOFTracksTrigRec_TRD{0};
  std::vector<o2::trd::TrackTRD>* mITSTPCTracks_TRD{0};
  std::vector<o2::trd::TrackTriggerRecord>* mITSTPCTracksTrigRec_TRD{0};
  std::vector<o2::trd::TrackTRD>* mTPCTracks_TRD{0};
  std::vector<o2::trd::TrackTriggerRecord>* mTPCTracksTrigRec_TRD{0};
  std::vector<TrackSegment> mITSTPCTracks_segments{0};
  std::vector<TrackSegment> mITSTPCTracks_TOF_segments{0};
  std::vector<TrackSegment> mITSTPCTRDTracks_TOF_segments{0};
  std::vector<TrackSegment> mITSTPCTOFTracks_TRD_segments{0}; 
  std::vector<TrackSegment> mITSTPCTracks_TRD_segments{0}; 
  std::vector<TrackSegment> mTPCTracks_TRD_segments{0}; 
  
  bool buildTrackSegments(bool onlydigits);
  // build the tracksegments from tracking 
  o2::trd::TrackletTransformer mTransformer;
  bool propagateToLayerX(o2::dataformats::TrackTPCITS& track, float xToGo, float e, float maxStep);
  void prepareTracking();//std::array<int32_t,540*mMaxTriggers>& trdTrackletIndexArray);
  int propagateTrack(o2::dataformats::TrackTPCITS& track, float e, float maxStep, float triggertime, int& glbTrkltIdxOffset,  int collisionId);
  int32_t getSector(float alpha);
  float getAlphaOfSector(const int32_t sec);
  int32_t getDetectorNumber(const float zPos, const float alpha, const int32_t layer);
  bool isGeoFindable(o2::dataformats::TrackTPCITS& track, const int32_t layer, const float alpha, const float zShiftTrk);
  void findChambersInRoad(o2::dataformats::TrackTPCITS& track, const float roadY, const float roadZ, const int32_t iLayer, std::array<int,18>& det, const float zMax, const float alpha, const float zShiftTrk);
  bool getYZAt(float xk, float b, float& y, float& z, o2::dataformats::TrackTPCITS& track);
  bool trackMatchesCollision(const float& triggertime, const float& tracktime);
  float getTriggerTime(const o2::trd::TriggerRecord& trig, const std::vector<o2::dataformats::TFIDInfo>& tfids, const int timeframeno);
  std::array<float,540> mR{0};
  std::array<float,540>  mRdriftstart{0};
  std::array<float,540> mRdriftend{0};
  std::array<float,540>  mRamplificationend{0};

  // access to Monte-Carlo events, tracks, hits
  TFile* mMCFile{0};
  TTree* mMCTree{0};
  // TTreeReader* mMCReader{0};
  std::vector<o2::dataformats::MCEventHeader>* mMCEventHeader{0};
  std::vector<o2::MCTrackT<Float_t>>* mMCTracks{0};
  std::vector<o2::trd::Hit>* mHits{0};

  // MC hits, converted to chamber coordinates
  std::vector<o2::trd::HitPoint> mHitPoints;

  // time frame information (for data only)
  std::vector<o2::dataformats::TFIDInfo>* mTFIDs{0};

  // collision context (for MC only)
  o2::steer::DigitizationContext* mCollisionContext{0};

  // current trigger record
  o2::trd::TriggerRecord mTriggerRecord;

  // time frame and event counters
  size_t mTimeFrameNo{0}, mEventNo{0};

  // matching parameters for tracks
  // float mMatchTimeMinTPC{-10.0}, mMatchTimeMaxTPC{20.0};

  // template <typename T>
  // void addReaderArray(TTreeReaderArray<T>*& array, std::filesystem::path file, std::string tree, std::string branch);
  //
  //
  o2::trd::Geometry* mGeo{0};     // TRD geometry
                                  // 
  o2::base::Propagator* mProp{0};  
  
  float mAverageRadii[6] = {300.2f, 312.8f, 325.4f, 338.0f, 350.6f, 363.2f}; // used as default value in case no transformation matrix can be obtained
//GPUTRDGeometry mGeogpu;     // TRD geometry
// only used for those times we write the track to tracklet matching to a root tree.
   TFile *mfile=nullptr;
   TTree *moutputtree=nullptr;
//
//
//const int kNStacks=5;
//const int kNLayers=6;
//const int kNChambers=540;
};

} // namespace o2::trd

#endif // ALICEO2_TRD_RawDataManager_H_
