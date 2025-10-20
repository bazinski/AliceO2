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

#ifndef ALICEO2_TRD_COORDINATE_TRANSFORMER_H_
#define ALICEO2_TRD_COORDINATE_TRANSFORMER_H_

///
/// \file   CoordinateTransformer.h
/// \author Thomas Dietel, tom@dietel.net
///

#include "DataFormatsTRD/Hit.h"

#include <array>

namespace o2::trd
{

class Geometry;
class CoordinateTransformer;

/// A position in spatial (x,y,z) and raw/digit coordinates (det,row,col,tb).
/// This class is intended to store the converted coordinates of a space point, to avoid
/// repeated transformations between spatial and row/col/tb coordinates.
class ChamberSpacePoint
{
 public:
  ChamberSpacePoint(int det = -999) : mDetector(det){};
  ChamberSpacePoint(int id, int detector, float x, float y, float z, std::array<float, 3> rct, bool inDrift)
    : mID(id), mDetector(detector), mX(x), mY(y), mZ(z), mPadrow(rct[0]), mPadcol(rct[1]), mTimebin(rct[2]), mInDrift(inDrift){};

  /// check if the space point has been initialized
  bool isValid() const { return mDetector >= 0; }

  /// check if the space point is in the drift region
  bool isFromDriftRegion() const { return mInDrift; }

  /// return the ID of the associated track (either MC or reconstructed)
  int getID() const { return mID; }

  /// spatial x coordinate of space point
  float getX() const { return mX; }

  /// spatial y coordinate of space point
  float getY() const { return mY; }

  /// spatial z coordinate of space point
  float getZ() const { return mZ; }

  /// detector number corresponding to space point
  int getDetector() const { return mDetector; }

  /// pad row within detector of space point
  int getPadRow() const { return int(floor(mPadrow)); }
  float getPadRowF() const { return mPadrow; }

  /// pad position (a.k.a. column) within pad row
  float getPadCol() const { return mPadcol; }

  /// time coordinate in drift direction
  float getTimeBin() const { return mTimebin; }

  /// calculate the channel number within the MCM. 0..21 if valid, -1 if not within this MCM
  float getMCMChannel(int mcmcol) const;

  bool isInMCM(int detector, int padrow, int mcmcol) const;

  void setCollisionId(int id) { mCollisionId = id; }

  /// calculate MCM corresponding to pad row/column
  // int getMCM() const { return o2::trd::HelperMethods::getMCMfromPad(mPadrow, mPadcol); }

  /// calculate readout board corresponding to pad row/column
  // int getROB() const { return o2::trd::HelperMethods::getROBfromPad(mPadrow, mPadcol); }

 protected:
  float mX, mY, mZ;
  bool mInDrift;
  int mID;
  int mDetector;
  float mPadrow, mPadcol, mTimebin;
  int mCollisionId{0};

  // static constexpr float xscale = 1.0 / (o2::trd::Geometry::cheight() + o2::trd::Geometry::cspace());
  // static constexpr float xoffset = o2::trd::Geometry::getTime0(0);
  // static constexpr float alphascale = 1.0 / o2::trd::Geometry::getAlpha();
};

std::ostream& operator<<(std::ostream& os, const ChamberSpacePoint& p);

/// A HitPoint is a ChamberSpacePoint with additional charge information, meant to to hold all information about
/// Monte-Carlo hits within a chamber.
class HitPoint : public ChamberSpacePoint
{
 public:
  HitPoint(ChamberSpacePoint position, float charge)
    : ChamberSpacePoint(position), mCharge(charge)
  {
  }

  HitPoint(){};

  float getCharge() { return mCharge; }

 private:
  float mCharge{0.0};
};

/// A track segment: a straight line connecting two points. The points are generally given in spatial coordinates,
/// which are then converted to pad row / column / timebin coordinates. These are then used to calculate the
/// position and slope of the track segment in a given pad row, which should correspond to the reconstructed
/// tracklet.
class TrackSegment
{
 public:
  TrackSegment(){}; // default constructor, will create invalid start and end points
  TrackSegment(ChamberSpacePoint start, ChamberSpacePoint end, int id)
    : mStartPoint(start), mEndPoint(end), mTrackID(id)
  {
    assert(start.getDetector() == end.getDetector());
  }

  /// check if the space point has been initialized
  bool isValid() const { return mStartPoint.isValid(); }

  /// position of track segment at timebin 0
  float getPadColAtTimeBin(float timebin = 0) const
  {
    return mStartPoint.getPadCol() - getSlope() * (mStartPoint.getTimeBin() - timebin);
  }

  float getSlope() const
  {
    return (mEndPoint.getPadCol() - mStartPoint.getPadCol()) / (mEndPoint.getTimeBin() - mStartPoint.getTimeBin());
  }

  /// detector number
  int getDetector() const { return mStartPoint.getDetector(); }
  float getPadRowF() const { return mStartPoint.getPadRowF();}
  int getPadRow() const { return mStartPoint.getPadRow();}
  int getPadCol() const { return mStartPoint.getPadCol();}
   
  int getCollisionId() const { return mCollisionId;}
  void setCollisionId(int id) { mCollisionId=id;}

  int getTriggerTime() const { return mTriggerTime;}
  void setTriggerTime(float time) { mTriggerTime=time;}
  int getTrackTime() const { return mTrackTime;}
  void setTrackTime(float time) { mTrackTime=time;}

  ChamberSpacePoint getStartPoint()const  { return mStartPoint; }
  ChamberSpacePoint getEndPoint() const { return mEndPoint; }
  void setStartPoint(const ChamberSpacePoint& startpoint)  { mStartPoint=startpoint; }
  void setEndPoint(const ChamberSpacePoint& endpoint) {mEndPoint=endpoint; }

  int getRefTPCId()const { return mRefTPCId;}
  int getRefITSId()const { return mRefITSId;}
  void setRefTPCId(int id) { mRefTPCId=id;}
  void setRefITSId(int id) { mRefITSId=id;}

  float getPt(){return mPt;}
  float getSnp(){return mSnp;}
  float getPhi(){return mPhi;}
  float getAlpha(){return mAlpha;}
  float getStartX(){return mStartX;}
  float getStartY(){return mStartY;}
  float getStartZ(){return mStartZ;}
  float getEventNo(){return mEventNo;}
  float getTimeFrame(){return mEventNo;}
  void setTimeFrame(int tf){mTimeFrame=tf;}
  void setEventNo(int eventno){mEventNo=eventno;}
  void setPt(float pt){ mPt=pt;}
  void setSnp(float snp){ mSnp=snp;}
  void setPhi(float phi){ mPhi=phi;}
  void setAlpha(float alpha){ mAlpha=alpha;}
  void setStartX(float x){ mStartX=x;}
  void setStartY(float y){ mStartY=y;}
  void setStartZ(float z){ mStartZ=z;}
  int getLayer(){ return mLayer;}
  void setLayer(int layer){mLayer=layer;}

 protected:
  ChamberSpacePoint mStartPoint, mEndPoint;
  std::array<ChamberSpacePoint,6> mDriftPoints, mAnodePoints; // [0] is the point of drift start and anode wires respectively
  int mTrackID;
  int mRefTPCId;
  int mRefITSId;
  int mCollisionId{0};
  float mTriggerTime{0.0};
  float mTrackTime{0.0};
  int mTimeFrame{0};
  int mEventNo{0};
  float mSagita; // this is findable from the TrackID
  int mDetector;
  int mDetector2; // for those instances where the track finishes a layer in a different detector to that which is starts in.
  float mPt{0.0};
  float mPhi{0.0};
  float mSnp{0.0};
  float mStartX{0.0};
  float mStartY{0.0};
  float mStartZ{0.0};
  float mAlpha{0.0};
  int mLayer; 


};

// std::ostream& operator<<(std::ostream& os, const ChamberSpacePoint& p);

/// CoordinateTransformer: translate between local spatial and pad/timebin coordinates
///
/// The intention of the CoordinateTransformer is to bundle calculations around coordinate transformations for the TRD and access to the calibration objects.
/// So far it translates local spatial (x/y/z) coordinates within a TRD chamber (detector) to a tuple of pad row, pad column and time bin.
/// In the future, it could translate e.g. track parameters to expected tracklet position and inclination.
/// At the time of writing, only constant drift velocity, Lorentz angle and T0 are considered, and must be set by hand.
class CoordinateTransformer
{
 public:
  static CoordinateTransformer* instance()
  {
    static CoordinateTransformer mCTrans;
    return &mCTrans;
  }

  /// Translate local spatial (x/y/z) coordinates within a TRD chamber (detector) to a tuple of pad row, column and time bin.
  ///
  /// Local x,y,z coordinates follow the convention for MC hist and are assumed to be corrrected for alignment.
  /// The x-coordinate points in time direction, y in pad and z in row direction.
  /// The result is an array of three floating point numbers for row, column and timebin.
  /// Time and column can directly be compared with digit or tracklet data.
  /// The pad row is returned as a floating point number that indicates also the position within the padrow.
  /// The fractional part of the pad row is not available for digits and tracklets, and only
  std::array<float, 3> Local2RCT(int det, float x, float y, float z);
  std::array<float, 3> Global2Local(int det, float x, float y, float z);

  /// Wrapper to conveniently calculate the row/column/time coordinate of a MC hit.
  std::array<float, 3> Local2RCT(o2::trd::Hit& hit)
  {
    return Local2RCT(hit.GetDetectorID(), hit.getLocalT(), hit.getLocalC(), hit.getLocalR());
  }

  o2::trd::ChamberSpacePoint MakeSpacePoint(o2::trd::Hit& hit);

  /// Legacy, less accurate method to convert local spatial to row/column/time coordinate.
  /// This method is only included for comparision, and should be removed in the future.
  std::array<float, 3> OrigLocal2RCT(int det, float x, float y, float z);

  // recalculate given a new t0, cdrift and exb
  std::array<float, 3> RecalculateRCT(int det, float x, float y, float z, float t0, float vdrift, float exb);
  std::array<float, 3> RecalculateRCT(int det, o2::trd::ChamberSpacePoint& point, float t0, float vdrift, float exb);
  std::array<float, 3> RecalculateRCT(int det,o2::trd::ChamberSpacePoint& point);
  std::array<float, 3> RecalculateRCT(int det, float x, float y, float z);

  /// Legacy, less accurate method to convert local spatial to row/column/time coordinate.
  /// This method is only included for comparision, and should be removed in the future.
  std::array<float, 3> OrigLocal2RCT(o2::trd::Hit& hit)
  {
    return OrigLocal2RCT(hit.GetDetectorID(), hit.getLocalT(), hit.getLocalC(), hit.getLocalR());
  }

  float GetVdrift() { return mVdrift; }
  void SetVdrift(float x) { mVdrift = x; }

  float GetT0() { return mT0; }
  void SetT0(float x) { mT0 = x; }

  float GetExB() { return mExB; }
  void SetExB(float x) { mExB = x; }

 protected:
  o2::trd::Geometry* mGeo;
  float mVdrift{1.5625}; ///< drift velocity in cm/us
  float mT0{4.0};        ///< time offset of start of drift region
  float mExB{0.140};     ///< tan(Lorentz angle): tan(8 deg) ~ 0.140

 private:
  CoordinateTransformer();
};

} // namespace o2::trd

#endif // ALICEO2_TRD_RAWDISPLAY_H_
