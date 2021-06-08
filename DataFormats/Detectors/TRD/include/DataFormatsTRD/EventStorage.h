// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// @file   CruRawReader.h
/// @author Sean Murray
/// @brief  Cru raw data reader, this is the part that parses the raw data
//          it runs on the flp(pre compression) or on the epn(pre tracklet64 array generation)
//          it hands off blocks of cru pay load to the parsers.

#ifndef O2_TRD_EVENTSTORAGE
#define O2_TRD_EVENTSTORAGE

#include <fstream>
#include <iostream>
#include <string>
#include <cstdint>
#include <array>
#include <vector>
#include <iterator>
#include "Headers/RAWDataHeader.h"
#include "Headers/RDHAny.h"
#include "CommonDataFormat/RangeReference.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsTRD/Constants.h"
#include "fairlogger/Logger.h"
//##include "DataFormatsTRD/CompressedDigit.h"

namespace o2::trd
{

//store indxes into the returned data structures from reading. Saves on a multitude of copies.

class EventStorage1
{

 using DataRange = o2::dataformats::RangeReference<int>;
 public:
  EventStorage1() = default;
  ~EventStorage1() = default;

  int addDigits(InteractionRecord& ir, int start, int count ){
    int index=0;for(auto irindex: mInteractionRecords){
                                                              if(irindex==ir){
                                                                mDigits[index].push_back(DataRange(start,count));
                                                              } index++; }
                                                            //no corresponding interaction record found.
                                                            mInteractionRecords.push_back(ir);
                                                            mDigits[mInteractionRecords.size()].push_back(DataRange(1,1));
                                                            return 1;
  }
  int addCompressedDigits(InteractionRecord& ir, int start, int count ){
    int index=0;for(auto irindex: mInteractionRecords){
                                                              if(irindex==ir){
                                                                mDigits[index].push_back(DataRange(start,count));
                                                              } index++; }
                                                            //no corresponding interaction record found.
                                                            mInteractionRecords.push_back(ir);
                                                            mCompressedDigits[mInteractionRecords.size()].push_back(DataRange(1,1));
                                                            return 1;
  }
  int addTracklets(InteractionRecord& ir, int start, int count ){
    int index=0;for(auto irindex: mInteractionRecords){
                                                              if(irindex==ir){
                                                                mDigits[index].push_back(DataRange(start,count));
                                                              } index++; }
                                                            //no corresponding interaction record found.
                                                            mInteractionRecords.push_back(ir);
                                                            mDigits[mInteractionRecords.size()].push_back(DataRange(1,1));
                                                            return 1;
  }
  //interaction records get added on no matching ones when adding the 3 above.
  //
  std::vector<DataRange>& getCompressedDigits(InteractionRecord& ir) {int count=0; for(auto intrec:mInteractionRecords) if(ir==intrec){return mCompressedDigits[count];}} ;
  std::vector<DataRange>& getDigits(InteractionRecord& ir) {return mDigits[0];};
  std::vector<InteractionRecord>& getItectionaRecords(InteractionRecord& ir) {return mInteractionRecords;};
  std::vector<DataRange>& getTracklets(InteractionRecord& ir) {if(mTracklets.size()==0) LOG(fatal) << "empty data range for tracklets";return mTracklets[0];};


  void clear()
  {
   for(auto trackletv : mTracklets){
     trackletv.clear();
   }
   for(auto digitv : mDigits){
     digitv.clear();
   }
   for(auto cdigitv : mCompressedDigits){
     cdigitv.clear();
   }
   for(auto irv : mInteractionRecords){
     irv.clear();
   }
  }

 protected:

  std::vector<std::vector<DataRange>> mTracklets; // when this runs properly it will only 6 for the flp its runnung on.
  std::vector<InteractionRecord> mInteractionRecords;
  std::vector<std::vector<DataRange>> mCompressedDigits;
  std::vector<std::vector<DataRange>> mDigits;
  /** summary data **/
};

} // namespace o2::trd

#endif
