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

// Cru raw data reader, this is the part that parses the raw data
// it runs on the flp(pre compression) or on the epn(pre tracklet64 array generation)
// it hands off blocks of cru pay load to the parsers.

#ifndef O2_TRD_RAWDATASTATS
#define O2_TRD_RAWDATASTATS

#include <iostream>
#include <string>
#include <cstdint>
#include <array>
#include <vector>
#include <gsl/span>
#include "DataFormatsTRD/Constants.h"

namespace o2::trd
{
enum ParsingErrors { TRDParsingNoError,
                     TRDParsingUnrecognisedFormat,
                     TRDParsingBadDigt,
                     TRDParsingBadTracklet,
                     TRDParsing };
class TRDDataCountersPerEvent
{ //thisis on a per event basis
 public:
  //TODO this should go into a dpl message for catching by qc ?? I think.
  std::array<uint32_t, 1080> mLinkWords;       //units of 256bits, read from the cru half chamber header
  std::array<uint32_t, 1080> mLinkWordsRead;   // units of 32 bits the data words read before dumping or finishing
  std::array<uint32_t, 1080> mLinkWordsDumped; // units of 32 bits the data dumped due to some or other error
  std::array<uint8_t, 1080> mLinkErrorFlag;    //status of the error flags for this event
  //from the above you can get the stats for supermodule and detector.
  std::array<bool, 1080> mLinkEmpty; // Link only has padding words only, probably not serious.
  //maybe change this to actual traps ?? but it will get large.
  std::array<int64_t, o2::trd::constants::MAXMCMCOUNT> mLinkMCMsWithData; // and its corresponding volume of data.
  std::array<std::array<uint16_t, 16>, 1080> mMCMsBeforeCorruption;       // count the mcms read before link was corrupt
  std::array<std::array<uint16_t, 16>, 1080> mMCMsOnlyClean;
  std::array<uint32_t, 1080> mLinkMCMCountBeforeCorruption;
  std::array<uint32_t, 1080> mLinkDigitCount;
  std::array<uint32_t, 1080> mLinkTrackletCount;
  std::array<uint32_t, constants::MAXMCMCOUNT> mTrapTrackletCount;
  std::array<uint16_t, constants::MAXMCMCOUNT> mMCMstats; // bit pattern for errors current event for a given mcm;
  std::vector<uint32_t> mEmptyTraps;                      // MCM indexes of traps that are empty ?? list might better
  double mTimeTaken;                                      // time take to process an event (summed trackletparsing and digitparsing) parts not accounted for.
  double mDataWordsRead;                                  // data words read in
  double mDataWordsSkipped;                               // data words skipped for various reasons.
};

class TRDDataCountersPerTimeFrame
{ //thisis on a per event basis
 public:
  double mTimeTaken;
  uint64_t mDigitsFound;
  uint64_t mTrackletsFound;
  uint64_t mDataWordsRead;
  uint64_t mDataWordsRejected;
};

//TODO not sure this class is needed
class TRDDataCountersRunning
{                                       //those counters that keep counting
  std::array<uint32_t, 1080> mLinkFreq; //units of 256bits "cru word"
  std::array<bool, 1080> mLinkEmpty;    // Link only has padding words only, probably not serious.
  std::array<uint64_t, 65535> mDataFormatRead; // 7bits.7bits major.minor version read from HCHeader.
};

} // namespace o2::trd

#endif
