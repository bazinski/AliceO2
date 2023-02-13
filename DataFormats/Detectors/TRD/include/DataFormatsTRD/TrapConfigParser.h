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

#ifndef O2_TRDTRAPCONFIGPARSER_H
#define O2_TRDTRAPCONFIGPARSER_H

#include "Rtypes.h"
#include "TH2F.h"
#include <fairlogger/Logger.h>
#include "DataFormatsTRD/Constants.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsTRD/TrapConfig3.h"
#include <chrono>  // chrono::system_clock
#include <ctime>   // localtime
#include <sstream> // stringstream
#include <iomanip> // put_time
#include <string>  // string
#include <vector>
#include <memory>
#include <map>
#include <tuple>
#include <bitset>

// Configuration of the TRD Tracklet Processor
// (TRD Front-End Electronics)
// There is a manual to describe all the internals of the electronics.
// A very detailed manual, will try put it in the source code.
// TRAP registers
// TRAP data memory (DMEM)
//

namespace o2::trd
{

/*struct regmap_addr {
  std::string name;
  int index;
  int bitwidth;
};

struct regmap {
  std::string name;
  int address;
  int bitwidth;
};
*/

enum TrapConfigErrorTypes {
  kTrapConfigAllGood,
  kTrapConfigNoEnd,
  kTrapConfigNoBegin,
  kTrapConfigNoGarbageEnd,
  kTrapConfigNoGarbageBegin,
  kTrapConfigLastError
}; // possible errors for parsing.

class TrapConfigParser
{
  // class used to unpack the trap config events, and figure out what to do to them.
 public:
  TrapConfigParser();
  ~TrapConfigParser();

  int getTrapReg(TrapReg reg, int det = -1, int rob = -1, int mcm = -1);

  unsigned int getDmemUnsigned(int addr, int det, int rob, int mcm);

  // helper methods
  std::string getConfigVersion() { return mTrapConfigVersion; }
  std::string getConfigName() { return mTrapConfigName; }
  void setConfigVersion(std::string version) { mTrapConfigVersion = version; } // these must some how be gotten from git or wingdb.
  void setConfigName(std::string name) { mTrapConfigName = name; }             // these must be gotten from git or wingdb.

  TrapReg getRegByAddress(int address);

  void printMCMRegisterCount(int hcid);
  void unpackBlockHeader(uint32_t& header, uint32_t& registerdata, uint32_t& step, uint32_t& bwidth, uint32_t& nwords, uint16_t& registeraddr, uint32_t& exit_flag);
  int parseLink(std::array<uint32_t, o2::trd::constants::HBFBUFFERMAX>& data, uint32_t start, uint32_t end, int hcid);

  bool checkRegister(uint16_t& registeraddr, int& registerbase, int& registeroffset, int& lastregisterindex,
                     int& currentregisterindex, int& registerscount, int& registererrorgap, uint32_t& lastregisterread,
                     int& registerwordscount, int& lastregidx, uint32_t& registerdata, int mcmid, int hcid, int mcm, int rob);

  void clear()
  {
    mMcmParsingStatus.fill(0);
    mHCIDhasConfig.reset();
  };                                 // clear out the mMcmParsingStatus, and hcid bitset for config event presence;
  void FillHistograms(int eventnum); //;TH2F *hists[6])
  void compareToTrackletsHCID(std::bitset<1080> trackletshcid);
  std::array<int, kLastReg>& getStartRegArray() { return mStartReg; }      // the number of time this register was read as the first register
  std::array<int, kLastReg>& getStopRegArray() { return mStopReg; }        // the number of time this register was read as the last register
  std::array<int, kLastReg>& getMissedRegArray() { return mMissedReg; }    // the number of times this register was not read
  std::array<int, kLastReg>& getRegisterCount() { return mRegisterCount; } // total count for each register

  // write the current config to a file
  void writeFile(int eventcount);

 private:
  //  std::map<int, std::tuple<std::string, int, int>> TrapRegisterMap_addr;
  //std::array<int, 0xe000> mTrapRegistersAddressIndex; // index by address into mTrapRegisters.
  std::bitset<o2::trd::constants::MAXMCMCOUNT> mMCMDataHasChanged;
  std::array<int, o2::trd::constants::MAXMCMCOUNT> mMcmParsingStatus{0};                        // status of what was found, errors types in the parsing
  std::array<uint32_t, o2::trd::constants::MAXMCMCOUNT * kLastReg> mCurrentMCMRegisters;        // store the registers for all the mcm registers.
  std::array<uint32_t, o2::trd::constants::MAXHALFCHAMBER> mHalfChamberLastSeen;                // timestamp
  std::array<uint32_t, o2::trd::constants::MAXHALFCHAMBER> mHalfChamberSeenSinceLastWritten;    // timestamp
  std::array<uint32_t, o2::trd::constants::MAXHALFCHAMBER> mHalfChamberFrequencyInAccumulation; // frequency in accumulation period.
  std::array<uint32_t, o2::trd::constants::MAXMCMCOUNT> mMCMLastSeen;                           // timestamp
  std::array<uint32_t, o2::trd::constants::MAXMCMCOUNT> mMCMSeenSinceLastWritten;               // timestamp
  std::array<uint32_t, o2::trd::constants::MAXMCMCOUNT> mMCMFrequencyInAccumulation;            // frequency in accumulation.
  std::array<int, 8 * 16> mcmSeen;                                                              // the mcm has been seen with or with out error, local to a link
  std::array<int, 8 * 16> mcmMCM;                                                               // the mcm has been seen with or with out error, local to a link
  std::array<int, 8 * 16> mcmROB;                                                               // the mcm has been seen with or with out error, local to a link
  std::array<int, 8 * 16> mcmSeenMissedRegister;                                                // the mcm does not have a complete set of registers, local to a link
  std::array<std::bitset<kLastReg>, 8 * 16> mcmMissedRegister;                                  // bitpattern of which registers were seen and not seen for a given mcm.
                                                                                                //  static bool mRegisterAddressMapInitialised;
  InteractionRecord mIR;
  std::time_t mConfigDate;
  std::bitset<1080> mHCIDhasConfig;            // to monitor when a link goes down, each hcid has a link
  std::array<int, kLastReg> mStartReg{0};      // count which register a config starts at
  std::array<int, kLastReg> mStopReg{0};       // count which register a config stops at.
  std::array<int, kLastReg> mMissedReg{0};     // count the number of missed reigsters
  std::array<int, kLastReg> mRegisterCount{0}; // register frequency, a count of how many times each register appears
  std::string mTrapConfigName;                 //TOOD figure how to pull this in or seperately put it in the CCDB
  std::string mTrapConfigVersion;
  TrapConfig3 mTrapConfig;
  std::array<std::map<uint32_t, uint32_t>, kLastReg> mTrapRegistersFrequencyMap;
  // TrapConfig3& operator=(const TrapConfig3& rhs); // not implemented
  // TrapConfig3(const TrapConfig3& cfg);            // not implemented

  ClassDefNV(TrapConfigParser, 1);
};

} //namespace o2::trd
#endif
