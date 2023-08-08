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

#ifndef O2_TRDTRAPCONFIGEVENT_H
#define O2_TRDTRAPCONFIGEVENT_H

#include "CommonDataFormat/InteractionRecord.h"
#include <fairlogger/Logger.h>

#include "DataFormatsTRD/Constants.h"
#include "DataFormatsTRD/RawData.h"
#include "DataFormatsTRD/MCMEvent.h"
#include "DataFormatsTRD/Digit.h"

#include <string>
#include <map>
#include <unordered_map>
#include <array>
#include <vector>
#include <bitset>
#include <gsl/span>

namespace o2::trd
{

struct fivebit {
  //             10987654321098765432109876543210
  // uint32_t:   00000000000000000000000000000000
  union {
    //             10987654321098765432109876543210
    // uint32_t:   00000000000000000000000000000000
    uint32_t word;
    struct {
      uint32_t z : 2;
      struct {
        uint64_t size : 5;
      } __attribute__((__packed__)) x[6]; // although this is 8 dont use index 0 as its part of reserved.
    } __attribute__((__packed__));
  };
};


struct sixbit {
  //             10987654321098765432109876543210
  // uint32_t:   00000000000000000000000000000000
  union {
    //             10987654321098765432109876543210
    // uint32_t:   00000000000000000000000000000000
    uint32_t word;
    struct {
      uint32_t z : 2;
      uint32_t a : 6;
      uint32_t b : 6;
      uint32_t c : 6;
      uint32_t d : 6;
      uint32_t e : 6;
    } __attribute__((__packed__));
  };
};

struct sevenbit {
  //             10987654321098765432109876543210
  // uint32_t:   00000000000000000000000000000000
  union {
    //             10987654321098765432109876543210
    // uint32_t:   00000000000000000000000000000000
    uint32_t word;
    struct {
      uint32_t z : 4;
      uint32_t a : 7;
      uint32_t b : 7;
      uint32_t c : 7;
      uint32_t d : 7;
    } __attribute__((__packed__));
  };
};

struct tenbit {
  //             10987654321098765432109876543210
  // uint32_t:   00000000000000000000000000000000
  union {
    //             10987654321098765432109876543210
    // uint32_t:   00000000000000000000000000000000
    uint32_t word;
    struct {
      uint32_t z : 2;
      uint32_t a : 10;
      uint32_t b : 10;
      uint32_t c : 10;
    } __attribute__((__packed__));
  };
};

struct fifteenbit {
  //             10987654321098765432109876543210
  // uint32_t:   00000000000000000000000000000000
  union {
    //             10987654321098765432109876543210
    // uint32_t:   00000000000000000000000000000000
    uint32_t word;
    struct {
      uint32_t z : 2;
      uint32_t a : 15;
      uint32_t b : 15;
    } __attribute__((__packed__));
  };
};
struct sixteenbit {
  //             10987654321098765432109876543210
  // uint32_t:   00000000000000000000000000000000
  union {
    //             10987654321098765432109876543210
    // uint32_t:   00000000000000000000000000000000
    uint32_t word;
    struct {
      uint32_t a : 16;
      uint32_t b : 16;
    } __attribute__((__packed__));
  };
};
struct thirtyonebit {
  //             10987654321098765432109876543210
  // uint32_t:   00000000000000000000000000000000
  union {
    //             10987654321098765432109876543210
    // uint32_t:   00000000000000000000000000000000
    uint32_t word;
    struct {
      uint32_t z : 1;
      uint32_t a : 31;
    } __attribute__((__packed__));
  };
};

struct mcmtrapevent {
  fivebit tpl[22];
  sixbit ng[8];
  uint32_t cpu[12];
//  uint32_t getRegister(const uint regidx) { return (int) tpl[0].x[1]; };
};

class TrapConfigEvent
{
  // class that is actually stored in ccdb.
  // it holds a compressed version of what comes in the config events

 public:
  TrapConfigEvent();
  ~TrapConfigEvent() = default;

  TrapConfigEvent(const TrapConfigEvent& A);

  // get a config register value by index, addr, and name, via mcm index
  uint32_t getRegisterValue(uint32_t regidx, int mcmidx);
  bool setRegisterValue(uint32_t data, uint32_t regidx, int mcmidx);
  uint32_t getRegisterValueByIdx(uint32_t regidx, int mcmidx);
  bool setRegisterValueByIdx(uint32_t data, uint32_t regidx, int mcmidx);
  uint32_t getRegisterValueByAddr(uint32_t addr, int mcmidx);
  bool setRegisterValueByAddr(uint32_t data, uint32_t addr, int mcmidx);
  uint32_t getRegisterValueByName(const std::string& name, int mcmidx);
  bool setRegisterValueByName(uint32_t data, const std::string& regname, int mcmidx);

  // get a registers name by addres and index;
  std::string getRegNameByAddr(uint16_t addr);
  std::string getRegNameByIdx(unsigned int regidx);

  // get a registers index (enum value) by address or name
  int32_t getRegIndexByAddr(unsigned int addr);
  int32_t getRegIndexByName(const std::string& name);

  // get a registers address by index or name
  int32_t getRegAddrByIdx(unsigned int regidx);
  int32_t getRegAddrByName(const std::string& name);

  bool isValidAddress(uint32_t addr);
  const std::string getRegisterName(unsigned int index) { return mTrapRegisters[index].getName(); }
  const uint32_t getRegisterMax(unsigned int index) { return mTrapRegisters[index].getMax(); }
  const uint32_t getRegisterNBits(unsigned int index) { return mTrapRegisters[index].getNbits(); }

  // get the identifying characteristics of a register given its address.
  void getRegisterByAddr(uint32_t registeraddr, std::string& regname, int32_t& newregidx, int32_t& numberbits);

  // return all the registers for a particular mcm
  void getAllRegisters(int mcmidx, std::array<uint32_t, kLastReg>& registers);

  // return all the mcm values for a particular register index
  void getAllMCMByIndex(int regindex, std::array<uint32_t, o2::trd::constants::MAXMCMCOUNT>& registers);

  // return all the mcm values for a particular register name
  void getAllMCMByName(std::string registername, std::array<uint32_t, o2::trd::constants::MAXMCMCOUNT>& mcms);

  // return all the mcm for a particular register address
  void getAllMCMByAddress(int registeraddress, std::array<uint32_t, o2::trd::constants::MAXMCMCOUNT>& mcms);

  // return all the config data in an unpacked array.
  void getAll(std::array<uint32_t, kLastReg * o2::trd::constants::MAXMCMCOUNT>& configdata);

  // return the full unpacked config event.
  bool printRegister(int regindex, int det, int rob, int mcm);
  // bool printRegister(TrapRegInfo* reg, int det, int rob, int mcm);

  const std::array<o2::trd::TrapRegInfo, kLastReg>& getTrapRegisters() { return mTrapRegisters; }

  uint32_t getRegisterMax(int idx) { return mTrapRegisters[idx].getMax(); }

  // population pending
  uint32_t getConfigVersion(const int mcmid=0) { return getRegisterValue(kC09CPU1,mcmid); }              // these must some how be gotten from git or wingdb.
  uint32_t getConfigName(const int mcmid=0) { return  getRegisterValue(kC09CPU1,mcmid); }                  // these must be gotten from git or wingdb.
  uint16_t getConfigSavedVersion(const int mcmid=0) { return getRegisterValue(kC09CPU1,mcmid); }    // the version that is saved, for the ability to later save the config differently.

  bool isConfigDifferent(const TrapConfigEvent& trapconfigevent)const;

  // for compliance with the same interface to o2::trd::TrapConfigEvent and its run1/2 inherited interface:
  // only these 2 are used in the simulations.
  uint32_t getDmemUnsigned(uint32_t address, int detector, int rob, int mcm);
  uint32_t getTrapReg(uint32_t index, int detector, int rob, int mcm);

  bool isHCIDPresent(int hcid)const  { return mHCIDPresent.test(hcid); }
  void HCIDIsPresent(int hcid) { mHCIDPresent.set(hcid); }
  uint32_t countHCIDPresent() { return mHCIDPresent.count(); }
  const std::bitset<constants::MAXHALFCHAMBER>& getHCIDPresent() const { return mHCIDPresent; }
  bool isMCMPresent(int mcmid) const { return mMCMPresent.test(mcmid); }
  void MCMIsPresent(int mcmid) { mMCMPresent.set(mcmid); }
  const std::bitset<constants::MAXMCMCOUNT>& getMCMPresent()const  { return mMCMPresent; }
  uint32_t countMCMPresent() { return mMCMPresent.count(); }
  void clearMCMPresent() { mMCMPresent.reset(); }
  void clearMCMEvent() { mConfigData.clear(); }
  void clear()
  {
    clearMCMPresent();
    clearMCMEvent();
  }
  bool ignoreWord(int offset) const { return mWordNumberIgnore.test(offset); }

  // required for a container for calibration
  void fill(const TrapConfigEvent& input);
  void fill(const gsl::span<const TrapConfigEvent> input); // dummy!
  void merge(const TrapConfigEvent* prev);
  void print();
  int getRegisterBase(int reg) { return mTrapRegisters[reg].getBase(); }
  const MCMEvent& getMCMEvent(int idx) { return mConfigData[idx]; }
  int getMCMEventSize() { return mConfigDataIndex.size(); }

  const uint32_t getRegBase(int idx){ return mTrapRegisters[idx].getBase();}
  const uint32_t getRegWordNumber(int idx){ return mTrapRegisters[idx].getWordNumber();}
  const uint32_t getRegMask(int idx){ return mTrapRegisters[idx].getMask();}
  const uint32_t getRegShift(int idx){ return mTrapRegisters[idx].getShift();}


 private:
  std::array<TrapRegInfo, kLastReg> mTrapRegisters;                 //!< store of layout of each block of mTrapRegisterSize, populated via initialiseRegisters
  std::bitset<constants::MAXMCMCOUNT> mMCMPresent{0};               //!< does the mcm actually receive data.
  std::bitset<constants::MAXHALFCHAMBER> mHCIDPresent{0};           //!< did the link actually receive data.
  std::vector<MCMEvent> mConfigData;                                //!< vector of register data blocks
  std::array<int32_t, constants::MAXMCMCOUNT> mConfigDataIndex{-1}; //!< one block of data per mcm, array as one wants to query if an mcm is present with having to walk the whole index.
  std::map<uint16_t, uint16_t> mTrapRegistersAddressIndexMap;       //!< map of address into mTrapRegisters, populated at the end of initialiseRegisters
  std::bitset<kTrapRegistersSize> mWordNumberIgnore;                //!< whether to ignore a register or not. Here to speed lookups up.
  void initialiseRegisters();

  //uint32_t mTrapConfigEventNumber;       //!< the version number coming from the config that is written to the traps, from config register C09CPU01
  //uint32_t mTrapConfigEventVersion;      //!< the version of the coming from the config that is written to the traps, from DigitHCHeader (svn version number)
  //uint16_t mTrapConfigEventSavedVersion; //!< the version that is saved, for the ability to later save the config differently.
  ClassDefNV(TrapConfigEvent, 1);
};

} // namespace o2::trd

#endif
