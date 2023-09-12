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
#include "DataFormatsTRD/TrapRegisters.h"
#include "DataFormatsTRD/TrapRegInfo.h"

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
    uint32_t word;
    struct {
      struct {
        uint8_t five : 5;
      } __attribute__((__packed__)) _five[6];
      uint32_t spare : 2;
    } __attribute__((__packed__));
  };
};

struct sixbit {
  //             10987654321098765432109876543210
  // uint32_t:   00000000000000000000000000000000
  union {
    uint32_t word;
    struct {
      struct {
        uint8_t six : 5;
      } __attribute__((__packed__)) _six[5];
      uint32_t spare : 2;
    } __attribute__((__packed__));
  };
};

struct sevenbit {
  //             10987654321098765432109876543210
  // uint32_t:   00000000000000000000000000000000
  union {
    uint32_t word;
    struct {
      struct {
        uint8_t seven : 7;
      } __attribute__((__packed__)) _seven[4];
      uint32_t spare : 4;
    } __attribute__((__packed__));
  };
};

struct tenbit {
  //             10987654321098765432109876543210
  // uint32_t:   00000000000000000000000000000000
  union {
    uint32_t word;
    struct {
      struct {
        uint8_t ten : 7;
      } __attribute__((__packed__)) _ten[3];
      uint32_t spare : 2;
    } __attribute__((__packed__));
  };
};

struct fifteenbit {
  //             10987654321098765432109876543210
  // uint32_t:   00000000000000000000000000000000
  union {
    uint32_t word;
    struct {
      struct {
        uint16_t fifteen : 15;
      } __attribute__((__packed__)) _fifteen[2];
      uint32_t spare : 2;
    } __attribute__((__packed__));
  };
};
struct sixteenbit {
  //             10987654321098765432109876543210
  // uint32_t:   00000000000000000000000000000000
  union {
    uint32_t word;
    struct {
      struct {
        uint16_t sixteen : 16;
      } __attribute__((__packed__)) _sixteen[2];
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
      uint32_t a : 31;
      uint32_t spare : 1;
    } __attribute__((__packed__));
  };
};

struct mcmtrapevent {
  short tpl[128];
  tenbit fgf[5];
  fivebit cpuclk[2];
  fifteenbit ebi[2];
  sixbit tp[2];
  fivebit ebin;
  sixbit fll[13];
  sevenbit tpp[4];
  tenbit fptc[2];
  fifteenbit fgta[2];
  tenbit ft[3];
  uint32_t adcmsk;
  fivebit adc;
  uint32_t adcpar;
  tenbit pasa[3];
  fivebit sadc;
  fifteenbit sml[2];
  sixteenbit smmode; // arbtim as well
  fifteenbit ni[2];
  fifteenbit irq[32];
  uint32_t ctgdini;
  sixteenbit ctgctrl; // half empty
  tenbit memrw[2];
  thirtyonebit nmod[5];
  sixteenbit nbnd; // half empty
  fifteenbit np[4];
  uint32_t cpu[19];
  thirtyonebit nes;
  uint32_t ncut;
  uint32_t pasachm;
  // still left with a mapping problem as in TrapRegisters, albeit different
  // map ktpl00 to mcmtrapevent.tpl[0]._five[0].five;
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
  uint32_t getRegisterValue(const uint32_t regidx, const int mcmidx);
  bool setRegisterValue(const uint32_t data, const uint32_t regidx, const int mcmidx);
  //  uint32_t getRegisterValueByIdx(const uint32_t regidx,const  int mcmidx);
  //  bool setRegisterValueByIdx(const uint32_t data, const uint32_t regidx, const int mcmidx);
  //  uint32_t getRegisterValueByAddr(const uint32_t addr, const int mcmidx);
  //  bool setRegisterValueByAddr(const uint32_t data, const uint32_t addr, const int mcmidx);
  //  uint32_t getRegisterValueByName(const std::string& name, const int mcmidx);
  //  bool setRegisterValueByName(const uint32_t data, const std::string& regname, const int mcmidx);

  // get a registers name by addres and index;
  std::string getRegNameByAddr(const uint16_t addr);
  const std::string getRegNameByIdx(const uint32_t regidx) { return mTrapRegisters[regidx].getName(); }

  // get a registers index (enum value) by address or name
  int32_t getRegIndexByAddr(const unsigned int addr);
  int32_t getRegIndexByName(const std::string& name);

  // get a registers address by index or name
  int32_t getRegAddrByIdx(const unsigned int regidx);
  int32_t getRegAddrByName(const std::string& name);

  bool isValidAddress(const uint32_t addr);
  const std::string getRegisterName(const unsigned int regidx) { return mTrapRegisters[regidx].getName(); }
  const uint32_t getRegisterMax(const unsigned int regidx) { return mTrapRegisters[regidx].getMax(); }
  const uint32_t getRegisterNBits(const unsigned int regidx) { return mTrapRegisters[regidx].getNbits(); }
  const TrapRegInfo& getRegisterInfo(const uint32_t regidx) { return mTrapRegisters[regidx]; }

  // get the identifying characteristics of a register given its address.
  void getRegisterByAddr(const uint32_t registeraddr, std::string& regname, int32_t& newregidx, int32_t& numberbits);

  // return all the registers for a particular mcm
  void getAllRegisters(const int mcmidx, std::array<uint32_t, o2::trd::TrapRegisters::kLastReg>& registers);

  // return all the mcm values for a particular register index
  void getAllMCMByIndex(const int regindex, std::array<uint32_t, o2::trd::constants::MAXMCMCOUNT>& registers);

  // return all the mcm values for a particular register name
  void getAllMCMByName(const std::string& registername, std::array<uint32_t, o2::trd::constants::MAXMCMCOUNT>& mcms);

  // return all the mcm for a particular register address
  void getAllMCMByAddress(const int registeraddress, std::array<uint32_t, o2::trd::constants::MAXMCMCOUNT>& mcms);

  // return all the config data in an unpacked array.
  void getAll(std::array<uint32_t, o2::trd::TrapRegisters::kLastReg * o2::trd::constants::MAXMCMCOUNT>& configdata);

  // population pending
  uint32_t getConfigVersion(const int mcmid = 0) { return getRegisterValue(TrapRegisters::kC09CPU1, mcmid); }      // these must some how be gotten from git or wingdb.
  uint32_t getConfigName(const int mcmid = 0) { return getRegisterValue(TrapRegisters::kC09CPU1, mcmid); }         // these must be gotten from git or wingdb.
  uint16_t getConfigSavedVersion(const int mcmid = 0) { return getRegisterValue(TrapRegisters::kC09CPU2, mcmid); } // the version that is saved, for the ability to later save the config differently.

  bool isConfigDifferent(const TrapConfigEvent& trapconfigevent) const;

  // for compliance with the same interface to o2::trd::TrapConfigEvent and its run1/2 inherited interface:
  // only these 2 are used in the simulations.
  uint32_t getDmemUnsigned(uint32_t address, int detector, int rob, int mcm);
  uint32_t getTrapReg(const uint32_t index, const int detector, const int rob, const int mcm);

  /*  bool isHCIDPresent(const int hcid) const { return mHCIDPresent.test(hcid); }
    void HCIDIsPresent(const int hcid) { mHCIDPresent.set(hcid); }
    uint32_t countHCIDPresent() { return mHCIDPresent.count(); }
    const std::bitset<constants::MAXHALFCHAMBER>& getHCIDPresent() const { return mHCIDPresent; }
    bool isMCMPresent(const int mcmid) const { return mMCMPresent.test(mcmid); }
    void MCMIsPresent(const int mcmid) { mMCMPresent.set(mcmid); }
    const std::bitset<constants::MAXMCMCOUNT>& getMCMPresent() const { return mMCMPresent; }
    uint32_t countMCMPresent() { return mMCMPresent.count(); }
    void clearMCMPresent() { mMCMPresent.reset(); }
    void clearMCMEvent() { mConfigData.clear(); }*/
  void clear()
  {
    //   clearMCMPresent();
    //  clearMCMEvent();
  }
  // bool ignoreWord(const int offset) const { return mWordNumberIgnore.test(offset); }

  // required for a container for calibration
  void fill(const TrapConfigEvent& input);
  void fill(const gsl::span<const TrapConfigEvent> input); // dummy!
  void merge(const TrapConfigEvent* prev);
  void print();
  int getRegisterBase(const int regidx) { return mTrapRegisters[regidx].getBase(); }
  const MCMEvent& getMCMEvent(const int mcmidx) { return mConfigData[mcmidx]; }
  int getMCMEventSize() { return mConfigDataIndex.size(); }
  // TODO put into QC

 private:
  TrapRegisters mTrapRegisters;
  // std::bitset<constants::MAXMCMCOUNT> mMCMPresent{0};               ///< does the mcm actually receive data.
  // std::bitset<constants::MAXHALFCHAMBER> mHCIDPresent{0};           ///< did the link actually receive data.
  std::vector<MCMEvent> mConfigData;                                ///< vector of register data blocks
  std::array<int32_t, constants::MAXMCMCOUNT> mConfigDataIndex{-1}; ///< one block of data per mcm, array as one wants to query if an mcm is present with having to walk the whole index.
  std::map<uint16_t, uint16_t> mTrapRegistersAddressIndexMap;       //!< map of address into mTrapRegisters, populated at the end of initialiseRegisters
  // std::bitset<kTrapRegistersSize> mWordNumberIgnore;                ///< whether to ignore a register or not. Here to speed lookups up.
  void initialiseRegisters(); // build some indexes
  ClassDefNV(TrapConfigEvent, 1);
};

} // namespace o2::trd

#endif
