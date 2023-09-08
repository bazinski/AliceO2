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

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  TRAP config as received from the trap chips in the form of            //
//  configuration events, configs are packed into digit events.           //
//  Header major version defines the config event payload follows as a    //
//  known block.                                                          //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "DataFormatsTRD/TrapConfigEvent.h"
#include "DataFormatsTRD/TrapRegInfo.h"
#include "DataFormatsTRD/TrapRegisters.h"
#include "DataFormatsTRD/Constants.h"
#include "DataFormatsTRD/HelperMethods.h"
#include "DataFormatsTRD/RawData.h"

#include <fairlogger/Logger.h>

#include <array>
#include <map>

using namespace o2::trd;

TrapConfigEvent::TrapConfigEvent()
{
  mConfigDataIndex.fill(-1); //!< one block of data per mcm.
  initialiseRegisters();
  //  setConfigSavedVersion(2);
}

TrapConfigEvent::TrapConfigEvent(const TrapConfigEvent& A)
{
  LOGP(info, "TrapConfigEvent Copy Constructor Called   .....");
  mConfigDataIndex.fill(-1); //!< one block of data per mcm.
}

void TrapConfigEvent::initialiseRegisters()
{
  // mTrapRegisters[kSMCMD].init("SMCMD", 0x0A04, 16, 0x0000);
  LOGP(info, " finished initing registers now to build the map : ");
  for (int reg = 0; reg < TrapRegisters::kLastReg; ++reg) {
    // reindex to speed things up, this time by address, map instead of a rather large lookup table.
    auto addr = mTrapRegisters[reg].getAddr();
    mTrapRegistersAddressIndexMap[addr] = reg;
    // build index of wordnumbers, to speed up comparisons.
    // this indexes ::  [raw block offset:ignorechange]
    auto wordnumber = mTrapRegisters[reg].getWordNumber();
    mWordNumberIgnore.set(wordnumber) = mTrapRegisters[reg].getIgnoreChange();
  }
  // build index of wordnumbers, to speed up comparisons.
}

uint32_t TrapConfigEvent::getRegisterValue(const uint32_t regidx, const int mcmidx)
{
  // get the register value based on the address of the register;
  // find register in mTrapRegisters.
  // calculate the offset from the base for the register and get the mask.
  if ((regidx < 0) || (regidx >= TrapRegisters::kLastReg) || (mcmidx < 0 || mcmidx >= o2::trd::constants::MAXMCMCOUNT)) {
    return 0; // TODO this could be a problem ?!
  }
  int mcmoffset = 0;                                                    // mcmidx * kTrapRegistersSize;                                 // get the start offset for this mcm
  int regbase = mcmoffset + mTrapRegisters[regidx].getBase();           // get the base of this register in the underlying storage block
  int regoffset = regbase + mTrapRegisters[regidx].getDataWordNumber(); // get the offset to the register in question

  uint32_t data = mConfigData[mConfigDataIndex[mcmidx]].getRegister(regidx, mTrapRegisters[regidx]);
  data &= mTrapRegisters[regidx].getMask(); // mask the data off as need be.
  return data;
}

bool TrapConfigEvent::setRegisterValue(uint32_t data, uint32_t regidx, int mcmidx)
{
  if (regidx < 0 || regidx >= TrapRegisters::kLastReg || mcmidx < 0 || mcmidx >= o2::trd::constants::MAXMCMCOUNT) {
    return false;
  }
  int mcmoffset = 0; // mcmidx * kTrapRegistersSize;                          // get the start offset for this mcm
  if (mConfigDataIndex[mcmidx] == -1) {
    // we dont have this mcm in the data store yet.
    mConfigData.emplace_back(mcmidx);
    mConfigDataIndex[mcmidx] = mConfigData.size() - 1;
  }
  uint32_t index = mConfigDataIndex[mcmidx];
  return mConfigData[mConfigDataIndex[mcmidx]].setRegister(data, regidx, mTrapRegisters[regidx]);
  // return mConfigData[mConfigDataIndex[mcmidx]].setRegister(data, regidx,  mTrapRegisters[regidx].getBase(),mTrapRegisters[regidx].getWordNumber(),mTrapRegisters[regidx].getShift(),mTrapRegisters[regidx].getMask());
}

// this is here to be able to use the on the fly index.
std::string TrapConfigEvent::getRegNameByAddr(uint16_t addr)
{
  std::string name = "";
  if (auto search = mTrapRegistersAddressIndexMap.find(addr); search != mTrapRegistersAddressIndexMap.end()) {
    name = mTrapRegisters[mTrapRegistersAddressIndexMap[addr]].getName();
  }
  return name;
}

// this is here to be able to use the on the fly index.
int32_t TrapConfigEvent::getRegIndexByAddr(unsigned int addr)
{
  if (isValidAddress(addr)) {
    return mTrapRegistersAddressIndexMap[addr];
  } else
    return -1;
}

// this is here to be able to use the on the fly index.
bool TrapConfigEvent::isValidAddress(uint32_t addr)
{
  auto search = mTrapRegistersAddressIndexMap.find(addr);
  return (search != mTrapRegistersAddressIndexMap.end());
}

int32_t TrapConfigEvent::getRegIndexByName(const std::string& name)
{
  // there is no index for this but its not used online
  return mTrapRegisters.getRegIndexByName(name);
}

int32_t TrapConfigEvent::getRegAddrByIdx(unsigned int regidx)
{
  if ((regidx < 0) && (regidx > TrapRegisters::kLastReg)) {
    return -1;
  }
  return mTrapRegisters[regidx].getAddr();
}

int32_t TrapConfigEvent::getRegAddrByName(const std::string& name)
{
  // there is no index for this but its not used online
  return mTrapRegisters.getRegAddrByName(name);
}

void TrapConfigEvent::getRegisterByAddr(uint32_t registeraddr, std::string& regname, int32_t& newregidx, int32_t& numberbits)
{
  int idx = -1;
  idx = getRegIndexByAddr(registeraddr);
  if (idx >= 0) {
    regname = mTrapRegisters[idx].getName();
    numberbits = mTrapRegisters[idx].getNbits();
    newregidx = idx;
  } else {
    regname = "";
    numberbits = -1;
    newregidx = -1;
  }
}

// get all the values for a given register values for a given mcm
void TrapConfigEvent::getAllRegisters(const int mcmidx, std::array<uint32_t, TrapRegisters::kLastReg>& mcmregisters)
{
  for (int reg = 0; reg < TrapRegisters::kLastReg; ++reg) {
    mcmregisters[reg] = getRegisterValue(reg, mcmidx);
  }
}

// get all the values for a given register, all mcms
void TrapConfigEvent::getAllMCMByIndex(const int regidx, std::array<uint32_t, o2::trd::constants::MAXMCMCOUNT>& mcms)
{
  if ((regidx >= 0) && (regidx < TrapRegisters::kLastReg)) {
  }
  for (int mcm = 0; mcm < o2::trd::constants::MAXMCMCOUNT; ++mcm) {
    mcms[mcm] = getRegisterValue(regidx, mcm);
  }
}

// get all the values for a given register name, all mcms
void TrapConfigEvent::getAllMCMByName(const std::string& registername, std::array<uint32_t, o2::trd::constants::MAXMCMCOUNT>& mcms) // return all the mcm values for a particular register name
{
  int regindex = getRegIndexByName(registername);
  for (int mcm = 0; mcm < o2::trd::constants::MAXMCMCOUNT; ++mcm) {
    mcms[mcm] = getRegisterValue(regindex, mcm);
  }
}

// get all the values for a given register address, all mcms
void TrapConfigEvent::getAllMCMByAddress(const int registeraddress, std::array<uint32_t, o2::trd::constants::MAXMCMCOUNT>& mcms)
{
  int regindex = getRegIndexByAddr(registeraddress);
  if (regindex >= 0) {
    for (int mcm = 0; mcm < o2::trd::constants::MAXMCMCOUNT; ++mcm) {
      mcms[mcm] = getRegisterValue(regindex, mcm);
    }
  }
}

// get all the values for a given register values for all the mcm
void TrapConfigEvent::getAll(std::array<uint32_t, TrapRegisters::kLastReg * o2::trd::constants::MAXMCMCOUNT>& configdata)
{
  for (int mcm = 0; mcm < constants::MAXMCMCOUNT; ++mcm) {
    for (int reg = 0; reg < TrapRegisters::kLastReg; ++reg) {
      configdata[mcm * TrapRegisters::kLastReg + reg] = getRegisterValue(reg, mcm);
    }
  }
}

uint32_t TrapConfigEvent::getDmemUnsigned(const uint32_t address, const int detector, const int rob, const int mcm)
{
  LOGP(error, "Dmem data is not provided in config events");
  return 0;
}

uint32_t TrapConfigEvent::getTrapReg(const uint32_t index, const int detector, const int rob, const int mcm)
{

  uint32_t data = 0;
  uint32_t mcmidx;
  if ((detector >= 0 && detector < o2::trd::constants::MAXCHAMBER) &&
      (rob >= 0 && rob < o2::trd::constants::NROBC1) &&
      (mcm >= 0 && mcm < o2::trd::constants::NMCMROB + 2)) {
    mcmidx = HelperMethods::getMCMId(detector, rob, mcm);
    data = getRegisterValue(index, mcmidx);
  }
  return data;
}

bool TrapConfigEvent::isConfigDifferent(const TrapConfigEvent& trapconfigevent) const
{
  // is other different from this.
  // v1 walk through all 32 bit ints and compare, ignore the internals, ones to ignore are 32 bit themselves.
  // TODO do we care where the difference is?
  uint32_t max = constants::MAXMCMCOUNT;
  uint32_t maxreg = TrapRegisters::kLastReg;
  // start with the biggest granularity and work down.
  // compare hcid in the 2.
  // this would be simpler with != but that is removed in c++20
  if (getHCIDPresent() == trapconfigevent.getHCIDPresent()) {
    // we can continue if the bitpattern of half chambers is the same.
  } else {
    for (int hcid = 0; hcid < constants::MAXHALFCHAMBER; ++hcid) {
      if (isHCIDPresent(hcid) == 0 && trapconfigevent.isHCIDPresent(hcid) == 1) {
        // it is not in the ccdb but now is present, this is a change and needs to be saved.
        LOGP(info, " hcid {} is present in new but not in ccdb version", hcid);
        return false;
      }
      // other cases we can continue.
      // 1. both present its ok.
      // 2. present in ccdb but not in current one.
    }
  }
  // could check the rob/side present/notpresent as an optimisation?
  if (getMCMPresent() == trapconfigevent.getMCMPresent()) {
    // we can continue the bit pattern of mcm in the config are the same
  } else {
    for (int mcmid = 0; mcmid < constants::MAXMCMCOUNT; ++mcmid) {
      if (isMCMPresent(mcmid) == 0 && trapconfigevent.isMCMPresent(mcmid) == 1) {
        // it is not in the ccdb but now is present, this is a change and needs to be saved.
        return false;
      }
      // other cases we can continue.
      // 1. both present its ok.
      // 2. present in ccdb but not in current one.
    }
  }
  // Now the long part, compare the register data
  // There is no need to unpack registers. We can simply store their underlying 32 bit stored value.
  for (int mcm = 0; mcm < max; ++mcm) {
    for (int rawoffset = 0; rawoffset < kTrapRegistersSize; ++rawoffset) {
      if (!ignoreWord(rawoffset)) { // we do indeed care if this register is different.
                                    //        if (mConfigData[mcm].mRegisterData[rawoffset] != trapconfigevent.mConfigData[mcm].mRegisterData[rawoffset]) {
        return false;
        //        }
      }
    }
  }
  return true;
}

TrapRegInfo::TrapRegInfo(const std::string& name, int addr, int nBits, int base, int wordoffset, bool ignorechange, uint32_t max)
{
  init(name, addr, nBits, base, wordoffset, ignorechange, max);
}

TrapRegInfo::~TrapRegInfo() = default;

void TrapRegInfo::init(const std::string& name, int addr, int nbits, int base, int wordnumber, bool ignorechange, uint32_t max)
{
  // initialise a TRAP register information
  // uint32_t mNbits,mBase,mWordNumber,mShift,mMask,mMax;
  // see headerfile for exact definition, WordNumber: the sequence of register in a block of data, like registers.
  int bitoffset;
  int gapbits = 0;
  int bitwordoffset32; // the beginning of the 32 bit word containing this reg
  int packedwordsize = 30;
  if (addr != 0) {
    mAddr = addr;
    mName = name;
    mNbits = nbits;
    mBase = base;
    mWordNumber = wordnumber;
    mMax = pow(2, max) - 1; // this can be different from the mask, not sure why.
    packedwordsize = 30;
    if (mNbits > 30 || mNbits == 16 || mNbits == 4) { // these are 32 bit aligned the rest are 30 bit aligned.
      packedwordsize = 32;
      if (mNbits == 31) {
        gapbits = 1;
      }
    }
    mDataWordNumber = mWordNumber * (nbits + gapbits) / packedwordsize; // the 32 bit word offset to the word containing the register in question
    bitoffset = mWordNumber * (nbits)-mDataWordNumber * packedwordsize; // the offset bit with in that word for the start of this register
    mMask = (1 << mNbits) - 1;                                          // e.g. 5 = 0x1f 10=0x3ff
    mShift = 32 - bitoffset - mNbits;                                   // 32-remainder= the right hand side of bits that need to be dropped.
    if (mNbits == 32) {
      mShift = 0;
      mMask = 0xffffffff;
    }
    if (mNbits == 31) {
      mShift = 1;
    }
  } else {
    LOGP(warn, "Initialising an TRAP register with address of {:08x} ", addr);
    mAddr = 0;
    mName = "";
    mNbits = 0;
    mBase = 0;
    mWordNumber = 0;
    mMax = 0;
    mShift = 0;
  }
}

void TrapRegInfo::logTrapRegInfo()
{
  LOGP(info, " TrapRegInfo : {} with nbits={} addr {:08x} mask {:04x} word number {} and baseword {} max {} ", getName(), getNbits(), getAddr(), getMask(), getWordNumber(), getBase(), getMax());
}

void TrapConfigEvent::print()
{
  // walk through MCMSeen, and print out the mcm seen for this config event.
  uint32_t startmcm = 0;
  uint32_t endmcm = 0;
  uint32_t lastseenmcm = 0;
  bool continuousmcm = false;
  std::string mcmList;
  std::string totalList;
  totalList += "MCMs seen:";
  LOGP(debug, "Which MCM were seen in this event so far:");

  uint32_t mcmposition = 0;
  for (int mcmcount = 0; mcmcount < constants::MAXCHAMBER; ++mcmcount) {
    auto seenmcm = 0; // TODO pull from map of maps isMCMPresent(mcmcount);

    LOGP(debug, "{} {} {} :: pos: {} : start: {} end: {} seen {} lastseen {} mcmList {}", __FILE__, __func__, __LINE__, mcmposition, startmcm, endmcm, seenmcm, lastseenmcm, mcmList);
    if (seenmcm == 0 && lastseenmcm == 1) {
      // dump string
      mcmList += fmt::format("-{},", mcmposition - 1);
      //    LOGP(info,"{}",mcmList);
      totalList += mcmList;
      mcmList = "";
    }
    if (lastseenmcm == 0 && seenmcm == 1) {
      startmcm = seenmcm;
      mcmList = fmt::format("{}", mcmposition);
    }
    if (lastseenmcm == 1 && seenmcm == 1) {
      // do nothing but add a - for 1 or more of these.
    }
    lastseenmcm = seenmcm;
    mcmposition++;
  }
  totalList += mcmList;
  LOGP(debug, "{}", totalList);
}

void TrapConfigEvent::merge(const TrapConfigEvent* prev)
{
  LOGP(info, " Merge called for TrapConfigEvent {} {} {}", __FILE__, __func__, __LINE__);
  // std::map<uint32_t, std::map<uint32_t, uint32_t>> mTrapValueFrequencyMap;                                      //!< count of different value in the registers for a mcm,register used to find most frequent value.
  // take the 2 slots (trapconfigeventslot) and merge the update the map of maps of value.
  // this will be collapsed in the finalise of the calibrator, for the object to be written to the ccdb.
}

void TrapConfigEvent::fill(const TrapConfigEvent& input)
{
  LOGP(info, " fill called for TrapConfigEvent {} {} {}", __FILE__, __func__, __LINE__);
}

void TrapConfigEvent::fill(const gsl::span<const TrapConfigEvent> input)
{
  LOGP(info, " fill called for TrapConfigEvent {} {} {}", __FILE__, __func__, __LINE__);
}
uint32_t o2::trd::mcmtrapevent::getRegister2(const uint regidx) { return (int)(tpl[0].x[0]); };
