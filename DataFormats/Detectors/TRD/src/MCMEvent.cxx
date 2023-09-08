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

#include "DataFormatsTRD/MCMEvent.h"

namespace o2::trd
{

const uint32_t MCMEvent::getRegister(const uint32_t regidx, const TrapRegInfo& trapreg) const
{
  // get the register value based on the address of the register;
  // find register in mTrapRegisters.
  // calculate the offset from the base for the register and get the mask.
  int mcmoffset = 0;                           // mcmidx * kTrapRegistersSize;                                 // get the start offset for this mcm
  int regbase = mcmoffset + trapreg.getBase(); // get the base of this register in the underlying storage block
  int regoffset = trapreg.getBase() + trapreg.getWordNumber();
  // int regoffset = trapreg.getBase() + trapreg.getWordNumber(); // get the offset to the register in question
  auto size = mRegisterData.size();
  uint32_t data1 = mRegisterData[0];
  uint32_t data = mRegisterData[regoffset];
  data = data >> trapreg.getShift();
  // data = data >> trapreg.getShift(); // get the data and shift it as needed
  data &= trapreg.getMask(); // mask the data off as need be.
                             // LOGP(info, " returning data of {:08x}", data);
  return data;
}

bool MCMEvent::setRegister(const uint32_t data, const uint32_t regidx, const TrapRegInfo& trapreg)
{
  uint32_t regvalue = data;
  int mcmoffset = 0;                                           // mcmidx * kTrapRegistersSize;                          // get the start offset for this mcm
  int regoffset = trapreg.getBase() + trapreg.getWordNumber(); // wordnumber; // get the offset to the register in question
  regvalue &= trapreg.getMask();                               // mask the data off as need be.
  uint32_t notdatamask = ~(trapreg.getMask() << trapreg.getShift());
  regvalue = regvalue << trapreg.getShift();
  auto trapregvalue = mRegisterData[regoffset];
  mRegisterData[regoffset] = mRegisterData[regoffset] & notdatamask;
  mRegisterData[regoffset] = mRegisterData[regoffset] | regvalue;
  return true;
}

} // namespace o2::trd
