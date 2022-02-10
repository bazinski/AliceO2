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

/// \file testTRDTrapConfig3
/// \brief This task tests the trap config
/// \author Sean Murray, murrays@cern.ch

#define BOOST_TEST_MODULE Test TRD_RawDataHeader
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include "DataFormatsTRD/TrapConfig3.h"

namespace o2
{
namespace trd
{

std::vector<int> mcmids = {1, 1, 127, 128, 12719, 69119, 69120, 69130};
//std::vector<int> mcmids={69119, 69120,69130};
std::vector<TrapRegInfo*> registersOfInterest;

void trapregcheck(TrapConfig3* trapconfig, uint32_t mcmidx)
{
  //loop over the mcm
  uint32_t newvalue;
  for (auto& mcm : mcmids) {
    //loop over the register
    for (int regcount = 0; regcount < registersOfInterest.size(); ++regcount) {
      uint32_t address = registersOfInterest[regcount]->getAddr();
      std::string name = registersOfInterest[regcount]->getName();
      uint32_t index = trapconfig->getRegIndexByAddr(address);
      for (int count = 0; count < 4; ++count) {
        // look at the zero value, first increment, mid point and max value.
        switch (count) {
          case 0:
            newvalue = 0; // zero value is ... zero
            break;
          case 1:
            newvalue = 1; // first non zero value is 1 bit, all registers are 1 or more bits.
            break;
          case 2:
            newvalue = registersOfInterest[regcount]->getMask() / 2; // midpoint value
            break;
          case 3:
            newvalue = registersOfInterest[regcount]->getMask(); // get max value, which is simply the mask;
            break;
        }
        uint32_t retval = trapconfig->setRegisterValueByIdx(newvalue, index, mcm);
        if (retval == -1) {
          newvalue = -1;
          // if the indices are wrong the set will return a -1 and so will the get methods.
          // This is to ensure they actually do.
        } else {
          BOOST_CHECK_EQUAL(trapconfig->getRegisterValueByIdx(index, mcm), newvalue);
          BOOST_CHECK_EQUAL(trapconfig->getRegisterValueByName(name, mcm), newvalue);
          BOOST_CHECK_EQUAL(trapconfig->getRegisterValueByAddr(address, mcm), newvalue);
        }
      }
    }
  }
}

/// \brief Test the trap register generation functions
//
BOOST_AUTO_TEST_CASE(TRDTrapConfig3Internals)
{
  // test trap register initialisation
  TrapRegInfo trapreg;
  trapreg.init("testreg", 0x3000, 6, 0, 1, false, 6);
  BOOST_CHECK_EQUAL(trapreg.getShift(), 20);
  BOOST_CHECK_EQUAL(trapreg.getMask(), 0x3f);
  BOOST_CHECK_EQUAL(trapreg.getDataWordNumber(), 0);
  trapreg.init("testreg", 0x3000, 6, 0, 5, false, 6);
  BOOST_CHECK_EQUAL(trapreg.getShift(), 26);
  BOOST_CHECK_EQUAL(trapreg.getMask(), 0x3f);
  BOOST_CHECK_EQUAL(trapreg.getDataWordNumber(), 1);
  trapreg.init("testreg", 0x3000, 5, 0, 1, false, 5);
  BOOST_CHECK_EQUAL(trapreg.getShift(), 22);
  BOOST_CHECK_EQUAL(trapreg.getMask(), 0x1f);
  BOOST_CHECK_EQUAL(trapreg.getDataWordNumber(), 0);
  trapreg.init("testreg", 0x3000, 5, 0, 4, false, 5);
  BOOST_CHECK_EQUAL(trapreg.getShift(), 7);
  BOOST_CHECK_EQUAL(trapreg.getMask(), 0x1f);
  BOOST_CHECK_EQUAL(trapreg.getDataWordNumber(), 0);
  trapreg.init("testreg", 0x3000, 5, 0, 5, false, 5);
  BOOST_CHECK_EQUAL(trapreg.getShift(), 2);
  BOOST_CHECK_EQUAL(trapreg.getMask(), 0x1f);
  BOOST_CHECK_EQUAL(trapreg.getDataWordNumber(), 0);
  trapreg.init("testreg", 0x3000, 5, 0, 6, false, 5);
  BOOST_CHECK_EQUAL(trapreg.getShift(), 27);
  BOOST_CHECK_EQUAL(trapreg.getMask(), 0x1f);
  BOOST_CHECK_EQUAL(trapreg.getDataWordNumber(), 1);
  trapreg.init("testreg", 0x3000, 31, 0, 1, false, 31);
  BOOST_CHECK_EQUAL(trapreg.getShift(), 1);
  BOOST_CHECK_EQUAL(trapreg.getMask(), 0x7fffffff);
  BOOST_CHECK_EQUAL(trapreg.getDataWordNumber(), 1);
  trapreg.init("testreg", 0x3000, 31, 0, 0, false, 31);
  BOOST_CHECK_EQUAL(trapreg.getShift(), 1);
  BOOST_CHECK_EQUAL(trapreg.getMask(), 0x7fffffff);
  BOOST_CHECK_EQUAL(trapreg.getDataWordNumber(), 0);
  trapreg.init("testreg", 0x3000, 32, 0, 0, false, 31);
  BOOST_CHECK_EQUAL(trapreg.getShift(), 0);
  BOOST_CHECK_EQUAL(trapreg.getMask(), 0xffffffff);
  BOOST_CHECK_EQUAL(trapreg.getDataWordNumber(), 0);
  trapreg.init("testreg", 0x3000, 32, 0, 1, false, 31);
  BOOST_CHECK_EQUAL(trapreg.getShift(), 0);
  BOOST_CHECK_EQUAL(trapreg.getMask(), 0xffffffff);
  BOOST_CHECK_EQUAL(trapreg.getDataWordNumber(), 1);
  /*  LOGP(info," trapreginfo done");
// now to systematically test the setting and reading back of trap register values by : 
  // a. size of register
  // b. method of access, address/name/index.
  LOGP(info," lets declare a trapconfig object");*/
}

BOOST_AUTO_TEST_CASE(TRDTrapConfig3GetSet)
{

  TrapConfig3* trapconfig = new trd::TrapConfig3();

  // setup the resgisters we will look at chosen for various reasons, size, on the edges, changes of register bit size
  trd::registersOfInterest.push_back(trapconfig->getTrapRegInfoByAddr(0x3180)); //TPL00 the first register
  trd::registersOfInterest.push_back(trapconfig->getTrapRegInfoByAddr(0x3185)); // TPL05 last register for the first 32 bit word
  trd::registersOfInterest.push_back(trapconfig->getTrapRegInfoByAddr(0x31ff)); // TPL7F
  trd::registersOfInterest.push_back(trapconfig->getTrapRegInfoByAddr(0x30A0)); // FGA0
  trd::registersOfInterest.push_back(trapconfig->getTrapRegInfoByAddr(0x30A4)); // last 6 bit reg in a 32bit word
  trd::registersOfInterest.push_back(trapconfig->getTrapRegInfoByAddr(0x30A5)); // the next 6 bit register first in the subsequent 32 bit word
  trd::registersOfInterest.push_back(trapconfig->getTrapRegInfoByAddr(0x30B4)); // last 6 bit register
  trd::registersOfInterest.push_back(trapconfig->getTrapRegInfoByAddr(0x3080)); // first 10 bit register
  trd::registersOfInterest.push_back(trapconfig->getTrapRegInfoByAddr(0x308C)); // next 2 are 10 bit registers spanning a 32 bit data word
  trd::registersOfInterest.push_back(trapconfig->getTrapRegInfoByAddr(0x308D));
  trd::registersOfInterest.push_back(trapconfig->getTrapRegInfoByAddr(0x313F)); //
  trd::registersOfInterest.push_back(trapconfig->getTrapRegInfoByAddr(0x3000));
  trd::registersOfInterest.push_back(trapconfig->getTrapRegInfoByAddr(0x3002));
  trd::registersOfInterest.push_back(trapconfig->getTrapRegInfoByAddr(0x3003));
  trd::registersOfInterest.push_back(trapconfig->getTrapRegInfoByAddr(0x300F));
  trd::registersOfInterest.push_back(trapconfig->getTrapRegInfoByAddr(0x3020));
  trd::registersOfInterest.push_back(trapconfig->getTrapRegInfoByAddr(0x3022)); // last 10 bit word
  trd::registersOfInterest.push_back(trapconfig->getTrapRegInfoByAddr(0x3028)); // first 15 bit register after the above 10bit
  trd::registersOfInterest.push_back(trapconfig->getTrapRegInfoByAddr(0x302A)); // last 15 bit
  trd::registersOfInterest.push_back(trapconfig->getTrapRegInfoByAddr(0x3030)); // first 10bit after the 15bitregisters
  trd::registersOfInterest.push_back(trapconfig->getTrapRegInfoByAddr(0x3050)); // a lone 32 bit register
  trd::registersOfInterest.push_back(trapconfig->getTrapRegInfoByAddr(0x3051)); // 2 5 bit registers
  trd::registersOfInterest.push_back(trapconfig->getTrapRegInfoByAddr(0x3052));
  trd::registersOfInterest.push_back(trapconfig->getTrapRegInfoByAddr(0x3053)); // then a 32 bit register again
  trd::registersOfInterest.push_back(trapconfig->getTrapRegInfoByAddr(0x0B6C)); // first 15 bitreg
  trd::registersOfInterest.push_back(trapconfig->getTrapRegInfoByAddr(0x0B6D)); // last 15 bit reg in the 32 bit word
  trd::registersOfInterest.push_back(trapconfig->getTrapRegInfoByAddr(0x0B6E)); // first
  trd::registersOfInterest.push_back(trapconfig->getTrapRegInfoByAddr(0x0B6F)); // second
  trd::registersOfInterest.push_back(trapconfig->getTrapRegInfoByAddr(0x0B80)); // lone 32 bit reg
  trd::registersOfInterest.push_back(trapconfig->getTrapRegInfoByAddr(0x0B81)); // lone 16 bit reg
  trd::registersOfInterest.push_back(trapconfig->getTrapRegInfoByAddr(0xD000)); // first 10 after a 16
  trd::registersOfInterest.push_back(trapconfig->getTrapRegInfoByAddr(0xD001)); // middle 10 in a 32 bit word
  trd::registersOfInterest.push_back(trapconfig->getTrapRegInfoByAddr(0xD002)); // last 10 bit in a 32 bit word
  trd::registersOfInterest.push_back(trapconfig->getTrapRegInfoByAddr(0xD003)); // 10 bit in the subsequent 32 bit word
  trd::registersOfInterest.push_back(trapconfig->getTrapRegInfoByAddr(0x0D40));
  trd::registersOfInterest.push_back(trapconfig->getTrapRegInfoByAddr(0x0D41));
  for (auto& mcm : mcmids) {
    trapregcheck(trapconfig, mcm);
  }
}
} // namespace trd
} // namespace o2
