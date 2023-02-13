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
//  TRAP config                                                           //
//                                                                        //
//  Author: S. Murray (murrays@cern.ch)                                   //
////////////////////////////////////////////////////////////////////////////

#include "DataFormatsTRD/TrapConfigParser.h"
#include "DataFormatsTRD/Constants.h"
#include <fairlogger/Logger.h>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <array>
#include <bitset>
#include "DataFormatsTRD/HelperMethods.h"
#include "DataFormatsTRD/RawData.h"

#include <TH2F.h>
#include <TFile.h>
#include <TCanvas.h>

using namespace o2::trd;

TrapConfigParser::TrapConfigParser()
{
  // default constructor
}

TrapConfigParser::~TrapConfigParser() = default;


bool TrapConfigParser::checkRegister(uint16_t& registeraddr, int& registerbase, int& registeroffset, int& lastregisterindex, int& currentregisterindex, int& registerscount,
                                     int& registererrorgap, uint32_t& lastregisteraddressread, int& registerwordscount, int& lastregidx, uint32_t& registerdata, int mcmid, int hcid, int mcm, int rob)
{
  lastregisteraddressread = registeraddr;
  int32_t numberbits = 0;
  std::string regname = "";
  int32_t newregidx = 0;
  bool badregister = true;
  mTrapConfig.getRegisterByAddr(registeraddr, regname, newregidx, numberbits);
  registeroffset = registerbase + newregidx;
  // Are the registers sequential
  if (newregidx != lastregidx + 1 && newregidx != 0) {
    mMcmParsingStatus[mcmid] = 5; // missing registers
    for (int miss = lastregidx; miss < newregidx; ++miss) {
      mMissedReg[miss]++; // count the gaps, we count the start and stop point of gaps elsewhere
    }
    LOGP(warn, " register step mismatch : went from {:06x} ({}) i:{} to {:06x} ({}) i:{} and mcmid: {} mcmid size : {} mMcmParsingStatus[{}]=5", lastregisteraddressread, mTrapConfig.getRegNameByAddr(lastregisteraddressread), lastregidx, registeraddr, mTrapConfig.getRegNameByAddr(registeraddr), newregidx, registerdata, mcmid, mMcmParsingStatus.size(), mcmid);
  }
  lastregidx = newregidx; //TODO rename lastregidx as its different from lastregisterindex :-(
  if (numberbits >= 0 || regname != "" || newregidx >= 0) {
    //this is a bogus or unknown register
    LOGP(debug, "good register : name:{} newregindex:{} numberofbits:{}, lastregindex:{} registeraddr:{:08x} ?= ", regname, newregidx, numberbits, lastregisterindex, registeraddr);
    currentregisterindex = newregidx;
    if (currentregisterindex < lastregisterindex) {
      mMcmParsingStatus[mcmid] = 4;     // no end
      mMcmParsingStatus[mcmid + 1] = 6; // no begining
      mStopReg[lastregisterindex]++;
      mStopReg[currentregisterindex]++;
      LOGP(warn, " current register index is less than previous, we have looped back mcmId[{}] = {} and mcmId[{}]= {}", mcmid, mMcmParsingStatus[mcmid], mcmid + 1, mMcmParsingStatus[mcmid + 1]);
      for (int miss = registerscount; miss < newregidx; ++miss) {
        mMissedReg[miss]++; // count the gaps, we count the start and stop point of gaps elsewhere
                            //   mcmMissedRegister[currentrob * constants::NMCMROB + currentmcm].set(miss);
      }
    } else {
      if (registeraddr != mTrapConfig.getRegAddrByIdx(registerscount)) {
        //if (registeraddr != std::get<1>(TrapRegisterMap[registerscount]))
        registererrorgap += abs((int)registerscount - (int)newregidx) + 1; // +1 as reg index is zero based.
        mMcmParsingStatus[mcmid] = 3;                                      // missing registers
        LOGP(debug, " get mTrapConfig.getRegNameByAddr( {:08x} ) at func:{} line:{}", registeraddr, __func__, __LINE__);
        auto tmpnamebyaddr = mTrapConfig.getRegNameByAddr(registeraddr);
        LOGP(debug, " got mTrapConfig.getRegNameByAddr( {:08x} ) at func:{} line:{}", registeraddr, __func__, __LINE__);
        std::string tmpnamebyidx = ""; //mTrapConfig.getRegNameByIdx(registerscount);
        auto tmpregidxbyaddr = mTrapConfig.getRegIndexByAddr(registeraddr);
        LOGP(warn, " current register is not sequential to previous register, we have a gap : {} from {} to {} registererrorgap now : {} name compare {} ?= {} at {} mMcmParsingStatus[{}]={}", abs((int)registerscount - (int)newregidx),
             registerscount, tmpregidxbyaddr, registererrorgap, tmpnamebyaddr, tmpnamebyidx, __LINE__, mcmid, mMcmParsingStatus[mcmid]);
        registerscount = newregidx + 1;
        for (int miss = registerscount; miss < newregidx; ++miss) {
          mMissedReg[miss]++; // count the gaps, we count the start and stop point of gaps elsewhere
        }
        LOGP(warn, " Registers are non-sequential : {} from {} to {} registererrorgap now : {} name compare {} ?= {} at {}", abs(registerscount - mTrapConfig.getRegIndexByAddr(registeraddr)), registerscount, mTrapConfig.getRegIndexByAddr(registeraddr), registererrorgap, mTrapConfig.getRegNameByAddr(registeraddr), "" /*mTrapConfig.getRegNameByIdx(registerscount)*/, __LINE__);
      } else {
        LOGP(debug, " Registers are sequential : {} from {} to {} registererrorgap now : {} name compare {} ?= {} at {}", abs(registerscount - mTrapConfig.getRegIndexByAddr(registeraddr)), registerscount, mTrapConfig.getRegIndexByAddr(registeraddr), registererrorgap, mTrapConfig.getRegNameByAddr(registeraddr), "" /*mTrapConfig.getRegNameByIdx(registerscount)*/, __LINE__);
        registerscount++;
        badregister = false;
      }
    }
    registerwordscount++;
  } else {
    LOGP(warn, "bad register : name:'{}' newregindex:{} numberofbits:{}, lastregindex:{} registeraddr:{:08x} ?= ", regname, newregidx, numberbits, lastregisterindex, registeraddr);
  }
  lastregisterindex = currentregisterindex;
  return badregister;
}

void TrapConfigParser::compareToTrackletsHCID(std::bitset<1080> trackletshcid)
{
  //loop over config hcid and if its not present check if the config event had tracklets.
  //this is of course not conclusive but a start.
  for (int i = 0; i < 1080; ++i) {
    if (mHCIDhasConfig.test(i)) {
      if (trackletshcid.test(i)) {
        LOGP(/*info*/debug, "Config event had tracklets for hcid {}", i);
      } else {
        LOGP(/*info*/debug, "Config event had no tracklets for hcid {}", i);
      }
    } else {
      if (trackletshcid.test(i) && !mHCIDhasConfig.test(i)) {
        LOGP(/*info*/debug, "No Config event but we have tracklets for HCID  {}", i);
      }
    }
  }
}

void TrapConfigParser::writeFile(int eventcount)
{
  //TFile *file = TFile(Form("ccdbtrapconfig_%d.root",eventcount);
  //mTrapConfig.Write();
  //file->Close();
}

void TrapConfigParser::FillHistograms(int eventnum)
{ //TODO move this to a script to run over a TrapConfig3
  // take the mMcmParsingStatus and fill in the per mcm per layer graph
  // 1 hist per layer, 6 layers
  TH2F* layerconfig[6];
  std::unique_ptr<TFile> file(TFile::Open(Form("Event_%d_Histograms.root", eventnum), "RECREATE"));
  LOGP(/*info*/debug, "Now to fill Event_%i_Histograms.root", eventnum);
  TH1F *missed, *start, *stop, *regcounts;
  missed = new TH1F(Form("MissedRegistersEvent_%d", eventnum), Form("Missed registers for mcm that read out in event %d;register;count", eventnum), 500, -0.5, 499.5);
  start = new TH1F(Form("StartRegistersEvent_%d", eventnum), Form("Starting registers for mcm that read out in event %d;register;count", eventnum), 500, -0.5, 499.5);
  stop = new TH1F(Form("StopRegistersEvent_%d", eventnum), Form("Stopping registers for mcm that read out in event %d;register;count", eventnum), 500, -0.5, 499.5);
  regcounts = new TH1F(Form("CounterRegistersEvent_%d", eventnum), Form("Counts of seen registers for mcm that read out in event %d;register;count", eventnum), 500, -0.5, 499.5);
  int loopcounter = 0;
  for (auto& counts : mStopReg) {
    LOG(info) << " in stopreg loop with count of " << loopcounter << " and counts of : " << counts;
    loopcounter++;
    stop->SetBinContent(loopcounter, counts);
  }
  loopcounter = 0;
  for (auto& counts : mStartReg) {
    loopcounter++;
    start->SetBinContent(loopcounter, counts);
  }
  loopcounter = 0;
  for (auto& counts : mMissedReg) {
    loopcounter++;
    missed->SetBinContent(loopcounter, counts);
  }
  loopcounter = 0;
  for (auto& counts : mRegisterCount) {
    loopcounter++;
    regcounts->SetBinContent(loopcounter, counts);
  }

  for (int i = 0; i < 6; i++) {
    layerconfig[i] = new TH2F(Form("ConfigErrorsPerLayerEvent_%d_layer_%d", eventnum, i), Form("Config Errors per mcm in layer %d event %d;stack;sector", i, eventnum), 76, -0.5, 75.5, 144, -0.5, 143.5);
  }
  for (int mcmid = 0; mcmid < o2::trd::constants::MAXMCMCOUNT; ++mcmid) {
    int detector = mcmid / 128; //hcid/2;
    int sm = detector / 30;
    int detLoc = detector % 30;
    int layer = detector % 6;
    int stack = detLoc / 6;
    int sec = detector / 30;
    int rem = mcmid - (mcmid / 128) * 128;
    int mcm = rem - (rem / 16) * 16;
    int rob = rem / 16;
    int row = rob + sm * 8;
    int col = stack * 16 + mcm;
    layerconfig[layer]->SetBinContent(col, row, mMcmParsingStatus[mcmid]);
    //       LOGP(/*info*/debug,"mcmindex : {} is non zero mcmid:{} col:{} row:{} event:{}",mMcmParsingStatus[mcmid],mcmid,col,row,eventnum);
    //    }
  }
  std::unique_ptr<TCanvas> canvas(new TCanvas("canvas", "Config Event parsing problems"));
  canvas->Divide(3, 2);
  for (int i = 0; i < 6; i++) {
    canvas->cd(i + 1);
    layerconfig[i]->SetStats(0);
    layerconfig[i]->Draw("col");
  }
  canvas->Modified(true);
  canvas->Update();
  canvas->Write();
  //canvas->Print("ConfigEvents.gif+");
  file->Write();
  mMcmParsingStatus.fill(0);
  LOG(info) << " Walk through frequency map : ";
  int regcount = 0;
  for (auto& regmap : mTrapRegistersFrequencyMap) {
    LOG(info) << "Register : " << regcount << " " << mTrapConfig.getRegNameByIdx(regcount) << " with " << regmap.size() << " entries";
    for (const auto& elem : regmap) {
      LOGP(/*info*/debug, "[{:08x}] = {}", elem.first, elem.second);
    }
    regcount++;
  }

  /*  this is huge, enable for debugging purposes at your peril
  for (int mcm = 0; mcm < o2::trd::constants::MAXMCMCOUNT; ++mcm) {
    LOGP(info, "MCM : {} ", mcm);
    for (int reg = 0; reg < kLastReg; ++reg) {
      LOGP(info, "register of {} = {:08x}", mTrapConfig.getRegNameByIdx(reg), mCurrentMCMRegisters[mcm * kLastReg + reg]);
  }
  }*/
}

void TrapConfigParser::printMCMRegisterCount(int hcid)
{
  int roboffset = 1;
  if (hcid % 2 == 0) {
    roboffset = 0;
  }
  std::stringstream errorMCM;
  LOGP(/*info*/debug, "bp rob for hcid : {}....", hcid);
  for (int robidx = roboffset; robidx < 8; robidx += 2) {
    std::stringstream display;
    display << "bp rob:" << robidx << " ";
    for (int mcmidx = 0; mcmidx < 16; ++mcmidx) {
      display << fmt::format("[{:04} ({:04})] ", mcmSeen[robidx * 16 + mcmidx], mcmSeenMissedRegister[robidx * 16 + mcmidx]);
      if (mcmSeen[robidx * 16 + mcmidx] != 433)
        errorMCM << fmt::format("[{:04} ({:04}) == hcid:{} mcm:{} rob:{} ] ", mcmSeen[robidx * 16 + mcmidx], mcmSeenMissedRegister[robidx * 16 + mcmidx], hcid, mcmMCM[robidx * 16 + mcmidx], mcmROB[robidx * 16 + mcmidx]);
    }
    LOG(info) << display.str();
  }
  LOG(info) << errorMCM.str();
  LOG(info) << "bp rob ....";
}

void TrapConfigParser::unpackBlockHeader(uint32_t& header, uint32_t& registerdata, uint32_t& step, uint32_t& bwidth, uint32_t& nwords, uint16_t& registeraddr, uint32_t& exit_flag)
{
  registerdata = (header >> 2) & 0xFFFF; // 16 bit data
  step = (header >> 1) & 0x0003;
  bwidth = ((header >> 3) & 0x001F) + 1;
  //check that the bit width matches what is should be
  nwords = (header >> 8) & 0x00FF;
  registeraddr = (header >> 16) & 0xFFFF;
  exit_flag = (step == 0) || (step == 3) || (nwords == 0);
}

int TrapConfigParser::parseLink(std::array<uint32_t, o2::trd::constants::HBFBUFFERMAX>& data, uint32_t start, uint32_t end, int hcid)
{
  mHCIDhasConfig.set(hcid);
  uint32_t step, bwidth, nwords, idx, err, exit_flag;
  int32_t bitcnt, werr;
  uint16_t registeraddr;
  uint32_t registerdata, msk, header, data_hi, rawdat;
  uint16_t currentmcm, currentrob;
  uint32_t previousregisterread = 0;
  int32_t mcmid = 0, registerbase = 0; // base offset for the first register of a particular mcm
  int32_t registeroffset = 0;          // offset of a register into the raw storage area
  int32_t previousregisterindex = 0;
  int32_t currentregisterindex = 0;
  DigitMCMHeader mcmheader;
  mcmSeen.fill(-1);
  mcmMCM.fill(-1);
  mcmROB.fill(-1);
  mcmSeenMissedRegister.fill(-1);
  //  for(auto regbit : mcmMissedRegister){ // bitpattern of which registers were seen and not seen for a given mcm.
  //    regbit.reset(); // set bitset to 0
  //  }
  // mcmheader, header, then data
  idx = 0; // index in the packed configuration
  err = 0;
  int lastregidx = 0;
  int previdx = 0;
  int datalength = data.size();

  int rawidx = 0;
  int registerwordscount = 0;
  int registerscount = 0;
  int registererrorgap = 0;
  bool endmarker = false;
  bool fastforward = false;
  // for (idx=start;idx<end/2;++idx) {
  //   LOGP(/*info*/debug, " data[{} = {:08x}]", idx, data[idx]);
  // }
  idx = start;
  auto whileloopstart = std::chrono::high_resolution_clock::now();
  bool firstfastforward = true;
  while (idx < end) {
    LOGP(debug, "************** Raw data : data[{}] = {:08x}", idx, data[idx]);
    if (fastforward) {
      //loop until the next digitmcmheader, or end.
      //search for :
      // a : end marker
      // b : mcmheader
      // c : end of data
      while (data[idx] && fastforward) {
        if (firstfastforward) {
          LOGP(/*info*/debug, "fastforwaring from idx:{} for mcm {}", idx, mcmid);
          firstfastforward = false;
        }
        //read until we find an end marker, ignoring the data coming in so as not to pollute the configs.
        if (data[idx] == 0x7fff00fe) {
          //LOGP(warn, " end marker found ");
          ++idx; // go past the end marker with the expectation of a DigitMCMHeader to come next.
          endmarker = true;
          fastforward = false;
          continue;
        } else {
          //LOGP(warn, " no end marker found ");
        }
        ++idx;
      }
      fastforward = false;
      continue;
    }
    if (idx == start || endmarker) {
      //mcm header
      mcmheader.word = data[idx];
      LOGP(debug, "header at data[{}] had : {:06} registers, last register read : {:08x} ({}) and a register gap of {}", idx, registerwordscount, previousregisterread, mTrapConfig.getRegNameByAddr(previousregisterread), registererrorgap);
      //printDigitMCMHeader(mcmheader);
      if (idx != start) {
        int index = currentrob * constants::NMCMROB + currentmcm;
        if (index < constants::NMCMROB * constants::NROBC1) {
          mcmSeen[currentrob * constants::NMCMROB + currentmcm] = registerwordscount;             // registers seen for this mcm
          mcmSeenMissedRegister[currentrob * constants::NMCMROB + currentmcm] = registererrorgap; // registers missed for this mcm
          mcmMCM[currentrob * constants::NMCMROB + currentmcm] = mcmheader.mcm;                   // registers seen for this mcm
          mcmROB[currentrob * constants::NMCMROB + currentmcm] = mcmheader.rob;                   // registers seen for this mcm
        }
      }
      registerwordscount = 0;    // count of register words read for this mcm
      registerscount = 0;        // count of registers read for this mcm
      registererrorgap = 0;      // reset the count as we are starting a fresh
      previousregisterindex = 0; // register index of the last register read.
      currentregisterindex = 0;  // index of the current register
      currentmcm = mcmheader.mcm;
      currentrob = mcmheader.rob;
      mcmid = (hcid / 2) * 128 + currentrob * constants::NMCMROB + currentmcm; // current rob is 0-7/0-5 and currentmcm is 0-16.
      mMcmParsingStatus[mcmid] = 1;
      if (data[idx] == 0) {
        LOG(warn) << "Breaking as a zero after endmarker";
        break;
      }
      mcmSeen[currentrob * 16 + currentmcm] = 0;
      endmarker = false;
      lastregidx = 0;
      registerbase = mcmid * 432;
      ++idx;
      continue; // dont parse and go back through the loop
    }
    if (data[idx] == 0x7fff00fe || (data[idx] == 0x7ffff && data[idx + 1] == 0x7fff00fe)) {
      //end marker next value which should be a header.
      LOGP(debug, "after end marker Header :  address:{0:08x}  words: {1:08x} width: {2:08x} astep: {3:08x} zero: {4:08x} at idx : {5:0d} gap: {6:0d}", (data[idx + 1] >> 16) & 0xffff, (data[idx + 1] >> 8) & 0xff, (data[idx + 1] >> 3) & 0x1f, (data[idx + 1] >> 1) & 0x3, data[idx + 1] & 0x1, idx, idx - previdx);
      previdx = idx;
      endmarker = true;
      ++idx;
      if (data[idx] == 0x7fff00fe) {
        ++idx; // handle the second case of the previous if statement.
      }
      continue;
    }
    header = data[idx];
    /*****************
    SINGLE DATA
    *****************/
    if (header & 0x01) {                      // single data
      registerdata = (header >> 2) & 0xFFFF;  // 16 bit data
      registeraddr = (header >> 18) & 0x3FFF; // 14 bit address
      LOGP(debug, "single data raw:{:08x} idx: addr : {:08x} data : {:08x}", data[idx], idx, registeraddr, registerdata);
      ++idx;
      if (registeraddr != 0x1FFF) {
        if (header & 0x02) { // check if > 16 bits
          data_hi = data[idx];
          LOGP(debug, "read {:08x} for data > 16 bits", data_hi);
          ++idx;
          err += ((data_hi ^ (registerdata | 1)) & 0xFFFF) != 0;
          registerdata = (data_hi & 0xFFFF0000) | registerdata;
        }
        auto badreg = checkRegister(registeraddr, registerbase, registeroffset, previousregisterindex, currentregisterindex, registerscount, registererrorgap, previousregisterread, registerwordscount, lastregidx, registerdata, mcmid, hcid, mcmheader.mcm, mcmheader.rob);
        if (badreg == true) {
          fastforward = true;
          continue;
        }
        if (mcmid < o2::trd::constants::MAXMCMCOUNT && registerscount - 1 < kLastReg) {
          LOGP(debug, "Adding single register {:08x} [{:08x}] name: {} for mcm {}, registerwordscount {} registerscount:{} regindex {} with badreg:{}", registeraddr, registerdata, mTrapConfig.getRegNameByAddr(registeraddr), mcmid, registerwordscount, registerscount, idx, badreg);
          //update frequency map:
          if (registerscount > 0) {
            mTrapConfig.setRegisterValueByIdx(registerdata, registerscount - 1, mcmid);
            auto regdata = registerdata; //mTrapConfig.getRegisterValueByIdx(registerscount-1,mcmid);
            auto regmax = mTrapConfig.getTrapRegInfoByIdx(registerscount - 1)->getMax();
            if (regdata > regmax) {
              LOGP(warn, "assumed corrupted data as register data is greater than the mask : {:08x} for max {:08x} registername : {} ", regdata, regmax, mTrapConfig.getRegNameByIdx(registerscount - 1));
            } else {
              mTrapRegistersFrequencyMap[registerscount - 1][regdata]++;
            }
            mCurrentMCMRegisters[mcmid * kLastReg + registerscount - 1] = regdata;
            //TODO this is not saving it to the CCDBConfig at all !
            //  if(mTrapConfig.getRegNameByIdx(registerscount-1) == "ADCMSK"){
            //  LOGP(/*info*/debug,"** just added {:08x} data to _ADCMSK for mcmid {}",registerdata,mcmid);
            //  }
          }
          mRegisterCount[mTrapConfig.getRegIndexByAddr(registeraddr)]++; // keep a count of seen and accepted registers
          if (mMcmParsingStatus[mcmid] < 2) {
            // this handles the gaps in registers, where it might be good (1) before and after the gap, but this should stay with status of gap.
            mMcmParsingStatus[mcmid] = 1;
          }
        }
        if (idx >= end) {
          LOGP(error, "(single-write): no more data, missing end marker\n");
          std::chrono::duration<double, std::micro> configparsingtime = std::chrono::high_resolution_clock::now() - whileloopstart;
          LOGP(warn, "Config while parsing took {} leaving parsing due to no more data at line {}", (double)std::chrono::duration_cast<std::chrono::microseconds>(configparsingtime).count(), __LINE__);
          return false;
        }

      } else {
        LOG(error) << "returning because of 1fff as register addr " << __func__ << " " << __LINE__;
        std::chrono::duration<double, std::micro> configparsingtime = std::chrono::high_resolution_clock::now() - whileloopstart;
        LOGP(warn, "Config while parsing took {} leaving parsing due to 1fff as addr", (double)std::chrono::duration_cast<std::chrono::microseconds>(configparsingtime).count());
        return err;
      }
    } else {
      /*****************
    BLOCK DATA
    *****************/
      LOGP(debug, "now to unpack");
      unpackBlockHeader(header, registerdata, step, bwidth, nwords, registeraddr, exit_flag);
      LOGP(debug, "unpacked, registeraddr:{:08x}", registeraddr);
      auto trapreg = mTrapConfig.getTrapRegInfoByAddr(registeraddr);
      if (trapreg != nullptr) {
        auto trapregindex = mTrapConfig.getRegIndexByAddr(registeraddr);
        if ((trapregindex != -1) && bwidth != trapreg->getNbits()) {
          //check that the bit width matches what it should be
          // something is corrupt. What we read does not match what we expect.
          // log to info for now until figured out. TODO take out info
          LOGP(warn, " probably corrupt data : bwidth of {} does not match expected bandwidth of {} for reg {} registeraddr of : {:08x} registerindex : {}", bwidth,
               trapreg->getNbits(), trapreg->getName(), registeraddr, trapregindex);
          trapreg->logTrapRegInfo();
          //TODO bail out but how far ? just mcm or whole link?
        }
      } else {
        LOGP(warn, "trapreg is nullptr");
      }
      if (exit_flag) {
        std::chrono::duration<double, std::micro> configparsingtime = std::chrono::high_resolution_clock::now() - whileloopstart;
        LOGP(warn, "Exit flag found. Config while parsing took {} leaving parsing due to the exit flag", (double)std::chrono::duration_cast<std::chrono::microseconds>(configparsingtime).count());
        return err;
      }
      if (bwidth == 31 || (bwidth > 4 && bwidth < 8) || bwidth == 10 || bwidth == 15) {
        //TODO the part after 31 is probably not required given the above if statement of bwidth, when its out mechanism is figured out.
        // only possible values for blocks of registers is 5, 6, 7, 10, 15, and 31
        msk = (1 << bwidth) - 1;
        bitcnt = 0;
        while (nwords > 0) {
          LOGP(debug, "bwidth {}: read {:08x}  ", bwidth, data[idx]);
          if (bwidth == 31)
            ++idx;
          --nwords;
          bitcnt -= bwidth;
          err += (data[idx] & 1);
          if (bwidth != 31 && bitcnt < 0) { // handle the settings for when there are multiple registers packed into 1 word
            LOGP(debug, "block next data [{}] {:08x} at line {} nwords {}", idx, data[idx], __LINE__, nwords);
            header = data[idx];
            idx++;
            err += (header & 1);
            header = header >> 1;
            bitcnt = 31 - bwidth;
            registerdata = header & (1 << bwidth) - 1;
          }
          auto badreg = checkRegister(registeraddr, registerbase, registeroffset, previousregisterindex, currentregisterindex, registerscount, registererrorgap, previousregisterread, registerwordscount, lastregidx, registerdata, mcmid, hcid, mcmheader.mcm, mcmheader.rob);
          if (badreg == true) {
            fastforward = true;
            continue;
          }
          if (mcmid < o2::trd::constants::MAXMCMCOUNT && registerscount < kLastReg) {
            LOGP(debug, "Adding block register {:08x} [{:08x}] name: {}  for mcm {} registerwordscount {} registerscount {} regindex {} header {:08x} with badreg:{}", registeraddr, registerdata, mTrapConfig.getRegNameByAddr(registeraddr), mcmid, registerwordscount, registerscount, idx, header, badreg);
            if (mCurrentMCMRegisters[mcmid * kLastReg + registerscount] != registerdata) {
              if (registerscount > 0) {
                mTrapConfig.setRegisterValueByIdx(registerdata, registerscount - 1, mcmid);
                auto regdata = registerdata; //mTrapConfig.getRegisterValueByIdx(registerscount-1,mcmid);
                auto regmax = mTrapConfig.getTrapRegInfoByIdx(registerscount - 1)->getMax();
                if (regdata > regmax) {
                  LOGP(warn, "assumed corrupted data as register data is greater than the mask : {:08x} for max {:08x} registername : {} ", regdata, regmax, "" /*mTrapConfig.getRegNameByIdx(registerscount-1)*/);
                } else {
                  mTrapRegistersFrequencyMap[registerscount - 1][regdata]++;
                  mCurrentMCMRegisters[mcmid * kLastReg + registerscount - 1] = registerdata;
                }
                mcmMissedRegister[currentrob * constants::NMCMROB + currentmcm].set(registerscount - 1);
              }
            } else {
              LOGP(debug, "if statement failed for currentmcmregister {:08x}  registerdata {:08x}", mCurrentMCMRegisters[mcmid * kLastReg + registerscount], registerdata);
            }

            mRegisterCount[mTrapConfig.getRegIndexByAddr(registeraddr)]++; // keep a count of seen and accepted registers
            // if( mTrapRegisters[].CanIgnore()==false){
            //   mMcMDataHasChanged[registerscount]=true;
            // }
            if (mMcmParsingStatus[mcmid] < 2) {
              // this handles the gaps in registers, where it might be good (1) before and after the gap, but this should stay with status of gap.
              mMcmParsingStatus[mcmid] = 1;
            }
          } else {
            LOGP(warn, "if (mCurrentMCMRegisters[mcmid * kLastReg + registerscount] != registerdata mcmid:{} kLastReg:{} registerscount:{} mCurrentMCMRegisters[mcmid*kLastReg+registerscount]=={} != {}", mcmid, kLastReg, registerscount, mCurrentMCMRegisters[mcmid * kLastReg + registerscount], registerdata);
          }
          registeraddr += step;
          header = header >> bwidth; // this is not used for the bwidth=31 case
          if (idx >= end) {
            LOGP(error, "no end markermore data, {} words read", idx);
            std::chrono::duration<double, std::micro> configparsingtime = std::chrono::high_resolution_clock::now() - whileloopstart;
            LOGP(warn, "Config while parsing took {} bombing out due to end of data at line : {}", (double)std::chrono::duration_cast<std::chrono::microseconds>(configparsingtime).count(), __LINE__);
            return false;
          }
        }
      } else {
        // must be corrupt data header.
        LOGP(warn, "Unknown bwidth from block header of {}", bwidth);
      }
      if (data[idx] != 0x7fff00fe) {
        LOGP(debug, "increment idx to {} due trailer of block data data[{}]={:08x} to data[{}]={:08x}", idx + 1, idx, data[idx], idx + 1, data[idx + 1]);
        ++idx; // jump over the block trailer
      }
      if (data[idx + 1] == 0x7fff00fe || data[idx] == 0x00007fff) {
        LOGP(debug, "end of block ignoring bogus : {:08x}?", data[idx]);
      }
    } // end block case
  }   // end while
  std::chrono::duration<double, std::micro> configparsingtime = std::chrono::high_resolution_clock::now() - whileloopstart;
  LOGP(/*info*/debug, "Config while parsing took {} end of loop at line : {}", (double)std::chrono::duration_cast<std::chrono::microseconds>(configparsingtime).count(), __LINE__);
  //loop over which mcms never sent data.
  printMCMRegisterCount(hcid);
  return false; // only if the max length of the block reached!
}
