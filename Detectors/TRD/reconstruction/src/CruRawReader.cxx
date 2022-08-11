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

/// @file   CruRawReader.h
/// @brief  TRD raw data translator

#include "DetectorsRaw/RDHUtils.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "Headers/RDHAny.h"
#include "TRDReconstruction/CruRawReader.h"
#include "TRDBase/FeeParam.h"
#include "DataFormatsTRD/RawData.h"
#include "DataFormatsTRD/Tracklet64.h"
#include "TRDReconstruction/DigitsParser.h"
#include "TRDReconstruction/TrackletsParser.h"
#include "DataFormatsTRD/Constants.h"
#include "DataFormatsTRD/HelperMethods.h"

#include "Framework/ControlService.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/RawDeviceService.h"
#include "Framework/DeviceSpec.h"
#include "Framework/DataSpecUtils.h"
#include "Framework/Output.h"
#include "Framework/InputRecordWalker.h"
#include "DataFormatsCTP/TriggerOffsetsParam.h"

#include <cstring>
#include <string>
#include <vector>
#include <array>
#include <numeric>
#include <iomanip>

namespace o2::trd
{

bool CruRawReader::skipRDH()
{
  // check rdh for being empty or only padding words.
  if (o2::raw::RDHUtils::getMemorySize(mOpenRDH) == o2::raw::RDHUtils::getHeaderSize(mOpenRDH)) {
    // empty rdh so we want to avoid parsing it for cru data.
    if (mOptions[TRDVerboseBit]) {
      LOG(info) << " skipping rdh (empty) packetcount of: " << std::hex << o2::raw::RDHUtils::getPacketCounter(mOpenRDH);
    }
    return true;
  } else {

    if (mHBFPayload[0] == o2::trd::constants::CRUPADDING32 && mHBFPayload[0] == o2::trd::constants::CRUPADDING32) {
      // event only contains paddings words.
      if (mOptions[TRDVerboseBit]) {
        LOG(info) << " skipping rdh (padding) with packetcounter of: " << std::hex << o2::raw::RDHUtils::getPacketCounter(mOpenRDH);
      }
      // mDataPointer+= o2::raw::RDHUtils::getOffsetToNext()/4;
      auto rdh = reinterpret_cast<const o2::header::RDHAny*>(mDataPointer);
      mDataPointer += o2::raw::RDHUtils::getOffsetToNext(rdh) / 4;
      // mDataPointer=reinterpret_cast<const uint32_t*>(reinterpret_cast<const char*>(rdh) + o2::raw::RDHUtils::getOffsetToNext(rdh));
      return true;
      return true;
    } else {
      return false;
    }
  }
}

void CruRawReader::OutputAllLinkRawData()
{
  for (int link = 0; link < 15; ++link) {
    OutputLinkRawData(link);
  }
}

void CruRawReader::OutputLinkRawData(int link)
{
  std::array<uint32_t, 1024>::iterator linkstart, linkend, bufferoffset;
  uint32_t lengthoflinksbefore = 0;
  int linkindex = 0;
  for (linkindex = 0; linkindex < link; ++linkindex) {
    lengthoflinksbefore += mCurrentHalfCRULinkLengths[linkindex];
  }
  lengthoflinksbefore *= 8;                                 // to get to 32 bit word units
  lengthoflinksbefore += sizeof(mCurrentHalfCRUHeader) / 4; // 8 bit words to 32 bit words
  linkstart = mHBFPayload.begin() + mHalfCRUStartOffset + lengthoflinksbefore;
  linkend = linkstart + mCurrentHalfCRULinkLengths[link] * 8;
  // mHalfCRUStartOffset=cruhbfstartoffset;
  bufferoffset = linkstart;
  uint32_t wordcount = 0;
  while (bufferoffset < linkend) {
    std::stringstream outputstring;
    uint32_t stringlength = 8;
    if (linkend - bufferoffset < 8) {
      stringlength = linkend - bufferoffset;
    }
    for (int z = 0; z < stringlength; ++z) {
      outputstring << "0x" << std::hex << std::setw(8) << std::setfill('0') << *(bufferoffset + z) << " ";
    }
    bufferoffset += 8;
    LOG(info) << wordcount++ << ":" << outputstring.str();
  }
}

void CruRawReader::OutputHalfCruRawData(int offset)
{
  LOG(info) << "Full 1/2 CRU dump for FEEID:0x" << std::hex << mFEEID.word;
  std::stringstream outputstring;
  for (int z = 0; z < 15; ++z) {
    outputstring << z << ":" << mCurrentHalfCRULinkLengths[z] << " ";
  }
  LOG(info) << outputstring.str();
  int linkcount = 0;
  uint32_t linkzsum = 0;
  int bufferoffset = 0;
  uint32_t totalhalfcrulength = std::accumulate(mCurrentHalfCRULinkLengths.begin(),
                                                mCurrentHalfCRULinkLengths.end(),
                                                decltype(mCurrentHalfCRULinkLengths)::value_type(0));
  totalhalfcrulength *= 8; // convert from 256 bits to 32 bits.
  outputstring.clear();
  for (int z = 0; z < 15; ++z) {
    outputstring << z << ":" << mHBFPayload[bufferoffset + z] << " ";
  }
  LOG(info) << "CRH0 " << outputstring.str();
  outputstring.clear();
  bufferoffset += 8;
  for (int z = 0; z < 15; ++z) {
    outputstring << z << ":" << mHBFPayload[bufferoffset + z] << " ";
  }
  LOG(info) << "CRH1 " << outputstring.str();
  int outputoffset = 2;
  bufferoffset += 8;
  int olengthoffset0 = 0;
  for (int link = 0; link < 15; ++link) {
    OutputLinkRawData(link);
  }
  LOG(info) << "CRU Next 16 words just for info";
  for (int bufferoffset = totalhalfcrulength; bufferoffset < totalhalfcrulength + 16; bufferoffset += 8) {
    LOGP(info, "0x{0:06x} :: {1:08x} {2:08x}  {3:08x} {4:08x} {5:08x} {6:08x} {7:08x} {8:08x} ", bufferoffset, mHBFPayload[bufferoffset], mHBFPayload[bufferoffset + 1], mHBFPayload[bufferoffset + 2], mHBFPayload[bufferoffset + 3], mHBFPayload[bufferoffset + 4], mHBFPayload[bufferoffset + 5], mHBFPayload[bufferoffset + 6], mHBFPayload[bufferoffset + 7]);
  }
  LOG(info) << "Full 1/2 CRU dump end ***************";
}

void CruRawReader::dumpRDHAndNextHeader(const o2::header::RDHAny* rdh)
{
  LOG(info) << "######################### Dumping RDH incoming buffer ##########################";
  //  o2::raw::RDHUtils::printRDH(rdh);
  LOG(info) << "Now for the buffer breakdown";
  auto offsetToNext = o2::raw::RDHUtils::getOffsetToNext(rdh);
  for (int i = 0; i < offsetToNext / 4; ++i) {
    LOG(info) << "ii:" << i << " 0x" << std::hex << *((uint32_t*)rdh + i); //*(uint32_t*)(reinterpret_cast<const uint32_t*>(rdh) + i);
  }
  LOG(info) << "Next 8....";
  for (int i = 0; i < 8; ++i) {
    LOG(info) << "iii:" << i << " 0x" << std::hex << *((uint32_t*)rdh + i + offsetToNext); //*(uint32_t*)(reinterpret_cast<const uint32_t*>(rdh) + i+offsetToNext);
  }
  LOG(info) << "######################### Finished RDH incoming buffer ##########################";
}

bool CruRawReader::checkRDH(const o2::header::RDHAny* rdh)
{
  // first check for FEEID from unconfigured CRU
  TRDFeeID feeid;
  feeid.word = o2::raw::RDHUtils::getFEEID(rdh);
  if (((feeid.word) >> 4) == 0xfff) { // error condition is 0xfff? as the end point is known to the cru, but the rest is configured.
    std::stringstream message;
    message << "failed due to 0xfff? : " << std::hex << feeid.word << " whole feeid : " << std::hex << (unsigned int)feeid.word;
    incrementErrors(TRDFEEIDIsFFFF, message, 0, 0, 0, 0, 0, 0);
    return false;
  }
  if (feeid.supermodule > 17) {
    std::stringstream message;
    message << "failed due to supermodule :  " << std::dec << (int)feeid.supermodule << " whole feeid : " << std::hex << (unsigned int)feeid.word;
    incrementErrors(TRDFEEIDBadSector, message, 0, true, 0, 0, 0, 0);
    return false;
  }
  return true;
}

bool CruRawReader::compareRDH(const o2::header::RDHAny* firstrdh, const o2::header::RDHAny* rdh)
{
  if (o2::raw::RDHUtils::getFEEID(firstrdh) != o2::raw::RDHUtils::getFEEID(rdh)) {
    if (mMaxWarnPrinted > 0) {
      LOG(warn) << "ERDH FEEID are not identical in rdh.";
      checkNoWarn();
    }
    std::stringstream message;
    message << "ERDH FEEID are not identical in rdh.";
    incrementErrors(TRDParsingBadRDHFEEID, message, 1, true, 0, 0, 0, 0);
    return false;
  }
  if (o2::raw::RDHUtils::getEndPointID(firstrdh) != o2::raw::RDHUtils::getEndPointID(rdh)) {
    if (mMaxWarnPrinted > 0) {
      checkNoWarn();
    }
    std::stringstream message;
    message << "ERDH  EndPointID are not identical in rdh.";
    incrementErrors(TRDParsingBadRDHEndPoint, message, 1, true, 0, 0, 0, 0);
    return false;
  }
  if (o2::raw::RDHUtils::getTriggerOrbit(firstrdh) != o2::raw::RDHUtils::getTriggerOrbit(rdh)) {
    std::stringstream message;
    message << "ERDH  Orbit are not identical in rdh.";
    incrementErrors(TRDParsingBadRDHOrbit, message, 1, true, 0, 0, 0, 0);
    return false;
  }
  if (o2::raw::RDHUtils::getCRUID(firstrdh) != o2::raw::RDHUtils::getCRUID(rdh)) {
    std::stringstream message;
    message << "ERDH  CRUID are not identical in rdh.";
    incrementErrors(TRDParsingBadRDHCRUID, message, 1, true, 0, 0, 0, 0);
    return false;
  }
  if (o2::raw::RDHUtils::getPacketCounter(firstrdh) == o2::raw::RDHUtils::getPacketCounter(rdh) + 1) {
    std::stringstream message;
    message << "ERDH  PacketCounters are not sequential in rdh.";
    incrementErrors(TRDParsingBadRDHPacketCounter, message, 1, true, 0, 0, 0, 0);
    return false;
  }
  return true;
}

bool CruRawReader::processHBFs(int datasizealreadyread, bool verbose)
{
  if (mOptions[TRDVerboseBit]) {
    LOG(info) << "PROCESS HBF starting at " << std::hex << (void*)mDataPointer << " already read in : " << datasizealreadyread;
  }
  mDataRDH = reinterpret_cast<const o2::header::RDHAny*>(mDataPointer);
  mOpenRDH = reinterpret_cast<const o2::header::RDHAny*>((const char*)mDataPointer);
  auto rdh = mDataRDH;
  uint32_t totaldataread = 0;

  mState = CRUStateHalfCRUHeader;
  uint32_t currentsaveddatacount = 0;
  mTotalHBFPayLoad = 0;
  int loopcount = 0;
  int counthalfcru = 0;
  // store the endpoint and cruid for later checking of subsequent rdh for integrity
  const o2::header::RDHAny* firstRDH;
  firstRDH = reinterpret_cast<const o2::header::RDHAny*>(reinterpret_cast<const char*>(rdh));
  mHBFoffset32 = 0;
  // loop until RDH stop header
  while (!o2::raw::RDHUtils::getStop(rdh)) { // carry on till the end of the event.
    if (!checkRDH(rdh)) {
      if (mMaxErrsPrinted > 0) {
        LOG(error) << "RDH validity failure : Call on call, an flp is not configured. Check the tracklets per halfchamber id for gaps of 30 or larger to figure out which one, feeid:" << std::hex << o2::raw::RDHUtils::getFEEID(rdh) << " cruid:" << o2::raw::RDHUtils::getCRUID(rdh) << " packetcounter:" << o2::raw::RDHUtils::getPacketCounter(rdh);
        checkNoErr();
      }
      return false; // dump and run.
    }
    if (!compareRDH(firstRDH, rdh)) { // compare previous rdh detector info to the current rdh, they must be the same.
      return false;                   // dump and run.
    }
    firstRDH = rdh;
    auto headerSize = o2::raw::RDHUtils::getHeaderSize(rdh);
    auto memorySize = o2::raw::RDHUtils::getMemorySize(rdh);
    if (memorySize == 0) {
      LOG(warn) << "rdh memory size is zero";
      return false; // get out of here if the rdh says it has nothing.
    } else if (memorySize < headerSize) {
      LOG(error) << "RDH payload memory is negative: memorySize(" << memorySize << ") - headerSize(" << headerSize << ")";
      return false; // this does not make sense and should probably not happen
    }
    auto offsetToNext = o2::raw::RDHUtils::getOffsetToNext(rdh);
    size_t rdhpayload = memorySize - headerSize;
    mFEEID.word = o2::raw::RDHUtils::getFEEID(rdh);       // TODO change this and just carry around the curreht RDH
    mCRUEndpoint = o2::raw::RDHUtils::getEndPointID(rdh); // the upper or lower half of the currently parsed cru 0-14 or 15-29
    mCRUID = o2::raw::RDHUtils::getCRUID(rdh);
    mIR = o2::raw::RDHUtils::getTriggerIR(rdh);
    int packetCount = o2::raw::RDHUtils::getPacketCounter(rdh);
    // mDataEndPointer = (uint32_t*)((char*)rdh + offsetToNext);
    if (mOptions[TRDVerboseWordBit]) {
      LOG(info) << "RDH : mFEEID:" << mFEEID.word << " mCRUEndpoint:" << mCRUEndpoint << " mCRUID:" << mCRUID << " packetCount:" << packetCount << "rdhpayload:" << rdhpayload << " offsettonext:" << offsetToNext << " dmDataEndPointer(after move to rdh+offsetToNext): 0x" << std::hex << (void*)mDataEndPointer << " rdh is currently at 0x" << std::hex << (void*)rdh;
    }
    // copy the contents of the current rdh into the buffer to be parsed
    std::memcpy((char*)&mHBFPayload[0] + currentsaveddatacount, ((char*)rdh) + headerSize, rdhpayload);
    mTotalHBFPayLoad += rdhpayload;
    currentsaveddatacount += rdhpayload;
    totaldataread += offsetToNext;
    // move to next rdh
    auto oldRDH = rdh;
    rdh = reinterpret_cast<const o2::header::RDHAny*>(reinterpret_cast<const char*>(rdh) + offsetToNext);
    // increment the data pointer by the size of the stop rdh.
    mDataPointer = reinterpret_cast<const uint32_t*>(reinterpret_cast<const char*>(rdh) + o2::raw::RDHUtils::getOffsetToNext(rdh));

    if (o2::raw::RDHUtils::getStop(rdh) || o2::raw::RDHUtils::getOffsetToNext(rdh) < mDataBufferSize - mHBFoffset32) {
      // we can still copy into this buffer.
    } else {
      if (mMaxWarnPrinted > 0) {
        LOG(warn) << "rdh bounds fail offsetToNext:" << offsetToNext << " rdh 0x" << (void*)rdh << " bufsize:" << mDataBufferSize << " payload start: 0x" << (void*)&mHBFPayload[0] << " mHBFoffset32 " << std::dec << mHBFoffset32;
        checkNoWarn();
      }
      if (mOptions[TRDVerboseBit]) {
        LOG(info) << "rdh bounds fail offsetToNext:" << offsetToNext << " rdh 0x" << (void*)rdh << " bufsize:" << mDataBufferSize << " payload start: 0x" << (void*)&mHBFPayload[0] << " mHBFoffset32 " << std::dec << mHBFoffset32;
      }
      return false;
    }
  }

  // at this point the entire HBF data payload is sitting in mHBFPayload and the total data count is mTotalHBFPayLoad
  while ((mHBFoffset32 < ((mTotalHBFPayLoad) / 4))) {
    if (mOptions[TRDVerboseBit]) {
      LOG(info) << "Looping over cruheaders in HBF, loop count " << counthalfcru << " current offset is" << mHBFoffset32 << " total payload is " << mTotalHBFPayLoad / 4 << "  raw :" << mTotalHBFPayLoad;
    }
    int halfcruprocess = processHalfCRU(mHBFoffset32, counthalfcru, currentsaveddatacount);
    if (mOptions[TRDVerboseBit]) {
      switch (halfcruprocess) {
        case -1:
          LOG(info) << "ignored rdh event ";
          break;
        case 0:
          LOG(warn) << "figure out what now";
          break;
        case 1:
          LOG(info) << "all good parsing half cru";
          LOG(info) << " mHBFoffset32:" << mHBFoffset32 << " and mTotalHBFPayload/4 : " << mTotalHBFPayLoad / 4;
          break;
        case 2:
          LOG(info) << "all good parsing half cru was blank double 0xe event";
          break;
        default:
          LOG(info) << "unhandled return from processhalfcru of : " << halfcruprocess;
          break;
      }
    }
    if (halfcruprocess == -2) {
      // dump rest of this rdh payload, something screwed up.
      mHBFoffset32 = mTotalHBFPayLoad / 4;
    }
    counthalfcru++;
  } // loop of halfcru's while there is still data in the heart beat frame.
  if (totaldataread > 0) {
    mDatareadfromhbf = totaldataread;
  }
  return true; // totaldataread;
}

int CruRawReader::checkTrackletHCHeader()
{
  if (!mOptions[TRDIgnoreTrackletHCHeaderBit]) { // we take half chamber header as authoritive
    return 0;
  }

  // check if HCHeader is ok.
  if (!sanityCheckTrackletHCHeader(mTrackletHCHeader, mOptions[TRDVerboseErrorsBit])) {
    return 1;
  }

  return 0;
}

int CruRawReader::checkDigitHCHeader()
{
  // index 0 is rdh data, index 1 is ori calculated data
  auto currentsector = mDigitHCHeader.supermodule;
  auto currentlayer = mDigitHCHeader.layer;
  auto currentstack = mDigitHCHeader.stack;
  auto currentside = mDigitHCHeader.side;
  // check rdh info vs half chamber header
  if (!mOptions[TRDIgnoreDigitHCHeaderBit]) { // we take half chamber header as authoritive
    // can use digithcheader for cross checking the sector/stack/layer
    /*if (currentstack != mStack[0]) { // || currentstack != mStack[1]) {
      //stack mismatch
      //count these
      incrementErrors(TRDParsingDigitStackMismatch, mFEEID.supermodule, mHalfChamberSide[0], mStack[0], mLayer[0]);
    }
    if (currentlayer != mLayer[0]) { //|| currentlayer != mLayer[1]) {
      //layer mismatch
      //count these
      incrementErrors(TRDParsingDigitLayerMismatch, mFEEID.supermodule, mHalfChamberSide[0], mStack[0], mLayer[0]);
    }*/
    // taken out as stack and layer changes with in a cru endoint, or halfcruheder of 15 links.
    // need to rather check that the subsequent ones are correct relative to the previous ones, but that can come in other rawreader.
    // sector does not change.
    if (currentsector != mSector[0]) { //} || currentsector != mSector[1]) {
      // sector mismatch, mDetector comes in from a construction via the feeid and ori.
      // count these
      std::stringstream message;
      message << " Sector mismatch in Digit HCHeader : currentsector : " << currentsector << " mSector:" << mSector[0];
      incrementErrors(TRDParsingDigitSectorMismatch, message, 1, true, mFEEID.supermodule, mHalfChamberSide[0], mStack[0], mLayer[0]);
    }
    mSector[2] = currentsector; // from hc header treating it as authoritative
    mLayer[2] = currentlayer;
    mStack[2] = currentstack;
    mDetector[2] = mLayer[2] + mStack[2] * constants::NLAYER + mSector[2] * constants::NLAYER * constants::NSTACK;
    return 2;
  } else { // ignore the halfcahmber headers contents so use the rdh
    // take mDetector, layer and stack from the rdh/cru, we have those already assigned on entry to here
    return 0; // the index in mSector/mStack etc.
  }
}

int CruRawReader::parseDigitHCHeader()
{
  // mHBFoffset is the current offset into the current buffer,
  //
  uint32_t dhcheader = mHBFPayload[mHBFoffset32++];
  std::array<uint32_t, 4> headers{0};
  if (mOptions[TRDByteSwapBit]) {
    // byte swap if needed.
    o2::trd::HelperMethods::swapByteOrder(dhcheader);
  }
  mDigitHCHeader.word = dhcheader;
  if (mDigitHCHeader.major == 0 && mDigitHCHeader.minor == 0 && mDigitHCHeader.numberHCW == 0) {
    // hack this data into something resembling usable.
    mDigitHCHeader.major = mHalfChamberMajor;
    mDigitHCHeader.minor = 42;
    mDigitHCHeader.numberHCW = mHalfChamberWords;
    if (mHalfChamberWords == 0 || mHalfChamberMajor == 0) {
      if (mMaxWarnPrinted > 0) {
        LOG(warn) << "halfchamber header is corrupted and you have only set the halfchamber command line option to zero, hex dump of data and revisit what it should be.";
        checkNoWarn();
      }
      // already in histograms
    }
  }

  int additionalHeaderWords = mDigitHCHeader.numberHCW;
  if (additionalHeaderWords >= 3) {
    std::stringstream message;
    message << "Error parsing DigitHCHeader, too many additional words count=" << additionalHeaderWords << " header:" << std::hex << mDigitHCHeader.word;
    incrementErrors(TRDParsingDigitHeaderCountGT3, message, 1, true, mFEEID.supermodule, mHalfChamberSide[0], mStack[0], mLayer[0]);
    return -1;
  }
  std::bitset<3> headersfound;

  for (int headerwordcount = 0; headerwordcount < additionalHeaderWords; ++headerwordcount) {
    headers[headerwordcount] = mHBFPayload[mHBFoffset32++];
    if (mOptions[TRDByteSwapBit]) {
      // byte swap if needed.
      o2::trd::HelperMethods::swapByteOrder(headers[headerwordcount]);
    }
    switch (getDigitHCHeaderWordType(headers[headerwordcount])) {
      case 1: // header header1;
        if (headersfound.test(0)) {
          // we have a problem, we already have a Digit HC Header1, we are hereby lost, so as Monty Python said, .... run away , run away, run away.
          std::stringstream message;
          message << "We have a >1 Digit HC Header 1  : " << std::hex << " raw: 0x" << headers[headerwordcount];
          incrementErrors(TRDParsingDigitHCHeader1, message, 1, true);
        }
        mDigitHCHeader1.word = headers[headerwordcount];
        headersfound.set(0);
        if (mDigitHCHeader1.res != 0x1) {
          std::stringstream message;
          message << "Digit HC Header 1 reserved : 0x" << std::hex << mDigitHCHeader1.res << " raw: 0x" << mDigitHCHeader1.word;
          incrementErrors(TRDParsingDigitHeaderWrong1, message, 1, true);
        }
        if ((mDigitHCHeader1.numtimebins > o2::trd::constants::TIMEBINS) || (mDigitHCHeader1.numtimebins < 3)) {
          // numtimebins is unsigned so no need to check for <1
          return -1;
        }
        mTimeBins = mDigitHCHeader1.numtimebins;
        if (mTimeBins < 1 && mTimeBins > o2::trd::constants::TIMEBINS) {
          // sanity check on the hcheader settings
          mTimeBins = o2::trd::constants::TIMEBINS;
          if (mMaxWarnPrinted > 0) {
            LOG(warn) << "Time bins in Digit HC Header 1 is " << mDigitHCHeader1.numtimebins << " this is absurd";
            checkNoWarn();
          }
        }
        break;
      case 2: // header header2;
        if (headersfound.test(1)) {
          // we have a problem, we already have a Digit HC Header2, we are hereby lost, so as Monty Python said, .... run away , run away, run away.
          std::stringstream message;
          message << "We have a >1 Digit HC Header 2  : " << std::hex << " raw: 0x" << headers[headerwordcount];
          incrementErrors(TRDParsingDigitHCHeader2, message, 1, true);
        }
        mDigitHCHeader2.word = headers[headerwordcount];
        headersfound.set(1);
        if (mDigitHCHeader2.res != 0b110001) {
          std::stringstream message;
          message << "Digit HC Header 2 reserved : " << std::hex << mDigitHCHeader2.res << " raw: 0x" << mDigitHCHeader2.word;
          incrementErrors(TRDParsingDigitHeaderWrong2, message, 1, true);
        }
        break;
      case 3: // header header3;
        if (headersfound.test(2)) {
          // we have a problem, we already have a Digit HC Header2, we are hereby lost, so as Monty Python said, .... run away , run away, run away.
          std::stringstream message;
          message << "We have a >1 Digit HC Header 2  : " << std::hex << " raw: 0x" << headers[headerwordcount];
          incrementErrors(TRDParsingDigitHCHeader3, message, 1, true);
        }
        mDigitHCHeader3.word = headers[headerwordcount];
        headersfound.set(2);
        if (mDigitHCHeader3.res != 0b110101) {
          std::stringstream message;
          message << "Digit HC Header 3 reserved : " << std::hex << mDigitHCHeader3.res << " raw: 0x" << mDigitHCHeader3.word;
          incrementErrors(TRDParsingDigitHeaderWrong3, message, 1, true);
        }
        if (mPreviousDigitHCHeadersvnver != 0xffffffff && mPreviousDigitHCHeadersvnrver != 0xffffffff) {
          if ((mDigitHCHeader3.svnver != mPreviousDigitHCHeadersvnver) && (mDigitHCHeader3.svnrver != mPreviousDigitHCHeadersvnrver)) {
            if (mMaxWarnPrinted > 0) {
              checkNoWarn();
            }
            std::stringstream message;
            message << "Digit HC Header 3 svn ver : " << std::hex << mDigitHCHeader3.svnver << " svn release ver : 0x" << mDigitHCHeader3.svnrver;
            incrementErrors(TRDParsingDigitHCHeaderSVNMismatch, message, 1, true);
            return -1;
          } else {
            // this is the first time seeing a DigitHCHeader3
            mPreviousDigitHCHeadersvnver = mDigitHCHeader3.svnver;
            mPreviousDigitHCHeadersvnrver = mDigitHCHeader3.svnrver;
          }
        }
        break;
      default:
        // LOG(warn) << "Error parsing DigitHCHeader at word:" << headerwordcount << " looking at 0x:" << std::hex << mHBFPayload[mHBFoffset32 - 1];
        std::stringstream message;
        message << " unknown error in switch staement for Digit HC Header";
        incrementErrors(TRDParsingDigitHeaderWrong4, message, 2, false);
    }
  }
  if (mOptions[TRDVerboseBit]) {
    printDigitHCHeader(mDigitHCHeader, &headers[0]);
  }

  return 1;
}

void CruRawReader::updateLinkErrorGraphs(int currentlinkindex, int supermodule_half, int stack_layer)
{
  mEventRecords.incLinkErrorFlags(mDetector[0], mHalfChamberSide[0], stack_layer, mCurrentHalfCRULinkErrorFlags[currentlinkindex]);
  if (mCurrentHalfCRULinkLengths[currentlinkindex] == 0) {
    mEventRecords.incLinkNoData(mDetector[0], mHalfChamberSide[0], stack_layer);
  }
}

int CruRawReader::processHalfCRU(uint32_t cruhbfstartoffset, int numberOfPreviousCRU, unsigned int maxdatawrittentobuffer)
{
  // It will clean this code up *alot*
  //  process a halfcru
  uint32_t currentlinkindex = 0;
  uint32_t currentlinkoffset = 0;
  uint32_t currentlinksize = 0;
  uint32_t currentlinksize32 = 0;
  uint32_t linksizeAccum32 = 0;
  uint32_t sumtrackletwords = 0;
  uint32_t sumdigitwords = 0;
  uint32_t sumlinklengths = 0;
  mDigitWordsRead = 0;
  mDigitWordsRejected = 0;
  mTrackletWordsRead = 0;
  mTrackletWordsRejected = 0;
  uint32_t cruwordsread = 9;
  // reject halfcru if it starts with padding words.
  // TODO put maxdatawrittentobuffer in a qc plot
  mHalfCRUStartOffset = cruhbfstartoffset;
  if (mHBFPayload.size() < cruhbfstartoffset || cruhbfstartoffset > maxdatawrittentobuffer) {
    mHBFoffset32++;
    std::stringstream message;
    message << "Error parsing HalfCRUHeader, HBFPayload size = " << mHBFPayload.size() << " payload offset:" << cruhbfstartoffset << " max data written to buffer : " << maxdatawrittentobuffer;
    incrementErrors(TRDProcessingBadPayloadOrOffset, message, 1, true);
    return -2;
  }
  // this should only hit that instance where the cru payload is a "blank event" of o2::trd::constants::CRUPADDING32
  if (mHBFPayload[cruhbfstartoffset] == o2::trd::constants::CRUPADDING32) { //} && mHBFPayload[cruhbfstartoffset + 1] == o2::trd::constants::CRUPADDING32) {
    if (mOptions[TRDVerboseBit]) {
      LOG(info) << "blank rdh payload data at " << cruhbfstartoffset << ": 0x " << std::hex << mHBFPayload[cruhbfstartoffset] << " and 0x" << mHBFPayload[cruhbfstartoffset + 1];
    }
    mHBFoffset32++; // increment past the word of the if statement and then any others that might be here.
    int loopcount = 0;
    while (mHBFPayload[mHBFoffset32] == o2::trd::constants::CRUPADDING32 && loopcount < 8) { // can only ever be an entire 256 bit word hence a limit of 8 here.
      mHBFoffset32++;
      mWordsAccepted++;
      loopcount++;
    }
    return 2;
  }
  if (mTotalHBFPayLoad == 0) {
    // empty payload
    return -1;
  }
  auto crustart = std::chrono::high_resolution_clock::now();
  // well then read the halfcruheader.
  memcpy((char*)&mCurrentHalfCRUHeader, (char*)(&mHBFPayload[cruhbfstartoffset]), sizeof(mCurrentHalfCRUHeader));
  mHBFoffset32 += sizeof(mCurrentHalfCRUHeader) / sizeof(mHBFoffset32); // advance past the header.
  if (mOptions[TRDVerboseWordBit]) {
    // output the cru half chamber header : raw/parsed
    //
    dumpHalfCRUHeader(mCurrentHalfCRUHeader);
  }

  o2::trd::getHalfCRULinkDataSizes(mCurrentHalfCRUHeader, mCurrentHalfCRULinkLengths);
  o2::trd::getHalfCRULinkErrorFlags(mCurrentHalfCRUHeader, mCurrentHalfCRULinkErrorFlags);
  mTotalHalfCRUDataLength256 = std::accumulate(mCurrentHalfCRULinkLengths.begin(),
                                               mCurrentHalfCRULinkLengths.end(),
                                               decltype(mCurrentHalfCRULinkLengths)::value_type(0));
  mTotalHalfCRUDataLength = mTotalHalfCRUDataLength256 * 32;      // convert to bytes.
  int mTotalHalfCRUDataLength32 = mTotalHalfCRUDataLength256 * 8; // convert to bytes.

  // in the interests of descerning real corrupt halfcruheaders from the sometimes garbage at the end of a half cru
  // if the first word is clearly garbage assume garbage and not a corrupt halfcruheader.
  if (numberOfPreviousCRU > 0) {
    if (mCurrentHalfCRUHeader.EndPoint != mPreviousHalfCRUHeader.EndPoint) {
      std::stringstream message;
      message << numberOfPreviousCRU << " current endpont : " << mCurrentHalfCRUHeader.EndPoint << " previous end point : " << mPreviousHalfCRUHeader.EndPoint;
      incrementErrors(TRDParsingHalfCRUCorrupt, message, 1, true);
      mWordsRejected += mTotalHalfCRUDataLength32;
      return -2;
    }
    // event type can change wit in a
    if (mCurrentHalfCRUHeader.StopBit != mPreviousHalfCRUHeader.StopBit) {
      std::stringstream message;
      message << numberOfPreviousCRU << " current stopbit: " << mCurrentHalfCRUHeader.StopBit << " previous stopbit: " << mPreviousHalfCRUHeader.StopBit;
      incrementErrors(TRDParsingHalfCRUCorrupt, message, 1, true);
      mWordsRejected += mTotalHalfCRUDataLength32;
      return -2;
    }
  }
  memcpy((char*)&mPreviousHalfCRUHeader, (char*)(&mHBFPayload[cruhbfstartoffset]), sizeof(mCurrentHalfCRUHeader));
  // can this half cru length fit into the available space of the rdh accumulated payload
  if (mTotalHalfCRUDataLength32 > mTotalHBFPayLoad - mHBFoffset32) {
    std::stringstream message;
    message << "Next HalfCRU header says it contains more data than in the rdh payloads! " << mTotalHalfCRUDataLength32 << " < " << mTotalHBFPayLoad << "-" << mHBFoffset32 << " sector:side:endpoint: " << (unsigned int)mFEEID.supermodule << ":" << (unsigned int)mFEEID.side << ":" << (unsigned int)mFEEID.endpoint;
    incrementErrors(TRDParsingHalfCRUSumLength, message, 1, true);
    mWordsRejected += mTotalHalfCRUDataLength32;
    mHBFoffset32 += mTotalHalfCRUDataLength32; // go to the end of this halfcruheader and payload.

    return -2;
  }
  if (!halfCRUHeaderSanityCheck(mCurrentHalfCRUHeader, mCurrentHalfCRULinkLengths, mCurrentHalfCRULinkErrorFlags)) {
    std::stringstream message;
    message << "HalfCRU header failed sanity check for FEEID with  sector:side:endpoint: " << (unsigned int)mFEEID.supermodule << ":" << (unsigned int)mFEEID.side << ":" << (unsigned int)mFEEID.endpoint;
    // let incrementErrors catch the undefined values of sector side stack and layer as if not set it will go so zero in the method, however if set, it means this is the second half cru header, and we have the values from the last one we read which
    // *SHOULD* be the same as this halfcruheader.
    incrementErrors(TRDParsingHalfCRUCorrupt, message, 1, true);
    mWordsRejected += mTotalHalfCRUDataLength32;
    mHBFoffset32 += mTotalHalfCRUDataLength32; // go to the end of this halfcruheader and payload.
    return -2;
  }

  // get eventrecord for event we are looking at
  mIR.bc = mCurrentHalfCRUHeader.BunchCrossing; // correct mIR to have the physics trigger bunchcrossing *NOT* the heartbeat trigger bunch crossing.
  // shift accordingly
  mIR.bc = mCurrentHalfCRUHeader.BunchCrossing - o2::ctp::TriggerOffsetsParam::Instance().LM_L0;
  if (mIR.bc < 0) {
    // dump to the end of this cruhalfchamberheader
    // data to dump is mTotalHalfCruDataLength32
    mHBFoffset32 += mTotalHalfCRUDataLength32;   // go to the end of this halfcruheader and payload.
    mWordsRejected += mTotalHalfCRUDataLength32; // add the rejected data to the accounting;
    std::stringstream message;
    message << "Bunchcrossing from previous orbit, is negative after shift and data is being rejected";
    incrementErrors(TRDParsingHalfCRUBadBC, message, 1, true);
    return 1; // nothing particularly wrong with the data, we just dont want it, as a trigger problem
  }
  if (mOptions[TRDOnlyCalibrationTriggerBit] && mCurrentHalfCRUHeader.EventType == o2::trd::constants::ETYPEPHYSICSTRIGGER) {
    mHBFoffset32 += mTotalHalfCRUDataLength32; // go to the end of this halfcruheader and payload.
    return 1;                                  // we dont want the physics triggers to polute the logging, possibly other reasons
  }
  InteractionRecord trdir(mIR);
  mCurrentEvent = &mEventRecords.getEventRecord(trdir);
  // check for cru errors :
  int linkerrorcounter = 0;
  for (auto& linkerror : mCurrentHalfCRULinkErrorFlags) {
    if (linkerror != 0) {
      if (mOptions[TRDVerboseBit]) {
        LOG(info) << "E link error FEEID:" << mFEEID.word << " CRUID:" << mCRUID << " Endpoint:" << mCRUEndpoint
                  << " on linkcount:" << linkerrorcounter++ << " errorval:0x" << std::hex << linkerror;
      }
    }
  }

  std::array<uint32_t, 1024>::iterator currentlinkstart = mHBFPayload.begin() + cruhbfstartoffset;

  if (mOptions[TRDVerboseHalfCruBit]) {
    OutputHalfCruRawData(cruhbfstartoffset);
  }
  std::array<uint32_t, 1024>::iterator linkstart, linkend;
  uint32_t dataoffsetstart32 = sizeof(mCurrentHalfCRUHeader) / sizeof(dataoffsetstart32) + cruhbfstartoffset; // in uint32
  // CHECK 1 does rdh endpoint match cru header end point.
  if (mCRUEndpoint != mCurrentHalfCRUHeader.EndPoint) {
    if (mMaxWarnPrinted > 0) {
      LOG(warn) << " Endpoint mismatch : CRU Half chamber header endpoint = " << mCurrentHalfCRUHeader.EndPoint << " rdh end point = " << mCRUEndpoint;
      checkNoWarn();
    }
    // disaster dump the rest of this hbf
    return -2;
  }

  // verify cru header vs rdh header
  // FEEID has supermodule/layer/stack/side in it.
  // CRU has

  linkstart = mHBFPayload.begin() + dataoffsetstart32;
  linkend = mHBFPayload.begin() + dataoffsetstart32;
  // loop over links
  for (currentlinkindex = 0; currentlinkindex < constants::NLINKSPERHALFCRU; currentlinkindex++) {
    auto linktimerstart = std::chrono::high_resolution_clock::now(); // measure total processing time
    mSector[0] = mFEEID.supermodule;
    mEndPoint[0] = mFEEID.endpoint;
    mSide[0] = mFEEID.side; // side of detector A/C
    uint32_t hbfoffsetatstartoflink = mHBFoffset32;
    // stack layer and side map to ori
    uint32_t oriindex = currentlinkindex + constants::NLINKSPERHALFCRU * mEndPoint[0]; // endpoint denotes the pci side, upper or lower for the pair of 15 fibres.
    FeeParam::unpackORI(oriindex, mSide[0], mStack[1], mLayer[1], mHalfChamberSide[1]);
    // sadly not all the data is redundant, probably a good thing, so stack and layer and halfchamber side is derived from the ori.
    mLayer[0] = mLayer[1];
    mStack[0] = mStack[1];
    mHalfChamberSide[0] = mHalfChamberSide[1];
    mSector[1] = oriindex / 30;
    mSide[1] = mSide[0];
    mDetector[0] = mStack[0] * constants::NLAYER + mLayer[0] + mSector[0] * constants::NLAYER * constants::NSTACK;
    mDetector[1] = mStack[1] * constants::NLAYER + mLayer[1] + mSector[1] * constants::NLAYER * constants::NSTACK;
    int supermodule_half = mSector[0] * 2 + mHalfChamberSide[0]; // will just go with the rdh one here its only for the hack graphing purposes.
    int stack_layer;
    stack_layer = mStack[0] * constants::NLAYER + mLayer[0]; // similarly this is also only for graphing so just use the rdh ones for now.
    updateLinkErrorGraphs(currentlinkindex, supermodule_half, stack_layer);

    mEventRecords.incLinkErrorFlags(mFEEID.supermodule, mHalfChamberSide[0], stack_layer, mCurrentHalfCRULinkErrorFlags[currentlinkindex]);
    currentlinksize = mCurrentHalfCRULinkLengths[currentlinkindex];
    // first parameter is the base of the link in the, so either the first or second part of the cru hence the *15
    mCurrentEvent->setDataPerLink((mFEEID.supermodule * 2 + mHalfChamberSide[0]) * 30 + currentlinkindex, currentlinksize);
    currentlinksize32 = currentlinksize * 8; // x8 to go from 256 bits to 32 bit;
    linkstart = mHBFPayload.begin() + dataoffsetstart32 + linksizeAccum32;
    linkend = linkstart + currentlinksize32;
    if (currentlinksize == 0) {
      mEventRecords.incLinkNoData(mDetector[0], mHalfChamberSide[0], stack_layer);
    }
    uint32_t linkzsum = 0;
    int dioffset = dataoffsetstart32 + linksizeAccum32;
    if (dioffset % 8 != 0) {
      if (mMaxErrsPrinted > 0) {
        LOG(warn) << " we are not 256 bit aligned ... this should never happen";
        checkNoErr();
      }
    }
    if (mHBFoffset32 != std::distance(mHBFPayload.begin(), linkstart)) {
      mHBFoffset32 = std::distance(mHBFPayload.begin(), linkstart);
    }
    if (mOptions[TRDVerboseBit]) {
      LOG(info) << "Cru link :" << currentlinkindex << " raw dump before processing begin linkstart:" << std::hex << linkstart << " to " << linkend << " mHBFoffset32=" << std::dec << mHBFoffset32 << " and distance from start is : " << std::distance(mHBFPayload.begin(), linkstart);
      for (int dumpoffset = dataoffsetstart32 + linksizeAccum32; dumpoffset < dataoffsetstart32 + linksizeAccum32 + currentlinksize32; dumpoffset += 8) {
        LOGP(info, "0x{0:06x} :: {1:08x} {2:08x}  {3:08x} {4:08x} {5:08x} {6:08x} {7:08x} {8:08x} ", dumpoffset, mHBFPayload[dumpoffset], mHBFPayload[dumpoffset + 1], mHBFPayload[dumpoffset + 2], mHBFPayload[dumpoffset + 3], mHBFPayload[dumpoffset + 4], mHBFPayload[dumpoffset + 5], mHBFPayload[dumpoffset + 6], mHBFPayload[dumpoffset + 7]);
      }
    }
    linksizeAccum32 += currentlinksize32;
    if (linkstart != linkend) { // if link is not empty
      bool cleardigits = false; // linkstart and linkend already have the multiple cruheaderoffsets built in
      auto trackletparsingstart = std::chrono::high_resolution_clock::now();
      if (mOptions[TRDVerboseBit]) {
        LOG(info) << "*** Tracklet Parser : starting at " << std::hex << linkstart << " at hbfoffset: " << std::dec << mHBFoffset32 << " linkhbf start pos:" << hbfoffsetatstartoflink;
      }
      if (std::distance(linkstart, linkend) > mCurrentHalfCRULinkLengths[currentlinkindex] * 8) { //*8 for lengths are stored in units of cru words(256bit) and iterators are 32bit.
        std::stringstream message;
        message << "linkstart - linkend for  LINK # " << currentlinkindex << " an FEEID:" << std::hex << mFEEID.word << " det:" << std::dec << mDetector[1] << " is > the lenght stored in the cruhalfchamber header : " << mCurrentHalfCRULinkLengths[currentlinkindex];
        incrementErrors(TRDParsingBadLinkstartend, message, 1, true, mFEEID.supermodule, mFEEID.side, mStack[0], mLayer[0]);
        // in the immortal words of ... run away ... run away ... run away
        return -2; // dump this buffer.
      }
      // for now we are using 0 i.e. from rdh FIXME figure out which is authoritative between rdh and ori tracklethcheader if we have it enabled.
      mTrackletWordsRead = mTrackletsParser.Parse(&mHBFPayload, linkstart, linkend, mFEEID, mHalfChamberSide[0], mDetector[0], mStack[0], mLayer[0], mCurrentEvent, &mEventRecords, mOptions, cleardigits, mTrackletHCHeaderState); // this will read up to the tracklet end marker.
      if (mTrackletsParser.dumpLink()) {
        // dump the link that cause the error
        //  the call to dumpLink resets the boolean to false;
        OutputLinkRawData(currentlinkindex);
      }
      if (mTrackletWordsRead == -1) {
        // something went wrong bailout of here.
        std::stringstream message;
        message << "TrackletParser returned -1 for  LINK # " << currentlinkindex << " an FEEID:" << std::hex << mFEEID.word << " det:" << std::dec << mDetector[1] << " is > the lenght stored in the cruhalfchamber header : " << mCurrentHalfCRULinkLengths[currentlinkindex];
        incrementErrors(TRDParsingTrackletsReturnedMinusOne, message, 1, true, mFEEID.supermodule, mFEEID.side, mStack[0], mLayer[0]);
        return -2;
      }
      mTrackletWordsRejected = mTrackletsParser.getDataWordsDumped();
      std::chrono::duration<double, std::micro> trackletparsingtime = std::chrono::high_resolution_clock::now() - trackletparsingstart;
      mCurrentEvent->incTrackletTime((double)std::chrono::duration_cast<std::chrono::microseconds>(trackletparsingtime).count());
      if (mOptions[TRDVerboseBit]) {
        LOG(info) << "trackletwordsread:" << mTrackletWordsRead << " trackletwordsrejected:" << mTrackletWordsRejected << "  mem copy with offset of : " << cruhbfstartoffset << " parsing with linkstart: " << linkstart << " ending at : " << linkend;
      }
      linkstart += mTrackletWordsRead + mTrackletWordsRejected;
      // now we have a tracklethcheader and a digithcheader.

      mHBFoffset32 += mTrackletWordsRead + mTrackletWordsRejected;
      mTotalTrackletsFound += mTrackletsParser.getTrackletsFound();
      mTotalTrackletWordsRejected += mTrackletWordsRejected;
      mTotalTrackletWordsRead += mTrackletWordsRead;
      mCurrentEvent->incWordsRead(mTrackletWordsRead);
      mCurrentEvent->incWordsRejected(mTrackletWordsRejected);
      mEventRecords.incLinkWordsRead(mFEEID.supermodule, mHalfChamberSide[0], stack_layer, mTrackletWordsRead);
      mEventRecords.incLinkWordsRejected(mFEEID.supermodule, mHalfChamberSide[0], stack_layer, mTrackletWordsRejected);
      if (mTrackletsParser.getTrackletParsingState()) {
        mHBFoffset32 += std::distance(linkstart, linkend);
        linkstart = linkend; // bail out as tracklet parsing bombed out. We are essentially lost.
      }
      if (mOptions[TRDVerboseBit]) {
        LOG(info) << "*** Tracklet Parser : trackletwordsread:" << mTrackletWordsRead << " ending " << std::hex << linkstart << " at hbfoffset: " << std::dec << mHBFoffset32;
      }

      /****************
      ** DIGITS NOW ***
      *****************/
      // Check if we have a calibration trigger ergo we do actually have digits data. check if we are now at the end of the data due to bugs, i.e. if trackletparsing read padding words.
      if (linkstart != linkend && (mCurrentHalfCRUHeader.EventType == o2::trd::constants::ETYPECALIBRATIONTRIGGER || mOptions[TRDIgnore2StageTrigger]) && (mHBFPayload[cruhbfstartoffset] != o2::trd::constants::CRUPADDING32)) { // calibration trigger and insure we dont come in here if we are on a padding word.
        if (mOptions[TRDVerboseBit]) {
          LOG(info) << "*** Digit Parsing : starting at " << std::hex << linkstart << " at hbfoffset: " << std::dec << mHBFoffset32 << " linkhbf start pos:" << hbfoffsetatstartoflink;
        }
        // linkstart advanced all the way to the end due to trackletparser parsing crupadding words (known bug or feature )
        auto hfboffsetbeforehcparse = mHBFoffset32;
        // now read the digit half chamber header
        auto hcparse = parseDigitHCHeader();
        if (hcparse != 1) {
          // disaster dump the rest of this hbf
          return -2;
        }
        mWhichData = checkDigitHCHeader();
        // move over the DigitHCHeader mHBFoffset32 has already been moved in the reading.
        if (mHBFoffset32 - hfboffsetbeforehcparse != 1 + mDigitHCHeader.numberHCW) {
          if (mMaxErrsPrinted > 0) {
            LOG(warn) << "Seems data offset is out of sync with number of HC Headers words " << mHBFoffset32 << "-" << hfboffsetbeforehcparse << "!=" << 1 << "+" << mDigitHCHeader.numberHCW;
            checkNoErr();
          }
        }
        if (hcparse == -1) {
          if (mMaxWarnPrinted > 0) {
            LOG(warn) << "Parsing Digit HCHeader returned a -1";
            checkNoWarn();
          }
        } else {
          linkstart += 1 + mDigitHCHeader.numberHCW;
        }
        mEventRecords.incMajorVersion(mDigitHCHeader.major); // 127 is max histogram goes to 256

        if (mDigitHCHeader.major == 0x47) {
          // config event so ignore for now and bail out of parsing.
          // advance data pointers to the end;
          linkstart = linkend;
          mHBFoffset32 = std::distance(mHBFPayload.begin(), linkend); // currentlinksize-mTrackletWordsRead-sizeof(digitHCHeader)/4; // advance to the end of the link
          mTotalDigitWordsRejected += std::distance(linkstart + mTrackletWordsRead + sizeof(DigitHCHeader) / 4, linkend);
          LOG(info) << "Configuration event  ";
        } else {
          mDigitWordsRead = 0;
          auto digitsparsingstart = std::chrono::high_resolution_clock::now();
          // linkstart and linkend already have the multiple cruheaderoffsets built in
          mDigitWordsRead = mDigitsParser.Parse(&mHBFPayload, linkstart, linkend, mDetector[mWhichData], mStack[mWhichData], mLayer[mWhichData], mHalfChamberSide[mWhichData], mDigitHCHeader, mTimeBins, mFEEID, currentlinkindex, mCurrentEvent, &mEventRecords, mOptions, cleardigits);
          std::chrono::duration<double, std::micro> digitsparsingtime = std::chrono::high_resolution_clock::now() - digitsparsingstart;
          mCurrentEvent->incDigitTime((double)std::chrono::duration_cast<std::chrono::microseconds>(digitsparsingtime).count());
          mDigitWordsRejected = mDigitsParser.getDumpedDataCount();
          mCurrentEvent->incWordsRead(mDigitWordsRead);
          mCurrentEvent->incWordsRejected(mDigitWordsRejected);
          mEventRecords.incLinkWordsRead(mFEEID.supermodule, mHalfChamberSide[0], stack_layer, mDigitWordsRead);
          mEventRecords.incLinkWordsRejected(mFEEID.supermodule, mHalfChamberSide[0], stack_layer, mDigitWordsRejected);
          if (mOptions[TRDVerboseBit]) {
            if (mDigitsParser.getDumpedDataCount() != 0) {
              LOG(info) << "FEEID: " << mFEEID.word << " LINK #" << oriindex << " bad datacount:" << mDigitsParser.getDataWordsParsed() << "::" << mDigitsParser.getDumpedDataCount();
            } else {
              LOG(info) << "FEEID: " << mFEEID.word << " LINK #" << oriindex << " good datacount:" << mDigitsParser.getDataWordsParsed() << "::" << mDigitsParser.getDumpedDataCount();
            }
          }
          mDigitWordsRejected = mDigitsParser.getDumpedDataCount();
          if (mDigitsParser.dumpLink()) {
            // dump the link that cause the error
            //  the call to dumpLink resets the boolean to false;
            OutputLinkRawData(currentlinkindex);
          }
          if (mOptions[TRDVerboseBit]) {
            if (mDigitsParser.getDumpedDataCount() != 0) {
              LOG(info) << "FEEID: " << mFEEID.word << " LINK #" << oriindex << " bad datacount:" << mDigitsParser.getDataWordsParsed() << "::" << mDigitsParser.getDumpedDataCount();
            } else {
              LOG(info) << "FEEID: " << mFEEID.word << " LINK #" << oriindex << " good datacount:" << mDigitsParser.getDataWordsParsed() << "::" << mDigitsParser.getDumpedDataCount();
            }
            std::stringstream message;
            message << "FEEID: " << mFEEID.word << " LINK #" << oriindex << " good datacount:" << mDigitsParser.getDataWordsParsed() << "::" << mDigitsParser.getDumpedDataCount();
            incrementErrors(TRDParsingDigitStackMismatch, message, 2, false, mFEEID.supermodule, mHalfChamberSide[0], mStack[0], mLayer[0]);
          }
          if (mDigitWordsRead + mDigitWordsRejected != std::distance(linkstart, linkend)) {
            // we have the data corruption problem of a pile of stuff at the end of a link, jump over it.
            if (mOptions[TRDFixDigitCorruptionBit]) {
              mDigitWordsRead = std::distance(linkstart, linkend);
            } else {
              std::stringstream message;
              message << "FEEID: " << mFEEID.word << " LINK #" << oriindex << " data still on link ";
              incrementErrors(TRDParsingDigitDataStillOnLink, message, 2, false, mFEEID.supermodule, mHalfChamberSide[0], mStack[0], mLayer[0]);
            }
          }
          mTotalDigitsFound += mDigitsParser.getDigitsFound();
          mHBFoffset32 += mDigitWordsRead + mDigitWordsRejected; // all 3 in 32bit units
          mTotalDigitWordsRead += mDigitWordsRead;
          mTotalDigitWordsRejected += mDigitWordsRejected;
          sumlinklengths += mCurrentHalfCRULinkLengths[currentlinkindex];
          sumtrackletwords += mTrackletWordsRead;
          sumdigitwords += mDigitWordsRead;
        }
      }
    } else {
      if (mCurrentHalfCRUHeader.EventType == o2::trd::constants::ETYPEPHYSICSTRIGGER) {
        mEventRecords.incMajorVersion(128); // 127 is max histogram goes to 256
      }
    }
  } // for loop over link index.
  // we have read in all the digits and tracklets for this event.
  // digits and tracklets are sitting inside the parsing classes.
  // extract the vectors and copy them to tracklets and digits here, building the indexing(triggerrecords)
  // as this is for a single cru half chamber header all the tracklets and digits are for the same trigger defined by the bc and orbit in the rdh which we hold in mIR

  int lasttrigger = 0, lastdigit = 0, lasttracklet = 0;
  std::chrono::duration<double, std::milli> cruparsingtime = std::chrono::high_resolution_clock::now() - crustart;
  mCurrentEvent->incTime(cruparsingtime.count());

  // if we get here all is ok.
  return 1;
}

bool CruRawReader::buildCRUPayLoad()
{
  // copy data for the current half cru, and when we eventually get to the end of the payload return 1
  // to say we are done.
  int cruid = 0;
  int additionalBytes = -1;
  int crudatasize = -1;
  return true;
}

void CruRawReader::resetCounters()
{
  // mStatCountersPerEvent.mLinkErrorFlag.fill(0);
  mEventCounter = 0;
  mFatalCounter = 0;
  mErrorCounter = 0;
}

void CruRawReader::checkSummary()
{
  char chname[2] = {'a', 'b'};

  LOG(info) << "--- SUMMARY COUNTERS: " << mEventCounter << " events "
            << " | " << mFatalCounter << " decode fatals "
            << " | " << mErrorCounter << " decode errors ";
}

bool CruRawReader::run()
{
  uint32_t dowhilecount = 0;
  uint32_t totaldataread = 0;
  rewind();
  mTotalDigitWordsRead = 0;
  mTotalDigitWordsRejected = 0;
  mTotalTrackletWordsRead = 0;
  mTotalTrackletWordsRejected = 0;
  do {
    mDatareadfromhbf = 0;
    auto goodprocessing = processHBFs(totaldataread);
    totaldataread += mDatareadfromhbf;
    if (!goodprocessing) {
      // processHBFs returned false, get out of here ...
      if (mMaxWarnPrinted > 0) {
        LOG(warn) << "Error processing heart beat frame ... dumping data, heart beat frame rejected";
        checkNoWarn();
      }
      break;
    }
    if (totaldataread == 0) {
      if (mMaxWarnPrinted > 0) {
        LOG(warn) << "EEE  we read zero data but bailing out of here for now.";
        checkNoWarn();
      }
      break;
    }
  } while (((char*)mDataPointer - mDataBuffer) < mDataBufferSize);

  return false;
};

void CruRawReader::getParsedObjects(std::vector<Tracklet64>& tracklets, std::vector<Digit>& digits, std::vector<TriggerRecord>& triggers)
{
  int digitcountsum = 0;
  int trackletcountsum = 0;
  mEventRecords.unpackData(triggers, tracklets, digits);
}

void CruRawReader::getParsedObjectsandClear(std::vector<Tracklet64>& tracklets, std::vector<Digit>& digits, std::vector<TriggerRecord>& triggers)
{
  getParsedObjects(tracklets, digits, triggers);
  clearall();
}

// write the output data directly to the given DataAllocator from the datareader task.
void CruRawReader::buildDPLOutputs(o2::framework::ProcessingContext& pc)
{
  mEventRecords.sendData(pc, mOptions[TRDGenerateStats]);
  clearall(); // having now written the messages clear for next.
}

void CruRawReader::checkNoWarn()
{
  if (!mOptions[TRDVerboseBit] && --mMaxWarnPrinted == 0) {
    LOG(alarm) << "Warnings limit reached, the following ones will be suppressed";
  }
}

void CruRawReader::checkNoErr()
{
  if (!mOptions[TRDVerboseBit] && --mMaxErrsPrinted == 0) {
    LOG(error) << "Errors limit reached, the following ones will be suppressed";
  }
}

} // namespace o2::trd
