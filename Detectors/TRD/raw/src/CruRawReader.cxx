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
/// @brief  TRD raw data translator

#include "DetectorsRaw/RDHUtils.h"
#include "Headers/RDHAny.h"
#include "TRDRaw/CruRawReader.h"
#include "DataFormatsTRD/RawData.h"
#include "DataFormatsTRD/Tracklet64.h"
#include "TRDRaw/DigitsParser.h"
#include "TRDRaw/TrackletsParser.h"
#include "DataFormatsTRD/Constants.h"

#include <cstring>
#include <string>
#include <vector>
#include <array>
#include <iostream>
#include <numeric>

namespace o2
{
namespace trd
{

uint32_t CruRawReader::processHBFs(int datasizealreadyread, bool verbose)
{

  LOG(info) << "PROCESS HBF starting at " << (void*)mDataPointer;

  mDataRDH = reinterpret_cast<const o2::header::RDHAny*>(mDataPointer);
  mOpenRDH = reinterpret_cast<o2::header::RDHAny*>((char*)mDataPointer);
  auto rdh = mDataRDH;
  uint64_t totaldataread = 0;
  /* loop until RDH stop header */
  mState = CRUStateHalfCRUHeader;
  while (!o2::raw::RDHUtils::getStop(rdh)) { // carry on till the end of the event.
    LOG(info) << "--- RDH open/continue detected";

    auto headerSize = o2::raw::RDHUtils::getHeaderSize(rdh);
    auto memorySize = o2::raw::RDHUtils::getMemorySize(rdh);
    auto offsetToNext = o2::raw::RDHUtils::getOffsetToNext(rdh);
    auto dataPayload = memorySize - headerSize;
    mFEEID = o2::raw::RDHUtils::getFEEID(rdh);            //TODO change this and just carry around the curreht RDH
    mCRUEndpoint = o2::raw::RDHUtils::getEndPointID(rdh); // the upper or lower half of the currently parsed cru 0-14 or 15-29
    mCRUID = o2::raw::RDHUtils::getCRUID(rdh);
    int packetCount = o2::raw::RDHUtils::getPacketCounter(rdh);
   // if(mVerbose)LOG(info) << "FEEID : " << mFEEID << " Packet: " << packetCount << " sizes : header" << headerSize << " memorysize:" << memorySize << " offsettonext:" << offsetToNext << " datapayload:" << dataPayload;
    mDataPointer = (uint32_t*)((char*)rdh + headerSize);
    mDataEndPointer = (const uint32_t*)((char*)rdh + offsetToNext);
    while ((void*)mDataPointer < (void*)mDataEndPointer) { // loop to handle the case where a halfcru ends/begins within the rdh data block
      mEventCounter++;
      if (memorySize == 96 && ((uint32_t*)mDataPointer)[0] == 0xeeeeeeee) {
        //header with only padding word (blank event), ignore and move on.
  //      if(mVerbose)LOG(info) << "MFEEID : " << std::hex << mFEEID << " rdh is at : " << (void*)rdh;
        mDataPointer = mDataEndPointer; //dump the rest
      } else {
        if (processHalfCRU()) {
            //TODO this should be done external to the getStop loop as getStop loop (line 46) will end up with a singular event buffer.
            //It will clean this code up *alot*
          // process a halfcru
          // or continue with the remainder of an rdh o2 payload if we got to the end of cru
          // or continue with a new rdh payload if we are not finished with the last cru payload.
          // TODO the above 2 lines are not possible.
          std::array<uint32_t, 1024>::iterator currentlinkstart = mCRUPayLoad.begin();
          std::array<uint32_t, 1024>::iterator linkstart, linkend;
          uint32_t currentlinkindex = 0;
          uint32_t currentlinkoffset = 0;
          uint32_t currentlinksize = 0;
          uint32_t currentlinksize32 = 0;
          uint32_t linksizeAccum32 = 0;
          linkstart = mCRUPayLoad.begin();
          linkend = mCRUPayLoad.begin();
          //loop over links
          for (currentlinkindex = 0; currentlinkindex < 15; currentlinkindex++) {
   //         if(mVerbose)LOG(info) << "******* LINK # " << currentlinkindex;
            currentlinksize = mCurrentHalfCRULinkLengths[currentlinkindex];
            currentlinksize32 = currentlinksize * 8; //x8 to go from 256 bits to 32 bit units;
   //         if(mVerbose)LOG(info) << " this link has size :" << currentlinksize << " in words its : " << currentlinksize32;
            linkstart = mCRUPayLoad.begin() + linksizeAccum32;
            linkend = linkstart + currentlinksize;
           // LOG(info) << "Parsing from " << std::distance(mCRUPayLoad.begin(), linkstart) << " --> " << std::distance(mCRUPayLoad.begin(), linkend);
          //  LOG(info) << " mState is : " << mState;
            linksizeAccum32 += currentlinksize32;
            int currentdetector = 1; // TODO fix this based on the above data.

            if (mState == CRUStateHalfCRUHeader) {
              //ergo mCRUPayLoad holds the whole links payload, so parse it.
              // tracklet first then digit ??
              // tracklets end with tracklet end marker(0x10001000 0x10001000), digits end with digit endmarker (0x0 0x0)
              if (linkstart != linkend) { // if link is not empty
              //  LOG(info) << "Now to parse for Tracklets with a buffer length of " << mHalfCRUPayLoadRead;
                mTrackletsParser.setVerbose(mVerbose);
                int trackletwordsread = mTrackletsParser.Parse(&mCRUPayLoad); // this will read up to the tracnklet end marker.
                linkstart += trackletwordsread;
              //  LOG(info) << "We read " << trackletwordsread << " words of tracklets, now to parse for Digits with a buffer length of " << std::distance(mCRUPayLoad.begin(), linkend) << " data read : " << mHalfCRUPayLoadRead;
                int digitwordsread = 0;
                mDigitsParser.setVerbose(mVerbose);
                digitwordsread = mDigitsParser.Parse(&mCRUPayLoad, linkstart, linkend, currentdetector);
              } else {
                LOG(info) << "link start and end are the same, link appears to be empty";
              }
            } else {
              LOG(info) << "mState not CRUStateHalfCRUHeader mState:" << mState;
            }
          } //for loop over link index.
          // we have read in all the digits and tracklets for this event.
          //digits and tracklets are sitting inside the parsing classes.
          //extract the vectors and copy them to tracklets and digits here, building the indexing(triggerrecords)
          //TODO version 2 remove the tracklet and digit class and write directly the binary format.
          mEventTracklets.insert(std::end(mEventTracklets),std::begin(mTrackletsParser.getTracklets()),std::end(mTrackletsParser.getTracklets()));
          mEventDigits.insert(std::end(mEventDigits),std::begin(mDigitsParser.getDigits()),std::end(mDigitsParser.getDigits()));
          if(mVerbose)LOG(info) << "Event digits after eventi # : " << mEventDigits.size() << " having added : " << mDigitsParser.getDigits().size();
          //break; // end of CRU
        } else {
          if(mVerbose)LOG(info) << "Processed part of a half cru, did not finish loop around again for next rdh mState:" << mState;
        }
      }
      mState = CRUStateHalfCRU;
    }
    // move to next HBF
    rdh = (o2::header::RDHAny*)((char*)(rdh) + offsetToNext);
    totaldataread += memorySize;
//    LOG(info) << "total data read so far is : " << totaldataread;
  }
  //at the end of an event we have an rdh stop event.
  //package digits and tracklets up into a block wrap in an rdh and send on their way.A
  uint32_t numberoftracklets=mEventTracklets.size(); // 64 bit tracklets.
  uint32_t numberofdigits=mEventDigits.size();
//  memcpy(); //opening rdh to outputbuffer.
 // memcpy(); // buffer.
  //construct closing rdh.
  //memcpy(); //copy rdh

  if (o2::raw::RDHUtils::getStop(rdh)) {
    if (mDataPointer != (const uint32_t*)((char*)rdh + o2::raw::RDHUtils::getOffsetToNext(rdh))) {
      if(mVerbose)LOG(warn) << " at end of parsing loop and mDataPointer is on next rdh";
    }
    mDataPointer = (const uint32_t*)((char*)rdh + o2::raw::RDHUtils::getOffsetToNext(rdh));
    // make sure mDataPointer is in the correct place.
  } else
    mDataPointer = (const uint32_t*)((char*)rdh);
//  LOGF(info, " at exiting processHBF after advancing to next rdh mDataPointer is %p  ?= %p data read in total is now : %d  0x%08x", (void*)mDataPointer, (void*)mDataEndPointer, totaldataread+datasizealreadyread, datasizealreadyread+ totaldataread);

 // if(mVerbose)LOG(info) << "--- END PROCESS HBF";

  /* move to next RDH */
  // mDataPointer = (uint32_t*)((char*)(rdh) + o2::raw::RDHUtils::getOffsetToNext(rdh));

  /* otherwise return */

  return totaldataread; //mDataEndPointer - mDataPointer;
}

bool CruRawReader::buildCRUPayLoad()
{
  // copy data for the current half cru, and when we eventually get to the end of the payload return 1
  // to say we are done.
  int cruid = 0;
  int additionalBytes = -1;
  int crudatasize = -1;
  LOG(info) << "--- Build CRU Payload, added " << additionalBytes << " bytes to CRU "
            << cruid << " with new size " << crudatasize;
  return true;
}

bool CruRawReader::processHalfCRU()
{

    //TODO this needto be reactored.
    //TODO wedo not need to trapse through each hbf, we cam simply dump from start to stop rdh into a buffer, read the halfcruheader and 
    //replace this function with that..... version 2.... not fiddling now.
    
  /* process a FeeID/halfcru, 15 links */
  bool returnstate = true;
  LOG(info) << "--- PROCESS HalfCRU FeeID:" << mFEEID << " and state is : " << mState;
  //advance pointer to next halfcruheader.
  if (mState == CRUStateHalfCRUHeader) {
    // well then read the halfcruheader.
    //   LOG(info) << "AAAAAAAAAAAAAAA";
    memcpy(&mCurrentHalfCRUHeader, (void*)(mDataPointer), sizeof(mCurrentHalfCRUHeader)); //TODO remove the copy just use pointer dereferencing, doubt it will improve the speed much though.
    //mEncoderRDH = reinterpret_cast<o2::header::RDHAny*>(mEncoderPointer);)
    mDataPointer += sizeof(mCurrentHalfCRUHeader) / 4; //mDataPointer is in units of 32 bit uint
    o2::trd::getlinkdatasizes(mCurrentHalfCRUHeader, mCurrentHalfCRULinkLengths);
    o2::trd::getlinkerrorflags(mCurrentHalfCRUHeader, mCurrentHalfCRULinkErrorFlags);
    mTotalHalfCRUDataLength256 = std::accumulate(mCurrentHalfCRULinkLengths.begin(),
                                                 mCurrentHalfCRULinkLengths.end(),
                                                 decltype(mCurrentHalfCRULinkLengths)::value_type(0));
    mTotalHalfCRUDataLength = mTotalHalfCRUDataLength256 * 32; //convert to bytes.
    //debugLinks();
    LOG(info) << " total length is " << mTotalHalfCRUDataLength;
    int linkcount = 0;
    uint32_t totallength = 0;
    for (auto length : mCurrentHalfCRULinkLengths) {
      LOG(info) << linkcount++ << " had length : " << length;
      totallength += length * 32;
    }
    if(mVerbose)LOG(info) << "Total length counted in loop is : " << totallength;
    if(mVerbose)LOG(info) << mCurrentHalfCRUHeader;
    if(mVerbose)LOG(info) << "Now to process";
    // we will always have at least a length of 1 fully padded for each link.
    //TODO Sanity check,1. each link is >=1, 2. link is < ?? what is the maximum link length.3. header values are sane. define sane?
    //size sanity check.
    if (mTotalHalfCRUDataLength > mMaxCRUBufferSize * 4) { // x4 is to convert to bytes, so same unit comparison
      LOG(fatal) << "Cru wont fit in the allocated buffer  " << mTotalHalfCRUDataLength << " > " << mMaxCRUBufferSize * 4;
    }
    if(mVerbose)LOG(info) << "Found  a HalfCRUHeader : ";
    if(mVerbose)LOG(info) << mCurrentHalfCRUHeader << " with payload total size of : " << mTotalHalfCRUDataLength;
    mState = CRUStateHalfCRU; // we expect a halfchamber header now
    //TODO maybe change name to something more generic, this is will have to change to handle other data types config/adcdata.
    int rdhpayloadleft = (char*)mDataEndPointer - (char*)mDataPointer;
    mHalfCRUPayLoadRead = 0;
    //now sort out if we copy the remainderof the rdh payload or only a portion of it (tillthe end off the current halfcruheader's body
    if (mTotalHalfCRUDataLength < rdhpayloadleft) {
      if(mVerbose)LOG(info) << "read a halfcruheader at the top of the rdh payload, and it fits with in the rdh payload : " << mTotalHalfCRUDataLength << " < " << rdhpayloadleft;
      if(mVerbose)LOG(info) << "copying from " << (char*)mDataPointer << " to " << (char*)mDataPointer + mTotalHalfCRUDataLength;
      memcpy(&mCRUPayLoad[0], (void*)(mDataPointer), sizeof(mTotalHalfCRUDataLength)); //0 as we have just read the header.
      //advance pointer to next halfcruheader.
      mDataPointer += mTotalHalfCRUDataLength; // this cru half chamber is contained with in the a single rdh payload.
      mHalfCRUPayLoadRead += mTotalHalfCRUDataLength;
      mState = CRUStateHalfCRUHeader; // now back on a halfcruheader with in the current rdh payload.
      returnstate = true;
    } else {
      //otherwise we copy till the end of the rdh payload, and place mDataPointer on the next rdh header.
      if(mVerbose)LOG(info) << "read a halfcruheader at the top of the rdh payload, and does not fit within the rdh payload : " << mTotalHalfCRUDataLength << " > " << rdhpayloadleft;
      memcpy(&mCRUPayLoad[0], (void*)(mDataPointer), (char*)mDataEndPointer - (char*)mDataPointer); //0 as we have just read the header.
      //advance pointer to next halfcruheader.
      mHalfCRUPayLoadRead += (char*)mDataEndPointer - (char*)mDataPointer;
      mDataPointer = mDataEndPointer;
      mState = CRUStateHalfCRU;
      returnstate = false;
    }
  } // end of cruhalfchamber header at the top of rdh.
  else {
    if (mState == CRUStateHalfCRU) {
      //we are still busy inside a halfcru
      //copy the remainder of the halfcru or the entire rdh payload, which ever is smaller.
      //mCurrentLinkDataPosition;
      //mCurrentHalfCRULinkHeaderPosition;

      int remainderofrdhpayload = (char*)mDataEndPointer - (char*)mDataPointer;
      int remainderofcrudatablock = mTotalHalfCRUDataLength - mHalfCRUPayLoadRead;
      if (remainderofcrudatablock > remainderofrdhpayload) {
        // the halfchamber block extends past the end of this rdh we are currently on.
        // copy the data from where we are to the end into the crupayload buffer.

        int remainderofrdhpayload = mTotalHalfCRUDataLength - mHalfCRUPayLoadRead;
        if(mVerbose)LOG(info) << "in state CRUStateHalfChamber with remainderofrdhpayload (" << remainderofrdhpayload << ") = " << mTotalHalfCRUDataLength << "-" << mHalfCRUPayLoadRead;
        memcpy(&mCRUPayLoad[mHalfCRUPayLoadRead / 4], (void*)(mDataPointer), remainderofrdhpayload);
        mHalfCRUPayLoadRead += (char*)mDataEndPointer - (char*)mDataPointer;
        mDataPointer = mDataEndPointer;
        // the halfchamber block extends past the end of this rdh we are currently on.
        mState = CRUStateHalfCRU; // state stays the same, written here to be explicit.
        returnstate = false;
      } else {
        // the current cru payload we are on finishes before the end of the current rdh block we are in.
        int remainderofrdhpayloadthatwewant = mTotalHalfCRUDataLength - mHalfCRUPayLoadRead;
        //sanity check :
        if (remainderofrdhpayloadthatwewant < ((char*)mDataEndPointer - (char*)mDataPointer)) {
          LOG(warn) << " something odd we are supposed to have the cruhalfchamber ending with in the current rdh however : remainder of the rdh payload we want is : " << remainderofrdhpayloadthatwewant << " yet the rdh block only has " << mDataEndPointer - mDataPointer << " data left";
        }
        memcpy(&mCRUPayLoad[mHalfCRUPayLoadRead / 4], (void*)(mDataPointer), remainderofrdhpayloadthatwewant);
        mDataPointer += remainderofrdhpayloadthatwewant;
        mHalfCRUPayLoadRead += remainderofrdhpayloadthatwewant;
        mState = CRUStateHalfCRUHeader;
        returnstate = true;
      }
    } else {
      LOG(warn) << "huh unknown CRUstate of " << mState;
    }
    }

    LOG(info) << "--- END PROCESS HalfCRU with state: " << mState;

    return returnstate;
}

bool CruRawReader::processCRULink()
{
  /* process a CRU Link 15 per half cru */
  //  checkFeeID(); // check the link we are working with corresponds with the FeeID we have in the current rdh.
  //  uint32_t slotId = GET_TRMDATAHEADER_SLOTID(*mDataPointer);
  return false;
}

void CruRawReader::resetCounters()
{
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
    LOG(info) << "And away we go, run method of Translator";
    uint32_t dowhilecount = 0;
    uint64_t totaldataread = 0;
    do {
        LOG(info) << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! start";
        LOG(info) << "do while loop count " << dowhilecount++;
        //      LOG(info) << " data readin : " << mDataReadIn;
        LOG(info) << " mDataBuffer :" << (void*)mDataBuffer;
        int datareadfromhbf = processHBFs(totaldataread,mVerbose);
        LOG(info) << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! end with " << datareadfromhbf;
        totaldataread += datareadfromhbf; //     LOG(info) << "mDataReadIn :" << mDataReadIn << " mDataBufferSize:" << mDataBufferSize;
        LOG(info) << " Total data read so far is : " << totaldataread;
    } while (mDataReadIn < mDataBufferSize);

    return false;
};

void checkSummary();

} // namespace trd
} // namespace o2