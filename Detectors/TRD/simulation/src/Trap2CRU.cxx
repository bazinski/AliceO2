// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD Trap2CRU class                                                       //
//  Class to take the trap output that arrives at the cru and produce        //
//  the cru output. A data mapping more than a cru simulator                 //
///////////////////////////////////////////////////////////////////////////////

#include <string>

#include "Headers/RAWDataHeader.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsTRD/TriggerRecord.h"
#include "DataFormatsTRD/LinkRecord.h"
#include "DataFormatsTRD/RawData.h"
#include "DataFormatsTRD/Tracklet64.h"
#include "DataFormatsTRD/Constants.h"
#include "DetectorsRaw/HBFUtils.h"
#include "DetectorsRaw/RawFileWriter.h"
#include "TRDSimulation/Trap2CRU.h"
#include "TRDSimulation/TrapSimulator.h"
#include "CommonUtils/StringUtils.h"
#include "TRDBase/CommonParam.h"
#include "TFile.h"
#include "TTree.h"
#include <TStopwatch.h>
#include <cassert>
#include <fstream>
#include <iostream>
#include <array>
#include <string>
#include <bitset>
#include <vector>
#include <gsl/span>

using namespace o2::raw;

namespace o2
{
namespace trd
{

Trap2CRU::Trap2CRU(const std::string& outputDir, const std::string& inputdigitsfilename, const std::string& inputtrackletsfilename)
{
    mOutputDir=outputDir;
    mSuperPageSizeInB=1024*1024;
//    mInputFileName=inputFilename;
    mInputDigitsFileName=inputdigitsfilename;
    mInputTrackletsFileName=inputtrackletsfilename;
  readTrapData();
}

void Trap2CRU::openInputFiles()
{
  mDigitsFile = TFile::Open(mInputDigitsFileName.data());
  if(mDigitsFile!=nullptr){
    mDigitsTree  = (TTree*)mDigitsFile->Get("o2sim");
    mDigitsTree->SetBranchAddress("TRDDigit",&mDigitsPtr); // the branch with the actual digits
    mDigitsTree->SetBranchAddress("TriggerRecord",&mDigitTriggerRecordsPtr); // branch with trigger records for digits
    LOG(info) << "Digit tree has " << mDigitsTree->GetEntries() << " entries";
  }
  else{
      LOG(fatal) << " cant open digit tree";
  }
  mTrackletsFile=TFile::Open(mInputTrackletsFileName.data());
  if(mTrackletsFile!=nullptr){
    mTrackletsTree  = (TTree*)mTrackletsFile->Get("o2sim");
    mTrackletsTree->SetBranchAddress("Tracklet",&mTrackletsPtr); // the branch with the actual tracklets.
    mTrackletsTree->SetBranchAddress("TrackTrg",&mTrackletTriggerRecordsPtr); // branch with trigger records for digits
    LOG(info) << "Tracklet tree has " << mTrackletsTree->GetEntries() << " entries";
  }
  else{
      LOG(fatal) << " cant open tracklet tree";
  }
}



bool trackletsindexcompare(const unsigned int A, const unsigned int B, const std::vector<o2::trd::Tracklet64>& originalTracklets)
{
  // sort into ORI:mcm
  // giving us a stream of data of link:mcm which is the structure of the raw data in the cru
  const o2::trd::Tracklet64 *a, *b;
  a = &originalTracklets[A];
  b = &originalTracklets[B];
  uint32_t oria=FeeParam::getORI(a->getDetector(),a->getROB());
  uint32_t orib=FeeParam::getORI(b->getDetector(),b->getROB());
  if (oria < orib) {
    return 1;
  }
  if (oria > orib) {
    return 0;
  }
  if (a->getMCM() < b->getMCM()) {
    return 0;
  }
  if (a->getMCM() > b->getMCM()) {
    return 1;
  }
  return 0;
}

bool digitindexcompare(const unsigned int A, const unsigned int B, const std::vector<o2::trd::Digit>& originalDigits)
{
  // sort into ORI:mcm
  // giving us a stream of data of link:mcm which is the structure of the raw data in the cru
  const o2::trd::Digit *a, *b;
  a = &originalDigits[A];
  b = &originalDigits[B];
  uint32_t oria=FeeParam::getORI(a->getDetector(),a->getROB());
  uint32_t orib=FeeParam::getORI(b->getDetector(),b->getROB());
  if (oria < orib) {
    return 1;
  }
  if (oria > orib) {
    return 0;
  }
  if (a->getMCM() < b->getMCM()) {
    return 0;
  }
  if (a->getMCM() > b->getMCM()) {
    return 1;
  }
  return 0;
}

void Trap2CRU::sortDataToLinks()
{
  //sort the digits array TODO refactor this intoa vector index sort and possibly generalise past merely digits.
  auto sortstart = std::chrono::high_resolution_clock::now();
  //build indexes 
  // digits first
  std::generate(mDigitsIndex.begin(), mDigitsIndex.end(), [n = 0]() mutable { return n++; });
  std::generate(mTrackletsIndex.begin(), mTrackletsIndex.end(), [n = 0]() mutable { return n++; });

  for (auto& trig : mDigitTriggerRecords) {
    std::stable_sort(std::begin(mDigitsIndex) + trig.getFirstEntry(), std::begin(mDigitsIndex) + trig.getNumberOfObjects() + trig.getFirstEntry(),
                     [this](auto&& PH1, auto&& PH2) { return digitindexcompare(PH1, PH2, mDigits); });
  }
  for (auto& trig : mTrackletTriggerRecords) {
    std::stable_sort(std::begin(mTrackletsIndex) + trig.getFirstEntry(), std::begin(mTrackletsIndex) + trig.getNumberOfObjects() + trig.getFirstEntry(),
                     [this](auto&& PH1, auto&& PH2) { return trackletsindexcompare(PH1, PH2, mTracklets); });
  }


  auto sortingtime = std::chrono::high_resolution_clock::now() - sortstart;
  LOG(debug) << "TRD Digit Sorting took " << sortingtime.count();


}

void Trap2CRU::readTrapData()
{
    //set things up, read the file and then deligate to convertTrapdata to do the conversion.
    //
    mRawData.reserve(1024 * 1024); //TODO take out the hardcoded 1MB its supposed to come in from the options
    LOG(debug) << "Trap2CRU::readTrapData";
    // data comes in index by event (triggerrecord) and link (linker record) both sequentially.
    // first 15 links go to cru0a, second 15 links go to cru0b, 3rd 15 links go to cru1a ... first 90 links to flp0 and then repeat for 12 flp
    // then do next event

    // lets register our links
    std::string prefix = mOutputDir;
    if (!prefix.empty() && prefix.back() != '/') {
        prefix += '/';
    }

    for (int link = 0; link < NumberOfHalfCRU; link++) {
        // FeeID *was* 0xFEED, now is indicates the cru Supermodule, side (A/C) and endpoint. See RawData.cxx for details.
        int supermodule = link / 4;
        int endpoint = link / 2;
        int side = link % 2 ? 1 : 0;
        mFeeID = buildTRDFeeID(supermodule, side, endpoint);
        mCruID = link / 2;
        mEndPointID = link % 2 ? 1 : 0; //TODO figure out a value ... endpoint needs a rebase to PR4106
        mLinkID = TRDLinkID;
        std::string trdside = link % 2 ? "c" : "a"; // the side of supermodule readout A or C, odd numbered CRU are A, even numbered CRU are C.
        // filenmae structure of trd_cru_[CRU#]_[upper/lower].raw
        std::string outputFilelink = o2::utils::concat_string(prefix, "trd_cru_", std::to_string(mCruID), "_", trdside, ".raw");
        mWriter.registerLink(mFeeID, mCruID, mLinkID, mEndPointID, outputFilelink);
    }
    openInputFiles();
    LOG(info) << "files opened now for some sorting";
    sortDataToLinks();
    LOG(info) << "digit and tracklets sorted ";
    LOG(info) << "lets look at tracklet tree for now";
    if(mDigitsTree->GetEntries() != 1) LOG(fatal) << "more than one entry in digits tree" ;
    if(mTrackletsTree->GetEntries() != 1) LOG(fatal) << "more than one entry in tracklets tree" ;

    for (int entry = 0; entry < mDigitsTree->GetEntries(); entry++) {
        mDigitsTree->GetEntry(entry);
        uint32_t linkcount = 0;
    }
    for (int entry = 0; entry < mTrackletsTree->GetEntries(); entry++) {
        mTrackletsTree->GetEntry(entry);
    }
    // everything is passed as the first entry in the tree. Ensure it *really* is only a single entry tree
    // and then work with the resulting vectors from the branches.
        uint32_t linkcount = 0;
        for (auto tracklettrigger : mTrackletTriggerRecords) {
            //get the event limits from TriggerRecord;
            uint32_t eventstart = tracklettrigger.getFirstEntry();
            uint32_t eventend = trackelttrigger.getFirstEntry() + tracklettrigger.getNumberOfObjects();
            for (auto digitstrigger : mDigitsTriggerRecords) { 
                //TODO come back and avoid the double looping.
                //get the event limits from TriggerRecord;
                uint32_t eventstart = digitstrigger.getFirstEntry();
                uint32_t eventend = digitstrigger.getFirstEntry() + digitstrigger.getNumberOfObjects();
                if(digitstrigger.BunchCrossing  == tracklettrigger.BunchCrossing){
                    convertTrapData(tracklettrigger,digitstrigger);
                }
                //TODO check if the triggerrecrods are *always* in sync?
            }
        }
}

int Trap2CRU::sortByORI()
{
    // data comes in sorted by padcolum, a row of 8 trap chips.
    // this is sadly not how the electronics is actually connected.
    // we therefore need to resort the data according to the ORI link.
    return 1;
}

void Trap2CRU::buildCRUPayLoad()
{
    // go through the data for the event in question, sort via above method, produce the raw stream for each cru.
    // i.e. 30 link per cru, 3cru per flp.
    // 30x [HalfCRUHeader, TrackletHCHeader0, [MCMHeader TrackletMCMData. .....] TrackletHCHeader1 ..... TrackletHCHeader30 ...]
    //
    // data must be padded into blocks of 256bit so on average 4 padding 32 bit words.
}

void Trap2CRU::linkSizePadding(uint32_t linksize, uint32_t& crudatasize, uint32_t& padding)
{
    // all data must be 256 bit aligned (8x64bit).
    // if zero the whole 256 bit must be padded (empty link)
    // crudatasize is the size to be stored in the cruheader, i.e. units of 256bits.
    // linksize is the incoming link size from the linkrecord,
    // padding is of course the amount of padding in 32bit words.
    uint32_t rem = 0;
    if (linksize != 0) {
        //data, so figure out padding cru word, the other case is simple, full padding if size=0
        rem = linksize % 8;
        if (rem != 0) {
            crudatasize = linksize / 8 + 1;
            padding = 8 - rem;
        } else {

            crudatasize = linksize / 8; // 32 bit word to 256 bit word.
            padding = 0;
        }
        LOG(debug) << "We have data with linkdatasize=" << linksize << " with size number in header of:" << crudatasize << " padded with " << padding << " 32bit words";
    } else {
        //linksize is zero so no data, pad fully.
        crudatasize = 1;
        padding = 8;
        LOG(debug) << "We have data with linkdatasize=" << linksize << " with size number in header of:" << crudatasize << " padded with " << padding << " 32bit words";
    }
    LOG(debug) << "linkSizePadding : CRUDATASIZE : " << crudatasize;
}



uint32_t Trap2CRU::countTrackletsizeForLink()
{

}

uint32_t Trap2CRU::countDigitsizeForLink()
{

}

uint32_t Trap2CRU::buildHalfCRUHeader(HalfCRUHeader& header, const uint32_t bc, const uint32_t halfcru)
{
    int bunchcrossing = bc;
    int stopbits = 0x01; // do we care about this and eventtype in simulations?
    int eventtype = 0x01;
    int crurdhversion = 6;
    int feeid = 0;
    int cruid = 0;
    uint32_t crudatasize = 0;           //link size in units of 256 bits.
    int endpoint = halfcru % 2 ? 1 : 0; 
    uint32_t padding = 0;
    setHalfCRUHeader(header, crurdhversion, bunchcrossing, stopbits, endpoint, eventtype, feeid, cruid); //TODO come back and pull this from somewhere.

    // halfcruheader from the relevant mLinkRecords.
    int totallinkdatasize = 0; //in units of 256bits
    for (int link = 0; link < NLinksPerHalfCRU; link++) {
        int linkid = link + halfcru * NLinksPerHalfCRU; // TODO this might have to change to a lut I dont think the mapping is linear.
        int errors = 0;
        int linksize = 0; // linkSizePadding will convert it to 1 for the no data case.
        //count tracklets
        //
        // # cru -- AorC 
        // 18 supermodules each with A or c side
        // 36 cru for 0A 0C 1A 1C 2A 2C ... 18C
        // each cru has 2 end points each wih 15 links
        
        cru%2 = 0=a side 1=c side
        then upper or lower given 15-29 or 0-14 inclusive
        int linktrackletsize = countTrackletsizeForLink();
        //count digits
        int linkdigitsize = countDigitsizeForLink();
        int linkdatasize=linktrackletsize+linkdigitsize;
        linkSizePadding(linksize, crudatasize, padding);
        setHalfCRUHeaderLinkData(header, link, crudatasize, errors); // write one padding block for empty links.
        totallinkdatasize += crudatasize;
    }
    return totallinkdatasize;
}

void Trap2CRU::convertTrapData(o2::trd::TriggerRecord const& trackletTrigRecord, o2::trd::TriggerRecord const& digitTriggerRecord)
{

    //build a HalfCRUHeader for this event/cru/endpoint
    //loop over cru's
    //  loop over all half chambers, thankfully they data is sorted.
    //    check if current chamber has a link
    //      if not blank, else fill in data from link records
    //  dump data to rawwriter
    //finished for event. this method is only called per event.
    int currentlinkrecord = 0;
    char* traprawdataptr = (char*)&mTrapRawData[0];

    for (int halfcru = 0; halfcru < o2::constants::NHalfCRU; halfcru++) {
        int supermodule = halfcru / 4;
        int endpoint = halfcru / 2;
        int side = halfcru % 2 ? 1 : 0;
        mFeeID = buildTRDFeeID(supermodule, side, endpoint);
        mCruID = halfcru / 2;
        mLinkID = TRDLinkID;
        mEndPointID = halfcru % 2 ? 1 : 0;
        //15 links per half cru or cru end point.
        int oristart=halfcru*15;
        int oriend=(halfcru+1)*15;
        memset(&mRawData[0], 0, sizeof(mRawData[0]) * mRawData.size()); //   zero the rawdata storage
        int numberofdetectors = o2::trd::constants::MAXCHAMBER;
        HalfCRUHeader halfcruheader;
        //now write the cruheader at the head of all the data for this halfcru.
        LOG(debug) << "cru before building cruheader for halfcru index : " << halfcru << " with contents \n"
            << halfcruheader;
        uint32_t totalhalfcrudatasize = buildHalfCRUHeader(halfcruheader, trackletTriggerRecord.getBCData().bc, halfcru);

        std::vector<char> rawdatavector(totalhalfcrudatasize * 32 + sizeof(halfcruheader)); // sum of link sizes + padding in units of bytes and some space for the header (512 bytes).
        char* rawdataptr = rawdatavector.data();

        //dumpHalfCRUHeader(halfcruheader);
        memcpy(rawdataptr, (char*)&halfcruheader, sizeof(halfcruheader));
        std::array<uint64_t, 8> raw{};
        memcpy((char*)&raw[0], rawdataptr, sizeof(halfcruheader));
        for (int i = 0; i < 2; i++) {
            int index = 4 * i;
            LOG(debug) << "[1/2rawdaptr " << i << " " << std::hex << raw[index + 3] << " " << raw[index + 2] << " " << raw[index + 1] << " " << raw[index + 0];
        }
        rawdataptr += sizeof(halfcruheader);
        int linkdatasize = 0; // in 32 bit words
        int link = halfcru / 2;
        for (int halfcrulink = 0; halfcrulink < NLinksPerHalfCRU; halfcrulink++) {
            //links run from 0 to 14, so linkid offset is halfcru*15;
            int linkid = halfcrulink + halfcru * NLinksPerHalfCRU;
            int ori=oristart+halfcrulink;
            LOG(debug) << "Currently checking for data on linkid : " << linkid << " from halfcru=" << halfcru << " and halfcrulink:" << halfcrulink << " ?? " << linkid << "==" << mLinkRecords[currentlinkrecord].getLinkId();
            int errors = 0;           // put no errors in for now.
            int size = 0;             // in 32 bit words
            int datastart = 0;        // in 32 bit words
            int dataend = 0;          // in 32 bit words
            uint32_t paddingsize = 0; // in 32 bit words
            uint32_t crudatasize = 0; // in 256 bit words.
            if (mLinkRecords[currentlinkrecord].getLinkId() == linkid) {
                //this link has data in the stream.
                linkdatasize = mLinkRecords[currentlinkrecord].getNumberOfObjects();
                datastart = mLinkRecords[currentlinkrecord].getFirstEntry();
                dataend = datastart + size;
                LOG(debug) << "We have data on linkid = " << linkid << " and linksize : " << linkdatasize << " so :" << linkdatasize / 8 << " 256 bit words";
                currentlinkrecord++;
            } else {
                assert(mLinkRecords[currentlinkrecord].getLinkId() < linkid);
                //blank data for this link
                // put in a 1 256 bit word of data for the link and padd with 0xeeee x 8
                linkdatasize = 0;
                paddingsize = 8;
            }
            // now copy data to rawdata, padding as and where needed.
            //
            linkSizePadding(linkdatasize, crudatasize, paddingsize); //TODO this can come out as we have already called it, but previously we have lost the #padding words, solve to remove.

            // now pad ....
            char* olddataptr = rawdataptr; // store the old pointer so we can do some sanity checks for how far we advance.
            //linkdatasize is the #of 32 bit words coming from the incoming tree.
            //paddingsize is the number of padding words to add 0xeeee
            uint32_t bytestocopy = linkdatasize * (sizeof(uint32_t));
            memcpy(rawdataptr, traprawdataptr, bytestocopy);
            //increment pointer
            rawdataptr += bytestocopy;
            traprawdataptr += bytestocopy;
            //now for padding
            uint16_t padbytes = paddingsize * sizeof(uint32_t);
            memset(rawdataptr, 0xee, padbytes);
            //increment pointer.
            rawdataptr += padbytes;
            if (padbytes + bytestocopy != crudatasize * 32) {
                LOG(debug) << "something wrong with data size writing padbytes:" << padbytes << " bytestocopy : " << bytestocopy << " crudatasize:" << crudatasize;
            }
            //sanity check for now:
            if (((char*)rawdataptr - (char*)olddataptr) != crudatasize * 32) { // cru words are 8 uint32 and comparison is in bytes.
            }
            if (crudatasize != o2::trd::getlinkdatasize(halfcruheader, halfcrulink)) {
                // we have written the wrong amount of data ....
                LOG(debug) << "crudata is ! = get link data size " << crudatasize << "!=" << o2::trd::getlinkdatasize(halfcruheader, halfcrulink);
            }
            LOG(debug) << "copied " << crudatasize * 32 << "bytes to halfcrurawdata which now has  size of " << rawdatavector.size() << " for " << link << ":" << endpoint;
        }
        mWriter.addData(mFeeID, mCruID, mLinkID, mEndPointID, trackletTriggerRecord.getBCData(), rawdatavector);
        if (DebugDataWriting) {
            std::ofstream out2("crutestdumprawdatavector");
            out2.write(rawdatavector.data(), rawdatavector.size());
            out2.flush();
            halfcru = NumberOfHalfCRU; // exit loop after 1 half cru for now.
        }
    }
}

} // end namespace trd
} // end namespace o2
