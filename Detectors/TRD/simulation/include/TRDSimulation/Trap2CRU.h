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
//  the cru output. I suppose a cru simulator                                //
///////////////////////////////////////////////////////////////////////////////

#ifndef ALICE_O2_TRD_TRAP2CRU_H
#define ALICE_O2_TRD_TRAP2CRU_H

#include <string>
#include "DataFormatsTRD/TriggerRecord.h"
#include "DataFormatsTRD/LinkRecord.h"
#include "DataFormatsTRD/RawData.h"
#include "DataFormatsTRD/Tracklet64.h"
#include "TRDBase/Digit.h"
//#include "DetectorsRaw/HBFUtils.h"
#include "DetectorsRaw/RawFileWriter.h"

namespace o2
{
namespace trd
{

class Trap2CRU
{
  static constexpr int NumberOfCRU = 36;
  static constexpr int NumberOfHalfCRU = 72;
  static constexpr int NumberOfFLP = 12;
  static constexpr int CRUperFLP = 3;
  static constexpr int WordSizeInBytes = 256;  // word size in bits, everything is done in 256 bit blocks.
  static constexpr int WordSize = 8;           // 8 standard 32bit words.
  static constexpr int NLinksPerHalfCRU = 15;  // number of fibers per each half cru
  static constexpr int TRDLinkID = 15;         // hard coded link id. TODO give reason?
  static constexpr uint32_t PaddWord = 0xeeee; // pad word to fill 256bit blocks or entire block for no data case.
  static constexpr bool DebugDataWriting = false; //dump put very first vector of data to see if the thing going into the rawwriter is correct.
                                                  //TODO come back and change the mapping of 1077 channels to a lut addnd probably configurable.
                                                  //
 public:
  Trap2CRU() = default;
  Trap2CRU(const std::string& outputDir, const std::string& inputDigitsFilename, const std::string& inputTrackletsFilename);
  //Trap2CRU(const std::string& outputDir, const std::string& inputFilename, const std::string& inputDigitsFilename, const std::string& inputTrackletsFilename);
  void readTrapData();
  void convertTrapData(o2::trd::TriggerRecord const& trackletTrigRecord, o2::trd::TriggerRecord const& digitTriggerRecord);
  // default for now will be file per half cru as per the files Guido did for us.
  // TODO come back and give users a choice.
  //       void setFilePerLink(){mfileGranularity = mgkFilesPerLink;};
  //       bool getFilePerLink(){return (mfileGranularity==mgkFilesPerLink);};
  //       void setFilePerHalfCRU(){mfileGranularity = mgkFilesPerHalfCRU;};
  //       bool getFilePerHalfCRU(){return (mfileGranularity==mgkFilesPerHalfCRU);};  //
  int getVerbosity() { return mVerbosity; }
  void setVerbosity(int verbosity) { mVerbosity = verbosity; }
  void buildCRUPayLoad();
  int sortByORI();
  void sortDataToLinks();
  o2::raw::RawFileWriter& getWriter() { return mWriter; }
  uint32_t buildHalfCRUHeader(HalfCRUHeader& header, const uint32_t bc, const uint32_t halfcru);
  void linkSizePadding(uint32_t linksize, uint32_t& crudatasize, uint32_t& padding);
  void openInputFiles();
  uint32_t countTrackletsizeForLink();
  uint32_t countDigitsizeForLink();
  bool isTrackletOnLink(int link, int trackletpos);
  bool isDigitOnLink(int link, int digitpos);
  int buildRawMCMData(const int trackletindex);
  int buildDigitRawData(const int digitindex, const std::array<Digit,21>&localParseDigits, char *dataptr);
  int buildTrackletRawData(const int trackletindex, char *dataptr);

 private:
  int mfileGranularity; /// per link or per half cru for each file
  uint32_t mLinkID;
  uint16_t mCruID;
  uint64_t mFeeID;
  uint32_t mEndPointID;
  std::string mInputFileName;
  std::string mOutputFileName;
  int mVerbosity = 0;
  std::string mOutputDir;
  uint32_t mSuperPageSizeInB;
  std::string mInputDigitsFileName;
  std::string mInputTrackletsFileName;
  HalfCRUHeader mHalfCRUHeader;
  TrackletMCMHeader mTrackletMCMHeader;
  TrackletMCMData mTrackletMCMData;

  std::vector<char> mRawData; // store for building data event for a single half cru
  uint32_t mRawDataPos = 0;
  TFile* mDigitsFile;
  TTree* mDigitsTree;
  TFile* mTrackletsFile;
  TTree* mTrackletsTree;
  // locations to store the incoming data branches
  std::vector<Digit> mDigits, *mDigitsPtr=&mDigits;
  std::vector<uint32_t> mDigitsIndex;
  std::vector<o2::trd::TriggerRecord> mDigitTriggerRecords, *mDigitTriggerRecordsPtr = &mDigitTriggerRecords;
  std::vector<Tracklet64> mTracklets, *mTrackletsPtr=&mTracklets;
  std::vector<uint32_t> mTrackletsIndex;
  std::vector<o2::trd::TriggerRecord> mTrackletTriggerRecords, *mTrackletTriggerRecordsPtr = &mTrackletTriggerRecords;
  std::vector<uint32_t> mTrapRawData, *mTrapRawDataPtr = &mTrapRawData;

  const o2::raw::HBFUtils& mSampler = o2::raw::HBFUtils::Instance();
  o2::raw::RawFileWriter mWriter{"TRD"};

  ClassDefNV(Trap2CRU, 2);
};

} // end namespace trd
} // end namespace o2
#endif
