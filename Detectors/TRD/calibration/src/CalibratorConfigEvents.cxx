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

/// \file CalibratorConfigEvents.cxx
/// \brief TimeSlot-based calibration of vDrift and ExB
/// \author Ole Schmidt

#include "Framework/ProcessingContext.h"
#include "Framework/TimingInfo.h"
#include "Framework/InputRecord.h"
#include "CCDB/CcdbApi.h"
#include "CCDB/BasicCCDBManager.h"
#include "CommonUtils/NameConf.h"
#include "CommonUtils/MemFileHelper.h"

#include "DataFormatsTRD/Constants.h"
#include "TRDCalibration/CalibratorConfigEvents.h"
#include "DataFormatsTRD/TrapConfigEvent.h"
#include "TRDBase/GeometryBase.h"

#include "TStopwatch.h"
#include <TFile.h>
#include <TTree.h>

#include <string>
#include <map>
#include <memory>
#include <ctime>

using namespace o2::trd::constants;

namespace o2::trd
{


void CalibratorConfigEvents::process(const gsl::span<const TrapConfigEvent>& trapconfigevents)
{
  if (mInitCompleted) {
    return;
  }

  // set tree addresses
  if (mEnableOutput) {
    // mOutTree->Branch("trapconfigevent", &mTrapConfigEvent);
    //  what do we want to save?
  }

  mInitCompleted = true;
  //take incoming trapconfigevents and populate the map of map of seen values for eventual
  //collapsing into a singular trapconfigevent at the end of the alloted accumulation time.
  //additionally fill in the mcmindex and the [mcm|hcid]present
  int count=0;
  for(auto& confevt : trapconfigevents){
//    for(auto& mcmidx : confevt.getMCMEvent()){
     LOGP(info,"processing config event of sequence #{}",count++); 
    //}
  }

}

bool CalibratorConfigEvents::hasEnoughData() const
{
  // determine if this slot has enough data ... normal calibration.
  // this method triggers a merge and a finaliseSlot.
  // if the slot does not have enough data this slot together with the next one are merged.

  // We therefore have 2 options :
  // a:keep all the difference and ergo, save to ccdb after each incoming config event.
  // b: the normal option, accumulate for 10,15 minutes and then compare and write to ccdb if new.
  bool timereached = 0;
  if (mSaveAllChanges) {
    return true;
  } else {
    // case b:
    //  we just care about the time taken since the beginning.
    //
    if (timereached) {
      return true;
    } else {
      return false;
    }
  }
}

/*void CalibratorConfigEvents::finalizeSlot()
{
  // do actual calibration for the data provided in the given slot
  TStopwatch timer;
  timer.Start();
 // initProcessing();

  timer.Stop();
  LOGF(info, "Done merging TrapConfigs. CPU time: %f, real time: %f", timer.CpuTime(), timer.RealTime());

  // Fill Tree and log to debug
  if (mEnableOutput) {
//    mOutTree->Fill();
  }

  // assemble CCDB object

  auto className = o2::utils::MemFileHelper::getClassName(mCCDBObject);
  auto fileName = o2::ccdb::CcdbApi::generateFileName(className);
  std::map<std::string, std::string> metadata; // TODO: do we want to store any meta data?
  long startValidity = 0;//slot.getStartTimeMS() - 10 * o2::ccdb::CcdbObjectInfo::SECOND;
  o2::ccdb::CcdbObjectInfo mInfoVector; ///< vector of CCDB infos; each element is filled with CCDB description of accompanying CCDB calibration object
  mInfoVector.setPath("TRD/TrapConfigEvent");
  mInfoVector.setObjectType(className);
  mInfoVector.setFileName(fileName);
  mInfoVector.setMetaData(metadata);
  mInfoVector.setStartValidityTimestamp(startValidity);
  mInfoVector.setEndValidityTimestamp(startValidity + o2::ccdb::CcdbObjectInfo::MINUTE * 15);
  // remove the 15 minute and change to update  on next start.
  // mObjectVector=.push_back(mCCDBObject);
  // mCCDBObjext is already filled.
}
*/
/*
void CalibratorConfigEvents::createFile()
{
  mEnableOutput = true;
  mOutFile = std::make_unique<TFile>("trd_configevents.root", "RECREATE");
  if (mOutFile->IsZombie()) {
    LOG(error) << "Failed to create output file!";
    mEnableOutput = false;
    return;
  }
  mOutTree = std::make_unique<TTree>("calib", "Config Events");
  LOG(info) << "Created output file trd_configevents.root";
}

void CalibratorConfigEvents::closeFile()
{
  if (!mEnableOutput) {
    return;
  }

  try {
    mOutFile->cd();
    mOutTree->Write();
    mOutTree.reset();
    mOutFile->Close();
    mOutFile.reset();
  } catch (std::exception const& e) {
    LOG(error) << "Failed to write calibration data file, reason: " << e.what();
  }
  // after closing, we won't open a new file
  mEnableOutput = false;
}
*/
} // namespace o2::trd
