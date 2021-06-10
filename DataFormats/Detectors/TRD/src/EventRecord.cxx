
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

#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsTRD/TriggerRecord.h"
#include "DataFormatsTRD/Tracklet64.h"
#include "DataFormatsTRD/Digit.h"
#include "DataFormatsTRD/EventRecord.h"
#include "DataFormatsTRD/CompressedDigit.h"
#include "DataFormatsTRD/Constants.h"
#include <cassert>
#include <array>
#include <string>
#include <bitset>
#include <vector>
#include <gsl/span>
#include <typeinfo>


namespace o2::trd
{

  //Digit information
  std::vector<Digit>& EventRecord::getDigits() { return mDigits; }
  std::vector<CompressedDigit>& EventRecord::getCompressedDigits() { return mCompressedDigits; }
  void EventRecord::addDigits(Digit& digit) { mDigits.push_back(digit); }
  void EventRecord::addCompressedDigits(CompressedDigit& digit) { mCompressedDigits.push_back(digit); }
  void EventRecord::addDigits(std::vector<Digit>::iterator& start, std::vector<Digit>::iterator& end) { mDigits.insert(std::end(mDigits), start, end); }
  void EventRecord::addCompressedDigits(std::vector<CompressedDigit>::iterator& start, std::vector<CompressedDigit>::iterator& end) { mCompressedDigits.insert(std::end(mCompressedDigits), start, end); }

  //tracklet information
  std::vector<Tracklet64>& EventRecord::getTracklets() { LOG(info)<<"mtracklets in get tracklets base ptr:"<< &mTracklets[0];return mTracklets; }
  void EventRecord::addTracklet(Tracklet64& tracklet) { mTracklets.push_back(tracklet); }
  void EventRecord::addTracklets(std::vector<Tracklet64>::iterator& start, std::vector<Tracklet64>::iterator& end) {
    LOG(info) << "inserting tracklets "<< std::distance(end,start);
    LOG(info) << "event record had tracklets of size:"<< mTracklets.size();
    LOG(info) << "!@!@!@!@!@!@!@!@!@!@!@!@!@@!@!@!@@!@ and we will crash here! ";
    LOG(info) << "insterting at : "<< &mTracklets[0] << " with data from :" << *start << " till " << *end;
    LOG(info) << "!@!@!@!@!@!@!@!@!@!@!@!@!@@!@!@!@@!@ and we will crash here! ";
    mTracklets.insert(std::end(mTracklets), start, end);
  }
  void EventRecord::addTracklets(std::vector<Tracklet64>& tracklets) 
  {
    LOG(info)<<"mtracklets in addtraclets base ptr:"<< &mTracklets[0];
    LOG(info) << "event record had tracklets of size:"<< mTracklets.size();
    LOG(info) << "!@!@!@!@!@!@!@!@!@!@!@!@!@@!@!@!@@!@ and we will crash here! ";
    for(auto tracklet : tracklets){
      mTracklets.push_back(tracklet);
    }
    //mTracklets.insert(mTracklets.end(), tracklets.begin(),tracklets.end()); 
  
  }

// now for event storage
  void EventStorage::addDigits(InteractionRecord& ir, Digit& digit)
  {
    bool added = false;
    for (auto eventrecord : mEventRecords) {
      if (ir == eventrecord.getBCData()) {
        //TODO replace this with a hash/map not a vector
        eventrecord.addDigits(digit);
        added = true;
      }
    }
    if (!added) {
      // unseen ir so add it
      mEventRecords.push_back(ir);
      mEventRecords.end()->addDigits(digit);
    }
  }
  void EventStorage::addCompressedDigits(InteractionRecord& ir, CompressedDigit& digit)
  {
    bool added = false;
    for (auto eventrecord : mEventRecords) {
      if (ir == eventrecord.getBCData()) {
        //TODO replace this with a hash/map not a vector
        eventrecord.addCompressedDigits(digit);
        added = true;
      }
    }
    if (!added) {
      // unseen ir so add it
      mEventRecords.push_back(ir);
      mEventRecords.end()->addCompressedDigits(digit);
    }
  }
  void EventStorage::addDigits(InteractionRecord& ir, std::vector<Digit>::iterator start, std::vector<Digit>::iterator end)
  {
    bool added = false;
    for (auto eventrecord : mEventRecords) {
      if (ir == eventrecord.getBCData()) {
        //TODO replace this with a hash/map not a vector
        eventrecord.addDigits(start, end);
        added = true;
      }
    }
    if (!added) {
      // unseen ir so add it
      mEventRecords.push_back(ir);
      mEventRecords.end()->addDigits(start, end);
    }
  }
  void EventStorage::addCompressedDigits(InteractionRecord& ir, std::vector<CompressedDigit>::iterator start, std::vector<CompressedDigit>::iterator end)
  {
    bool added = false;
    for (auto eventrecord : mEventRecords) {
      if (ir == eventrecord.getBCData()) {
        //TODO replace this with a hash/map not a vector
        eventrecord.addCompressedDigits(start, end);
        added = true;
      }
    }
    if (!added) {
      // unseen ir so add it
      mEventRecords.push_back(ir);
      mEventRecords.end()->addCompressedDigits(start, end);
    }
  }
  void EventStorage::addTracklet(InteractionRecord& ir, Tracklet64& tracklet)
  {
    bool added = false;
    for (auto eventrecord : mEventRecords) {
      if (ir == eventrecord.getBCData()) {
        //TODO replace this with a hash/map not a vector
        eventrecord.addTracklet(tracklet);
        added = true;
      }
    }
    if (!added) {
      // unseen ir so add it
      mEventRecords.push_back(ir);
      mEventRecords.end()->addTracklet(tracklet);
    }
  }
  void EventStorage::addTracklets(InteractionRecord& ir, std::vector<Tracklet64>& tracklets)
  {
    bool added = false;
    for (auto eventrecord : mEventRecords) {
      if (ir == eventrecord.getBCData()) {
        //TODO replace this with a hash/map not a vector
        eventrecord.addTracklets(tracklets); //mTracklets.insert(mTracklets.end(),start,end);
        added = true;
      }
    }
    if (!added) {
      // unseen ir so add it
      mEventRecords.push_back(ir);
      LOG(info) << "added a  new IR:" << ir;
      LOG(info) << "incoming tracklets have size:"<< tracklets.size();
      LOG(info)<<"mtracklets in addtraclets base ptr:"<< &(mEventRecords.back().getTracklets()[0]);
      LOG(info) << "event record had tracklets of size:"<< mEventRecords.back().getTracklets().size();
      LOG(info) << "now to actually add the tracklet";
      mEventRecords.end()->addTracklets(tracklets);
    }
  }
  void EventStorage::addTracklets(InteractionRecord& ir, std::vector<Tracklet64>::iterator& start, std::vector<Tracklet64>::iterator& end)
  {
    bool added = false;
    for (auto eventrecord : mEventRecords) {
      if (ir == eventrecord.getBCData()) {
        //TODO replace this with a hash/map not a vector
        eventrecord.addTracklets(start, end); //mTracklets.insert(mTracklets.end(),start,end);
        added = true;
      }
    }
    if (!added) {
      // unseen ir so add it
      LOG(info)<<"adding tracklets";
      mEventRecords.push_back(ir);
      mEventRecords.end()->addTracklets(start, end);
    }
  }
  void EventStorage::unpackDataForSending(std::vector<TriggerRecord>& triggers, std::vector<Tracklet64>& tracklets, std::vector<Digit>& digits)
  {
    int digitcount = 0;
    int trackletcount = 0;
    for (auto event : mEventRecords) {
      tracklets.insert(std::end(tracklets), std::begin(event.getTracklets()), std::end(event.getTracklets()));
      digits.insert(std::end(digits), std::begin(event.getDigits()), std::end(event.getDigits()));
      triggers.emplace_back(event.getBCData(), digitcount, event.getDigits().size(), trackletcount, event.getTracklets().size());
      digitcount += event.getDigits().size();
      trackletcount += event.getTracklets().size();
    }
  }
  void EventStorage::unpackDataForSending(std::vector<TriggerRecord>& triggers, std::vector<Tracklet64>& tracklets, std::vector<CompressedDigit>& digits)
  {
    int digitcount = 0;
    int trackletcount = 0;
    for (auto event : mEventRecords) {
      tracklets.insert(std::end(tracklets), std::begin(event.getTracklets()), std::end(event.getTracklets()));
      digits.insert(std::end(digits), std::begin(event.getCompressedDigits()), std::end(event.getCompressedDigits()));
      triggers.emplace_back(event.getBCData(), digitcount, event.getDigits().size(), trackletcount, event.getTracklets().size());
      digitcount += event.getDigits().size();
      trackletcount += event.getTracklets().size();
    }
  }
  int EventStorage::sumTracklets()
  {
    int sum = 0;
    for (auto event : mEventRecords) {
      sum += event.getTracklets().size();
    }
    return sum;
  }
  int EventStorage::sumDigits()
  {
    int sum = 0;
    for (auto event : mEventRecords) {
      sum += event.getDigits().size();
    }
    return sum;
  }
  std::vector<Tracklet64>& EventStorage::getTracklets(InteractionRecord& ir)
  {
    bool found = false;
    for (int i=0;i<mEventRecords.size();++i) {
      if (ir == mEventRecords[i].getBCData()) {
        found = true;
        return mEventRecords[i].getTracklets();
      }
    }
    LOG(warn) << "attempted to get tracklets from IR: " << ir << " total tracklets of:" << sumTracklets();
    printIR();
    return mDummyTracklets;
  }
  std::vector<Digit>& EventStorage::getDigits(InteractionRecord& ir)
  {
    bool found = false;
    for (auto event : mEventRecords) {
      if (ir == event.getBCData()) {
        found = true;
        return event.getDigits();
      }
    }
    LOG(warn) << "attempted to get digits from IR: " << ir << " total digits of:" << sumDigits();
    printIR();
    return mDummyDigits;
  }

  std::vector<CompressedDigit>& EventStorage::getCompressedDigits(InteractionRecord& ir)
  {
    bool found = false;
    for (auto event : mEventRecords) {
      if (ir == event.getBCData()) {
        found = true;
        return event.getCompressedDigits();
      }
    }
    LOG(warn) << "attempted to get digits from IR: " << ir << " total digits of:" << sumDigits();
    printIR();
    return mDummyCompressedDigits;
  }

  void EventStorage::printIR()
  {
    std::string records;
    int count = 0;
    for (auto event : mEventRecords) {
      LOG(info) << "[" << count << "]" << event.getBCData() << " ";
      count++;
    }
  }

} // namespace o2::trd

