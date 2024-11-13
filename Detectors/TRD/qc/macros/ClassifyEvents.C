// Copyright 2019-2023 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright
// holders. All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <array>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <numeric>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include <TBox.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TPolyLine.h>
#include <TString.h>
#include <TTree.h>

#include "CommonUtils/NameConf.h"
#include "DataFormatsTRD/CalVdriftExB.h"
#include "DataFormatsTRD/CalibratedTracklet.h"
#include "DataFormatsTRD/Constants.h"
#include "DataFormatsTRD/Digit.h"
#include "DataFormatsTRD/HelperMethods.h"
#include "DataFormatsTRD/Hit.h"
#include "DataFormatsTRD/Tracklet64.h"
#include "DataFormatsTRD/TriggerRecord.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "MathUtils/Cartesian.h"
#include "SimulationDataFormat/TrackReference.h"
#include "TFile.h"
#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TH1F.h"
#include "TRDBase/Calibrations.h"
#include "TRDBase/Geometry.h"
#include "TRDBase/TrackletTransformer.h"
#include "TRDQC/CoordinateTransformer.h"
#include "TRDSimulation/Detector.h"
#include "TTree.h"
#include "TVector3.h"
#include "VirtualMC.h"
#include <fairlogger/Logger.h>

#include "DetectorsCommonDataFormats/DetID.h"
#include "MathUtils/Cartesian.h"
#include "SimulationDataFormat/MCTrack.h"

#endif

// using namespace o2::detectors;
// using namespace o2::trd;

bool isDigitSaturated(o2::trd::Digit &digit) {
  auto adc = digit.getADC();
  int min = 10000;
  int max = -1000;
  for (int tb = 0; tb < 30; ++tb) {
    if (min > adc[tb])
      min = adc[tb];
    if (max < adc[tb])
      max = adc[tb];
  }
  if (max - min < 10) {
    return true;
  }
  return false;
}
bool isDigitSaturatedVal(o2::trd::Digit &digit, int min, int max) {
  auto adc = digit.getADC();
  int outofgap = 0;
  int ingap = 0;
  for (int tb = 0; tb < 30; ++tb) {
    if (adc[tb] <= max && adc[tb] >= min) {
      ingap++;
    } else {
      outofgap++;
    }
  }
  if (ingap == 30)
    return true;
  return false;
}
bool isFullMCMData(mcm) {
  for (auto &digit : mcm.digits) {
    auto adc = digit.getADC();
    for (int tb = 0; tb < 30; ++tb) {
      if (adc[tb] <= max && adc[tb] >= min) {
        ingap++;
      } else {
        outofgap++;
      }
    }

    return false;
  }
  /// This macro demonstrates how to use the o2::trd::RawDataManager and
  /// o2::trd:MCMDisplay to visualize TRD digits and tracklets.
  /// You probably want to copy it and adjust it for your use case.
  void ClassifyEvents(std::string dirname = ".") {
    // instantiate the class that handles all the data access
    auto dman = o2::trd::RawDataManager(dirname);
    cout << dman.describeFiles() << endl;

    TH1F *hSaturatedMCM = new TH1F("saturagemcm", "Saturation count per mcm.",
                                   70000, -0.5, 69998.5);
    TH1F *hTotalMCM =
        new TH1F("totalmcm", "Total per mcm.", 70000, -0.5, 69998.5);
    TH1F *hDifferenceMCM = new TH1F("saturagedifferencemcm",
                                    "Total counts - Saturation count per mcm.",
                                    70000, -0.5, 69998.5);
    // Set a number of plots you want for any MCM, not specified above
    int ndrawany = 50;
    int ndrawoffset = 50;
    int countdrawnmcm = 0;
    int drawstart = 50;
    int drawend = 75;
    int timeframecount = 0;
    int runningeventcount = -1;
    // --------------------------------------------------------------------
    // loop over timeframes
    while (dman.nextTimeFrame()) {
      // if(timeframecount>0 ) continue;
      timeframecount++;
      runningeventcount++;
      //  if(timeframecount<3 ) continue;
      //  if(timeframecount>4) break;
      cout << dman.describeTimeFrame() << endl;

      int eventcount = 0;
      // loop over events
      while (dman.nextEvent()) {
        //     if(eventcount>0) continue;
        eventcount++;
        auto ev = dman.getEvent();

        // skip events without digits immediately
        if (ev.digits.size() == 0) {
          continue;
        }
        // if (ev.digits.length() > 800000) { continue; }
        int saturationcount = 0;
        cout << dman.describeEvent() << endl;
        for (auto &mcm : dman.getEvent().iterateByMCM()) {
          bool isFullMCM = isFullMCMData(mcm);
          // skip MCMs without digits
          if (mcm.digits.size() == 0) {
            continue;
          }
          if (mcm.tracklets.size() == 0) {
            continue;
          }
          //        if (mcm.simtracklets.size() ==0){
          //          continue;
          //        }
          //        if(mcm.tracklets.size()!=1){
          //          continue;
          //        }
          //        if(mcm.simtracklets.size()!=1){
          //          continue;
          //        }
          // we skipped MCMs without digits, so we can use the first digit to
          // find out where we are
          auto firstdigit = *mcm.digits.begin();
          //        if(firstdigit.getDetector()==50) continue;
          array<int, 3> key = {firstdigit.getDetector(), firstdigit.getROB(),
                               firstdigit.getMCM()};

          bool isSaturated = false;
          for (auto &digit : mcm.digits) {
            bool isSaturated266 = isDigitSaturatedVal(digit, 264, 267);
            bool isSaturated522 = isDigitSaturatedVal(digit, 521, 524);
          }
          /*    if (isSaturated) {
                // LOGP(info,"Saturated {}: {} : {}
                //
             ",firstdigit.getDetector(),firstdigit.getROB(),firstdigit.getMCM())
                // ;
                hSaturatedMCM->Fill(o2::trd::HelperMethods::getMCMId(
                    firstdigit.getDetector(), firstdigit.getROB(),
                    firstdigit.getMCM()));
                hTotalMCM->Fill(o2::trd::HelperMethods::getMCMId(
                    firstdigit.getDetector(), firstdigit.getROB(),
                    firstdigit.getMCM()));
                saturationcount++;
              } else {
                hTotalMCM->Fill(o2::trd::HelperMethods::getMCMId(
                    firstdigit.getDetector(), firstdigit.getROB(),
                    firstdigit.getMCM()));
                hDifferenceMCM->Fill(o2::trd::HelperMethods::getMCMId(
                    firstdigit.getDetector(), firstdigit.getROB(),
                    firstdigit.getMCM()));
              }*/

          /*        if (ndrawany > 0) {
                    --ndrawany;
                  } else {
                    continue;
                  }*/

          //  cout <<
          //  "=============================================================="
          //       << endl;

          // the actual drawing
          if (/*isSaturated) { //} && */ firstdigit.getDetector() == 52 &&
              firstdigit.getROB() == 4 && firstdigit.getMCM() == 3) {
            countdrawnmcm++;
            if (countdrawnmcm > drawstart && countdrawnmcm < drawend) {
              LOGP(info,
                   "Drawing for countdrawmcm : {} with drawstart of {} drawend "
                   "{}",
                   countdrawnmcm, drawstart, drawend);
              std::string text = " ";
              o2::trd::MCMDisplay disp(mcm, runningeventcount, text);
              // disp.Draw();
              disp.drawDigits("colz");
              disp.drawDigits("text,same");
              disp.drawClusters();
              disp.drawTracklets();
              disp.drawSimTracklets();
            }
          }
        }

        runningeventcount++;
        LOGP(info, "Saturation count : {}", saturationcount);

      } // event/trigger record loop
    } // time frame loop
    TCanvas *c1 = new TCanvas();
    c1->cd();
    hSaturatedMCM->Draw();
    TCanvas *c2 = new TCanvas();
    c2->cd();
    hTotalMCM->Draw();
    TCanvas *c3 = new TCanvas();
    c3->cd();
    hDifferenceMCM->Draw();
  }
