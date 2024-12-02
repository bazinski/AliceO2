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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD Standalone trap simulator                                            //
//  Take a single MCM ADCvstimebin                                           //
//  Now also digits.
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "fairlogger/Logger.h"

#include "CommonUtils/StringUtils.h"
#include "CommonUtils/ConfigurableParam.h"
#include "CommonUtils/NameConf.h"
#include "DetectorsRaw/HBFUtils.h"
#include "TRDSimulation/TrapSimulator.h"
#include "DataFormatsTRD/Constants.h"
#include "DataFormatsTRD/Digit.h"
#include "DataFormatsTRD/RawData.h"
#include "DataFormatsTRD/Tracklet64.h"
#include "DataFormatsTRD/TrapConfigEvent.h"

#include <boost/program_options.hpp>
#include <filesystem>
#include <TStopwatch.h>
#include <string>
#include <iomanip>
#include <iostream>
#include <iostream>
#include <fstream>
#include <string>
#include <TFile.h>

using namespace std;

namespace bpo = boost::program_options;

void runsim(const std::string& inpDigitsName, int verbosity );

int main(int argc, char** argv)
{
  bpo::variables_map vm;
  bpo::options_description opt_general("Usage:\n  " + std::string(argv[0]) +
                                       "Trap simulator on a given mcm adcvstimebin array");
  bpo::options_description opt_hidden("");
  bpo::options_description opt_all;
  bpo::positional_options_description opt_pos;

  try {
    auto add_option = opt_general.add_options();
    add_option("help,h", "Print this help message");
    add_option("verbosity,v", bpo::value<int>()->default_value(0), "verbosity level");
    add_option("input-file-digits,d", bpo::value<std::string>()->default_value("rawmcm.dat"), "input digits");
    add_option("output-dir,o", bpo::value<std::string>()->default_value("./"), "output directory for raw data for corresponding tracklets");

    opt_all.add(opt_general).add(opt_hidden);
    bpo::store(bpo::command_line_parser(argc, argv).options(opt_all).positional(opt_pos).run(), vm);

    if (vm.count("help")) {
      std::cout << opt_general << std::endl;
      exit(0);
    }

    bpo::notify(vm);
  } catch (bpo::error& e) {
    std::cerr << "ERROR: " << e.what() << std::endl
              << std::endl;
    std::cerr << opt_general << std::endl;
    exit(1);
  } catch (std::exception& e) {
    std::cerr << e.what() << ", application will now exit" << std::endl;
    exit(2);
  }

  //runsim(vm["input-file-digits"].as<std::string>(), vm["input-file-tracklets"].as<std::string>(), vm["output-dir"].as<std::string>(), vm["verbosity"].as<int>(), vm["file-for"].as<std::string>(), vm["rdh-version"].as<uint32_t>(), vm["no-empty-hbf"].as<bool>(), vm["tracklethcheader"].as<int>(), 1024 * 1024, hbfu.startTime);
  runsim(vm["input-file-digits"].as<std::string>(), vm["verbosity"].as<int>());

  return 0;
}

void runsim(const std::string& inpDigitsName, int verbosity )
{
  TStopwatch swTot;
  swTot.Start();
  LOG(info) << "timer started";
  std::array<o2::trd::ArrayADC,30> mcmdata; 
  o2::trd::TrapSimulator trapsim;
  std::unique_ptr<TFile> file(TFile::Open("trdconfigevents.root"));
  if (!file || file->IsZombie()) {
    std::cerr << "Error opening trdconfigevent file" << std::endl;
    exit(-1);
  }
  o2::trd::TrapConfigEvent *mTrapConfigEvent = file->Get<o2::trd::TrapConfigEvent>("ccdb_object");
  std::unique_ptr<o2::trd::TrapConfigEvent> configevent(file->Get<o2::trd::TrapConfigEvent>("ccdb_object"));
  if (mTrapConfigEvent) {
    LOGP(debug, " we have a valid trapconfigevent object");
  } else {
    LOGP(debug, " we have a invalid trapconfigevent object");
  }
  trapsim.init(configevent.get(),1,1,1);

//void TrapSimulator::init(TrapConfigEvent* trapconfigevent, int det, int robPos, int mcmPos)
   ifstream inpfile(inpDigitsName.c_str());
    string content;
   int adccounter=0;
   int tb=0;
    while(inpfile >> content) {
        cout << content << ' ';
        mcmdata[adccounter][tb]=std::stoi(content);
        if(adccounter>19) {
          cout << std::endl;
        adccounter=0;
        tb++;
        }
        else adccounter++;
    }
    for(int adc=0;adc<21;adc++){
    trapsim.setData(adc,mcmdata[adc],adc);
    }
    std::cout << trapsim;
  
    if (!trapsim.isDataSet()) {
      LOGP(error,"trap sim data is not set, leaving");
      return;
    }
    std::cout << trapsim << std::endl;
    trapsim.setBaselines();
    trapsim.filter();
    trapsim.tracklet();
    auto trackletsOut = trapsim.getTrackletArray64();
    LOGP(info ,"Received {} tracklets from simlator",trackletsOut.size());
    for(auto &trklt : trackletsOut){
      std::cout << trklt << std::endl;
    }
    LOGP(info ,"end of Received {} tracklets from simlator",trackletsOut.size());
  trapsim.draw(0xf,1);  
  swTot.Stop();
  swTot.Print();
}
