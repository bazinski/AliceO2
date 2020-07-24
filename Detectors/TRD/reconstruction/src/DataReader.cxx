// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// @file   raw2tracklet.cxx
/// @author Sean Murray
/// @brief  Basic DPL workflow for TRD CRU output(raw) to tracklet data.
///         There may or may not be some compression in this at some point.

#include "TRDReconstruction/DataReaderTask.h"
#include "Framework/WorkflowSpec.h"
#include "Framework/ConfigParamSpec.h"
#include "Framework/ConcreteDataMatcher.h"
#include "Framework/Logger.h"
#include "DetectorsRaw/RDHUtils.h"

using namespace o2::framework;

// add workflow options, note that customization needs to be declared before
// including Framework/runDataProcessing
void customize(std::vector<ConfigParamSpec>& workflowOptions)
{
  auto config = ConfigParamSpec{"trd-datareader-config", VariantType::String, "A:TRD/RAWDATA", {"TRD raw data config"}};
  auto outputDesc = ConfigParamSpec{"trd-datareader-output-desc", VariantType::String, "TRDTLT", {"Output specs description string"}};
  auto verbosity = ConfigParamSpec{"trd-datareader-verbose", VariantType::Bool, false, {"Enable verbose epn data reading"}};
  auto headerverbosity = ConfigParamSpec{"trd-datareader-headerverbose", VariantType::Bool, false, {"Enable verbose header info"}};
  auto dataverbosity = ConfigParamSpec{"trd-datareader-dataverbose", VariantType::Bool, false, {"Enable verbose data info"}};
  auto compresseddata = ConfigParamSpec{"trd-datareader-compresseddata", VariantType::Bool, false, {"The incoming data is compressed or not"}};

  workflowOptions.push_back(config);
  workflowOptions.push_back(outputDesc);
  workflowOptions.push_back(verbosity);
  workflowOptions.push_back(headerverbosity);
  workflowOptions.push_back(dataverbosity);
}

#include "Framework/runDataProcessing.h" // the main driver

/// This function hooks up the the workflow specifications into the DPL driver.
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{

  auto config = cfgc.options().get<std::string>("trd-datareader-config");
  auto inputspec = cfgc.options().get<std::string>("trd-datareader-inputspec");
  auto verbosity = cfgc.options().get<bool>("trd-datareader-verbose");
  std::vector<OutputSpec> outputs;
  outputs.emplace_back(OutputSpec(ConcreteDataTypeMatcher{"TRD", "TRDTRACKLET"}));
  outputs.emplace_back(OutputSpec(ConcreteDataTypeMatcher{"TRD", "TRDDIGIT"}));
  //outputs.emplace_back(OutputSpec(ConcreteDataTypeMatcher{"TRD", "FLPSTAT"}));

  AlgorithmSpec algoSpec;
  algoSpec = AlgorithmSpec{adaptFromTask<o2::trd::DataReaderTask>()};

  WorkflowSpec workflow;

  /*
   * This is replicated from TOF
     We define at run time the number of devices to be attached
     to the workflow and the input matching string of the device.
     This is is done with a configuration string like the following
     one, where the input matching for each device is provide in
     comma-separated strings. For instance
  */

  std::stringstream ssconfig(config);
  std::string iconfig;
  std::string inputDescription;
  int idevice = 0;
  LOG(info) << " config string is : " << config;
  // this is probably never going to be used but would to nice to know hence here.
  workflow.emplace_back(DataProcessorSpec{
    std::string("trd-datareader"), // left as a string cast incase we append stuff to the string
    select(std::string("x:TRD/" + inputspec).c_str()),
    outputs,
    algoSpec,
    Options{
      {"trd-datareader-verbose", VariantType::Bool, false, {"verbose flag"}}}});

  return workflow;
}
