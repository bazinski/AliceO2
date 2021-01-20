// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// @file   CruCompressorTask.cxx
/// @author Sean Murray
/// @brief  TRD cru output to tracklet task

#include "TRDRaw/CruCompressorTask.h"
#include "TRDRaw/CruRawReader.h"
#include "Framework/ControlService.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/RawDeviceService.h"
#include "Framework/DeviceSpec.h"
#include "Framework/DataSpecUtils.h"
#include "DataFormatsTRD/Constants.h"
#include "DataFormatsTRD/CompressedHeader.h"
#include <fairmq/FairMQDevice.h>

using namespace o2::framework;

namespace o2
{
namespace trd
{

void CruCompressorTask::init(InitContext& ic)
{
  LOG(INFO) << "FLP Compressore Task init";

  auto finishFunction = [this]() {
    mReader.checkSummary();
  };

  ic.services().get<CallbackService>().set(CallbackService::Id::Stop, finishFunction);
}

uint64_t CruCompressorTask::buildBlobOutput()
{
    //mReader holds the vectors of tracklets and digits.
    // tracklets are 64 bit
    // digits are DigitMCMHeader and DigitMCMData
    
    uint64_t currentpos=0;
    uint64_t trailer = 0xeeeeeeeeeeeeeeeeLL;
    CompressedRawHeader *trackletheader=(CompressedRawHeader*)&mOutBuffer[0];
    trackletheader->format = 1;
    trackletheader->eventtime=99;
    currentpos=1;
    trackletheader->size = mReader.getTracklets().size()*8; // to get to bytes.
    for(auto tracklet : mReader.getTracklets()){
        //convert tracklet to 64 bit and add to blob
        mOutBuffer[currentpos++]=tracklet.getTrackletWord();
    }
    CompressedRawTrailer *tracklettrailer=(CompressedRawTrailer*)&mOutBuffer[currentpos];
    tracklettrailer->word=trailer;
    currentpos++;
    CompressedRawHeader *digitheader=(CompressedRawHeader*)&mOutBuffer[currentpos];
    currentpos++;
    digitheader->format = 2;
    digitheader->eventtime=99;
 
    for(auto digit : mReader.getDigits()){
        uint64_t *word;
        word= &mOutBuffer[currentpos];
        DigitMCMHeader mcmheader;
        mcmheader.eventcount=1;
        mcmheader.mcm=digit.getMCM();
        mcmheader.rob=digit.getROB();
        mcmheader.yearflag=1;
        //unpack the digit.
        mcmheader.word=(*word)>>32;
        currentpos++;
       
    }
    CompressedRawTrailer *digittrailer=(CompressedRawTrailer*)&mOutBuffer[currentpos];
    digittrailer->word=trailer;
    currentpos++;
    //as far as I can tell this is almost always going to be blank.
    CompressedRawHeader *configheader=(CompressedRawHeader*)&mOutBuffer[currentpos];
    currentpos++;
    configheader->size=2;
    configheader->format=3;
    configheader->eventtime=99;
    CompressedRawTrailer *configtrailer=(CompressedRawTrailer*)&mOutBuffer[currentpos];
    configtrailer->word=trailer;
return currentpos;
}

void CruCompressorTask::run(ProcessingContext& pc)
{
  LOG(info) << "TRD Compression Task run method";

  auto device = pc.services().get<o2::framework::RawDeviceService>().device();
  auto outputRoutes = pc.services().get<o2::framework::RawDeviceService>().spec().outputs;
  auto fairMQChannel = outputRoutes.at(0).channel;
  int inputcount = 0;
  /* loop over inputs routes */
  for (auto iit = pc.inputs().begin(), iend = pc.inputs().end(); iit != iend; ++iit) {
    if (!iit.isValid())
      continue;
    //LOG(info) << "iit.mInputs  " << iit.mInputs.
    /* prepare output parts */
    FairMQParts parts;

    /* loop over input parts */
    for (auto const& ref : iit) {

      auto headerIn = DataRefUtils::getHeader<o2::header::DataHeader*>(ref);
      auto dataProcessingHeaderIn = DataRefUtils::getHeader<o2::framework::DataProcessingHeader*>(ref);
      auto payloadIn = ref.payload;
      auto payloadInSize = headerIn->payloadSize;
      mReader.setDataBuffer(payloadIn);
      mReader.setDataBufferSize(payloadInSize);
//      mReader.setVerbosity(mVerbose);
      /* run */
      mReader.run();
      
      auto payloadOutSize = mReader.buildBlobOutput();
      auto payloadMessage = device->NewMessage(payloadOutSize);
      std::memcpy(payloadMessage->GetData(), (char*)bufferOut, payloadOutSize*8);

      /* output */
      auto headerOut = *headerIn;
      auto dataProcessingHeaderOut = *dataProcessingHeaderIn;
      headerOut.dataDescription = "TRDRAW";
      headerOut.payloadSize = payloadOutSize;
      o2::header::Stack headerStack{headerOut, dataProcessingHeaderOut};
      auto headerMessage = device->NewMessage(headerStack.size());
      std::memcpy(headerMessage->GetData(), headerStack.data(), headerStack.size());

      /* add parts */
      parts.AddPart(std::move(headerMessage));
      parts.AddPart(std::move(payloadMessage));
    }

    /* send message */
    device->Send(parts, fairMQChannel);
  }
}

} // namespace trd
} // namespace o2
