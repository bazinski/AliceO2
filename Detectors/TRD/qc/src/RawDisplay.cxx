// Copyright 2019-2023 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "TRDQC/RawDisplay.h"
#include "TRDQC/RawDataManager.h"
#include "DataFormatsTRD/Constants.h"
#include "DataFormatsTRD/HelperMethods.h"
#include "TRDBase/TrackletTransformer.h"

#include <TColor.h>
#include <TVirtualPad.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TLine.h>
#include <TMarker.h>

using namespace o2::trd;

namespace o2::trd
{

/// Modified version of o2::trd::Tracklet64::getPadCol returning a float
float PadColF(o2::trd::Tracklet64& tracklet)
{
  // obtain pad number relative to MCM center
  float padLocal = tracklet.getPositionBinSigned() * constants::GRANULARITYTRKLPOS;
  // MCM number in column direction (0..7)
  int mcmCol = (tracklet.getMCM() % constants::NMCMROBINCOL) + constants::NMCMROBINCOL * (tracklet.getROB() % 2);

  // original calculation
  // FIXME: understand why the offset seems to be 6 pads and not nChannels / 2 = 10.5
  // return CAMath::Round(6.f + mcmCol * ((float)constants::NCOLMCM) + padLocal);

  // my calculation
  return float((mcmCol + 1) * constants::NCOLMCM) + padLocal - 10.0;
}

}; // namespace o2::trd

RawDisplay::RawDisplay(RawDataSpan& dataspan, TVirtualPad* pad)
  : mDataSpan(dataspan), mPad(pad)
{
  /*
  std::string calvdexbfile = "~/alice/thesis/o2-trd-CalVdriftExB.root";
  TFile *vdriftexbf = TFile::Open(calvdexbfile.data());
  o2::trd::CalVdriftExB *calvdriftexb;
  vdriftexbf->GetObject("ccdb_object", calvdriftexb);
  mTransformer.init();
  mTransformer.setCalVdriftExB(calvdriftexb);
  LOGP(info,"VDrift and ExB");
  for(uint32_t i=0;i<constants::MAXCHAMBER;++i){
    LOGP(info,"{} : {} {} ", i, vdriftexbf->getVdrift(i), vdriftexbf->getExB(i));
  }
  */
}

MCMDisplay::MCMDisplay(RawDataSpan& mcmdata, int event, std::string text, TVirtualPad* pad)
  : RawDisplay(mcmdata, pad) // initializes mDataSpan, mPad
{
  int det = -1, rob = -1, mcm = -1;

  if (std::distance(mDataSpan.digits.begin(), mDataSpan.digits.end())) {
    auto x = *mDataSpan.digits.begin();
    det = x.getDetector();
    rob = x.getROB();
    mcm = x.getMCM();
  } else if (std::distance(mDataSpan.tracklets.begin(), mDataSpan.tracklets.end())) {
    auto x = *mDataSpan.tracklets.begin();
    det = x.getDetector();
    rob = x.getROB();
    mcm = x.getMCM();
  } else {
    O2ERROR("found neither digits nor tracklets in MCM");
    assert(false);
  }

  if (event != -1) {
    mName = Form("det%03d_rob%d_mcm%02d_e%d", det, rob, mcm, event);
    mDesc = Form("Detector %02d_%d_%d (%03d) - MCM %d:%02d e:%d %s", det / 30, (det % 30) / 6, det % 6, det, rob, mcm, event, text.c_str());
  } else {
    mName = Form("det%03d_rob%d_mcm%02d", det, rob, mcm);
    mDesc = Form("Detector %02d_%d_%d (%03d) - MCM %d:%02d %s", det / 30, (det % 30) / 6, det % 6, det, rob, mcm, event, text.c_str());
  }

  // MCM column number on ROC [0..7]
  int mcmcol = mcm % constants::NMCMROBINCOL + HelperMethods::getROBSide(rob) * constants::NMCMROBINCOL;

  mFirstPad = mcmcol * constants::NCOLMCM - 1;
  mLastPad = (mcmcol + 1) * constants::NCOLMCM + 2;

  if (pad == nullptr) {
    mPad = new TCanvas(mName.c_str(), mDesc.c_str(), 800, 600);
  } else {
    mPad = pad;
    mPad->SetName(mName.c_str());
    mPad->SetTitle(mDesc.c_str());
  }

  mDigitsHisto = new TH2F(mName.c_str(), (mDesc + ";pad;time bin").c_str(), (mLastPad - mFirstPad), mFirstPad, mLastPad, 30, 0., 30.);

  for (auto digit : mDataSpan.digits) {
    auto adc = digit.getADC();
    for (int tb = 0; tb < 30; ++tb) {
      mDigitsHisto->Fill(digit.getPadCol(), tb, adc[tb]);
    }
  }
  mDigitsHisto->SetStats(0);
}

void RawDisplay::drawDigits(std::string opt)
{
  mPad->cd();
  mDigitsHisto->Draw(opt.c_str());
}

void RawDisplay::drawTracklets()
{
  mPad->cd();

  TLine trkl;
  trkl.SetLineColor(kRed);
  trkl.SetLineWidth(3);

  for (auto tracklet : mDataSpan.tracklets) {
    auto pos = PadColF(tracklet);
    auto slope = -tracklet.getSlopeBinSigned() * constants::GRANULARITYTRKLSLOPE / constants::ADDBITSHIFTSLOPE;
    trkl.DrawLine(pos, 0, pos + 30 * slope, 30);
  }
}

void RawDisplay::drawShiftedTracklets()
{
  mPad->cd();

  TLine trkl;
  trkl.SetLineColor(kGreen);
  trkl.SetLineWidth(3);

  for (auto tracklet : mDataSpan.tracklets) {
    auto pos = PadColF(tracklet);
    auto slope = -tracklet.getSlopeBinSigned() * constants::GRANULARITYTRKLSLOPE / constants::ADDBITSHIFTSLOPE;
    trkl.DrawLine(pos + 1, 0, pos + 1 + 30 * slope, 30);
  }
}

void RawDisplay::drawCalibratedTracklets() // CalibratedTracklet& caltracklet)
{
  mPad->cd();

  TLine trkl;
  trkl.SetLineColor(kPink);
  trkl.SetLineWidth(3);

  for (auto tracklet : mDataSpan.tracklets) {
    auto pos = PadColF(tracklet);
    auto slope = -tracklet.getSlopeBinSigned() * constants::GRANULARITYTRKLSLOPE / constants::ADDBITSHIFTSLOPE;
    trkl.DrawLine(pos, 0, pos + 1 + 30 * slope, 30);
  }
}

void RawDisplay::drawSimTracklets()
{
  mPad->cd();

  TLine trkl;
  trkl.SetLineColor(kGreenRedViolet);
  trkl.SetLineWidth(3);

  for (auto tracklet : mDataSpan.simtracklets) {
    auto pos = PadColF(tracklet);
    auto slope = -tracklet.getSlopeBinSigned() * constants::GRANULARITYTRKLSLOPE / constants::ADDBITSHIFTSLOPE;
    trkl.DrawLine(pos, 0, pos + 30 * slope, 30);
  }
}

void RawDisplay::drawClusters()
{
  mPad->cd();

  drawDigits();

  TMarker clustermarker;
  clustermarker.SetMarkerColor(kRed);
  clustermarker.SetMarkerStyle(2);
  clustermarker.SetMarkerSize(1.5);

  TMarker cogmarker;
  cogmarker.SetMarkerColor(kGreen);
  cogmarker.SetMarkerStyle(3);
  cogmarker.SetMarkerSize(1.5);
  for (int t = 1; t <= mDigitsHisto->GetNbinsY(); ++t) {
    for (int p = 2; p <= mDigitsHisto->GetNbinsX() - 1; ++p) {
      // cout << p << "/" << t << " -> " << mDigitsHisto->GetBinContent(i,j) << endl;
      double baseline = 9.5;
      double left = mDigitsHisto->GetBinContent(p - 1, t) - baseline;
      double centre = mDigitsHisto->GetBinContent(p, t) - baseline;
      double right = mDigitsHisto->GetBinContent(p + 1, t) - baseline;
      if (centre > left && centre > right && (centre + left + right) > mClusterThreshold) {
        double pos = 0.5 * log(right / left) / log(centre * centre / left / right);
        double clpos = mDigitsHisto->GetXaxis()->GetBinCenter(p) + pos;
        double cog = (right - left) / (right + centre + left) + mDigitsHisto->GetXaxis()->GetBinCenter(p);
        // cout << "t=" << t << " p=" << p
        //      << ":   ADCs = " << left << " / " << centre << " / " << right
        //      << "   pos = " << pos << " ~ " << clpos
        //      << endl;
        clustermarker.DrawMarker(clpos, t - 0.5);
        // cogmarker.DrawMarker(cog, t - 0.5);
      }
    }
  }
}


void RawDisplay::drawSimHits()
{
  mPad->cd();

  TMarker clustermarker;
  clustermarker.SetMarkerColor(kGreen);
  clustermarker.SetMarkerStyle(2);
  clustermarker.SetMarkerSize(1.5);
  int TPVBY=0;
  int TPVT=313;
  int TPHT=313;
  int TPFP=25;
  std::array<int,128> TPL00;
  TMarker cogmarker;
  cogmarker.SetMarkerColor(kGreen);
  cogmarker.SetMarkerStyle(3);
  cogmarker.SetMarkerSize(1.5);
  for (unsigned int timebin = 1; timebin < 24; ++timebin) {
    // first find the hit candidates and store the total cluster charge in qTotal array
    // in case of not hit store 0 there.
    std::array<unsigned short, 20> qTotal{}; //[19 + 1]; // the last is dummy
    int adcLeft, adcCentral, adcRight;
    for (int adcch = 2; adcch < 21- 2; adcch++) {
      adcLeft =((int)mDigitsHisto->GetBinContent(adcch - 1, timebin))<<2;
      adcCentral = ((int)mDigitsHisto->GetBinContent(adcch,timebin))<<2;
      adcRight = ((int)mDigitsHisto->GetBinContent(adcch + 1,timebin))<<2;
      bool hitQual = false;
      if (TPVBY == 0) {
        // bypass the cluster verification
        hitQual = true;
      } else {
        hitQual = ((adcLeft * adcRight) <
                   ((TPVT * adcCentral * adcCentral) >> 10));
        if (hitQual) {
          LOG(info) << "clus q cut " << adcLeft << ", " << adcCentral << ", "
                    << adcRight << " - th:" << TPVT
                    << " -> " << TPVT;
        }
      }
      hitQual=true;
      // The accumulated charge is with the pedestal!!!
      int qtotTemp = adcLeft + adcCentral + adcRight;
      //LOGP(info,"qtotTemp {} L {} C {} R {} ",qtotTemp,adcLeft,adcCentral,adcRight);
      if ((hitQual) &&
          (qtotTemp >= TPHT) &&
          (adcLeft <= adcCentral) &&
          (adcCentral > adcRight)) {
        qTotal[adcch] = qtotTemp;
       LOGP(info," adding qtotTemp {} to new qTotal[{}]={}",qtotTemp,adcch,qTotal[adcch]);
      } else {
        qTotal[adcch] = 0;
       // LOGP(info,"!hitQual, setting qTotal[{}]=0",adcch);
      }
    }

    short fromLeft = -1;
    short adcch = 0;
    short found = 0;
    std::array<unsigned short, 6> marked{};
    marked[4] = 19; // invalid channel
    marked[5] = 19; // invalid channel
    qTotal[19] = 0;
    while ((adcch < 16) && (found < 3)) {
      if (qTotal[adcch] > 0) {
        fromLeft = adcch;
        marked[2 * found + 1] = adcch;
        found++;
      }
      adcch++;
    }

    short fromRight = -1;
    adcch = 18;
    found = 0;
    while ((adcch > 2) && (found < 3)) {
      if (qTotal[adcch] > 0) {
        marked[2 * found] = adcch;
        found++;
        fromRight = adcch;
      }
      adcch--;
    }

    // here mask the hit candidates in the middle, if any
    LOGP(info,"Fromleft={}, Fromright={}\n",fromLeft, fromRight);
    if ((fromLeft >= 0) && (fromRight >= 0) && (fromLeft < fromRight)) {
      for (adcch = fromLeft + 1; adcch < fromRight; adcch++) {
        qTotal[adcch] = 0;
      }
    }

    found = 0;
    for (adcch = 0; adcch < 19; adcch++) {
      if (qTotal[adcch] > 0) {
        found++;
        // NOT READY
      }
    }
    if (found > 4) // sorting like in the TRAP in case of 5 or 6 candidates!
    {
      if (marked[4] == marked[5]) {
        marked[5] = 19;
      }
      std::array<unsigned short, 6> qMarked;
      for (found = 0; found < 6; found++) {
        qMarked[found] = qTotal[marked[found]] >> 4;
        LOGP(info,"ch_{} Qtot {} Qtots {} |",marked[found],qTotal[marked[found]],qMarked[found]);
      }

      unsigned short worse1, worse2;
      sort6To2Worst(marked[0], marked[3], marked[4], marked[1], marked[2], marked[5],
                    qMarked[0],
                    qMarked[3],
                    qMarked[4],
                    qMarked[1],
                    qMarked[2],
                    qMarked[5],
                    &worse1, &worse2);
      // Now mask the two channels with the smallest charge
      if (worse1 < 19) {
        qTotal[worse1] = 0;
      }
      if (worse2 < 19) {
        qTotal[worse2] = 0;
      }
    }

    for (adcch = 0; adcch < 19; adcch++) {
      if (qTotal[adcch] > 0) // the channel is marked for processing
      {
      adcLeft =((int)mDigitsHisto->GetBinContent(adcch , timebin))<<2;
      adcCentral = ((int)mDigitsHisto->GetBinContent(adcch+1,timebin))<<2;
      adcRight = ((int)mDigitsHisto->GetBinContent(adcch + 2,timebin))<<2;
        LOGF(info, "ch(%i): left(%i), central(%i), right(%i)", adcch, adcLeft, adcCentral, adcRight);
        //  hit detected, in TRAP we have 4 units and a hit-selection, here we proceed all channels!
        //  subtract the pedestal TPFP, clipping instead of wrapping

        int regTPFP = TPFP; // TODO put this together with the others as members of trapsim, which is initiliased by det,rob,mcm.
        LOG(info) << "YY Hit found, time=" << timebin << ", adcch=" << adcch << "/" << adcch + 1 << "/"
                  << adcch + 2 << ", adc values=" << adcLeft << "/" << adcCentral << "/"
                  << adcRight << ", regTPFP=" << regTPFP << ", TPHT=" << TPHT;
        // regTPFP >>= 2; // OS: this line should be commented out when checking real data. It's only needed for comparison with Venelin's simulation if in addition mgkAddDigits == 0
        if (adcLeft < regTPFP) {
          adcLeft = 0;
        } else {
          adcLeft -= regTPFP;
        }
        if (adcCentral < regTPFP) {
          adcCentral = 0;
        } else {
          adcCentral -= regTPFP;
        }
        if (adcRight < regTPFP) {
          adcRight = 0;
        } else {
          adcRight -= regTPFP;
        }

        // Calculate the center of gravity
        // checking for adcCentral != 0 (in case of "bad" configuration)
        if (adcCentral == 0) {
          LOGP(info, " bad configuration detected adcCentral={} adcRight:{} adcLeft:{}  regTPFP:{}", adcCentral, adcRight, adcLeft, TPFP);
          continue;
        }
        int ypos = (adcRight - adcLeft) / adcCentral;
        if (ypos < 0) {
          ypos = -ypos;
        }
        // ypos element of [0:128]
        //  make the correction using the position LUT
        // LOG(info) << "ypos raw is " << ypos << "  adcrigh-adcleft/adccentral " << adcRight << "-" << adcLeft << "/" << adcCentral << "==" << (adcRight - adcLeft) / adcCentral << " 128 * numerator : " << 128 * (adcRight - adcLeft) / adcCentral;
        LOGP(info,"ypos before lut correction : {}   adcright-adcleft : {}-{}={}, adcentral : {}", ypos, adcRight,adcLeft,adcRight-adcLeft,adcCentral);
        //ypos = ypos + TPL00[(ypos & 0x7F)];
        // ypos += LUT_POS[ypos & 0x7f]; // FIXME use this LUT to obtain the same results as Venelin
        LOGP(info,"ypos after lut correction : {}",ypos);
        //   LOG(info) << "ypos after lut correction : " << ypos;
        if (adcLeft > adcRight) {
          ypos = -ypos;
        }
        LOGP(info,"ypos after leftright correction : {}", ypos);
        LOGP(info, "XX Add hit ch({}): left({}), central({}), right({}), ypos({}) qtot= {}, timebin:{}", adcch, adcLeft, adcCentral, adcRight, ypos, qTotal[adcch],timebin);
        //now draw the hit as if we would have added it to the fitreg in the trapsimualtor 
        //double clpos = ypos;
        //double clpos = mDigitsHisto->GetXaxis()->GetBinCenter(adcch) + ypos;
        //double cog = (right - left) / (right + centre + left) + mDigitsHisto->GetXaxis()->GetBinCenter(p);
        double pos=0;
        if(adcRight!=0 && adcLeft !=0 && adcCentral!=0) pos = (adcRight- adcLeft) / log(adcCentral + adcLeft + adcRight);
        //if(adcRight!=0 && adcLeft !=0 && adcCentral!=0) pos = 0.5 * log(adcRight / adcLeft) / log(adcCentral * adcCentral / adcLeft / adcRight);
        double clpos = mDigitsHisto->GetXaxis()->GetBinCenter(adcch) + pos;
        //double cog = (right - left) / (right + centre + left) + mDigitsHisto->GetXaxis()->GetBinCenter(p);
        LOGP(info, "XX Drawing {} {}", clpos,timebin-0.5);
        clustermarker.DrawMarker(adcch+clpos, timebin - 0.5); //0.5 to put it in the middle of the timebin
      }
    }
  }
}


void RawDisplay::drawHits()
{
  TMarker hitmarker;
  hitmarker.SetMarkerColor(kBlue);
  hitmarker.SetMarkerStyle(38);
  for (auto hit : mDataSpan.hits) {
    if (hit.getCharge() > 0.0) {
      hitmarker.SetMarkerSize(log10(hit.getCharge()));
      hitmarker.DrawMarker(hit.getPadCol(), hit.getTimeBin());
    }
  }
}

void RawDisplay::drawMCTrackSegments()
{
  TLine line;
  line.SetLineColor(kBlue);
  line.SetLineWidth(2.0);

  for (auto& trkl : mDataSpan.makeMCTrackSegments()) {
    line.DrawLine(trkl.getStartPoint().getPadCol(), trkl.getStartPoint().getTimeBin(), trkl.getEndPoint().getPadCol(), trkl.getEndPoint().getTimeBin());
  }
}


void RawDisplay::sort2(uint16_t idx1i, uint16_t idx2i,
                          uint16_t val1i, uint16_t val2i,
                          uint16_t* idx1o, uint16_t* idx2o,
                          uint16_t* val1o, uint16_t* val2o) const
{
  // sorting for tracklet selection

  if (val1i > val2i) {
    *idx1o = idx1i;
    *idx2o = idx2i;
    *val1o = val1i;
    *val2o = val2i;
  } else {
    *idx1o = idx2i;
    *idx2o = idx1i;
    *val1o = val2i;
    *val2o = val1i;
  }
}

void RawDisplay::sort3(uint16_t idx1i, uint16_t idx2i, uint16_t idx3i,
                          uint16_t val1i, uint16_t val2i, uint16_t val3i,
                          uint16_t* idx1o, uint16_t* idx2o, uint16_t* idx3o,
                          uint16_t* val1o, uint16_t* val2o, uint16_t* val3o) const
{
  // sorting for tracklet selection

  int sel;

  if (val1i > val2i) {
    sel = 4;
  } else {
    sel = 0;
  }
  if (val2i > val3i) {
    sel = sel + 2;
  }
  if (val3i > val1i) {
    sel = sel + 1;
  }
  switch (sel) {
    case 6: // 1 >  2  >  3            => 1 2 3
    case 0: // 1 =  2  =  3            => 1 2 3 : in this case doesn't matter, but so is in hardware!
      *idx1o = idx1i;
      *idx2o = idx2i;
      *idx3o = idx3i;
      *val1o = val1i;
      *val2o = val2i;
      *val3o = val3i;
      break;

    case 4: // 1 >  2, 2 <= 3, 3 <= 1  => 1 3 2
      *idx1o = idx1i;
      *idx2o = idx3i;
      *idx3o = idx2i;
      *val1o = val1i;
      *val2o = val3i;
      *val3o = val2i;
      break;

    case 2: // 1 <= 2, 2 > 3, 3 <= 1   => 2 1 3
      *idx1o = idx2i;
      *idx2o = idx1i;
      *idx3o = idx3i;
      *val1o = val2i;
      *val2o = val1i;
      *val3o = val3i;
      break;

    case 3: // 1 <= 2, 2 > 3, 3  > 1   => 2 3 1
      *idx1o = idx2i;
      *idx2o = idx3i;
      *idx3o = idx1i;
      *val1o = val2i;
      *val2o = val3i;
      *val3o = val1i;
      break;

    case 1: // 1 <= 2, 2 <= 3, 3 > 1   => 3 2 1
      *idx1o = idx3i;
      *idx2o = idx2i;
      *idx3o = idx1i;
      *val1o = val3i;
      *val2o = val2i;
      *val3o = val1i;
      break;

    case 5: // 1 > 2, 2 <= 3, 3 >  1   => 3 1 2
      *idx1o = idx3i;
      *idx2o = idx1i;
      *idx3o = idx2i;
      *val1o = val3i;
      *val2o = val1i;
      *val3o = val2i;
      break;

    default: // the rest should NEVER happen!
      LOG(error) << "ERROR in Sort3!!!";
      break;
  }
}

void RawDisplay::sort6To4(uint16_t idx1i, uint16_t idx2i, uint16_t idx3i, uint16_t idx4i, uint16_t idx5i, uint16_t idx6i,
                             uint16_t val1i, uint16_t val2i, uint16_t val3i, uint16_t val4i, uint16_t val5i, uint16_t val6i,
                             uint16_t* idx1o, uint16_t* idx2o, uint16_t* idx3o, uint16_t* idx4o,
                             uint16_t* val1o, uint16_t* val2o, uint16_t* val3o, uint16_t* val4o) const
{
  // sorting for tracklet selection

  uint16_t idx21s, idx22s, idx23s, dummy;
  uint16_t val21s, val22s, val23s;
  uint16_t idx23as, idx23bs;
  uint16_t val23as, val23bs;

  sort3(idx1i, idx2i, idx3i, val1i, val2i, val3i,
        idx1o, &idx21s, &idx23as,
        val1o, &val21s, &val23as);

  sort3(idx4i, idx5i, idx6i, val4i, val5i, val6i,
        idx2o, &idx22s, &idx23bs,
        val2o, &val22s, &val23bs);

  sort2(idx23as, idx23bs, val23as, val23bs, &idx23s, &dummy, &val23s, &dummy);

  sort3(idx21s, idx22s, idx23s, val21s, val22s, val23s,
        idx3o, idx4o, &dummy,
        val3o, val4o, &dummy);
}

void RawDisplay::sort6To2Worst(uint16_t idx1i, uint16_t idx2i, uint16_t idx3i, uint16_t idx4i, uint16_t idx5i, uint16_t idx6i,
                                  uint16_t val1i, uint16_t val2i, uint16_t val3i, uint16_t val4i, uint16_t val5i, uint16_t val6i,
                                  uint16_t* idx5o, uint16_t* idx6o) const
{
  // sorting for tracklet selection

  uint16_t idx21s, idx22s, idx23s, dummy1, dummy2, dummy3, dummy4, dummy5;
  uint16_t val21s, val22s, val23s;
  uint16_t idx23as, idx23bs;
  uint16_t val23as, val23bs;

  sort3(idx1i, idx2i, idx3i, val1i, val2i, val3i,
        &dummy1, &idx21s, &idx23as,
        &dummy2, &val21s, &val23as);

  sort3(idx4i, idx5i, idx6i, val4i, val5i, val6i,
        &dummy1, &idx22s, &idx23bs,
        &dummy2, &val22s, &val23bs);

  sort2(idx23as, idx23bs, val23as, val23bs, &idx23s, idx5o, &val23s, &dummy1);

  sort3(idx21s, idx22s, idx23s, val21s, val22s, val23s,
        &dummy1, &dummy2, idx6o,
        &dummy3, &dummy4, &dummy5);
}
