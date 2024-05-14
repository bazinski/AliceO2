/// \file FindTrackRefIntersection.C
/// \brief Macro to find the intersection of a trackref with a trd chamber

#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <array>
#include <cmath>
#include <fstream>
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
#include "TRDBase/Calibrations.h"
#include "TRDBase/Geometry.h"
#include "TRDBase/TrackletTransformer.h"
#include "TRDQC/CoordinateTransformer.h"
#include "TRDSimulation/Detector.h"
#include "VirtualMC.h"
#include "TFile.h"
#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TH1F.h"
#include "TTree.h"
#include "TVector3.h"
#include <fairlogger/Logger.h>

#include "DetectorsCommonDataFormats/DetID.h"
#include "MathUtils/Cartesian.h"
#include "SimulationDataFormat/MCTrack.h"

#endif

using namespace o2::detectors;

/*std::vector<float> transformtracklet(o2::trd::Tracklet64 tracklet) {
  // this is ported kind of from tracklettransformer but that needs ccdb
  // connections and other stuff.
}*/
class MCTracklet {
public:
  MCTracklet(o2::TrackReference &in, o2::TrackReference &out, int trackid,
             int det)
      : mEnter(in), mExit(out), mMCTrackId(trackid), det(det) {
    double xdiff = mEnter.X() - mExit.X();
    double ydiff = mEnter.Y() - mExit.Y();
    double zdiff = mEnter.Z() - mExit.Z();
    mDistanceXY = std::sqrt(xdiff * xdiff + ydiff * ydiff);
    mDistanceXYZ =
        std::sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff);
    mDistanceXY = std::sqrt(xdiff * xdiff + ydiff * ydiff);
    mDistanceX = std::sqrt(xdiff * xdiff);
    mDistanceY = std::sqrt(ydiff * ydiff);
    mDistanceZ = std::sqrt(zdiff * zdiff);
  }
  o2::TrackReference mEnter;
  o2::TrackReference mExit;
  int mMCTrackId;
  int det;
  float mDistanceX, mDistanceY, mDistanceXY, mDistanceXYZ, mDistanceZ;
};

class MCPoint {
public:
  MCPoint(o2::math_utils::Point3D<float> pos, int id) : mId(id), mPos(pos){};
  int mId;
  o2::math_utils::Point3D<float> mPos;
};


class RefTracklet{
public : 
  RefTracklet(){};
  o2::trd::Tracklet64 mTracklet;
  o2::trd::CalibratedTracklet mCalibratedTracklet;
  float mdY;
  float mdCalY;
  uint32_t mHCID;
  std::array<float,3> mXYZ;
  std::array<float,3> mTransformedXYZ;
  std::array<double,3> mRotatedXYZ;
  void transformL2T(std::array<float,3> &transformedxyz){mTransformedXYZ=transformedxyz;}
  void rotateL2G(std::array<double,3> &rotatedxyz){mRotatedXYZ=rotatedxyz;}
};


std::array<double, 3> rotateL2G(std::array<double, 3> trdpoints, int sector) {

  float phi = 2.0 * TMath::Pi() / (float)18 * ((float)sector + 0.5);
  std::array<double, 3> alicepoints;
  alicepoints[0] = trdpoints[0] * std::cos(phi) - trdpoints[1] * std::sin(phi);
  alicepoints[1] = trdpoints[0] * std::sin(phi) + trdpoints[1] * std::cos(phi);
  alicepoints[2] = trdpoints[2];
  // std::cout << "point changed from : " << trdpoints[0] << ":" << trdpoints[1]
  //          << ":" << trdpoints[2] << " :->: " << alicepoints[0] << ":"
  //         << alicepoints[1] << ":" << alicepoints[2] << std::endl;
  return alicepoints;
}

std::array<double, 3> rotateL2G(double x, double y, double z, int sector) {
  std::array<double, 3> point;
  point[0] = x;
  point[1] = y;
  point[2] = z;
  return rotateL2G(point, sector);
}

std::array<double, 3> rotateL2G(std::array<float, 3> trdpoints, int sector) {
  return rotateL2G(trdpoints[0], trdpoints[1], trdpoints[2], sector);
}

bool comparetrackrefs(o2::TrackReference &i, o2::TrackReference &j) {

//  LOGP(info, "{} refi:{} refj:{} itime {} jtime{}", __LINE__, i.getTrackID(),
 //      j.getTrackID(), i.getTime(), j.getTime());
  /* if (i.getTime() < j.getTime()) {
     return true;
   }*/
  /* if (i.getTrackID() < j.getTrackID()) {
     LOGP(info, "{} trackidi < trackidj", __LINE__);
     return true;
   }*/
  int deti = i.getUserId() >> 2;
  int detj = j.getUserId() >> 2;
  if (deti < detj)
    return true;
  //LOGP(info, "{} trackidi < trackidj", __LINE__);
  //LOGP(info, "{} trackidi {} < trackidj {}", __LINE__, i.getTrackID(),
    //   j.getTrackID());
  if (i.getTrackID() < j.getTrackID()) {
    LOGP(info, "{} trackidi < trackidj", __LINE__);
    return true;
  }
  return false;
}

void PlotTrackRef(
    const int trackid = 0,
    const int detector = 114, // 340, // 29, // 247,//, 242, 50
    std::string digifile = "trddigits.root",
    std::string trackletfile = "trdtracklets.root",
    std::string hitfile = "o2sim_HitsTRD.root",
    std::string inputGeom = "o2sim_geometry.root",
    std::string paramfile = "o2sim_par.root",
    std::string calvdexbfile = "~/alice/thesis/o2-trd-CalVdriftExB.root") {
  // o2::trd::Calibrations calib;
  // calib.getCCDBObjects(297595);
  // autoo2::trd::Geometry *mgeo;
  auto mgeo = o2::trd::Geometry::instance();
  // mgeo->instance();

  // const std::string& volMapFile = "MCStepLoggerVolMap.dat",
  TGeoManager::Import("o2sim_geometry.root"); // inputGeom.c_str());
  auto nDet = o2::detectors::DetID::getNDetectors();
  std::vector<std::string> detectorNames(nDet);

  TFile *fin = TFile::Open(hitfile.data());
  TTree *hitTree = (TTree *)fin->Get("o2sim");
  std::vector<o2::trd::Hit> *hits = nullptr;
  hitTree->SetBranchAddress("TRDHit", &hits);
  int nev = hitTree->GetEntries();

  std::vector<o2::TrackReference> *trackrefs = nullptr;
  TFile *ftrackref = TFile::Open("o2sim_Kine.root");
  TTree *trackrefTree = (TTree *)ftrackref->Get("o2sim");
  trackrefTree->SetBranchAddress("TrackRefs", &trackrefs);
  int nevtrackref = trackrefTree->GetEntries();

  TFile *dfin = TFile::Open(digifile.data());
  TTree *digitTree = (TTree *)dfin->Get("o2sim");
  std::vector<o2::trd::Digit> *digits = nullptr;
  digitTree->SetBranchAddress("TRDDigit", &digits);
  int ndigitev = digitTree->GetEntries();

  TFile *trackletfin = TFile::Open(trackletfile.data());
  TTree *trackletTree = (TTree *)trackletfin->Get("o2sim");
  std::vector<o2::trd::Tracklet64> *tracklets = nullptr;
  trackletTree->SetBranchAddress("Tracklet", &tracklets);
  std::vector<o2::trd::TriggerRecord> *trigRecs = nullptr;
  trackletTree->SetBranchAddress("TrackTrg", &trigRecs);
  int ntrackletev = trackletTree->GetEntries();

  TFile *vdriftexbf = TFile::Open(calvdexbfile.data());
  o2::trd::CalVdriftExB *calvdriftexb;
  vdriftexbf->GetObject("ccdb_object", calvdriftexb);
  o2::trd::TrackletTransformer transformer;
  transformer.init();
  transformer.setCalVdriftExB(calvdriftexb);

  TH1F *hdistances[4];
  hdistances[0] = new TH1F(
      "distances_0",
      "/distance from trackref to tracklet radial;radial distance;count", 200,
      -0.5, 19.5);
  hdistances[1] = new TH1F(
      "distances_1",
      "distance from trackref to tracklet x-distance; x distance;count", 200,
      -0.5, 19.5);
  hdistances[2] = new TH1F(
      "distances_2",
      "distance from trackref to tracklet y-distance; y distance;count", 200,
      -0.5, 19.5);
  hdistances[3] = new TH1F(
      "distances_3",
      "distance from trackref to tracklet z-distance; z distance;count", 200,
      -0.5, 19.5);
  TH1F *trackletangles = new TH1F(
      "trackletangles", "Tracklet Deflection; dY; count", 10000, -5, 5);
  TH1F *hdistance[5];
  for (int hd = 0; hd < 5; ++hd) {
    hdistance[hd] = new TH1F(Form("distance%d", hd),
                             Form("stack %d distance from trackref to "
                                  "tracklet;radial distance;counts",
                                  hd),
                             2000, -0.5, 19.5);
  }

  LOG(info) << nev << " hit entries found and " << nevtrackref
            << " mctree events found, digit events : " << ndigitev
            << " trackletevents : " << ntrackletev;
  // std::vector<float> xenter, xexit, yenter, yexit, xUJ, yUJ,
  // xUK, yUK;
  std::vector<float> trackletx, tracklety;
  std::vector<o2::math_utils::Point3D<float>> trackletf, enterunmatched, exitunmatched;
  std::vector<MCTracklet> mctracklets;
  std::vector<RefTracklet> simulationtracklets; 

  std::vector<o2::TrackReference> trackrefsvectorwithextra;
  std::vector<o2::TrackReference> trackrefsvector; 
  // track references indexed by trackid and detector.

  // we can define the above here as the hits tree has only 1 event, all of
  // them

  for (int iev = 0; iev < nev; ++iev) {
    hitTree->GetEvent(iev);
    trackrefTree->GetEvent(iev);
    //  sort track ref by time and enter before exit.
    for (const auto &trackref : *trackrefs) {
 //     LOGP(info, "trackrefs  det id : {} ", trackref.getDetectorId());
      if (trackref.getDetectorId() == 2) {
        // only copy the trd trackrefs.
        // std::map<int, o2::TrackReference> inner;
        // inner.insert(std::make_pair(det, trackref));
        // trackreferences.insert(trackref.getTrackID(), inner);
        trackrefsvectorwithextra.push_back(trackref);
        switch (trackref.getUserId() & 0x3) {
        case 0x1:
          // direction entering
          enterunmatched.push_back(o2::math_utils::Point3D<float>(
              trackref.X(), trackref.Y(), trackref.Z()));
          break;
        case 0x2:
          exitunmatched.push_back(o2::math_utils::Point3D<float>(
              trackref.X(), trackref.Y(), trackref.Z()));
          // direction exiting
          break;
        }
      }
    }

    for (auto &trf : trackrefsvectorwithextra) {
//      LOGP(info, "TRFa : {} s:s:l {}:{}:{} uid {:08x} det:{}",
           // LOGP(info, "TRF : {} t:{} s:s:l {}:{}:{} uid {:08x} det:{}",
//           trf.getTrackID(),
 //          o2::trd::HelperMethods::getSector(trf.getUserId() >> 2),
  //         o2::trd::HelperMethods::getStack(trf.getUserId() >> 2),
   //        o2::trd::HelperMethods::getLayer(trf.getUserId() >> 2),
    //       trf.getUserId(), trf.getDetectorId());
      // LOG(info) << "trackref : " << trf.getTrackID();
    }
   /* std::sort(std::begin(trackrefsvectorwithextra),
              std::end(trackrefsvectorwithextra),
              [&trackrefsvectorwithextra](o2::TrackReference &i,
                                          o2::TrackReference &j) {
                return comparetrackrefs(i, j);
              });*/
    // now go through the sorted list and remove those that are the trackid and
    // same detector but extra enter and exit points. Mostly extra entry points.
    std::cout << " sorting done ... " << std::endl;
    // check thesorting
    int countrfb=0;
    for (auto &trf : trackrefsvectorwithextra) {
 /*     if(countrfb!=0)LOGP(info, "TRFb : {} s:s:l {}:{}:{} uid {:08x} det:{}",
           // LOGP(info, "TRF : {} t:{} s:s:l {}:{}:{} uid {:08x} det:{}",
           trf.getTrackID(),
           o2::trd::HelperMethods::getSector(trf.getUserId() >> 2),
           o2::trd::HelperMethods::getStack(trf.getUserId() >> 2),
           o2::trd::HelperMethods::getLayer(trf.getUserId() >> 2),
           trf.getUserId(), trf.getDetectorId());
      LOG(info) << "trackref : " << trf.getTrackID() << trf;*/
      countrfb++;
    }
    // LOGP(info, "trackrefswithextra size : : {}  ",
    //     trackrefsvectorwithextra.size());
    for (int ref = 0; ref < trackrefsvectorwithextra.size(); ++ref) {
      int trackid = trackrefsvectorwithextra[ref].getTrackID();
      //LOGP(info, "Trackid : {} trackid {} ",
//           trackrefsvectorwithextra[ref].getTrackID(), trackid);
      while (trackrefsvectorwithextra[ref].getTrackID() == trackid) {
        // while we are on the same trackid.
        int det = trackrefsvectorwithextra[ref].getUserId() >> 2;
        std::vector<o2::TrackReference> references;
        //LOGP(info, "trackid : {} det {} ", trackid, det);
        while (((trackrefsvectorwithextra[ref].getUserId() >> 2) == det) &&
               (trackrefsvectorwithextra[ref].getTrackID() == trackid)) {
          // we are on the same track and detector.
          references.push_back(trackrefsvectorwithextra[ref]);
       //   LOG(info) << "adding to trackrefsvectorwithextra "
        //            << trackrefsvectorwithextra[ref];
          ++ref;
        }
        // now strip out that which is not the first entry and last exit.
        // base it on the distance from the interaction point
        std::pair<double, int> min = std::make_pair(999999, 0);
        std::pair<double, int> max = std::make_pair(0, 0);
        int idx = 0;
        //LOGP(info,
        //     "------ start references distances idx:{} size of references:{}",
         //    idx, references.size());
        if (references.size() > 1) {
          for (auto &reftrk : references) {
            double distance =
                sqrt(reftrk.X() * reftrk.X() + reftrk.Y() * reftrk.Y());
         //   LOGP(info, "distance : {} idx:{} x:y {}:{}", distance, idx,
//                 reftrk.X(), reftrk.Y());
            if (distance < std::abs(min.first)) {
              min = std::make_pair(distance, idx);
          //    LOGP(info, "try get depth of  references");
              int depth1 = std::distance(references.end(), references.begin());
           //   LOGP(info, "depth of references {}", depth1);
           //   LOGP(info, "try get distance into references");
              int depth = 0; // std::distance(reftrk, references.begin());
            //  LOGP(info, "depth {}", depth);
              // LOGP(info, "MIN: distance : {} mindist:{} idx:{} x:y {}:{}",
              //    distance, min.first, idx, reftrk.X(), reftrk.Y());
            } else {
              if (distance > std::abs(max.first)) {
                max = std::make_pair(distance, idx);
                // LOGP(info, "MAX: distance : {} maxdist:{} idx:{} x:y {}:{}",
                //    distance, max.first, idx, reftrk.X(), reftrk.Y());
              }
            }
            ++idx;
          }
          //LOGP(info, "mctracklet from [{}] {}:{} to [{}] {}:{}", min.second,
 //j              references[min.second].X(), references[min.second].Y(),
   //            max.second, references[max.second].X(),
    //           references[max.second].Y());
          mctracklets.push_back(MCTracklet(
              references[min.second], references[max.second], trackid, det));
        }
        //LOGP(info, "--------------- end references distances");
      }
    }
  }
  // now we have a vector of mctracklets defining the enter and exit points
  // together with the detector and trackid.
  //

  for (int iev = 0; iev < ntrackletev; ++iev) {
    // if (iev > 1)
    //  continue;
    digitTree->GetEvent(iev);
    trackletTree->GetEvent(iev);
    //   std::cout << " # Digits : " << digits->size() << " # tracklets " <<
    //   tracklets->size() << std::endl;
    int counttriggers = 0;
    for (auto &trigrec : *trigRecs) {
      counttriggers++;
      // if (counttriggers != 2)
      //  continue;
      //	    cout << " time : " << trigrec.getBCData().bc2ns() << " ns "
      //<< std::endl;
      int digitcounter = 0;
      for (int digit = trigrec.getFirstDigit();
           digit < trigrec.getFirstDigit() + trigrec.getNumberOfDigits();
           ++digit) {
        if ((*digits)[digit].getDetector() == detector) {
          digitcounter++;
        }
      }
      int trkltcounter = 0;
      int trackrefscounter = 0;
      for (int trklt = trigrec.getFirstTracklet();
           trklt < trigrec.getFirstTracklet() + trigrec.getNumberOfTracklets();
           ++trklt) {
        trkltcounter++;
        auto hcid = (*tracklets)[trklt].getHCID();
        auto trackletdet = hcid / 2;
        int trackletsector = o2::trd::HelperMethods::getSector(hcid / 2);
        int trackletstack = o2::trd::HelperMethods::getStack(hcid / 2);
        int trackletlayer = o2::trd::HelperMethods::getLayer(hcid / 2);
        o2::trd::Tracklet64 tracklet;
        RefTracklet reftracklet;
        reftracklet.mTracklet.setTrackletWord((*tracklets)[trklt].getTrackletWord());

        tracklet.setTrackletWord((*tracklets)[trklt].getTrackletWord());
        trackletangles->Fill(tracklet.getSlopeFloat());
        auto caltracklet = transformer.transformTracklet(tracklet, false);
        reftracklet.mCalibratedTracklet = transformer.transformTracklet(tracklet, false);
         
        std::array<double, 3> xyz;
        std::array<double, 3> txyzfd;
        std::array<float, 3> txyzf;
        xyz[0] = caltracklet.getX();
        xyz[1] = caltracklet.getY();
        xyz[2] = caltracklet.getZ();
        std::array<float, 3> txyz =
            transformer.transformL2T(hcid / 2, xyz);
        reftracklet.transformL2T(txyz);
        txyzfd = rotateL2G(txyz, trackletsector);
        reftracklet.rotateL2G(txyzfd);
         
        txyzf[0] = txyzfd[0];
        txyzf[1] = txyzfd[1];
        txyzf[2] = txyzfd[2];
        trackletf.push_back(
            o2::math_utils::Point3D<float>(txyzf[0], txyzf[1], txyzf[2]));
        simulationtracklets.push_back(reftracklet); 
        // find closest trackref :
        double mindistance = 9999999; // something way outside the trd
        std::array<double, 3> mindist = {0, 0, 0};
        for (auto &trf : mctracklets) {
          int det = trf.det;
          int refstack = o2::trd::HelperMethods::getStack(det);
          if (refstack != trackletstack)
            continue;
          if (det != trackletdet)
            continue;

          double x = trf.mExit.X() - txyzfd[0];
          double y = trf.mExit.Y() - txyzfd[1];
          double z = trf.mExit.Z() - txyzfd[2];
          double distance = std::sqrt(x * x + y * y + z * z);
          if (distance < mindistance) {
            mindistance = distance;
            mindist[0] = trf.mExit.X() - txyzfd[0];
            mindist[1] = trf.mExit.Y() - txyzfd[1];
            mindist[2] = trf.mExit.Z() - txyzfd[2];
          }
          trackrefscounter++;
        }
        hdistance[trackletstack]->Fill(mindistance);
        hdistances[0]->Fill(std::abs(mindistance));
        hdistances[1]->Fill(std::abs(mindist[0]));
        hdistances[2]->Fill(std::abs(mindist[1]));
        hdistances[3]->Fill(std::abs(mindist[2]));
      }
    } // end of triggerrecord loop
  }
  auto numtrackletf = trackletf.size();
  auto numsimulationtracklets = simulationtracklets.size();
  LOGP(info, "sizes :  trackletf{}", numtrackletf);
  LOGP(info, "sizes :  simulationtracklets{}", numsimulationtracklets);
  std::vector<float> x, y;

  x.clear();
  y.clear();
  for (auto &xy : enterunmatched) {
    x.push_back(xy.X());
    y.push_back(xy.Y());
  }
  TGraph *htrackrefxyenter =
      new TGraph(enterunmatched.size(), x.data(), y.data());
  htrackrefxyenter->SetTitle("TrackReference x vs y entering;x(m);y(m)");
  htrackrefxyenter->SetMarkerStyle(20);
  x.clear();
  y.clear();
  for (auto &xy : exitunmatched) {
    x.push_back(xy.X());
    y.push_back(xy.Y());
  }
  TGraph *htrackrefxyexit =
      new TGraph(exitunmatched.size(), x.data(), y.data());
  htrackrefxyexit->SetTitle("TrackReference x vs y exiting;x(m);y(m)");
  htrackrefxyexit->SetMarkerStyle(24);

  x.clear();
  y.clear();
  for (auto &xy : trackletf) {
    x.push_back(xy.X());
    y.push_back(xy.Y());
  }
  TGraph *htrackletsf = new TGraph(numtrackletf, x.data(), y.data());
  htrackletsf->SetTitle("Tracklets x vs y rotated correctly;x(m);y(m)");

  TCanvas *c3 = new TCanvas("c3", "trd x vs y for trackref entering", 800, 600);
  c3->cd();
  htrackrefxyenter->Draw("PA");

  TCanvas *c4 = new TCanvas("c4", "trd x vs y for trackref exiting", 800, 600);
  c4->cd();
  htrackrefxyexit->Draw("PA");

  TCanvas *c7 =
      new TCanvas("c7", "trd x vs y tracklets rotated correctly", 800, 600);
  c7->cd();
  htrackletsf->Draw("A*");

  TCanvas *c6 =
      new TCanvas("c6", "trd x vs y for tracklets and trackref", 800, 600);
  c6->cd();
  TMultiGraph *mg = new TMultiGraph();
  htrackrefxyexit->SetMarkerColor(kRed);
  htrackrefxyenter->SetMarkerColor(kBlue);
  htrackletsf->SetMarkerColor(kMagenta);
  int index = 0;
  //  for (auto &label : identer) {
  //    TText *a;
  //    a = new TText(xenter[index], yenter[index], Format("%d", label));
  //    index++;
  //  }
  //  mg->Add(htrackrefxyUJ, "*");
  htrackrefxyexit->SetMarkerStyle(24);
  mg->Add(htrackrefxyexit, "P");
  htrackrefxyenter->SetMarkerStyle(20);
  mg->Add(htrackrefxyenter, "P");
  mg->Add(htrackletsf, "*");
  mg->Draw("a");
  LOGP(info, "we have {} mctracklet", mctracklets.size());
  for (auto &mctracklet : mctracklets) {
    // for each mctracklet draw a line
    TLine *a = new TLine(mctracklet.mEnter.X(), mctracklet.mEnter.Y(),
                         mctracklet.mExit.X(), mctracklet.mExit.Y());
    a->Draw();
    TText *t =
        new TText(mctracklet.mExit.X(), mctracklet.mExit.Y(),
                  Form("%d z-%f", mctracklet.mMCTrackId, mctracklet.mExit.Z()));
    t->SetTextSize(0.010);
    t->Draw();
    // LOGP(info, "line from {}:{} to {}:{}", mctracklet.mEnter.X(),
    //    mctracklet.mEnter.Y(), mctracklet.mExit.X(), mctracklet.mExit.Y());
  }
/* for (auto &mctracklet : simulationtracklets) {
    // for each mctracklet draw a line
    TLine *a = new TLine(mctracklet.mEnter.X(), mctracklet.mEnter.Y(),
                         mctracklet.mExit.X(),  mctracklet.mExit.Y());
    a->Draw();
    TText *t =
        new TText(mctracklet.mExit.X(), mctracklet.mExit.Y(),
                  Form("%d z-%f", mctracklet.mMCTrackId, mctracklet.mExit.Z()));
    t->SetTextSize(0.010);
    t->Draw();
    // LOGP(info, "line from {}:{} to {}:{}", mctracklet.mEnter.X(),
    //    mctracklet.mEnter.Y(), mctracklet.mExit.X(), mctracklet.mExit.Y());
  }*/
  // ki  mgeo->Draw();
  for (int sector = 0; sector < 18; ++sector) {
    for (int stack = 0; stack < 1 ; ++stack) {
    for (int layer = 0; layer < 6; ++layer) {
      int idx = stack + layer * 5;
      double rmin = mgeo->getTime0(layer);
      double rmax = mgeo->getTime0(layer) - mgeo->drThick() - mgeo->amThick();
      double ymin = -mgeo->getChamberWidth(layer) / 2;
      double ymax = mgeo->getChamberWidth(layer) / 2;
      double zmin = -mgeo->getChamberLength(layer, stack) / 2;
      double zmax = mgeo->getChamberLength(layer, stack) / 2;
      std::array<double, 3> tl, br, tr, bl;
      std::array<double, 3> tlx, brx, trx, blx;
      //  ymin,ymax   ... zmin zmax
      //  r is along the x-axis in the TGraph
      tl[0] = rmin; // std::sin(ymax / rmax);
      tl[1] = ymax;
      tl[2] = 0;
      tr[0] = rmax;
      tr[1] = ymax;
      tr[2] = 0;
      bl[0] = rmin;
      bl[1] = ymin;
      bl[2] = 0;
      br[0] = rmax;
      br[1] = ymin;
      br[2] = 0;
      tlx = rotateL2G(tl, sector);
      trx = rotateL2G(tr, sector);
      blx = rotateL2G(bl, sector);
      brx = rotateL2G(br, sector);
      double x[] = {tlx[0], trx[0], brx[0], blx[0], tlx[0]};
      double y[] = {tlx[1], trx[1], brx[1], blx[1], tlx[1]};
      TPolyLine *pl = new TPolyLine(5, x, y);
      pl->SetLineColor(kBlue);
      pl->SetFillColorAlpha(0, 0);
      // pl->Draw("f");
      pl->Draw();
      rmin = rmax;
      rmax = rmax - mgeo->craHght();
      tl[0] = rmin; // std::sin(ymax / rmax);
      tl[1] = ymax;
      tl[2] = 0;
      tr[0] = rmax;
      tr[1] = ymax;
      tr[2] = 0;
      bl[0] = rmin;
      bl[1] = ymin;
      bl[2] = 0;
      br[0] = rmax;
      br[1] = ymin;
      br[2] = 0;
      tlx = rotateL2G(tl, sector);
      trx = rotateL2G(tr, sector);
      blx = rotateL2G(bl, sector);
      brx = rotateL2G(br, sector);
      double x1[] = {tlx[0], trx[0], brx[0], blx[0], tlx[0]};
      double y1[] = {tlx[1], trx[1], brx[1], blx[1], tlx[1]};
      TPolyLine *pla = new TPolyLine(5, x1, y1);
      pla->SetLineColor(kMagenta);
      pla->SetFillColorAlpha(0, 0);
      // pl->Draw("f");
      pla->Draw();
        //TODO
        //take the angle calculate as deviation from a line from center to tracklet dot.
        //draw a line from dot at that angle intersecting start of drift region.
        //probably put this in a function returning the "start point" and draw.
    }
  }
}
/*TCanvas *c9 = new TCanvas("c9", "trd trackref distance ", 800, 600);
c9->Divide(3, 2);
c9->cd(1);
hdistance[0]->Draw();
c9->cd(2);
hdistance[1]->Draw();
c9->cd(3);
hdistance[2]->Draw();
c9->cd(4);
hdistance[3]->Draw();
c9->cd(5);
hdistance[4]->Draw();
TCanvas *ca = new TCanvas(
    "ca", "trd trackref distance to tracklet radial distance", 800, 600);
ca->cd();
hdistances[0]->Draw();
TCanvas *cb = new TCanvas(
    "cb", "trd trackref distance to tracklet x distance", 800, 600);
cb->cd();
hdistances[1]->Draw();
TCanvas *cc = new TCanvas(
    "cc", "trd trackref distance to tracklet y distance", 800, 600);
cc->cd();
hdistances[2]->Draw();
TCanvas *cd = new TCanvas(
    "cd", "trd trackref distance to tracklet z distance", 800, 600);
cd->cd();
hdistances[3]->Draw();
TCanvas *ce = new TCanvas("ce", "trd angle (deflection", 800, 600);
ce->cd();
trackletangles->Draw();*/
}
