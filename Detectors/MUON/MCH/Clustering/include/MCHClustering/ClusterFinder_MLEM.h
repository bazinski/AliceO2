


#ifndef O2_MCH_CLUSTERFINDERMLEM_H
#define O2_MCH_CLUSTERFINDERMLEM_H

/// \ingroup rec
/// \class ClusterFinder_MLEM
/// \brief Cluster finder in MUON arm of ALICE
///
//  Author Alexander Zinchenko, JINR Dubna; Laurent Aphecetche, SUBATECH
//

class TH2D;
class TMinuit;


class Pad;

#include "MUONVClusterFinder.h"

class ClusterSplitter_MLEM;

class ClusterFinder_MLEM : public AliMUONVClusterFinder
{
public:
  ClusterFinder_MLEM(bool plot, AliMUONVClusterFinder* clusterFinder); // Constructor
  virtual ~ClusterFinder_MLEM(); // Destructor

  virtual bool needSegmentation() const { return kTRUE; }
  
  using AliMUONVClusterFinder::Prepare;

  virtual bool prepare(int detElemId,
                         TObjArray* pads[2],
                         const AliMpArea& area,
                         const AliMpVSegmentation* segmentations[2]);
  
  virtual Cluster* nextCluster();
  
  virtual void setChargeHints(double lowestPadCharge, double lowestClusterCharge);
  
  // Status flags for pads

               /// Return pad "basic" state flag
  static int getZeroFlag()       { return fgkZero; }
               /// Return do not kill flag
  static int getMustKeepFlag()   { return fgkMustKeep; }
               /// Return should be used for fit flag
  static int getUseForFitFlag()  { return fgkUseForFit; }
               /// Return processing is over flag
  static int getOverFlag()       { return fgkOver; }
               /// Return modified pad charge flag
  static int getModifiedFlag()   { return fgkModified; }
               /// Return coupled pad flag
  static int getCoupledFlag()    { return fgkCoupled; }
  
private:
  /// Not implemented
  ClusterFinder_MLEM(const ClusterFinder_MLEM& rhs);
  /// Not implemented
  ClusterFinder_MLEM& operator=(const ClusterFinder_MLEM& rhs);

  bool workOnPreCluster();

  /// Check precluster to simplify it (if possible), and return the simplified cluster
  Cluster* checkPrecluster(const Cluster& cluster); 
  Cluster* checkPreclusterTwoCathodes(Cluster* cluster); 
  
  /// Checks whether a pad and a pixel have an overlapping area.
  bool overlap(const Pad& pad, const Pad& pixel); 
  
  /// build array of pixels
  void buildPixArray(Cluster& cluster); 
  void buildPixArrayOneCathode(Cluster& cluster); 
  void padOverHist(int idir, int ix0, int iy0, Pad *pad,
		   TH2D *hist1, TH2D *hist2);

  void removePixel(int i);
  
  Pad* pixel(int i) const;
  
  bool mainLoop(Cluster& cluster, int iSimple); // repeat MLEM algorithm until pixels become sufficiently small
  
  void   mlem(Cluster& cluster, const double *coef, double *probi, int nIter); // use MLEM for cluster finding
  
  void   findCOG(double *xyc); // find COG position around maximum bin
  int  findNearest(const Pad *pixPtr0); // find nearest neighbouring pixel to the given one

  int findLocalMaxima(TObjArray *pixArray, int *localMax, double *maxVal); // find local maxima 
  void  flagLocalMax(TH2D *hist, int i, int j, int *isLocalMax); // flag local max
  void  findCluster(Cluster& cluster, const int *localMax, int iMax); // find cluster around local max
  void  addVirtualPad(Cluster& cluster); // add virtual pads for some clusters (if necessary)
  
  void  padsInXandY(Cluster& cluster, int &nInX, int &nInY) const; // get number of pads in X and Y

  /// Process simple cluster
  void simple(Cluster& cluster); 
  
  void plot(const char* outputfile);
    
  void computeCoefficients(Cluster& cluster, 
                           double* coef, double* probi);
  
  void checkOverlaps();
  void addBinSimple(TH2D *mlem, int ic, int jc);
  void maskPeaks(int mask);

private:
  // Status flags for pads
  static const int mgkZero; ///< pad "basic" state
  static const int mgkMustKeep; ///< do not kill (for pixels)
  static const int mgkUseForFit; ///< should be used for fit
  static const int mgkOver; ///< processing is over
  static const int mgkModified; ///< modified pad charge 
  static const int mgkCoupled; ///< coupled pad  
      
  // Some constants
  static const double mgkDistancePrecision; ///< used to check overlaps and so on
  static const TVector2 mgkIncreaseSize; ///< idem
  static const TVector2 mgkDecreaseSize; ///< idem
  
  AliMUONVClusterFinder* mPreClusterFinder; //!<! the pre-clustering worker
  Cluster* mPreCluster; //!<! current pre-cluster
  TObjArray mClusterList; //!<! clusters corresponding to the current pre-cluster
  
  int mEventNumber; //!<! current event being processed
  int mDetElemId; //!<! current DE being processed
  int mClusterNumber; //!<! current cluster number
  
  const AliMpVSegmentation *mkSegmentation[2]; //!<! new segmentation
  
  //int fCathBeg;               //!<! starting cathode (for combined cluster / track reco)
  //int fPadBeg[2];             //!<! starting pads (for combined cluster / track reco)
  
  //static     TMinuit* fgMinuit; //!<! Fitter
  TH2D *mHistMlem; //!<! histogram for MLEM procedure
  TH2D *mHistAnode; //!<! histogram for local maxima search
  
  TObjArray* fPixArray; //!<! collection of pixels
  int mDebug; //!<! debug level
  bool mPlot; //!<! whether we should plot thing (for debug only, quite slow!)
  
  ClusterSplitter_MLEM* mSplitter; //!<! helper class to go from pixel arrays to clusters
  int mNClusters; //!<! total number of clusters
  int mNAddVirtualPads; //!<! number of clusters for which we added virtual pads
  
  double mLowestPixelCharge; //!<! see AliMUONRecoParam
  double mLowestPadCharge; //!<! see AliMUONRecoParam
  double mLowestClusterCharge; //!<! see AliMUONRecoParam
  
  ClassDef(ClusterFinder_MLEM,0) // cluster finder in MUON arm of ALICE
};

#endif
