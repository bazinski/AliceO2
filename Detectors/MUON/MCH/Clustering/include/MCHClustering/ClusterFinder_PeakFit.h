#ifndef O2_MCH_CLUSTERFINDER_PEAKFIT_H
#define O2_MCH_CLUSTERFINDER_PEAKFIT_H

/// \ingroup rec
/// \class ClusterFinder_PeakFit
/// \brief Cluster finder in MUON arm of ALICE
///
//  Author Alexander Zinchenko, JINR Dubna; Laurent Aphecetche, SUBATECH
//

class AliMUONMathieson;
class VirtualClusterFinder;
class TH2D;
class Pad;
class Mathieson;

class ClusterFinder_PeakFit : public VirtualClusterFinder
{
public:
  ClusterFinder_PeakFit(bool plot, VirtualClusterFinder* clusterFinder); // Constructor
  virtual ~ClusterFinder_PeakFit(); // Destructor

  /// It needs segmentation
  virtual bool needSegmentation() const { return kTRUE; }

  using VirtualClusterFinder::prepare;

  virtual bool prepare(int detElemId, TObjArray* pads[2],
                         const AliMpArea& area, const AliMpVSegmentation* seg[2]);

  virtual Cluster* nextCluster();
  
  /// Return the number of local maxima
  int getNMax() const { return fNMax; }
  
  virtual void print(Option_t* opt="") const;

private:
  /// Not implemented
  ClusterFinder_PeakFit(const ClusterFinder_PeakFit& rhs);
  /// Not implemented
  ClusterFinder_PeakFit& operator=(const ClusterFinder_PeakFit& rhs);

  bool workOnPreCluster();

  /// Check precluster to simplify it (if possible), and return the simplified cluster
  Cluster* checkPrecluster(const Cluster& cluster); 
  Cluster* checkPreclusterTwoCathodes(Cluster* cluster); 
  
  /// Checks whether a pad and a pixel have an overlapping area.
  bool overlap(const AliMUONPad& pad, const AliMUONPad& pixel); 
  
  /// build array of pixels
  void buildPixArray(Cluster& cluster); 
  void buildPixArrayOneCathode(Cluster& cluster); 
  void padOverHist(int idir, int ix0, int iy0, AliMUONPad *pad, TH2D *h1, TH2D *h2);

  void removePixel(int i);
  
  Pad* pixel(int i) const;
  
  int findLocalMaxima(TObjArray *pixArray, int *localMax, double *maxVal); // find local maxima 
  void  flagLocalMax(TH2D *hist, int i, int j, int *isLocalMax); // flag local max
  void  findClusterCOG(Cluster& cluster, const int *localMax, int iMax); // find cluster around local max with COG
  void  findClusterFit(Cluster& cluster, const int *localMax, const int *maxPos, int nMax); // find cluster around local max with fitting
  void  padsInXandY(Cluster& cluster, int &nInX, int &nInY) const; // get number of pads in X and Y

  void checkOverlaps();

private:
  // Status flags for pads
  static const int mgkZero; ///< pad "basic" state
  static const int mgkMustKeep; ///< do not kill (for pixels)
  static const int mgkUseForFit; ///< should be used for fit
  static const int mgkOver; ///< processing is over
  static const int mgkModified; ///< modified pad charge 
  static const int mgkCoupled; ///< coupled pad    
    
  // Some constants
  static const double mgkZeroSuppression; ///< average zero suppression value
  static const double mgkDistancePrecision; ///< used to check overlaps and so on
  static const TVector2 mgkIncreaseSize; ///< idem
  static const TVector2 mgkDecreaseSize; ///< idem
  
  VirtualClusterFinder* mPreClusterFinder; //!<! the pre-clustering worker
  Cluster* mPreCluster; //!<! current pre-cluster
  TObjArray mClusterList; //!<! clusters corresponding to the current pre-cluster
  Mathieson* mMathieson; //!<! Mathieson to compute the charge repartition
  
  int mEventNumber; //!<! current event being processed
  int mDetElemId; //!<! current DE being processed
  int mClusterNumber; //!<! current cluster number
  int mNMax; //!<! number of local maxima
  TH2D *mHistAnode; //!<! histo for local maxima search
  
  const AliMpVSegmentation *mkSegmentation[2]; //!<! new segmentation
  
  TObjArray* mPixArray; //!<! collection of pixels
  int mDebug; //!<! debug level
  bool mPlot; //!<! whether we should plot thing (for debug only, quite slow!)
  
  int mNClusters; //!<! total number of clusters
  int mNAddVirtualPads; //!<! number of clusters for which we added virtual pads
  
  ClassDef(ClusterFinder_PeakFit,0) // cluster finder in MUON arm of ALICE
};

#endif
