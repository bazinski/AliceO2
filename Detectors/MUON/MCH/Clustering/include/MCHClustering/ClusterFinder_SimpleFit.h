#ifndef O2_MCH_CLUSTERFINDER_SIMPLEFIT_H
#define O2_MCH_CLUSTERFINDER_SIMPLEFIT_H


/// \ingroup rec
/// \class ClusterFinder_SimpleFit
/// \brief Basic cluster finder 
/// 
// Author Laurent Aphecetche, Subatech

class VirtualClusterFinder;
class AliMUONMathieson;
class MpArea;

class ClusterFinder_SimpleFit : public VirtualClusterFinder
{
public:
  ClusterFinder_SimpleFit(VirtualClusterFinder* clusterFinder);
  virtual ~ClusterFinder_SimpleFit();
  
  using VirtualClusterFinder::prepare;

  virtual Bool_t Prepare(Int_t detElemId,
                         TObjArray* pads[2],
                         const MpArea& area);
  
  virtual AliMUONCluster* nextCluster();
  
  virtual void setChargeHints(Double_t /*lowestPadCharge*/, Double_t lowestClusterCharge) { 
    mLowestClusterCharge=lowestClusterCharge; 
  }
  
private:
  /// Not implemented
  ClusterFinder_SimpleFit(const ClusterFinder_SimpleFit& rhs);
  /// Not implemented
  ClusterFinder_SimpleFit& operator=(const ClusterFinder_SimpleFit& rhs);

  void computePosition(Cluster& cluster);

private:
  VirtualClusterFinder* mClusterFinder; //!<! the preclustering we use
  Mathieson* mMathieson; //!<! Mathieson to compute the charge repartition
  Double_t mLowestClusterCharge; //!<! minimum cluster charge we allow
  
  ClassDef(ClusterFinder_SimpleFit,0) // Basic cluster finder
};

#endif
