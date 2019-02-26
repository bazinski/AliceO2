




#ifndef O2_MCH_VIRTUALCLUSTERFINDER_H
#define O2_MCH_VIRTUALCLUSTERFINDER_H


/// \ingroup rec
/// \class VirtualClusterFinder
/// \brief Interface of a cluster finder.
/// 
//  Author Laurent Aphecetche

class Cluster;
class AliMUONRecoParam;
class AliMpVSegmentation;
class NPad;
class MpArea;

class VirtualClusterFinder 
{
public:
  VirtualClusterFinder();
  virtual ~VirtualClusterFinder();
  
  /// \todo add comment

  virtual bool NeedSegmentation() const { return kFALSE; }
  
  virtual bool Prepare(int detElemId,
                         TObjArray* pads[2],
                         const MpArea& area);

  virtual bool Prepare(int detElemId,
                         TObjArray* pads[2],
                         const MpArea& area,
                         const AliMpVSegmentation* segmentations[2]);  
  
  /// \todo add comment
  virtual Cluster* NextCluster() = 0;
  
  /** Add a pad to the list of pads to be considered for clustering.
    Typical usage is to "put-back-in-business" a pad that was part 
    of a previous cluster (returned by NextCluster) but was externally
    identified of not being part of that cluster, so it must be reuseable.
    Might not be implemented by all cluster finders...
    (in which case it must returns kFALSE)
    */
  virtual bool UsePad(const Pad& pad);
  
  /** Specify a couple of charge hints. We call them hints because some
   clustering need them and use them directly, other cook them before
   using them, and some others yet simply don't care about them.
   */
  virtual void SetChargeHints(double /*lowestPadCharge*/, double /*lowestClusterCharge*/) { }
  
  ClassDef(VirtualClusterFinder,0) // Interface of a MUON cluster finder.
};

#endif
