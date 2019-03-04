#ifndef O2_MUON_CLUSTERSTORE_H
#define O2_MUON_CLUSTERSTORE_H

/// \ingroup rec
/// \class ClusterStore
/// \brief Implementation of VClusterStore
/// 
// Author Philippe Pillot, Subatech

class VirtualClusterStore;
class TClonesArray;

class ClusterStore : public VirtualClusterStore
{
  friend class ClusterStoreIterator;
  
public:
  ClusterStore();
  ClusterStore(TRootIOCtor* dummy);
  ClusterStore(const ClusterStore& store);
  ClusterStore& operator=(const ClusterStore& store);  
  virtual ~ClusterStore();
  
  virtual void Clear(Option_t* opt="");
  
  /// Whether the Connect(TTree&) method is implemented
  virtual bool CanConnect() const { return kTRUE; }
  virtual bool Connect(TTree& tree, bool alone=kTRUE) const;
  
  /// Create an empty copy of this
  virtual ClusterStore* Create() const { return new ClusterStore; }
  
  virtual VirtualCluster* CreateCluster(int chamberId, int detElemId, int clusterIndex) const;
  
  using VirtualClusterStore::Add;
  
  virtual VirtualCluster* Add(const VirtualCluster& Cluster);
  virtual VirtualCluster* Add(int chamberId, int detElemId, int clusterIndex);

  virtual VirtualCluster* Remove(VirtualCluster& cluster);

  using VirtualClusterStore::GetSize;
  
  /// Return the number of clusters we hold
  virtual int GetSize() const {return fClusters->GetLast()+1;}
  
  using AliMUONVStore::FindObject;
  
  VirtualCluster* FindObject(const TObject* object) const;
  VirtualCluster* FindObject(unsigned int uniqueID) const;
  
  virtual TIterator* CreateIterator() const;
  virtual TIterator* CreateChamberIterator(int firstChamberId, int lastChamberId) const;
  
private:
  void ReMap();
  void UpdateMap(VirtualCluster& cluster);
  
private:
  TClonesArray* fClusters; ///< collection of clusters
  TClonesArray* fMap;      //!<! index map for fast cluster retrieval
  bool        fMapped;   //!<! whether our internal indices are uptodate
  
  ClassDef(ClusterStore,1) // Implementation of VClusterStore
};

#endif
