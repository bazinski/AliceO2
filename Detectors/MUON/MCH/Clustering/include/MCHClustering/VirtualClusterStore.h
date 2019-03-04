#ifndef O2_MUON_VIRTUALCLUSTERSTORE_H
#define O2_MUON_VIRTUALCLUSTERSTORE_H

/// \ingroup rec
/// \class VirtualClusterStore
/// \brief Interface of a cluster container
/// 
// Author Laurent Aphecetche, Subatech

class VirtualCluster;

class VirtualClusterStore : public AliMUONVStore
{
public:
  VirtualClusterStore();
  virtual ~VirtualClusterStore();
  
  virtual bool Add(TObject* object);

  /// Add a cluster object to the store
  virtual VirtualCluster* Add(const VirtualCluster& Cluster) = 0;
  /// Create a new cluster with an unique ID and add it to the store
  virtual VirtualCluster* Add(int chamberId, int detElemId, int clusterIndex) = 0;

  using AliMUONVStore::Create;
  
  static VirtualClusterStore* Create(TTree& tree);
  
  /// Create a cluster
  virtual VirtualCluster* CreateCluster(int chamberId, int detElemId, int clusterIndex) const = 0;
  
  /// Return an iterator to loop over the whole store
  virtual TIterator* CreateIterator() const = 0;

  /// Return an iterator to loop over the store in the given chamber range
  virtual TIterator* CreateChamberIterator(int firstChamberId, int lastChamberId) const = 0;
  
  /// Clear container
  virtual void Clear(Option_t* opt="") = 0;
  
  /// Remove a cluster object to the store
  virtual VirtualCluster* Remove(VirtualCluster& cluster) = 0;
    
  using AliMUONVStore::FindObject;

  // Find an object (default is the same as in AliMUONVStore)
  virtual VirtualCluster* FindObject(const TObject* object) const;
  
  // Find an object by its uniqueID (default is the same as in AliMUONVStore)
  virtual VirtualCluster* FindObject(unsigned int uniqueID) const;
  
  ClassDef(VirtualClusterStore,1) // Cluster container interface
};

#endif

