#ifndef O2_MUON_CLUSTERSTORE_ITERATOR_H
#define O2_MUON_CLUSTERSTORE_ITERATOR_H

/// \ingroup base
/// \class AliMUONClusterStoreV2Iterator
/// \brief Base implementation of TIterator for AliMUONClusterStoreV2
/// 
//  Author Philippe Pillot, Subatech

class Iterator;
class ClusterStore;

class ClusterStoreIterator : public TIterator
{
public:
  AliMUONClusterStoreV2Iterator(const ClusterStore* store,
                                int firstChamberId, int lastChamberId);
  
  virtual ~AliMUONClusterStoreV2Iterator();
  
  TObject* Next();
  
  void Reset();
  
  /// Return 0 as we're not dealing with TCollection objects really
  virtual const TCollection* GetCollection() const { return 0x0; }
  
private:
  TObject* NextInCurrentChamber() const;
  
private:
  /// Not implemented
  ClusterStoreIterator(const ClusterStoreIterator& rhs);
  /// Not implemented
  ClusterStoreIterator& operator=(const ClusterStoreIterator& rhs);
  /// Overriden TIterator virtual operator=
  ClusterStoreIterator& operator=(const TIterator& rhs);

  const ClusterStore* fkStore; ///< store to iterate upon
  int fFirstChamberId; ///< first chamber
  int fLastChamberId; ///< last chamber
  int fCurrentChamberId; ///< current chamber
  TIterator* fChamberIterator; ///< helper iterator
  
  ClassDef(ClusterStoreIterator,0) // Implementation of TIterator
};

#endif
