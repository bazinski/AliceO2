
//-----------------------------------------------------------------------------
/// \class ClusterStoreIterator
///
/// Implementation of TIterator for AliMUONClusterStoreV2
///
/// \author Philippe Pillot, Subatech
///
//-----------------------------------------------------------------------------

#include "ClusterStoreIterator.h"
#include "AliMUONClusterStoreV2.h"

#include "AliMpExMapIterator.h"
#include "AliMpExMap.h"

#include "fairlogger/Logger.h"

//_____________________________________________________________________________
ClusterStoreIterator::ClusterStoreIterator(const AliMUONClusterStoreV2* store,
                                                             int firstChamberId, int lastChamberId)
: TIterator(),
  fkStore(store),
  fFirstChamberId(firstChamberId),
  fLastChamberId(lastChamberId),
  fCurrentChamberId(-1),
  fChamberIterator(0x0)
{
  /// Constructor for partial iteration
  if (fFirstChamberId > fLastChamberId) {
    fLastChamberId = fFirstChamberId;
    fFirstChamberId = lastChamberId;
  }
  Reset();
}

//_____________________________________________________________________________
ClusterStoreIterator& 
ClusterStoreIterator::operator=(const TIterator& /*iter*/)
{
  // Overriden operator= (imposed by Root's definition of TIterator::operator= ?)

  LOG (fatal) << "ClusterStoreIterator::operator= reimplement me";
  return *this;
}

//_____________________________________________________________________________
ClusterStoreIterator::~ClusterStoreIterator()
{
  /// Destructor
  delete fChamberIterator;
}

//_____________________________________________________________________________
TObject* ClusterStoreIterator::NextInCurrentChamber() const
{
  /// Return the value corresponding to theKey in iterator iter
  
  return fChamberIterator->Next();
}

//_____________________________________________________________________________
TObject* ClusterStoreIterator::Next()
{
  /// Return next cluster in store
  TObject* o = NextInCurrentChamber();
  
  while (!o) {
    // fChamberIterator exhausted, try to get the next ones
    if (fCurrentChamberId == fLastChamberId) return 0x0; // we reached the end
    
    fCurrentChamberId++;
    delete fChamberIterator;
    fChamberIterator = static_cast<AliMpExMap*>(fkStore->fMap->UncheckedAt(fCurrentChamberId))->CreateIterator();
    
    o = NextInCurrentChamber();
  }
  
  return o;
}

//_____________________________________________________________________________
void ClusterStoreIterator::Reset()
{
  /// Reset the iterator
  fCurrentChamberId = fFirstChamberId;
  delete fChamberIterator;
  fChamberIterator = static_cast<AliMpExMap*>(fkStore->fMap->UncheckedAt(fCurrentChamberId))->CreateIterator();
}
