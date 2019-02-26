#ifndef ALIMUONVCLUSTERSERVER_H
#define ALIMUONVCLUSTERSERVER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONVClusterServer
/// \brief Interface of a cluster finder for combined tracking
/// 
// Author Laurent Aphecetche, Subatech


class AliMUONVDigitStore;
class VirtualClusterStore;
class AliMUONVTriggerTrackStore;
class AliMUONRecoParam;
class MpArea;
class TIter;

class VirtualClusterServer 
{
public:
  AliMUONVClusterServer();
  virtual ~AliMUONVClusterServer();
  
  /// Find and add clusters from a given region of a given chamber to the store.
  virtual int clusterize(int chamberId, 
                           VirtualClusterStore& clusterStore,
                           const AliMpArea& area,
			   const AliMUONRecoParam* recoParam = 0x0) = 0;
  
  /// Specify an iterator to loop over the digits needed to perform our job.
  virtual void useDigits(TIter& next, AliMUONVDigitStore* digitStore = 0x0) = 0;
  
  /// Use trigger tracks. Return kFALSE if not used.
  virtual bool useTriggerTrackStore(AliMUONVTriggerTrackStore* /*trackStore*/) { return kFALSE; }
  
  ClassDef(VirtualClusterServer,1) // Cluster server interface
};

#endif
