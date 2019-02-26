#ifndef ALIMUONCLUSTERFINDERCOG_H
#define ALIMUONCLUSTERFINDERCOG_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class ClusterFinder_COG
/// \brief A very basic (and mostly useless, probably) cluster finder
/// 
// Author Laurent Aphecetche, Subatech

#ifndef AliMUONVCLUSTERFINDER_H
#  include "AliMUONVClusterFinder.h"
#endif

class ClusterFinder_COG : public AliMUONVClusterFinder
{
public:
  ClusterFinder_COG(AliMUONVClusterFinder* clusterFinder);
  virtual ~ClusterFinder_COG();

  using AliMUONVClusterFinder::Prepare;
  
  virtual Bool_t Prepare(Int_t detElemId,
                         TObjArray* pads[2],
                         const AliMpArea& area);
  
  virtual AliMUONCluster* NextCluster();
  
private:
  /// Not implemented
  ClusterFinder_COG(const ClusterFinder_COG& rhs);
  /// Not implemented
  ClusterFinder_COG& operator=(const ClusterFinder_COG& rhs);

  void ComputePosition(AliMUONCluster& cluster);

private:
    AliMUONVClusterFinder* fPreClusterFinder; ///< the preclustering we use

  ClassDef(ClusterFinder_COG,1) // A very basic (and mostly useless, probably) cluster finder
};

#endif
