



#ifndef O2_MCH_CLUSTER_H
#define O2_MCH_CLUSTER_H


/// \ingroup rec
/// \class Cluster
/// \brief A group of adjacent pads
/// 
// Author Laurent Aphecetche


class AliMUONPad;
class MpArea;

class Cluster 
{
public:
  Cluster();
  Cluster(const Cluster& rhs);
  Cluster& operator=(const Cluster& rhs);
  
  ~Cluster();
  
  bool Contains(const AliMUONPad& pad) const;
  
  std::string asString() const;
  
  static bool areOverlapping(const Cluster& c1, const Cluster& c2);
  
  AliMUONPad* addPad(const AliMUONPad& pad);

  /// Area that contains all the pads (whatever the cathode)
  AliMpArea area() const;

  /// Area that contains all the pads of a given cathode
  AliMpArea area(int cathode) const;

  float charge() const;
  float charge(int cathode) const;

  /// Return the cathode's charges asymmetry
  float chargeAsymmetry() const;

  /// Return chi2 of the RawCharge fit (if any)
  float chi2() const { return fChi2; }

  /// Return false for pre-cluster
  bool hasPosition() const { return fHasPosition; }

  /// Whether we have at least one saturated pad in a given cathode 
  bool isSaturated(int cathode) const { return fIsSaturated[cathode]; }
  
  /// Whether we have one saturated pad on *each* cathode
  bool isSaturated() const { return IsSaturated(0) && IsSaturated(1); }
  
  /// Return the max charge on the chathod
  int maxChargeCathode() const { return Charge(0) > Charge(1) ? 0:1; }

  /// Return the max raw charge on the chathod
  int maxRawChargeCathode() const { return RawCharge(0) > RawCharge(1) ? 0:1; }

  /// Return the biggest pad dimensions for a given cathode
  TVector2 maxPadDimensions(int cathode, int statusMask, bool matchMask) const;
  
  /// Return the biggest pad dimensions
  TVector2 maxPadDimensions(int statusMask, bool matchMask) const;
  
  /// Return the smallest pad dimensions for a given cathode
  TVector2 minPadDimensions(int cathode, int statusMask, bool matchMask) const;
  
  /// Return the smallest pad dimensions
  TVector2 minPadDimensions(int statusMask, bool matchMask) const;
  
  int multiplicity() const;
  int multiplicity(int cathode) const;

  /// Compute number of pads in X and Y direction for a given cathode.  
  Long_t nofPads(int cathode, int statusMask, bool matchMask) const;
  
  /// Number of pads in (X,Y) direction, whatever the cathode.
  Long_t nofPads(int statusMask, bool matchMask=kTRUE) const;
  
  /// Return true as the function Compare is implemented
  bool isSortable() const { return kTRUE; }
  
  AliMUONPad* pad(int index) const;
  
  /// Return (x,y) of that cluster
  TVector2 position() const { return fPosition; }
  /// Return errors on (x,y)
  TVector2 positionError() const { return fPositionError; }

  virtual void print(Option_t* opt="") const;
      
  /// By default, return the average of both cathode RawCharges.
  float rawCharge() const;
  
  /// Returns the RawCharge on the given cathode.
  float rawCharge(int cathode) const;
    
  /// Return the cathode's raw charges asymmetry
  float rawChargeAsymmetry() const;
  
  void removePad(AliMUONPad* pad);

  /// Set cathode (re)computed charges  
  void setCharge(float chargeCath0, float chargeCath1)
  { fHasCharge = kTRUE; fCharge[0]=chargeCath0; fCharge[1]=chargeCath1; }

  /// Set chi2 of the RawCharge fit   
  void setChi2(float chi2) { fChi2 = chi2; }

  /// Set (x,y) of that cluster and errors   
  void setPosition(const TVector2& pos, const TVector2& errorOnPos) 
  { fHasPosition = kTRUE; fPosition = pos; fPositionError = errorOnPos; }
  
  int cathode() const;
  
  void addCluster(const Cluster& cluster);

  void clear(Option_t* opt="");
  
  bool isMonoCathode() const;

//private:
    void dumpMe() const;
  
private:
  std::vector<MUONPads> mPads; ///< AliMUONPad(s) composing this cluster
  bool mHasPosition; ///< false for pre-cluster (i.e. not yet computed)
  TVector2 mPosition; ///< (x,y) of that cluster (only valid if fHasPosition is kTRUE)
  TVector2 mPositionError; ///< errors on (x,y)
  int mMultiplicity[2]; ///< number of pads in each cathode
  float mRawCharge[2]; ///< cathode RawCharges
  bool mHasCharge; ///< false if SetCharge has not been called
  float mCharge[2]; ///< cathode (re)computed charges
  float mChi2; ///< chi2 of the RawCharge fit (if any)
  bool mIsSaturated[2]; ///< saturation status of cathodes
  
  ClassDef(Cluster,3) // A cluster of AliMUONPad
};

#endif
