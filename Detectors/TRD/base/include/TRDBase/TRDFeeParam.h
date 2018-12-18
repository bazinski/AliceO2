
// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef O2_TRDFEEPARAM_H
#define O2_TRDFEEPARAM_H

//Forwards to standard header with protection for GPU compilation
#include "AliTPCCommonRtypes.h" // for ClassDef

namespace o2
{
namespace trd
{

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  TRD front end electronics parameters class                            //
//  Contains all FEE (MCM, TRAP, PASA) related                            //
//  parameters, constants, and mapping.                                   //
//                                                                        //
//  Author:                                                               //
//    Ken Oyama (oyama@physi.uni-heidelberg.de)                           //
//                                                                        //
//  many things now configured by AliTRDtrapConfig reflecting             //
//  the real memory structure of the TRAP (Jochen)                        //
//                                                                        //
////////////////////////////////////////////////////////////////////////////


class TRootIoCtor;

class TRDCommonParam;
class TRDpadPlane;
class TRDgeometry;

//_____________________________________________________________________________
class TRDFeeParam 
{

 public:

  TRDFeeParam(TRootIoCtor *);
  TRDFeeParam(const TRDFeeParam &p);
  virtual           ~TRDFeeParam();
  TRDFeeParam    &operator=(const TRDFeeParam &p);
  virtual void       Copy(TRDFeeParam  &p) const;

  static TRDFeeParam *Instance();  // Singleton
  static void            Terminate();

  // Translation from MCM to Pad and vice versa
  virtual Int_t    getPadRowFromMCM(Int_t irob, Int_t imcm) const;
  virtual Int_t    getPadColFromADC(Int_t irob, Int_t imcm, Int_t iadc) const;
  virtual Int_t    getExtendedPadColFromADC(Int_t irob, Int_t imcm, Int_t iadc) const;
  virtual Int_t    getMCMfromPad(Int_t irow, Int_t icol) const;
  virtual Int_t    getMCMfromSharedPad(Int_t irow, Int_t icol) const;
  virtual Int_t    getROBfromPad(Int_t irow, Int_t icol) const;
  virtual Int_t    getROBfromSharedPad(Int_t irow, Int_t icol) const;
  virtual Int_t    getRobSide(Int_t irob) const;
  virtual Int_t    getColSide(Int_t icol) const;

  // SCSN-related
  static  UInt_t   AliToExtAli(Int_t rob, Int_t aliid);  // Converts the MCM-ROB combination to the extended MCM ALICE ID (used to address MCMs on the SCSN Bus)
  static  Int_t    ExtAliToAli( UInt_t dest, UShort_t linkpair, UShort_t rocType, Int_t *list, Int_t listSize);  // translates an extended MCM ALICE ID to a list of MCMs
  static  Short_t  ChipmaskToMCMlist( UInt_t cmA, UInt_t cmB, UShort_t linkpair, Int_t *mcmList, Int_t listSize );
  static  Short_t  getRobAB( UShort_t robsel, UShort_t linkpair );  // Returns the chamber side (A=0, B=0) of a ROB

  // geometry
  static  Float_t  getSamplingFrequency() { return (Float_t)fgkLHCfrequency / 4000000.0; }
  static  Int_t    getNmcmRob()           { return fgkNmcmRob;      }
  static  Int_t    getNmcmRobInRow()      { return fgkNmcmRobInRow; }
  static  Int_t    getNmcmRobInCol()      { return fgkNmcmRobInCol; }
  static  Int_t    getNrobC0()            { return fgkNrobC0;       }
  static  Int_t    getNrobC1()            { return fgkNrobC1;       }
  static  Int_t    getNadcMcm()           { return fgkNadcMcm;      }
  static  Int_t    getNcol()              { return fgkNcol;         }
  static  Int_t    getNcolMcm()           { return fgkNcolMcm;      }
  static  Int_t    getNrowC0()            { return fgkNrowC0;       }
  static  Int_t    getNrowC1()            { return fgkNrowC1;       }

  // tracklet simulation
	  Bool_t   getTracklet()         const { return fgTracklet; } 
  static  void     setTracklet(Bool_t trackletSim = kTRUE) { fgTracklet = trackletSim; }
          Bool_t   getRejectMultipleTracklets() const { return fgRejectMultipleTracklets; }
  static  void     setRejectMultipleTracklets(Bool_t rej = kTRUE) { fgRejectMultipleTracklets = rej; }
          Bool_t   getUseMisalignCorr()            const { return fgUseMisalignCorr; }
  static  void     setUseMisalignCorr(Bool_t misalign = kTRUE) { fgUseMisalignCorr = misalign; }
          Bool_t   getUseTimeOffset()              const { return fgUseTimeOffset; }
  static  void     setUseTimeOffset(Bool_t timeOffset = kTRUE) { fgUseTimeOffset = timeOffset; }

  // Concerning raw data format
          Int_t    getRAWversion() const                    { return fRAWversion;        }
          void     setRAWversion( Int_t rawver );

 protected:

  static TRDFeeParam *fgInstance;         // Singleton instance
  static Bool_t          fgTerminated;       // Defines if this class has already been terminated

  TRDCommonParam     *fCP;                // TRD common parameters class

  // Remark: ISO C++ allows initialization of static const values only for integer.

  // Basic Geometrical numbers
  static const Int_t    fgkLHCfrequency      = 40079000 ; // [Hz] LHC clock (should be moved to STEER?)
  static const Int_t    fgkNmcmRob           = 16;        // Number of MCMs per ROB         (old fgkMCMmax)
  static const Int_t    fgkNmcmRobInRow      = 4;         // Number of MCMs per ROB in row dir. (old fgkMCMrow)
  static const Int_t    fgkNmcmRobInCol      = 4;         // Number of MCMs per ROB in col dir. (old fgkMCMrow)
  static const Int_t    fgkNrobC0            = 6;         // Number of ROBs per C0 chamber  (old fgkROBmaxC0)
  static const Int_t    fgkNrobC1            = 8;         // Number of ROBs per C1 chamber  (old fgkROBmaxC1)
  static const Int_t    fgkNadcMcm           = 21;        // Number of ADC channels per MCM (old fgkADCmax)
  static const Int_t    fgkNcol              = 144;       // Number of pads per padplane row(old fgkColmax)
  static const Int_t    fgkNcolMcm           = 18;        // Number of pads per MCM         (old fgkPadmax)
  static const Int_t    fgkNrowC0            = 12;        // Number of Rows per C0 chamber  (old fgkRowmaxC0)
  static const Int_t    fgkNrowC1            = 16;        // Number of Rows per C1 chamber  (old fgkRowmaxC1)

 // Tracklet  processing on/off 
  static       Bool_t   fgTracklet; // tracklet processing
  static       Bool_t   fgRejectMultipleTracklets; // only accept best tracklet if found more than once
  static       Bool_t   fgUseMisalignCorr; // add correction for mis-alignment in y
  static       Bool_t   fgUseTimeOffset; // add time offset in calculation of fit sums

  // For raw production
               Int_t    fRAWversion;                      // Raw data production version
  static const Int_t    fgkMaxRAWversion      = 3;        // Maximum raw version number supported

 private:

  TRDFeeParam();

  ClassDef(TRDFeeParam,4)                              // The TRD front end electronics parameter

};

}//namespace trd
}//namespace o2
#endif

