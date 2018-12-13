
// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef O2_TRDCALSINGLECHAMBERSTATUS_H
#define O2_TRDCALSINGLECHAMBERSTATUS_H

namespace o2
{
namespace trd
{


////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  TRD calibration base class containing status values for one ROC       //
//                                                                        //
////////////////////////////////////////////////////////////////////////////


//_____________________________________________________________________________
class TRDCalSingleChamberStatus

 public:

  enum class { kMasked          = 2
             , kPadBridgedLeft  = 4
             , kPadBridgedRight = 8
             , kReadSecond      = 16 
             , kNotConnected    = 32};

  TRDCalSingleChamberStatus();
  TRDCalSingleChamberStatus(Int_t p, Int_t c, Int_t cols);
  TRDCalSingleChamberStatus(const TRDCalSingleChamberStatus &c);
  virtual                      ~TRDCalSingleChamberStatus();
  TRDCalSingleChamberStatus &operator=(const TRDCalSingleChamberStatus &c);


          Bool_t  IsMasked(Int_t col, Int_t row) const       { return ((GetStatus(col,row) & kMasked) 
                                                                       ? kTRUE 
                                                                       : kFALSE);                 };
	  Bool_t  IsBridgedLeft(Int_t col, Int_t row) const  { return ((GetStatus(col, row) & kPadBridgedLeft)                                            ? kTRUE                                                                         : kFALSE);                 };
	  Bool_t  IsBridgedRight(Int_t col, Int_t row) const { return ((GetStatus(col, row) & kPadBridgedRight)                                           ? kTRUE                                                                         : kFALSE);                 };
	  Bool_t  IsNotConnected(Int_t col, Int_t row) const { return ((GetStatus(col, row) & kNotConnected)                                           ? kTRUE                                                                         : kFALSE);                 };
          Int_t   GetNrows() const                           { return fNrows;                     };
          Int_t   GetNcols() const                           { return fNcols;                     };

          Int_t   GetChannel(Int_t col, Int_t row) const     { return row+col*fNrows;             };
          Int_t   GetNchannels() const                       { return fNchannels;                 };
          Char_t  GetStatus(Int_t ich) const                 { return fData[ich];                 };
          Char_t  GetStatus(Int_t col, Int_t row) const      { return fData[GetChannel(col,row)]; };

          void    SetStatus(Int_t ich, Char_t vd)            { fData[ich] = vd;                   };
          void    SetStatus(Int_t col, Int_t row, Char_t vd) { fData[GetChannel(col,row)] = vd;   };

 protected:

          Int_t   fPla;                    //  Plane number
          Int_t   fCha;                    //  Chamber number

          Int_t   fNrows;                  //  Number of rows
          Int_t   fNcols;                  //  Number of columns

          Int_t   fNchannels;              //  Number of channels
          Char_t *fData;                   //[fNchannels] Data

  ClassDef(TRDCalSingleChamberStatus,1) //  TRD ROC calibration class

};

}// end trd namespace
}//end o2 namespace
#endif
