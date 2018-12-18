// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef O2_TRDARRAYADC_H
#define O2_TRDARRAYADC_H


namespace o2
{
namespace trd
{


///////////////////////////////////////////////
//                                           //
// Container class for ADC values            //
//                                           // 
///////////////////////////////////////////////

// ported from cblume original code in AliRoot, possibly some modifications
//
#include "Rtypes.h"
#include "TObject.h"

class TRDSignalIndex;

class TRDArrayADC : public TObject
{
 public:

  enum {
    kDataInvalid = 14
  };

  TRDArrayADC();
  TRDArrayADC(Int_t nrow, Int_t ncol, Int_t ntime);
  TRDArrayADC(const TRDArrayADC &b);
  ~TRDArrayADC();
  TRDArrayADC &operator=(const TRDArrayADC &b);

  void    Allocate(Int_t nrow, Int_t ncol, Int_t ntime);
  void    setNdet(Int_t ndet) {fNdet=ndet;};  
  Int_t   getNdet()  const {return fNdet;};
  void    setDataByAdcCol(Int_t nrow, Int_t ncol, Int_t ntime, Short_t value)
                         {fADC[(nrow*fNumberOfChannels+ncol)*fNtime+ntime]=value;}
  Bool_t  HasData() const {return fNtime ? 1 : 0;};
  Short_t getDataByAdcCol(Int_t nrow, Int_t ncol, Int_t ntime) const
                         {return fADC[(nrow*fNumberOfChannels+ncol)*fNtime+ntime];};
  inline  void getData(Int_t r, Int_t c, Int_t t, Int_t n, Short_t *vals) const;
  Short_t getDataBits(Int_t nrow, Int_t ncol, Int_t ntime) const;
  UChar_t getPadStatus(Int_t nrow, Int_t ncol, Int_t ntime) const;
  void    setPadStatus(Int_t nrow, Int_t ncol, Int_t ntime, UChar_t status);
  Bool_t  IsPadCorrupted(Int_t nrow, Int_t ncol, Int_t ntime);
  void    Compress();
  void    Expand();
  Int_t   getNtime() const {return fNtime;};
  Int_t   getNrow() const {return fNrow;};
  Int_t   getNcol() const {return fNcol;};
  Int_t   getDim() const {return fNAdim;};
  void    DeleteNegatives();
  void    Reset();
  void    ConditionalReset(TRDSignalIndex* idx);
  inline  Short_t* getDataAddress(Int_t nrow, Int_t ncol, Int_t ntime = 0) const;
  inline  Short_t getData(Int_t nrow, Int_t ncol, Int_t ntime) const;
  inline  void    setData(Int_t nrow, Int_t ncol, Int_t ntime, Short_t value);
  static  void    CreateLut(); 
  Bool_t  IsValid() { return !TestBit(kDataInvalid); }
  void    SetDataInvalid() { SetBit(kDataInvalid); }
  void    SetDataValid() { ResetBit(kDataInvalid); }

 protected:
  
  Int_t fNdet;    //ID number of the chamber
  Int_t fNrow;    //Number of rows
  Int_t fNcol;    //Number of columns(pads)
  Int_t fNumberOfChannels;  //  Number of MCM channels per row
  Int_t fNtime;   //Number of time bins
  Int_t fNAdim;   //Dimension of the ADC array
//  std::vector<Short_t> fADC;  //[fNAdim]   //Pointer to adc values
  Short_t *fADC;  //[fNAdim]   //Pointer to adc values
  //static std::vector<Short_t> fgLutPadNumbering;   //  [fNcol] Look Up Table
  static Short_t *fgLutPadNumbering;   //  [fNcol] Look Up Table
  ClassDef(TRDArrayADC,1) //ADC container class
    
};

//________________________________________________________________________________
Short_t* TRDArrayADC::getDataAddress(Int_t nrow, Int_t ncol, Int_t ntime) const
{
  //
  // get the address of the given pad
  //

  Int_t corrcolumn = fgLutPadNumbering[ncol];

  return &fADC[(nrow*fNumberOfChannels+corrcolumn)*fNtime+ntime];

}
//________________________________________________________________________________
Short_t TRDArrayADC::getData(Int_t nrow, Int_t ncol, Int_t ntime) const
{
  //
  // get the data using the pad numbering.
  // To access data using the mcm scheme use instead
  // the method getDataByAdcCol
  //

  Int_t corrcolumn = fgLutPadNumbering[ncol];

  return fADC[(nrow*fNumberOfChannels+corrcolumn)*fNtime+ntime];

}
//________________________________________________________________________________
void TRDArrayADC::setData(Int_t nrow, Int_t ncol, Int_t ntime, Short_t value)
{
  //
  // set the data using the pad numbering.
  // To write data using the mcm scheme use instead
  // the method setDataByAdcCol
  //

  Int_t colnumb = fgLutPadNumbering[ncol];

  fADC[(nrow*fNumberOfChannels+colnumb)*fNtime+ntime] = value;

}

void TRDArrayADC::getData(Int_t r, Int_t c, Int_t t, Int_t n, Short_t *vals) const
{
  Int_t colNum = fgLutPadNumbering[c];
  for(Int_t ic=n, idx = (r*fNumberOfChannels+colNum)*fNtime+t; ic--; idx+=fNtime) vals[ic] = fADC[idx];
 }

} // end of namespace trd
} // end of namespace o2
#endif 

