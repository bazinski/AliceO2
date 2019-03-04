// Class MpArea
// ----------------
// Class that defines a rectangle area positioned in plane..
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay
//-----------------------------------------------------------------------------

#include "MCHClustering/MpArea.h"

#include "AliMpConstants.h"

#include <Riostream.h>

using namespace std;
using namespace o2::mch;

//_____________________________________________________________________________
MpArea::MpArea(double x, double y, 
                     double dx, double dy)
  : mPositionX(x),
    mPositionY(y),
    mDimensionX(dx),
    mDimensionY(dy),
    mValidity(true) 
{
/// Standard constructor

  // Check dimensions
  if (  mDimensionX < - AliMpConstants::lengthTolerance() || 
        mDimensionY < - AliMpConstants::lengthTolerance() || 
      ( mDimensionX < AliMpConstants::lengthTolerance() && 
        mDimensionY < AliMpConstants::lengthTolerance() ) )
  {
    mDimensionX = 0.;
    mDimensionY = 0.;
    mValidity = false;
  }  
}

//_____________________________________________________________________________
MpArea::MpArea() = default;

//_____________________________________________________________________________
MpArea::MpArea(const MpArea& rhs):
  mPositionX(rhs.fPositionX),
  mPositionY(rhs.fPositionY),
  mDimensionX(rhs.fDimensionX), 
  mDimensionY(rhs.fDimensionY), 
  mValidity(rhs.fValidity) 
{
/// Copy constructor
}

//_____________________________________________________________________________
MpArea::~MpArea() = default;

//
// operators
//

//______________________________________________________________________________
MpArea& MpArea::operator = (const MpArea& right)
{
/// Assignment operator

  // check assignment to self
  if (this == &right) return *this;

  mPositionX = right.mPositionX;
  mPositionY = right.mPositionY;
  mDimensionX = right.mDimensionX;
  mDimensionY = right.mDimensionY;
  mValidity = right.mValidity;

  return *this;
} 

//
// public methods
//

//_____________________________________________________________________________
double MpArea::leftBorder() const
{
/// Return the position of the left edge.

  return mPositionX - mDimensionX;
}

//_____________________________________________________________________________
double MpArea::rightBorder() const
{
/// Return the position of right edge.

  return mPositionX + mDimensionX;
}

//_____________________________________________________________________________
double MpArea::upBorder() const
{
/// Return the position of the up edge.

  return mPositionY + mDimensionY;
}

//_____________________________________________________________________________
double MpArea::downBorder() const
{
/// Return the position of the down edge.

  return mPositionY - mDimensionY;
}

//_____________________________________________________________________________
void MpArea::leftDownCorner(double& x, double& y) const
{
/// Return position of the left down corner.

  x = leftBorder();
  y = downBorder();
}  

//_____________________________________________________________________________
void MpArea::leftUpCorner(double& x, double& y) const
{
/// Return position of the left up corner.

  x = leftBorder();
  y = upBorder();
}  

//_____________________________________________________________________________
void MpArea::rightDownCorner(double& x, double& y) const
{
/// Return position of the right down corner.

  x = rightBorder();
  y = downBorder();
}  


//_____________________________________________________________________________
void MpArea::rightUpCorner(double& x, double& y) const
{
/// Return position of the right up corner.

  x = rightBorder();
  y = upBorder();
}  

//_____________________________________________________________________________
bool MpArea::contains(const MpArea& area) const
{
/// Whether area is contained within this
  
  if ( area.leftBorder() < leftBorder() ||
       area.rightBorder() > rightBorder() ||
       area.downBorder() < downBorder() ||
       area.upBorder() > upBorder() ) 
  {
    return kFALSE;
  }
  else
  {
    return kTRUE;
  }
}

//_____________________________________________________________________________
MpArea MpArea::intersect(const MpArea& area) const
{ 
/// Return the common part of area and this

  double xmin = TMath::Max(area.leftBorder(),leftBorder());
  double xmax = TMath::Min(area.rightBorder(),rightBorder());
  double ymin = TMath::Max(area.downBorder(),downBorder());
  double ymax = TMath::Min(area.upBorder(),upBorder());

  return MpArea( (xmin+xmax)/2.0, (ymin+ymax)/2.0 ,
                    (xmax-xmin)/2.0, (ymax-ymin)/2.0 );
}

//_____________________________________________________________________________
bool MpArea::overlap(const MpArea& area) const
{
/// Return true if this overlaps with given area

  if ( leftBorder() > area.rightBorder() - AliMpConstants::lengthTolerance() ||
       rightBorder() < area.leftBorder() + AliMpConstants::lengthTolerance() )
  {
    return kFALSE;
  }

  if ( downBorder() > area.upBorder() - AliMpConstants::lengthTolerance() ||
       upBorder() < area.downBorder() + AliMpConstants::lengthTolerance() )
  {
    return kFALSE;
  }
  return kTRUE;
  
}

//_____________________________________________________________________________
void
MpArea::print(Option_t* opt) const
{
/// Printing
/// When option is set to B (borders), the area boreders will be printed 
/// instead of default parameters

  
  if ( opt[0] == 'B' ) {
    cout << "Area x-borders: (" 
         << leftBorder() << ", " << rightBorder() << ") " 
	 << " y-borders: (" 
         << downBorder() << ", " << upBorder() << ") " 
	 << endl;
    return;

  }       

  cout << (*this) << endl;
}

//_____________________________________________________________________________
void      
MpArea::getParameters(double& x, double& y,
                         double& dx, double& dy) const
{
/// Fill the parameters: x, y position and dimensions
                         
  x = fPositionX;
  y = fPositionY;
  dx = fDimensionX;
  dy = fDimensionY;
}  

//_____________________________________________________________________________
ostream& operator<< (ostream &stream,const MpArea& area)
{
/// Output streaming

  stream << "Area: position: (" 
         << area.GetPositionX() << ", " << area.GetPositionY() << ") " 
	 << " dimensions: (" 
         << area.GetDimensionX() << ", " << area.GetDimensionY() << ") " 
  << " valid: " << (area.IsValid()==true ? "YES":"NO")
	 << endl;
  return stream;
}

