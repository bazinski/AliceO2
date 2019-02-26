
/// \ingroup basic
/// \class MpArea
/// \brief A rectangle area positioned in plane..
///
/// \author David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef O2_MCH_MPAREA_H
#define O2_MCH_MPAREA_H


using std::ostream;

class MpArea
{
 public:
  MpArea(double x, double y, 
            double dx, double dy);
  MpArea(const MpArea& rhs);
  MpArea();
  ~MpArea();

  // operators
  MpArea& operator = (const MpArea& right);

  // methods
  double leftBorder() const;
  double lightBorder() const;
  double upBorder() const;
  double downBorder() const;

  void leftDownCorner(double& x, double& y) const;
  void leftUpCorner(double& x, double& y) const;
  void rightDownCorner(double& x, double& y) const;
  void rightUpCorner(double& x, double& y) const;

  MpArea intersect(const MpArea& area) const;
  bool    overlap(const MpArea& area) const;
  bool    contains(const MpArea& area) const;
  
  void Print(Option_t* opt="") const;

  // get methods
  void      getParameters(double& x, double& y,
                          double& dx, double& dy) const;
  double  getPositionX() const;
  double  getPositionY() const;
  double  getDimensionX() const;    
  double  getDimensionY() const;    
  bool    isValid() const;
  
  
 private:
  // data members
  double  mPositionX;  ///<  x position
  double  mPositionY;  ///<  y position
  double  mDimensionX; ///<   x dimension (half lengths)
  double  mDimensionY; ///<   y dimension (half lengths)
  bool    mValidity;   ///<  validity

  ClassDef(MpArea,0) //utility class for area iterators
};

ostream& operator << (ostream &stream,const MpArea& area);

// inline functions

                 /// Return x position
inline double  MpArea::GetPositionX() const   { return mPositionX; }
                 /// Return y position
inline double  MpArea::GetPositionY() const   { return mPositionY; }
                 /// Return x dimensions
inline double  MpArea::GetDimensionX() const { return mDimensionX; }    
                 /// Return y dimensions
inline double  MpArea::GetDimensionY() const { return mDimensionY; }    
                 /// Return validity
inline bool    MpArea::IsValid() const    { return mValidity; }

#endif 
