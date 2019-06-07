/* -*- c++ -*- ----------------------------------------------------------*/


#ifdef REGION_CLASS

RegionStyle(cylinder,RegCylinder)

#else

#ifndef MPM_REGION_CYLINDER_H
#define MPM_REGION_CYLINDER_H

#include "region.h"

class RegCylinder : public Region {

 public:
  RegCylinder(class MPM *, vector<string>);
  ~RegCylinder();
  int inside(double, double, double);
  vector<double> limits();

 protected:
  double c1, c2, R, RSq, lo, hi;
  char axis;
  double xlo, xhi, ylo, yhi, zlo, zhi;
};

#endif
#endif
