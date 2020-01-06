/* -*- c++ -*- ----------------------------------------------------------*/


#ifdef REGION_CLASS

RegionStyle(sphere,RegSphere)

#else

#ifndef MPM_REGION_SPHERE_H
#define MPM_REGION_SPHERE_H

#include "region.h"

class RegSphere : public Region {

 public:
  RegSphere(class MPM *, vector<string>);
  ~RegSphere();
  int inside(double, double, double);
  vector<double> limits();

 protected:
  double c1, c2, c3, R, RSq;
  double xlo, xhi, ylo, yhi, zlo, zhi;
  string usage[3] = {"Usage in 1D: region(region-ID, sphere, center_x, R)\n",
		     "Usage in 3D: region(region-ID, sphere, center_x, center_y, R)\n",
		     "Usage in 3D: region(region-ID, sphere, center_x, center_y, center_z, R)\n"};
  int Nargs[3] = {4, 5, 6};
};

#endif
#endif
