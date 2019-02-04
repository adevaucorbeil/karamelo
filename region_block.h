/* -*- c++ -*- ----------------------------------------------------------*/


#ifdef REGION_CLASS

RegionStyle(block,RegBlock)

#else

#ifndef MPM_REGION_BLOCK_H
#define MPM_REGION_BLOCK_H

#include "region.h"

class RegBlock : public Region {

 public:
  RegBlock(class MPM *, vector<string>);
  ~RegBlock();
  int inside(double, double, double);
  vector<double> limits();

 protected:
  double xlo,xhi,ylo,yhi,zlo,zhi;
};

#endif
#endif
