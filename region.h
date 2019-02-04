/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_REGION_H
#define MPM_REGION_H

#include "pointers.h"
#include <vector>

class Region : protected Pointers {
 public:
  string id;

  Region(class MPM *, vector<string>);
  virtual ~Region();
  virtual void init();
  void options(vector<string> *, vector<string>::iterator);

  // implemented by each region

  virtual vector<double> limits() {return vector<double>();};
  virtual int inside(double, double, double) = 0;
  //protected:
};

#endif
