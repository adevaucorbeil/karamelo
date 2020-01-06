/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_TEMPERATURE_H
#define MPM_TEMPERATURE_H

#include "pointers.h"
#include <vector>
#include <Eigen/Eigen>

class Temperature : protected Pointers {
 public:
  string id;

  Temperature(MPM *, vector<string>);
  virtual ~Temperature();
  virtual void init();
  void options(vector<string> *, vector<string>::iterator);

  // implemented by each Temperature
  virtual void compute_temperature(double &, const double &, const double &) = 0;
};

#endif
