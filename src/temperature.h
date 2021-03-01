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
  virtual double cp() = 0;
  virtual double kappa() = 0;
  virtual void compute_heat_source(double, double &, const double &, const double &) = 0;
  virtual double compute_thermal_pressure(double) = 0;
};

#endif
