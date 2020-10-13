/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_TEMPERATURE_H
#define MPM_TEMPERATURE_H

#include "pointers.h"
#include <vector>
#include <Eigen/Eigen>

class Temperature : protected Pointers {
 public:
  string id;                        ///< Temperature identification string
  string style;                     ///< Temperature style

  Temperature(MPM *, vector<string>);
  virtual ~Temperature();
  virtual void init();
  void options(vector<string> *, vector<string>::iterator);

  virtual void write_restart(ofstream*) = 0;
  virtual void read_restart(ifstream*) = 0;

  // implemented by each Temperature
  virtual void compute_temperature(double &, const double &, const double &) = 0;
};

#endif
