/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_TEMPERATURE_H
#define MPM_TEMPERATURE_H

#include <pointers.h>
#include <vector>
#include <matrix.h>

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
  virtual double cp() = 0;
  virtual double kappa() = 0;
  virtual double alpha() = 0;
  virtual void compute_heat_source(double, double &, const double &, const double &) = 0;

protected:
  double cp_, kappa_, alpha_ = 0;
};

#endif
