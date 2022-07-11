/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_TEMPERATURE_H
#define MPM_TEMPERATURE_H

#include <pointers.h>
#include <vector>
#include <matrix.h>
#include <grid.h>

class Solid;

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
  virtual float cp() = 0;
  virtual float kappa() = 0;
  virtual float alpha() = 0;
  virtual void compute_heat_source(Solid &solid,
                                   Kokkos::View<Matrix3d*> &sigma_dev) const = 0;

protected:
  float cp_, kappa_, alpha_ = 0;
};

#endif
