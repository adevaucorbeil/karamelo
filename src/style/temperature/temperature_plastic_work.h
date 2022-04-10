/* -*- c++ -*- ----------------------------------------------------------*/


#ifdef TEMPERATURE_CLASS

TemperatureStyle(plastic_work,TemperaturePlasticWork)

#else

#ifndef MPM_TEMPERATURE_PLASTIC_WORK_H
#define MPM_TEMPERATURE_PLASTIC_WORK_H

#include <temperature.h>
#include <matrix.h>

class TemperaturePlasticWork : public Temperature {

public:
  TemperaturePlasticWork(class MPM *, vector<string>);
  ~TemperaturePlasticWork() {};

  inline double cp() { return cp_; }
  inline double kappa() { return kappa_; }
  inline double alpha() { return alpha_; }
  void compute_heat_source(Solid &solid,
                           Kokkos::View<Matrix3d*, MemorySpace> &sigma_dev) const override;

  void write_restart(ofstream *);
  void read_restart(ifstream *);

protected:
  double chi, T0, Tm;
  string usage = "Usage: temperature(temp-ID, plastic_work, chi, cp, kappa, alpha, T0, Tm)\n";
  int Nargs = 8;
};

#endif
#endif
