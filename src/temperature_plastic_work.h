/* -*- c++ -*- ----------------------------------------------------------*/


#ifdef TEMPERATURE_CLASS

TemperatureStyle(plastic_work,TemperaturePlasticWork)

#else

#ifndef MPM_TEMPERATURE_PLASTIC_WORK_H
#define MPM_TEMPERATURE_PLASTIC_WORK_H

#include "temperature.h"
#include <Eigen/Eigen>

class TemperaturePlasticWork : public Temperature {

public:
  TemperaturePlasticWork(class MPM *, vector<string>);
  ~TemperaturePlasticWork() {};

  void compute_temperature(double &, const double &, const double &);

protected:
  double chi, rho, cp, alpha;
  string usage = "Usage: temperature(temp-ID, plastic_work, chi, rho, cp)\n";
  int Nargs = 5;
};

#endif
#endif
