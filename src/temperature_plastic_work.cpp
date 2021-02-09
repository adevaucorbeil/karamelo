#include "temperature_plastic_work.h"
#include "domain.h"
#include "error.h"
#include "input.h"
#include "mpm_math.h"
#include "update.h"
#include "var.h"
#include <Eigen/Eigen>
#include <iostream>

using namespace std;
using namespace Eigen;
using namespace MPM_Math;

#define SQRT_3_OVER_2 1.224744871 // sqrt(3.0/2.0)

TemperaturePlasticWork::TemperaturePlasticWork(MPM *mpm, vector<string> args) : Temperature(mpm, args)
{
  cout << "Initiate TemperaturePlasticWork" << endl;

  if (args.size() < Nargs) {
    error->all(FLERR, "Error: not enough arguments.\n" + usage);
  }
  if (args.size() > Nargs) {
    error->all(FLERR, "Error: too many arguments.\n" + usage);
  }

  //options(&args, args.begin()+3);
  chi = input->parsev(args[2]);
  cp_ = input->parsev(args[3]);
  kappa_ = input->parsev(args[4]);
  alpha = input->parsev(args[5]);
  T0 = input->parsev(args[6]);

  cout << "Plastic work material temperature model:\n";
  cout << "\tTaylor-Quinney coefficient chi:" << chi << endl;
  cout << "\tSpecific heat at constant pressure cp:" << cp_ << endl;
  cout << "\tThermal conductivity kappa:" << kappa_ << endl;
  cout << "\tCoefficient of thermal expansion alpha:" << alpha << endl;
  cout << "\tInitial temperature T0:" << T0 << endl;  
}

void TemperaturePlasticWork::compute_heat_source(double &gamma, const double &flow_stress, const double &eff_plastic_strain_rate)
{
  gamma = chi * flow_stress * eff_plastic_strain_rate;
}

double TemperaturePlasticWork::compute_thermal_pressure(double T) {
  return alpha * (T0 - T);
}
