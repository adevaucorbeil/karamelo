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

  if (args.size() < 3) {
    error->all(FLERR, "Error: too few arguments for the strength command.\n");
  }

  if (args[2].compare("restart") ==
      0) { // If the keyword restart, we are expecting to have read_restart()
           // launched right after.
    chi, cp_, kappa_, alpha, T0, Tm = 0;
    return;
  }

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
  Tm = input->parsev(args[7]);

  cout << "Plastic work material temperature model:\n";
  cout << "\tTaylor-Quinney coefficient chi:" << chi << endl;
  cout << "\tSpecific heat at constant pressure cp:" << cp_ << endl;
  cout << "\tThermal conductivity kappa:" << kappa_ << endl;
  cout << "\tCoefficient of thermal expansion alpha:" << alpha << endl;
  cout << "\tInitial temperature T0:" << T0 << endl;
  cout << "\tMelting temperature Tm:" << Tm << endl;
}

void TemperaturePlasticWork::compute_heat_source(double T, double &gamma, const double &flow_stress, const double &eff_plastic_strain_rate)
{
  if (T < Tm)
    gamma = chi * flow_stress * eff_plastic_strain_rate;
  else
    gamma = 0;
}

double TemperaturePlasticWork::compute_thermal_pressure(double T) {
  return alpha * (T0 - T);
}

void TemperaturePlasticWork::write_restart(ofstream *of) {
  of->write(reinterpret_cast<const char *>(&chi), sizeof(double));
  of->write(reinterpret_cast<const char *>(&kappa_), sizeof(double));
  of->write(reinterpret_cast<const char *>(&cp_), sizeof(double));
  of->write(reinterpret_cast<const char *>(&alpha), sizeof(double));
  of->write(reinterpret_cast<const char *>(&T0), sizeof(double));
  of->write(reinterpret_cast<const char *>(&Tm), sizeof(double));
}

void TemperaturePlasticWork::read_restart(ifstream *ifr) {
  cout << "Restart TemperaturePlasticWork" << endl;
  ifr->read(reinterpret_cast<char *>(&chi), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&kappa_), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&cp_), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&alpha), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&T0), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&Tm), sizeof(double));
}
