#include <iostream>
#include "temperature_plastic_work.h"
#include "input.h"
#include "domain.h"
#include "update.h"
#include "mpm_math.h"
#include <Eigen/Eigen>
#include "var.h"

using namespace std;
using namespace Eigen;
using namespace MPM_Math;

#define SQRT_3_OVER_2 1.224744871 // sqrt(3.0/2.0)

TemperaturePlasticWork::TemperaturePlasticWork(MPM *mpm, vector<string> args) : Temperature(mpm, args)
{
  cout << "Initiate TemperaturePlasticWork" << endl;

  if (args.size() < Nargs) {
    cout << "Error: too few arguments for the temperature command" << endl;
    cout << usage;
    exit(1);
  }
  if (args.size() > Nargs) {
    cout << "Error: too many arguments for the temperature command" << endl;
    cout << usage;
    exit(1);
  }

  //options(&args, args.begin()+3);
  chi = input->parsev(args[2]);
  rho = input->parsev(args[3]);
  cp = input->parsev(args[4]);

  cout << "Plastic work material temperature model:\n";
  cout << "\tTaylor-Quinney coefficient chi:" << chi << endl;
  cout << "\tdensity rho:" << rho << endl;
  cout << "\tSpecific heat at constant pressure cp:" << cp << endl;
  alpha = chi/(rho*cp);
}

void TemperaturePlasticWork::compute_temperature(double &T, const double &flow_stress, const double &plastic_strain_increment)
{
  T += alpha*flow_stress*plastic_strain_increment;
}
