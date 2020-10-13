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
    chi, rho, cp, alpha = 0;
    return;
  }

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

void TemperaturePlasticWork::write_restart(ofstream *of) {
  of->write(reinterpret_cast<const char *>(&chi), sizeof(double));
  of->write(reinterpret_cast<const char *>(&rho), sizeof(double));
  of->write(reinterpret_cast<const char *>(&cp), sizeof(double));
}

void TemperaturePlasticWork::read_restart(ifstream *ifr) {
  cout << "Restart TemperaturePlasticWork" << endl;
  ifr->read(reinterpret_cast<char *>(&chi), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&rho), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&cp), sizeof(double));
  alpha = chi/(rho*cp);
}
