#include <strength_fluid.h>
#include <domain.h>
#include <error.h>
#include <input.h>
#include <mpm_math.h>
#include <universe.h>
#include <update.h>
#include <var.h>
#include <matrix.h>
#include <iostream>

using namespace std;

using namespace MPM_Math;


StrengthFluid::StrengthFluid(MPM *mpm, vector<string> args) : Strength(mpm, args)
{
  if (universe->me == 0) {
    cout << "Initiate StrengthFluid" << endl;
  }

  if (args.size() < 3) {
    error->all(FLERR, "Error: too few arguments for the strength command.\n");
  }

  if (args[2].compare("restart") ==
      0) { // If the keyword restart, we are expecting to have read_restart()
           // launched right after.
    G_ = 0;
    return;
  }

  //options(&args, args.begin()+3);
  G_ = input->parsev(args[2]);
  if (universe->me == 0) {
    cout << "Fluid strength model:\n";
    cout << "\tmu: shear viscosity " << G_ << endl;
  }
}

double StrengthFluid::G(){
  return G_;
}

Matrix3d StrengthFluid::update_deviatoric_stress(const Matrix3d& sigma,
						 const Matrix3d& D,
						 double &plastic_strain_increment,
						 const double eff_plastic_strain,
						 const double epsdot,
						 const double damage,
						 const double temperature)
{
  Matrix3d dev_rate;

  dev_rate = 2.0 * G_ * Deviator(D);
  return dev_rate;
}

void StrengthFluid::write_restart(ofstream *of) {
  of->write(reinterpret_cast<const char *>(&G_), sizeof(double));
}

void StrengthFluid::read_restart(ifstream *ifr) {
  if (universe->me == 0) {
    cout << "Restart StrengthFluid" << endl;
  }
  ifr->read(reinterpret_cast<char *>(&G_), sizeof(double));
}
