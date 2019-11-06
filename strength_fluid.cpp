#include <iostream>
#include "strength_fluid.h"
#include "input.h"
#include "domain.h"
#include "update.h"
#include "mpm_math.h"
#include <Eigen/Eigen>
#include "var.h"

using namespace std;
using namespace Eigen;
using namespace MPM_Math;


StrengthFluid::StrengthFluid(MPM *mpm, vector<string> args) : Strength(mpm, args)
{
  cout << "Initiate StrengthFluid" << endl;

  if (args.size()<2) {
    cout << "Error: too few arguments for the strength command" << endl;
    exit(1);
  }
  //options(&args, args.begin()+3);
  G_ = input->parsev(args[2]);
  cout << "Fluid strength model:\n";
  cout << "\tmu: shear viscosity " << G_ << endl;
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
