#include <iostream>
#include "strength_linear.h"
#include "input.h"
#include "domain.h"
#include "update.h"
#include "mpm_math.h"
#include <Eigen/Eigen>
#include "var.h"

using namespace std;
using namespace Eigen;
using namespace MPM_Math;


StrengthLinear::StrengthLinear(MPM *mpm, vector<string> args) : Strength(mpm, args)
{
  cout << "Initiate StrengthLinear" << endl;

  if (args.size()<2) {
    cout << "Error: too few arguments for the strength command" << endl;
    exit(1);
  }
  //options(&args, args.begin()+3);
  G_ = input->parsev(args[2]);
  cout << "Linear strength model:\n";
  cout << "\tG: shear modulus " << G_ << endl;
}

double StrengthLinear::G(){
  return G_;
}

Matrix3d StrengthLinear::update_deviatoric_stress(const Eigen::Matrix3d& sigma,
							const Eigen::Matrix3d& D,
							double &               plastic_strain_increment,
							const double           eff_plastic_strain,
							const double           epsdot,
							const double           damage,
							const double           temperature)
{
  Matrix3d dev_rate;

  dev_rate = 2.0 * G_* (1-damage) * Deviator(D);
  return Deviator(sigma) + update->dt * dev_rate;
}
