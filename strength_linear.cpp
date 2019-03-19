#include <iostream>
#include "strength_linear.h"
#include "input.h"
#include "domain.h"
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
    cout << "Error: region command not enough arguments" << endl;
    exit(1);
  }
  //options(&args, args.begin()+3);
  G_ = input->parsev(args[2]);
  cout << "Set G to " << G_ << endl;
}

double StrengthLinear::G(){
  return G_;
}

Matrix3d StrengthLinear::update_deviatoric_stress(const Matrix3d strain_increment, const Matrix3d sigma)
{
  return Deviator(sigma) + 2 * G_ * Deviator(strain_increment);
}
