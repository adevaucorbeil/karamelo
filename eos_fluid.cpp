#include <iostream>
#include "eos_fluid.h"
#include "input.h"
#include "domain.h"
#include "mpm_math.h"
#include "math_special.h"
#include <Eigen/Eigen>
#include "var.h"
#include <math.h>

using namespace std;
using namespace Eigen;
using namespace MPM_Math;
using namespace MathSpecial;


EOSFluid::EOSFluid(MPM *mpm, vector<string> args) : EOS(mpm, args)
{
  cout << "Initiate EOSFluid" << endl;

  if (args.size()<55) {
    cout << "Error: region command not enough arguments" << endl;
    exit(1);
  }
  //options(&args, args.begin()+3);
  rho0_ = input->parsev(args[2]);
  cout << "Set rho0 to " << rho0_ << endl;

  K_ = input->parsev(args[3]);
  cout << "Set K to " << K_ << endl;

  Gamma = input->parsev(args[4]);
  cout << "Set gamma to " << Gamma << endl;
}


EOSFluid::~EOSFluid()
{

}

double EOSFluid::rho0(){
  return rho0_;
}

double EOSFluid::K(){
  return K_;
}

double EOSFluid::compute_pressure(const double J, const double rho, const double e, const double damage){
  double mu = rho / rho0_;
  double pH = K_ * (pow(mu, Gamma) - 1.0);

  return pH;
}

