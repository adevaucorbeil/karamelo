#include <iostream>
#include "eos_linear.h"
#include "input.h"
#include "domain.h"
#include "mpm_math.h"
#include <Eigen/Eigen>
#include "var.h"

using namespace std;
using namespace Eigen;
using namespace MPM_Math;


EOSLinear::EOSLinear(MPM *mpm, vector<string> args) : EOS(mpm, args)
{
  cout << "Initiate EOSLinear" << endl;

  if (args.size()<3) {
    cout << "Error: region command not enough arguments" << endl;
    exit(1);
  }
  //options(&args, args.begin()+3);
  rho0_ = input->parsev(args[2]);
  cout << "Set rho0 to " << rho0_ << endl;
  K_ = input->parsev(args[3]);
  cout << "Set K to " << K_ << endl;
}


EOSLinear::~EOSLinear()
{

}

double EOSLinear::rho0(){
  return rho0_;
}

double EOSLinear::K(){
  return K_;
}

double EOSLinear::compute_pressure(const double J, const double rho, const double e, const double damage){
  return K_*(1-J)*(1-damage);
}


