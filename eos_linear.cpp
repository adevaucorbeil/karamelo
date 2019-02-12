#include <iostream>
#include "eos_linear.h"
#include "input.h"
#include "domain.h"
#include "mpm_math.h"
#include <Eigen/Eigen>

using namespace std;
using namespace Eigen;
using namespace MPM_Math;


EOSLinear::EOSLinear(MPM *mpm, vector<string> args) : EOS(mpm, args)
{
  cout << "Initiate EOSLinear" << endl;

  if (args.size()<4) {
    cout << "Error: region command not enough arguments" << endl;
    exit(1);
  }
  //options(&args, args.begin()+3);
  rho0_ = input->parse(args[2]);
  cout << "Set rho0 to " << rho0_ << endl;
  K_ = input->parse(args[3]);
  cout << "Set K to " << K_ << endl;
  G_ = input->parse(args[4]);
  cout << "Set G to " << G_ << endl;
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

double EOSLinear::G(){
  return G_;
}

double EOSLinear::compute_pressure(double J){
  return K_*(1-J);
}

Matrix3d EOSLinear::update_deviatoric_stress(Matrix3d strain_increment, Matrix3d sigma)
{
  return Deviator(sigma) + 2 * G_ * Deviator(strain_increment);
}

void EOSLinear::update_stress(Matrix3d& sigma, Matrix3d strain_increment, double J)
{
  double p = compute_pressure(J);
  Matrix3d sigma_dev = update_deviatoric_stress(strain_increment, sigma);
  
  Matrix3d eye;
  eye.setIdentity();
  sigma = -p*eye + sigma_dev;

  return;
}

