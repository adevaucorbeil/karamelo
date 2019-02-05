#include <iostream>
#include "eos_linear.h"
#include "input.h"

using namespace std;


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

double EOSLinear::compute_pressure(double mu){
  return mu;
}
