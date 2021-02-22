/* ----------------------------------------------------------------------
 *
 *                    ***       Karamelo       ***
 *               Parallel Material Point Method Simulator
 * 
 * Copyright (2019) Alban de Vaucorbeil, alban.devaucorbeil@monash.edu
 * Materials Science and Engineering, Monash University
 * Clayton VIC 3800, Australia

 * This software is distributed under the GNU General Public License.
 *
 * ----------------------------------------------------------------------- */

#include <iostream>
#include "eos_linear.h"
#include "input.h"
#include "domain.h"
#include "mpm_math.h"
#include <Eigen/Eigen>
#include "var.h"
#include "error.h"

using namespace std;
using namespace Eigen;
using namespace MPM_Math;


EOSLinear::EOSLinear(MPM *mpm, vector<string> args) : EOS(mpm, args)
{
  cout << "Initiate EOSLinear" << endl;

  if (args.size()<3) {
    error->all(FLERR, "Error: not enough arguments.\n");
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

void EOSLinear::compute_pressure(double &pFinal, double &e, const double J, const double rho, const double damage, const Eigen::Matrix3d D, const double cellsize, const double T){
  e = 0;
  pFinal = K_*(1-J)*(1-damage);
}


