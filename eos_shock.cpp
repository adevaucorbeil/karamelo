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
#include "eos_shock.h"
#include "input.h"
#include "domain.h"
#include "mpm_math.h"
#include "math_special.h"
#include <Eigen/Eigen>
#include "var.h"
#include "error.h"

using namespace std;
using namespace Eigen;
using namespace MPM_Math;
using namespace MathSpecial;


EOSShock::EOSShock(MPM *mpm, vector<string> args) : EOS(mpm, args)
{
  cout << "Initiate EOSShock" << endl;

  if (args.size() < Nargs) {
    error->all(FLERR, "Error: not enough arguments.\n" + usage);
  }
  //options(&args, args.begin()+3);
  rho0_ = input->parsev(args[2]);
  cout << "Set rho0 to " << rho0_ << endl;

  K_ = input->parsev(args[3]);
  cout << "Set K to " << K_ << endl;

  c0 = input->parsev(args[4]);
  cout << "Set c0 to " << c0 << endl;

  S = input->parsev(args[5]);
  cout << "Set S to " << S << endl;

  Gamma = input->parsev(args[6]);
  cout << "Set gamma to " << Gamma << endl;

  cv = input->parsev(args[7]);
  cout << "Set cv to " << cv << endl;

  Tr = input->parsev(args[8]);
  cout << "Set Tr to " << Tr << endl;

  alpha = cv*rho0_;
}


EOSShock::~EOSShock()
{

}

double EOSShock::rho0(){
  return rho0_;
}

double EOSShock::K(){
  return K_;
}

void EOSShock::compute_pressure(double &pFinal, double &e, const double J, const double rho, const double T, const double damage){
  double mu = rho / rho0_ - 1.0;
  double pH = rho0_ * square(c0) * mu * (1.0 + mu) / square(1.0 - (S - 1.0) * mu);
  pFinal = (pH + rho * Gamma * (e - e0));

  if ( damage > 0.0 ) {
    if ( pFinal < 0.0 ) {
      if ( damage >= 1.0) {
	pFinal = rho0_ * Gamma * (e - e0);
      } else {
	double mu_damaged = (1.0 - damage) * mu;
	double pH_damaged = rho0_ * (1.0 - damage) * square(c0) * mu_damaged * (1.0 + mu_damaged) / square(1.0 - (S - 1.0) * mu_damaged);
	pFinal = (pH_damaged + rho0_ * (1 + mu_damaged) * Gamma * (e - e0));;
      }
    }
  }

  e = alpha*(T - Tr);
}

