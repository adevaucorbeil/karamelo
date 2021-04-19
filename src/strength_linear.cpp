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
#include <Eigen/Eigen>
#include "strength_linear.h"
#include "input.h"
#include "domain.h"
#include "update.h"
#include "mpm_math.h"
#include "var.h"
#include "error.h"

using namespace std;
using namespace Eigen;
using namespace MPM_Math;


StrengthLinear::StrengthLinear(MPM *mpm, vector<string> args) : Strength(mpm, args)
{
  cout << "Initiate StrengthLinear" << endl;

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

void StrengthLinear::write_restart(ofstream *of) {
  of->write(reinterpret_cast<const char *>(&G_), sizeof(double));
}

void StrengthLinear::read_restart(ifstream *ifr) {
  cout << "Restart StrengthLinear" << endl;
  ifr->read(reinterpret_cast<char *>(&G_), sizeof(double));
}

