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

#include <strength_plastic.h>
#include <domain.h>
#include <error.h>
#include <input.h>
#include <mpm_math.h>
#include <universe.h>
#include <update.h>
#include <var.h>
#include <matrix.h>
#include <iostream>

using namespace std;

using namespace MPM_Math;


StrengthPlastic::StrengthPlastic(MPM *mpm, vector<string> args) : Strength(mpm, args)
{
  if (universe->me == 0) {
    cout << "Initiate StrengthPlastic" << endl;
  }

  if (args.size() < 3) {
    error->all(FLERR, "Error: too few arguments for the strength command.\n");
  }

  if (args[2].compare("restart") ==
      0) { // If the keyword restart, we are expecting to have read_restart()
           // launched right after.
    G_, yieldStress = 0;
    return;
  }

  //options(&args, args.begin()+3);
  G_ = input->parsev(args[2]);
  yieldStress = input->parsev(args[3]);
  if (universe->me == 0) {
    cout << "Linear plastic strength model:\n";
    cout << "\tG: shear modulus " << G_ << endl;

    cout << "\t yield stress: " << yieldStress << endl;
  }
}

double StrengthPlastic::G(){
  return G_;
}

Matrix3d StrengthPlastic::update_deviatoric_stress(const Matrix3d& sigma,
							const Matrix3d& D,
							double &               plastic_strain_increment,
							const double           eff_plastic_strain,
							const double           epsdot,
							const double           damage,
							const double           temperature)
{
  Matrix3d sigmaInitial_dev, sigmaFinal_dev, sigmaTrial_dev, dev_rate;
  double J2, Gd, yieldStressD;

  /*
   * deviatoric rate of unrotated stress
   */
  Gd = G_*(1-damage);
  yieldStressD = yieldStress*(1-damage);
  dev_rate = 2.0 * Gd * Deviator(D);
  sigmaInitial_dev = Deviator(sigma);

  /*
   * perform a trial elastic update to the deviatoric stress
   */
  sigmaTrial_dev = sigmaInitial_dev + update->dt * dev_rate; // increment stress deviator using deviatoric rate

  /*
   * check yield condition
   */
  J2 = sqrt(3. / 2.) * sigmaTrial_dev.norm();
  sigmaFinal_dev = sigmaTrial_dev;

  if (J2 < yieldStressD) {
    /*
     * no yielding has occured.
     * final deviatoric stress is trial deviatoric stress
     */
    plastic_strain_increment = 0.0;

  } else {
    //printf("yiedl\n");
    /*
     * yielding has occured
     */

    plastic_strain_increment = (J2 - yieldStressD) / (3.0 * Gd);
    /*
     * new deviatoric stress:
     * obtain by scaling the trial stress deviator
     */
    sigmaFinal_dev *= (yieldStressD / J2);

  }

  return sigmaFinal_dev;
}

void StrengthPlastic::write_restart(ofstream *of) {
  of->write(reinterpret_cast<const char *>(&G_), sizeof(double));
  of->write(reinterpret_cast<const char *>(&yieldStress), sizeof(double));
}

void StrengthPlastic::read_restart(ifstream *ifr) {
  if (universe->me == 0) {
    cout << "Restart StrengthPlastic" << endl;
  }
  ifr->read(reinterpret_cast<char *>(&G_), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&yieldStress), sizeof(double));
}
