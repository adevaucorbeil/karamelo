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

#include "strength_jc.h"
#include "domain.h"
#include "error.h"
#include "input.h"
#include "mpm_math.h"
#include "universe.h"
#include "update.h"
#include "var.h"
#include <Eigen/Eigen>
#include <iostream>

using namespace std;
using namespace Eigen;
using namespace MPM_Math;

#define SQRT_3_OVER_2 1.224744871 // sqrt(3.0/2.0)

StrengthJohnsonCook::StrengthJohnsonCook(MPM *mpm, vector<string> args)
    : Strength(mpm, args)
{
  if (universe->me == 0) {
    cout << "Initiate StrengthJohnsonCook" << endl;
  }

  if (args.size() < 3) {
    error->all(FLERR, "Error: too few arguments for the strength command.\n");
  }

  if (args[2].compare("restart") ==
      0) { // If the keyword restart, we are expecting to have read_restart()
           // launched right after.
    G_, A, B, n, m, epsdot0, C, Tr, Tm, Tmr = 0;
    return;
  }

  if (args.size() < Nargs) {
    error->all(FLERR,
               "Error: too few arguments for the strength command.\n" + usage);
  }

  // options(&args, args.begin()+3);
  G_      = input->parsev(args[2]);
  A       = input->parsev(args[3]);
  B       = input->parsev(args[4]);
  n       = input->parsev(args[5]);
  epsdot0 = input->parsev(args[6]);
  C       = input->parsev(args[7]);
  m       = input->parsev(args[8]);
  Tr      = input->parsev(args[9]);
  Tm      = input->parsev(args[10]);

  if (universe->me == 0) {
    cout << "Johnson Cook material strength model:\n";
    cout << "\tG: shear modulus " << G_ << endl;
    cout << "\tA: initial yield stress " << A << endl;
    cout << "\tB: proportionality factor for plastic strain dependency " << B
         << endl;
    cout << "\tn: exponent for plastic strain dependency " << n << endl;
    cout << "\tepsdot0: reference strain rate " << epsdot0 << endl;
    cout << "\tC: proportionality factor for logarithmic plastic strain rate "
            "dependency "
         << C << endl;
    cout << "\tm: factor for thermal softening " << m << endl;
    cout << "\tTr: reference temperature " << Tr << endl;
    cout << "\tTm: melting temperature " << Tm << endl;
    cout << "\tFlow stress: sigma_f = [" << A << " + " << B << " * (eps_p)^"
         << n << "]";
    if (C != 0)
      cout << "[1 + " << C << " * ln(eps_p_dot/" << epsdot0 << ")]";
    if (m != 0)
      cout << "[1 - ((T - " << Tr << ")/(" << Tm << " - " << Tr << "))^" << m
           << "]";
    cout << "(1 - D)\n";
  }

  if (Tr == Tm) {
    error->all(FLERR, "Error: reference temperature Tr=" + to_string(Tr) +
                          " equals melting temperature Tm=" + to_string(Tm) +
                          ".\n");
  }

  Tmr = Tm - Tr;
}

double StrengthJohnsonCook::G() { return G_; }

Matrix3d StrengthJohnsonCook::update_deviatoric_stress(
    const Eigen::Matrix3d &sigma, const Eigen::Matrix3d &D,
    double &plastic_strain_increment, const double eff_plastic_strain,
    const double epsdot, const double damage, const double T)
{
  if (damage >= 1.0)
  {
    Matrix3d sigmaFinal_dev;
    sigmaFinal_dev.setZero();
    return sigmaFinal_dev;
  }

  Matrix3d sigmaFinal_dev, sigmaTrial, sigmaTrial_dev;
  double J2, Gd, yieldStress;

  double epsdot_ratio = epsdot / epsdot0;
  epsdot_ratio        = MAX(epsdot_ratio, 1.0);

  if (eff_plastic_strain < 1.0e-10) {
    yieldStress = A;
  } else {
    yieldStress = A + B * pow(eff_plastic_strain, n);
  }
  if (C != 0) {
    yieldStress *= pow(1.0 + epsdot_ratio, C);
  }

  if (T < Tm) {
    if (m != 0 && T >= Tr) {
      yieldStress *= 1.0 - pow((T - Tr) / Tmr, m);
    }
  } else {
    yieldStress = 0;
  }

  /*
   * deviatoric rate of unrotated stress
   */
  Gd               = G_;

  if (damage > 0)
  {
    Gd *= (1 - damage);
    yieldStress *= (1 - damage);
  }

  // sigmaInitial_dev = Deviator(sigma);

  /*
   * perform a trial elastic update to the deviatoric stress
   */
  sigmaTrial = sigma + update->dt * 2.0 * Gd * D;
  sigmaTrial_dev = Deviator(sigmaTrial); // increment stress deviator using deviatoric rate

  /*
   * check yield condition
   */
  J2             = SQRT_3_OVER_2 * sigmaTrial_dev.norm();
  sigmaFinal_dev = sigmaTrial_dev;

  if (J2 < yieldStress)
  {
    /*
     * no yielding has occured.
     * final deviatoric stress is trial deviatoric stress
     */
    plastic_strain_increment = 0.0;
  }
  else
  {
    // printf("yiedl\n");
    /*
     * yielding has occured
     */

    plastic_strain_increment = (J2 - yieldStress) / (3.0 * Gd);
    /*
     * new deviatoric stress:
     * obtain by scaling the trial stress deviator
     */
    sigmaFinal_dev *= (yieldStress / J2);
  }

  return sigmaFinal_dev;
}

void StrengthJohnsonCook::write_restart(ofstream *of) {
  of->write(reinterpret_cast<const char *>(&G_), sizeof(double));
  of->write(reinterpret_cast<const char *>(&A), sizeof(double));
  of->write(reinterpret_cast<const char *>(&B), sizeof(double));
  of->write(reinterpret_cast<const char *>(&n), sizeof(double));
  of->write(reinterpret_cast<const char *>(&m), sizeof(double));
  of->write(reinterpret_cast<const char *>(&epsdot0), sizeof(double));
  of->write(reinterpret_cast<const char *>(&C), sizeof(double));
  of->write(reinterpret_cast<const char *>(&Tr), sizeof(double));
  of->write(reinterpret_cast<const char *>(&Tm), sizeof(double));
}

void StrengthJohnsonCook::read_restart(ifstream *ifr) {
  if (universe->me == 0) {
    cout << "Restart StrengthJohnsonCook" << endl;
  }
  ifr->read(reinterpret_cast<char *>(&G_), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&A), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&B), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&n), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&m), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&epsdot0), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&C), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&Tr), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&Tm), sizeof(double));
  Tmr = Tm - Tr;
}
