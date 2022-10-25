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

#include <damage_jc.h>
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

#define SQRT_3_OVER_2 1.224744871 // sqrt(3.0/2.0)

DamageJohnsonCook::DamageJohnsonCook(MPM *mpm, vector<string> args)
    : Damage(mpm, args)
{
  if (args.size() < 3) {
    error->all(FLERR, "Error: too few arguments for the strength command.\n");
  }

  if (args[2].compare("restart") ==
      0) { // If the keyword restart, we are expecting to have read_restart()
           // launched right after.
    d1 = d2 = d3 = d4 = d5 = epsdot0 = Tr = Tm = Tmr = 0;
    return;
  }

  if (universe->me == 0)
    cout << "Initiate DamageJohnsonCook" << endl;

  if (args.size() < Nargs) {
    error->all(FLERR,
               "Error: too few arguments for the strength command\n" + usage);
  }

  // options(&args, args.begin()+3);
  d1      = input->parsev(args[2]);
  d2      = input->parsev(args[3]);
  d3      = input->parsev(args[4]);
  d4      = input->parsev(args[5]);
  d5      = input->parsev(args[6]);
  epsdot0 = input->parsev(args[7]);
  Tr      = input->parsev(args[8]);
  Tm      = input->parsev(args[9]);

  if (universe->me == 0) {
    cout << "Johnson Cook material damage model:\n";
    cout << "\tparameter d1:" << d1 << endl;
    cout << "\tparameter d2:" << d2 << endl;
    cout << "\tparameter d3:" << d3 << endl;
    cout << "\tparameter d4:" << d4 << endl;
    cout << "\tparameter d5:" << d5 << endl;
    cout << "\tepsdot0: reference strain rate " << epsdot0 << endl;
    cout << "\tTr: reference temperature " << Tr << endl;
    cout << "\tTm: melting temperature " << Tm << endl;
    cout << "\tFailure strain equation: eps_f = [" << d1 << " + " << d2
         << " * exp(" << d3 << "*sigma*)]";
    if (d4 != 0)
      cout << "[ 1 + " << d4 << " * ln(epsdot/" << epsdot0 << ")]";
    if (d5 != 0)
      cout << "[1 + " << d5 << " * (T - " << Tr << ")/(" << Tm << " - " << Tr
           << ")]";
    cout << endl;
  }

  if (Tr == Tm) {
    error->all(FLERR, "Error: reference temperature Tr=" + to_string(Tr) +
                          " equals melting temperature Tm=" + to_string(Tm) +
                          ".\n");
  }
  Tmr = Tm - Tr;
}

//void DamageJohnsonCook::compute_damage(float &damage_init, float &damage,
//                                       const float pH,
//                                       const Matrix3d Sdev,
//                                       const float epsdot,
//                                       const float plastic_strain_increment,
//                                       const float T) {


void
DamageJohnsonCook::compute_damage(Solid &solid,
                                  Kokkos::View<float*> &pH,
                                  Kokkos::View<Matrix3d*> &sigma_dev,
                                  Kokkos::View<float*> &plastic_strain_increment) const
{
  float d1 = this->d1;
  float d2 = this->d2;
  float d3 = this->d3;
  float d4 = this->d4;
  float d5 = this->d5;
  float epsdot0 = this->epsdot0;
  float Tr = this->Tr;
  float Tmr = this->Tmr;

  float cp = solid.mat->cp;

  Kokkos::View<float*> sdamage = solid.damage;
  Kokkos::View<float*> sdamage_init = solid.damage_init;
  Kokkos::View<float*> sT = solid.T;
  Kokkos::View<float*> seff_plastic_strain_rate = solid.eff_plastic_strain_rate;

  Kokkos::parallel_for("DamageJohnsonCook::compute_damage", solid.np_local,
  KOKKOS_LAMBDA (const int &ip)
  {
    if (plastic_strain_increment[ip] == 0 && sdamage[ip] >= 1.0)
      return;

    float vm = SQRT_3_OVER_2*sigma_dev[ip].norm(); // von-Mises equivalent stress

    //if (vm < 0.0)
    //{
    //  // cout << "this is sdev " << endl << Sdev << endl;
    //  char str[128];
    //  sprintf(str,"vm=%f < 0.0, surely must be an error\n", vm);
    //  error->all(FLERR, str);
    //}

    // determine stress triaxiality
    float triax = 0.0;
    if (pH[ip] != 0.0 && vm != 0.0)
    {
      triax = -pH[ip] / (vm + 0.001 * Kokkos::Experimental::fabs(pH[ip])); // have softening in denominator to
                                             // avoid divison by zero
    }

    if (triax <= -3) {
      sdamage_init[ip] = 0;
      return ;
    }

    // Johnson-Cook failure strain, dependence on stress triaxiality
    float jc_failure_strain = d1 + d2 * Kokkos::Experimental::exp(d3 * triax);

    // include strain rate dependency if parameter d4 is defined and current
    // plastic strain rate exceeds reference strain rate
    if (d4 > 0.0) {
      if (seff_plastic_strain_rate[ip] > epsdot0) {
        float epdot_ratio = seff_plastic_strain_rate[ip]/epsdot0;
        jc_failure_strain *= (1.0 + d4 * Kokkos::Experimental::log(epdot_ratio));
      }
    }

    if (d5 > 0.0 && cp && sT[ip] >= Tr)
      jc_failure_strain *= 1 + d5 * (sT[ip] - Tr) / Tmr;

    sdamage_init[ip] += plastic_strain_increment[ip]/jc_failure_strain;

    if (sdamage_init[ip] >= 1.0)
      sdamage[ip] = MIN((sdamage_init[ip] - 1.0) * 10, 1.0);
  });
}

void DamageJohnsonCook::write_restart(ofstream *of) {
  of->write(reinterpret_cast<const char *>(&d1), sizeof(float));
  of->write(reinterpret_cast<const char *>(&d2), sizeof(float));
  of->write(reinterpret_cast<const char *>(&d3), sizeof(float));
  of->write(reinterpret_cast<const char *>(&d4), sizeof(float));
  of->write(reinterpret_cast<const char *>(&d5), sizeof(float));
  of->write(reinterpret_cast<const char *>(&epsdot0), sizeof(float));
  of->write(reinterpret_cast<const char *>(&Tr), sizeof(float));
  of->write(reinterpret_cast<const char *>(&Tm), sizeof(float));
}

void DamageJohnsonCook::read_restart(ifstream *ifr) {
  // cout << "Restart DamageJohnsonCook" << endl;
  ifr->read(reinterpret_cast<char *>(&d1), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&d2), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&d3), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&d4), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&d5), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&epsdot0), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&Tr), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&Tm), sizeof(float));
  Tmr = Tm - Tr;
}