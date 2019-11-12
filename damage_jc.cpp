#include <iostream>
#include "damage_jc.h"
#include "input.h"
#include "domain.h"
#include "update.h"
#include "mpm_math.h"
#include <Eigen/Eigen>
#include "var.h"

using namespace std;
using namespace Eigen;
using namespace MPM_Math;

#define SQRT_3_OVER_2 1.224744871 // sqrt(3.0/2.0)

DamageJohnsonCook::DamageJohnsonCook(MPM *mpm, vector<string> args) : Damage(mpm, args)
{
  cout << "Initiate DamageJohnsonCook" << endl;

  if (args.size() < Nargs) {
    cout << "Error: too few arguments for the strength command" << endl;
    cout << usage;
    exit(1);
  }

  //options(&args, args.begin()+3);
  d1 = input->parsev(args[2]);
  d2 = input->parsev(args[3]);
  d3 = input->parsev(args[4]);
  d4 = input->parsev(args[5]);
  d5 = input->parsev(args[6]);
  epsdot0 = input->parsev(args[6]);
  Tr      = input->parsev(args[7]);
  Tm      = input->parsev(args[8]);

  cout << "Johnson Cook material damage model:\n";
  cout << "\tparameter d1:" << d1 << endl;
  cout << "\tparameter d2:" << d2 << endl;
  cout << "\tparameter d3:" << d3 << endl;
  cout << "\tparameter d4:" << d4 << endl;
  cout << "\tparameter d5:" << d5 << endl;
  cout << "\tepsdot0: reference strain rate " << epsdot0 << endl;
  cout << "\tTr: reference temperature " << Tr << endl;
  cout << "\tTm: melting temperature " << Tm << endl;
  
  if (Tr == Tm) {
    cout << "Error: reference temperature Tr=" << Tr << " equals melting temperature Tm=" << Tm << endl;
    exit(1);
  }
  Tmr = Tm - Tr;
}

void DamageJohnsonCook::compute_damage(double &damage_init, double &damage, const double pH, const Eigen::Matrix3d Sdev, const double epsdot, const double plastic_strain_increment, const double T)
{
  double vm = SQRT_3_OVER_2 * Sdev.norm(); // von-Mises equivalent stress
  if (vm < 0.0) {
    cout << "this is sdev " << endl << Sdev << endl;
    printf("vm=%f < 0.0, surely must be an error\n", vm);
    exit(1);
  }

  // determine stress triaxiality
  double triax = 0.0;
  if (pH != 0.0 && vm != 0.0) {
    triax = -pH / (vm + 0.01 * fabs(pH)); // have softening in denominator to avoid divison by zero
  }
  if (triax > 3.0) {                                                                                                                                                                                
    triax = 3.0;
  }
  // Johnson-Cook failure strain, dependence on stress triaxiality
  double jc_failure_strain = d1 + d2 * exp(d3 * triax);

  // include strain rate dependency if parameter d4 is defined and current plastic strain rate exceeds reference strain rate
  if (d4 > 0.0) { //
    if (epsdot > epsdot0) {
      double epdot_ratio = epsdot / epsdot0;
      jc_failure_strain *= (1.0 + d4 * log(epdot_ratio));
    }
  }

  if (d5 > 0.0 && T >= Tr) jc_failure_strain *= 1 + d5*(T-Tr)/Tmr;

  damage_init += plastic_strain_increment/jc_failure_strain;

  if (damage_init >= 1.0) damage = MIN((damage_init-1.0)*10, 1.0);
}
