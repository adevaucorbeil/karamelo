/* -*- c++ -*- ----------------------------------------------------------
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


#ifdef EOS_CLASS

EOSStyle(shock,EOSShock)

#else

#ifndef MPM_EOS_SHOCK_H
#define MPM_EOS_SHOCK_H

#include "eos.h"
#include <Eigen/Eigen>

class EOSShock : public EOS {

public:
  EOSShock(class MPM *, vector<string>);
  ~EOSShock();

  double rho0();
  double K();
  double G();
  void compute_pressure(double &, double &, const double, const double, const double, const double);

protected:
  double rho0_, K_, e0, c0, S, Gamma, Tr, cv, alpha;
  string usage = "Usage: eos(eos-ID, shock, rho, K, c0, S, Gamma, cv, Tr)\n";
  int Nargs = 9;
};

#endif
#endif