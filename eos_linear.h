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

EOSStyle(linear,EOSLinear)

#else

#ifndef MPM_EOS_LINEAR_H
#define MPM_EOS_LINEAR_H

#include "eos.h"
#include <Eigen/Eigen>

class EOSLinear : public EOS {

public:
  EOSLinear(class MPM *, vector<string>);
  ~EOSLinear();

  double rho0();
  double K();
  double G();
  void compute_pressure(double &, double &, const double, const double, const double, const double);

protected:
  double rho0_, K_;
};

#endif
#endif
