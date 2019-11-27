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

#ifndef MPM_EOS_H
#define MPM_EOS_H

#include "pointers.h"
#include <vector>
#include <Eigen/Eigen>

class EOS : protected Pointers {
 public:
  string id;

  EOS(class MPM *, vector<string>);
  virtual ~EOS();
  virtual void init();
  void options(vector<string> *, vector<string>::iterator);

  // implemented by each EOS
  //virtual compute_pressure()
  virtual double rho0() = 0;
  virtual double K() = 0;
  virtual void compute_pressure(double &, double&, const double, const double, const double, const double) = 0;
  //protected:
};

#endif
