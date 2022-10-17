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

#ifdef FIX_CLASS

FixStyle(contact/minimize_penetration, FixContactMinPenetration)

#else

#ifndef MPM_FIX_CONTACT_MIN_PENETRATION_H
#define MPM_FIX_CONTACT_MIN_PENETRATION_H

#include "fix.h"
#include "var.h"
#include "fix_contact.h"
#include <vector>

class FixContactMinPenetration : public FixContact
{
public:
  FixContactMinPenetration(class MPM *, vector<string>);
  ~FixContactMinPenetration();

  void initial_integrate();
  void force_increment(Eigen::Vector3d &dx, Eigen::Vector3d &f, Eigen::Vector3d &ftot,
                       Solid *s1, Solid *s2,
                       const int ip1, const int ip2, 
                       const double r, const double Rp1, const double Rp2);

  void write_restart(ofstream *);
  void read_restart(ifstream *);

private:
  const string USAGE = "Usage: fix(fix-ID, contact/minimize_penetration, solid1, solid2, mu)\n";
  ;
  const int NARGS = 5;
  int solid1, solid2;
  double alpha;
  double mu; // Friction coefficient
};

#endif
#endif
