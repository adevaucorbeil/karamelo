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

FixStyle(contact/hertz, FixContactHertz)

#else

#ifndef MPM_FIX_CONTACT_HERTZ_H
#define MPM_FIX_CONTACT_HERTZ_H

#include "fix.h"
#include "var.h"
#include "fix_contact.h"
#include <vector>

class FixContactHertz : public FixContact
{
public:
  FixContactHertz(class MPM *, vector<string>);
  ~FixContactHertz();

  void initial_integrate();
  void init();
  void force_increment(Eigen::Vector3d &dx, Eigen::Vector3d &ftot,
                       Solid *s1, Solid *s2,
                       const int ip1, const int ip2, 
                       const double r, const double Rp1, const double Rp2);

private:
  const string USAGE = "Usage: fix(fix-ID, contact/hertz, solid1, solid2)\n";
  const int NARGS = 4;
  int solid1, solid2;
  double Estar;
};

#endif
#endif
