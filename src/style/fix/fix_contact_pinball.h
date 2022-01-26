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

FixStyle(contact/pinball, FixContactPinball)

#else

#ifndef MPM_FIX_CONTACT_PINBALL_H
#define MPM_FIX_CONTACT_PINBALL_H

#include <fix.h>
#include <var.h>
#include <vector>

class FixContactPinball : public Fix {
public:
  FixContactPinball(MPM *, vector<string>);

  void initial_integrate();

  void write_restart(ofstream *);
  void read_restart(ifstream *);

private:
  string usage = "Usage: fix(fix-ID, contact/pinball, solid1, solid2, K)\n";
  int Nargs = 5;
  int solid1, solid2;
  double K;
};

#endif
#endif

