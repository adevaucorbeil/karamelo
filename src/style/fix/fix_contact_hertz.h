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

#include <fix.h>
#include <matrix.h>

class FixContactHertz : public Fix {
public:
  FixContactHertz(MPM *, vector<string>);

  void initial_integrate();

  void prepare();
  void reduce();

  void write_restart(ofstream *);
  void read_restart(ifstream *);

private:
  string usage = "Usage: fix(fix-ID, contact/hertz, solid1, solid2)\n";
  int Nargs = 4;
  int solid1, solid2;
  Vector3d ftot;
};

#endif
#endif

