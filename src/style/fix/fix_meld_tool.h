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

FixStyle(meldtool, FixMeldTool)

#else

#ifndef MPM_FIX_MELD_TOOL_H
#define MPM_FIX_MELD_TOOL_H

#include <fix.h>
#include <var.h>
#include <matrix.h>

class Expression;

class FixMeldTool : public Fix {
public:
  FixMeldTool(MPM *, vector<string>);

  void prepare();
  void reduce();

  void initial_integrate(Solid &solid);

  void write_restart(ofstream *);
  void read_restart(ifstream *);

private:
  string usage = "Usage: fix(fix-ID, meldtool, group, K, dim, half-width, c1, c2, theta, lo, hi, Rmax)\n";
  int Nargs = 12;
  int dim, axis0, axis1;
  float K, w, lo, hi, Rmax, RmaxSq;

  Expression *c1, *c2, *theta;         //< Position and angle of the tool
  Vector3d ftot;

  enum Axis { X, Y, Z, };
};

#endif
#endif

