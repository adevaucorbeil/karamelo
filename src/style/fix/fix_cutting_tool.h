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

FixStyle(cuttingtool, FixCuttingTool)

#else

#ifndef MPM_FIX_CUTTING_TOOL_H
#define MPM_FIX_CUTTING_TOOL_H

#include <fix.h>
#include <var.h>
#include <matrix.h>

class FixCuttingTool : public Fix {
public:
  FixCuttingTool(MPM *, vector<string>);

  void prepare();
  void reduce();

  void initial_integrate();

  void write_restart(ofstream *);
  void read_restart(ifstream *);

private:
  string usage = "Usage: fix(fix-ID, cuttingtool, group, K, x_tip, y_tip, "
                 "z_tip, vx_tip, vy_tip, vz_tip, xA, yA, xB, yB)\n";
  int Nargs = 14;
  double K;
  Var xtvalue, ytvalue, ztvalue;
  Var vtxvalue, vtyvalue, vtzvalue;
  Var xAvalue, yAvalue;
  Var xBvalue, yBvalue;
  Vector3d ftot;
};

#endif
#endif

