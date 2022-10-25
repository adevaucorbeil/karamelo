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

FixStyle(indent/minimize_penetration, FixIndentMinimizePenetration)

#else

#ifndef MPM_FIX_INDENT_MINIMIZE_PENETRATION_H
#define MPM_FIX_INDENT_MINIMIZE_PENETRATION_H

#include <fix.h>
#include <var.h>
#include <matrix.h>

class Expression;

class FixIndentMinimizePenetration : public Fix {
public:
  FixIndentMinimizePenetration(MPM *, vector<string>);

  void prepare();
  void reduce();

  void initial_integrate(Solid &solid);

  void write_restart(ofstream *);
  void read_restart(ifstream *);

private:
  string type; // sphere
  string usage = "Usage: fix(fix-ID, indent/minimize_penetration, group, "
                 "sphere, R, x_center, y_center, z_center, vx_center, "
                 "vy_center, vz_center, mu)\n";
  int Nargs = 12;
  Expression *xvalue, *yvalue, *zvalue, *vxvalue, *vyvalue, *vzvalue;
  float R;  //< Sphere radius
  float mu; //< Friction coefficient
  float A;  //< Contact area?
  Vector3d ftot;
};

#endif
#endif
