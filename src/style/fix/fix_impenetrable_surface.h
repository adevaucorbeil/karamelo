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

FixStyle(impenetrablesurface, FixImpenetrableSurface)

#else

#ifndef MPM_FIX_IMPENETRABLE_SURFACE_H
#define MPM_FIX_IMPENETRABLE_SURFACE_H

#include <fix.h>
#include <var.h>
#include <matrix.h>

class FixImpenetrableSurface : public Fix {
public:
  FixImpenetrableSurface(MPM *, vector<string>);

  void prepare();
  void reduce();

  void initial_integrate();

  void write_restart(ofstream *);
  void read_restart(ifstream *);

private:
  string usage = "Usage: fix(fix-ID, impenetrablesurface, group, K, xs, ys, "
                 "zs, nx, ny, nz)\n";
  int Nargs = 10;

  Var xs_x, xs_y, xs_z;                  //< Position of a point on the surface
  Var nx, ny, nz;                        //< Normal to the plane
  double K;                                    //< Contact stiffness
  Vector3d ftot;
};

#endif
#endif

