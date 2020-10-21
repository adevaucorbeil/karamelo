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

#include "fix.h"
#include "var.h"
#include <vector>

class FixMeldTool : public Fix {
public:
  FixMeldTool(class MPM *, vector<string>);
  ~FixMeldTool();
  void setmask();
  void init();
  void setup();

  void initial_integrate();
  void post_particles_to_grid(){};
  void post_update_grid_state(){};
  void post_grid_to_point(){};
  void post_advance_particles(){};
  void post_velocities_to_grid(){};
  void final_integrate(){};

private:
  string usage = "Usage: fix(fix-ID, meldtool, group, K, dim, half-width, c1, c2, theta, lo, hi, Rmax)\n";
  int Nargs = 12;
  int dim, axis0, axis1, Rmax, RmaxSq;
  double K, w, lo, hi;
  int c1pos = 6;
  int c2pos = 7;
  int thetapos = 8; // Positions the arguments

  enum Axis { X, Y, Z, };
};

#endif
#endif

