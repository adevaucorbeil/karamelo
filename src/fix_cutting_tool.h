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

#include "fix.h"
#include "var.h"
#include <vector>

class FixCuttingTool : public Fix {
public:
  FixCuttingTool(class MPM *, vector<string>);
  ~FixCuttingTool();
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
  string usage = "Usage: fix(fix-ID, cuttingtool, group, K, x_tip, y_tip, "
                 "z_tip, vx_tip, vy_tip, vz_tip, xA, yA, xB, yB)\n";
  int Nargs = 14;
  int Kpos, xtpos, ytpos, ztpos, vtxpos, vtypos, vtzpos, xApos, yApos, xBpos,
      yBpos; // Positions the arguments
};

#endif
#endif

