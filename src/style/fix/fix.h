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

#ifndef MPM_FIX_H
#define MPM_FIX_H

#include <pointers.h>
#include <vector>

/*! Parent class of all the different kinds of fixes that can be used.
 *
 * Stores the fix id, the information about the group involved,  as well
 * as the mask that select the nodes or particles  within this group.
 */

class Fix : protected Pointers {
 public:
  string id, style;
  int igroup, groupbit;
  int mask;

  bool requires_ghost_particles = false;

  Fix(MPM *mpm, const vector<string> &args, int mask);

  virtual void init() {}
  virtual void setup() {}

  virtual void prepare() {}
  virtual void reduce() {}

  virtual void initial_integrate() {}
  virtual void post_particles_to_grid() {}
  virtual void post_update_grid_state() {}
  virtual void post_grid_to_point() {}
  virtual void post_advance_particles() {}
  virtual void post_velocities_to_grid() {}
  virtual void final_integrate() {}

  virtual void write_restart(ofstream *of) = 0;
  virtual void read_restart(ifstream *ifr) = 0;
};

namespace FixConst {
  static const int INITIAL_INTEGRATE =       1<<0;
  static const int POST_PARTICLES_TO_GRID =  1<<1;
  static const int POST_UPDATE_GRID_STATE =  1<<2;
  static const int POST_GRID_TO_POINT =      1<<3;
  static const int POST_ADVANCE_PARTICLES =  1<<4;
  static const int POST_VELOCITIES_TO_GRID = 1<<5;
  static const int FINAL_INTEGRATE =         1<<6;
}
#endif


/*! \defgroup fix fix

\section Syntax Syntax
\code
fix(fix-ID, fix_type, group-ID, fix_specific_arguments)
\endcode

<ul>
<li>fix-ID: name of the fix to be created.</li>
<li>fix-type: velocity_nodes, force_nodes, kinetic_energy, ...</li>
<li>group-ID: name of the group onto which the fix will be applied. 
If 'all' is used, all particles or nodes will be selected.</li>
<li>fix_specific_arguments: list of arguments specific to the fix type used. </li>
</ul>

\section Examples Examples
\code
group(sym2n, nodes, region, sym2, solid, solid1)
fix(BC_sym2, velocity_nodes, sym2n, 0, 0, 0)
\endcode
Defines a fix 'BC_sym2' which sets the velocity of the nodes of the background grid linked to solid
'solid1' and lying in region 'sym2' at the beginning of the simulation to 0.

\code
fix(Ek, kinetic_energy, all)
\endcode
Define a fix called 'Ek' which calculates the cumulative kinetic energy of all the particles in the simulation box.

\section Description Description

This command defines a fix. Fixes are used to change the behaviour of nodes or particles, 
or to calculate some properties. Using fixes, users can implement many operations that can 
alter the system such as: changing particles or node attributes (position, velocity, forces, etc.), 
implement boundary conditions, reading/writing data, or even saving information about particles for 
future use (previous positions, for instance).

This command also generates three variables: fix-ID_s, fix-ID_x, fix-ID_y and fix-ID_z. They correspond to one scalar and the x, y, and z components of a vector. Each fix is free to allocate the desired quantities to these variables.

\section Fix_Subclasses Fixes types

To access the different fix types supported, and how to use the fix() command with them, refer to the corresponding class.
*/
