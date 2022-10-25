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

FixStyle(body_force,FixBodyForce)

#else

#ifndef MPM_FIX_BODY_FORCE_H
#define MPM_FIX_BODY_FORCE_H

#include <fix.h>
#include <matrix.h>

class Expression;

/*! \ingroup fix fixbodyforce fix_body_force

\section Syntax Syntax
\code
fix(fix-ID, body_force, group-ID, fx, fy, fz)
\endcode

<ul>
<li>fix-ID: name of the fix to be created.</li>
<li>fix-type: velocity_nodes, force_nodes, kinetic_energy, ...</li>
<li>group-ID: name of the group onto which the fix will be applied. 
If 'all' is used, all particles will be selected.</li>
<li>fx, fy and fz: x, y and z components of the body force to be applied, respectively. </li>
</ul>

\section Examples Examples
\code
fix(fbody, body_force, all, 0, -9.81, 0)
\endcode
Applies a body force \f$\boldsymbol{f}_b = [0, -9.81, 0]\f$ onto all particles in the simulation box.

\section Description Description

This command defines a body force \f$\boldsymbol{f}_b\f$ that is applied to all the particles in the group 'group-ID'non-linear equation of state:
\f{equation}{
\boldsymbol{f}_b = \begin{bmatrix}
fx\\
fy\\
fz
\end{bmatrix}
\f}

\section Class Class description
*/
class FixBodyForce : public Fix {
 public:
  FixBodyForce(MPM *, vector<string>);
  
  void prepare();
  void reduce();

  void post_particles_to_grid(Grid &grid);

  void write_restart(ofstream *);
  void read_restart(ifstream *);

  
private:
  const map<int, string> usage = {
      {1, "Usage: fix(fix-ID, body_force, group, fbx)\n"},
      {2, "Usage: fix(fix-ID, body_force, group, fbx, fby)\n"},
      {3, "Usage: fix(fix-ID, body_force, group, fbx, fby, fbz)\n"}};
  const map<int, int> Nargs = {{1, 4}, {2, 5}, {3, 6}};

  bool xset, yset, zset;                             //< Does the fix set the x, y, and z velocities of the group?
  Vector3d ftot;

  Expression *fb[3];
};

#endif
#endif
