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

FixStyle(neumann_bc_mech,FixNeumannBCMech)

#else

#ifndef MPM_FIX_NEUMANN_BC_MECH_H
#define MPM_FIX_NEUMANN_BC_MECH_H

#include <fix.h>
#include <var.h>
#include <matrix.h>

class FixNeumannBCMech : public Fix {
 public:
  FixNeumannBCMech(MPM *, vector<string>);

  void prepare();
  void reduce();
  
  void initial_integrate(Solid &solid);

  void write_restart(ofstream *);
  void read_restart(ifstream *);

private:
  const map<int, string> usage = {
      {1, "Usage: fix(fix-ID, neumann_bc_mech, group, tx)\n"},
      {2, "Usage: fix(fix-ID, neumann_bc_mech, group, tx, ty)\n"},
      {3, "Usage: fix(fix-ID, neumann_bc_mech, group, tx, ty, tz)\n"}};
  const map<int, int> Nargs = {{1, 4}, {2, 5}, {3, 6}};

  Var t[3];                     //< Surface traction/compression pressure
  Vector3d ftot;
};

#endif
#endif

