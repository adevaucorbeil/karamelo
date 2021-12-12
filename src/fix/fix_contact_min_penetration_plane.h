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

FixStyle(contact/minimize_penetration/plane, FixContactMinPenetrationPlane)

#else

#ifndef MPM_FIX_CONTACT_MIN_PENETRATION__PLANE_H
#define MPM_FIX_CONTACT_MIN_PENETRATION_PLANE_H

#include <fix.h>
#include <var.h>
#include <Eigen/Eigen>
#include <vector>

class FixContactMinPenetrationPlane : public Fix {
public:
  FixContactMinPenetrationPlane(class MPM *, vector<string>);
  ~FixContactMinPenetrationPlane();
  void setmask();
  void init();
  void setup();

  void initial_integrate();
  void post_particles_to_grid(){};
  void post_update_grid_state(){};
  void post_grid_to_point(){};
  void post_advance_particles() {};
  void post_velocities_to_grid(){};
  void final_integrate(){};

  void write_restart(ofstream *);
  void read_restart(ifstream *);

private:
  string usage = "Usage: fix(fix-ID, contact/minimize_penetration/plane, solid, x, y, z, nx, ny, nz, mu)\n";
  int Nargs = 10;
  int solid;
  double D; // Coordinates of the plane.
  Eigen::Vector3d xq; // Point in the plane.
  Eigen::Vector3d n;  // Normal to the plane.
  double mu;          // Friction coefficient.
};

#endif
#endif

