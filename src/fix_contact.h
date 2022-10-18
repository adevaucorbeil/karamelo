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

#ifndef MPM_FIX_CONTACT_H
#define MPM_FIX_CONTACT_H

#include "fix.h"
#include "var.h"
#include "solid.h"
#include <vector>
#include <Eigen/Eigen>

class FixContact : public Fix
{
public:
  FixContact(class MPM *, vector<string>);
  virtual ~FixContact(){};
  void setmask();
  void init();
  void setup();

  void initial_integrate();
  virtual void force_increment(Eigen::Vector3d &dx, Eigen::Vector3d &ftot,
                               Solid *s1, Solid *s2,
                               const int ip1, const int ip2,
                               const double r, const double Rp1, const double Rp2) = 0;
  void post_particles_to_grid(){};
  void post_update_grid_state(){};
  void post_grid_to_point(){};
  void post_advance_particles(){};
  void post_velocities_to_grid(){};
  void final_integrate(){};

  void write_restart(ofstream *);
  void read_restart(ifstream *);

private:
  const int NARGS = 4;
  const string USAGE;
  int solid1, solid2;
};

#endif // FIX_CONTACT_H
