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

#ifndef MPM_MODIFY_H
#define MPM_MODIFY_H

#include <pointers.h>
#include <map>
#include <string>
#include <vector>

class Solid;
class Grid;

using namespace std;

class Modify : protected Pointers
{
public:
  int restart_pbc_any;      // 1 if any fix sets restart_pbc
  int nfix_restart_global;  // stored fix global info from restart file
  int nfix_restart_peratom; // stored fix peratom info from restart file

  vector<class Fix *> fix;  // list of fixes

  vector<class Compute *> compute; // list of computes

  Modify(class MPM *);
  virtual ~Modify();
  virtual void init();

  void add_fix(vector<string>);
  void delete_fix(string);
  void delete_fix(int);
  int find_fix(string);
  // int check_package(const char *);

  void add_compute(vector<string>);
  void delete_compute(string);
  void delete_compute(int);
  int find_compute(string);

  void initial_integrate(Solid &solid);
  void post_particles_to_grid(Grid &grid);
  void post_update_grid_state(Grid &grid);
  void post_grid_to_point(Solid &solid);
  void post_advance_particles(Solid &solid);
  void post_velocities_to_grid(Grid &grid);
  void final_integrate(Solid &solid);
  void prepare();
  void reduce();

  void run_computes(Solid &solid);

  void write_restart(ofstream*);      ///< Write fixes and computes to restart file
  void read_restart(ifstream*);       ///< Read fixes and computes from restart file


protected:
  // lists of fixes to apply at different stages of timestep
  vector<int> list_initial_integrate;
  vector<int> list_post_particles_to_grid;
  vector<int> list_post_update_grid_state;
  vector<int> list_post_grid_to_point;
  vector<int> list_post_advance_particles;
  vector<int> list_post_velocities_to_grid;
  vector<int> list_final_integrate;

  int *end_of_step_every;

  int n_timeflag; // list of computes that store time invocation
  int *list_timeflag;

  char **id_restart_global;    // stored fix global info
  char **style_restart_global; // from read-in restart file
  char **state_restart_global;

  char **id_restart_peratom;    // stored fix peratom info
  char **style_restart_peratom; // from read-in restart file
  int *index_restart_peratom;

  int index_permanent; // fix/compute index returned to library call

  void list_init(int, vector<int> &);

private:
  typedef Compute *(*ComputeCreator)(MPM *, vector<string>);
  typedef map<string, ComputeCreator> ComputeCreatorMap;
  ComputeCreatorMap *compute_map;

  typedef Fix *(*FixCreator)(MPM *, vector<string>);
  typedef map<string, FixCreator> FixCreatorMap;
  FixCreatorMap *fix_map;

  template <typename T> static Compute *compute_creator(MPM *, vector<string>);
  template <typename T> static Fix *fix_creator(MPM *, vector<string>);
};

#endif
