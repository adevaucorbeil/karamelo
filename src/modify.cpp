/* ----------------------------------------------------------------------
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

#include "modify.h"
#include "compute.h"
#include "error.h"
#include "fix.h"
#include "style_compute.h"
#include "style_fix.h"
#include "update.h"

using namespace FixConst;

Modify::Modify(MPM *mpm) : Pointers(mpm)
{
  end_of_step_every = NULL;

  list_timeflag = NULL;

  nfix_restart_global = 0;
  id_restart_global = style_restart_global = state_restart_global = NULL;
  nfix_restart_peratom = 0;
  id_restart_peratom = style_restart_peratom = NULL;
  index_restart_peratom = NULL;


  // fill map with fixes listed in style_fix.h

  fix_map = new FixCreatorMap();

#define FIX_CLASS
#define FixStyle(key,Class) \
  (*fix_map)[#key] = &fix_creator<Class>;
#include "style_fix.h"
#undef FixStyle
#undef FIX_CLASS

  // fill map with computes listed in style_compute.h

  compute_map = new ComputeCreatorMap();

#define COMPUTE_CLASS
#define ComputeStyle(key,Class) \
  (*compute_map)[#key] = &compute_creator<Class>;
#include "style_compute.h"
#undef ComputeStyle
#undef COMPUTE_CLASS
}

/* ---------------------------------------------------------------------- */

Modify::~Modify()
{
  delete [] end_of_step_every;
  delete [] list_timeflag;

  delete compute_map;
  delete fix_map;

  while (fix.size()) delete_fix(0);

  while (compute.size()) delete_compute(0);
}

/* ----------------------------------------------------------------------
   initialize all fixes and computes
------------------------------------------------------------------------- */

void Modify::init()
{
  for (int i = 0; i < fix.size(); i++) fix[i]->init();
  for (int i = 0; i < compute.size(); i++) compute[i]->init();

  list_init(INITIAL_INTEGRATE, list_initial_integrate);
  list_init(POST_PARTICLES_TO_GRID, list_post_particles_to_grid);
  list_init(POST_UPDATE_GRID_STATE, list_post_update_grid_state);
  list_init(POST_GRID_TO_POINT, list_post_grid_to_point);
  list_init(POST_ADVANCE_PARTICLES, list_post_advance_particles);
  list_init(POST_VELOCITIES_TO_GRID, list_post_velocities_to_grid);
  list_init(FINAL_INTEGRATE, list_final_integrate);
}

/* ----------------------------------------------------------------------
   setup for run, calls setup() of all fixes and computes
   called from scheme
------------------------------------------------------------------------- */

void Modify::setup()
{
  for (int i = 0; i < fix.size(); i++) fix[i]->setup();
  for (int i = 0; i < compute.size(); i++) compute[i]->setup();
}

/* ----------------------------------------------------------------------
   create a new region
------------------------------------------------------------------------- */

void Modify::add_fix(vector<string> args){
  cout << "In add_fix" << endl;

  int ifix = find_fix(args[0]);

  if (ifix >= 0) { // Fix already exists: modify it
    delete fix[ifix];
    FixCreator fix_creator = (*fix_map)[args[1]];
    fix[ifix] = fix_creator(mpm, args);
    fix[ifix]->init();
    if (fix[ifix]->requires_ghost_particles)
      ghost_particles_any = true;
  }
  else if (fix_map->find(args[1]) != fix_map->end()) {
    ifix = fix.size();
    FixCreator fix_creator = (*fix_map)[args[1]];
    fix.push_back(fix_creator(mpm, args));
    fix.back()->init();

    if (update->method_type.compare("tlmpm") == 0 &&
        fix.back()->requires_ghost_particles) {
      ghost_particles_any = true;
      cout << "Fix " << args[1] << " requires ghost particles.\n";
    }
  }
  else {
    error->all(FLERR, "Unknown fix style " + args[1] + ".\n");
  }

  fix[ifix]->setmask();

}

int Modify::find_fix(string name)
{
  for (int ifix = 0; ifix < fix.size(); ifix++) {
    //cout << "fix["<< ifix <<"]->id=" << fix[ifix]->id << endl;
    if (name.compare(fix[ifix]->id) == 0) return ifix;
  }
  return -1;
}

/* ----------------------------------------------------------------------
   one instance per fix style in style_fix.h
------------------------------------------------------------------------- */

template <typename T>
Fix *Modify::fix_creator(MPM *mpm, vector<string> args)
{
  return new T(mpm, args);
}

/* ----------------------------------------------------------------------
   delete a Fix from list of Fixes
------------------------------------------------------------------------- */

void Modify::delete_fix(string id)
{
  int ifix = find_fix(id);
  if (ifix < 0) {
    error->all(FLERR, "Could not find fix ID to delete.\n");
  }
  delete_fix(ifix);
}

void Modify::delete_fix(int ifix)
{
  if (fix[ifix]) delete fix[ifix];
  fix.erase(fix.begin()+ifix);
}

/* ----------------------------------------------------------------------
   create a new compute
------------------------------------------------------------------------- */

void Modify::add_compute(vector<string> args){
  cout << "In add_compute" << endl;

  if (find_compute(args[0]) >= 0) {
    error->all(FLERR, "Error: reuse of compute ID.\n");
  }

    // create the Compute

  if (compute_map->find(args[1]) != compute_map->end()) {
    ComputeCreator compute_creator = (*compute_map)[args[1]];
    compute.push_back(compute_creator(mpm, args));
    compute.back()->init();
  }
  else {
    error->all(FLERR, "Unknown compute style " + args[1] + ".\n");
  }
  
}

int Modify::find_compute(string name)
{
  for (int icompute = 0; icompute < compute.size(); icompute++) {
    cout << "compute["<< icompute <<"]->id=" << compute[icompute]->id << endl;
    if (name.compare(compute[icompute]->id) == 0) return icompute;
  }
  return -1;
}

/* ----------------------------------------------------------------------
   one instance per compute style in style_compute.h
------------------------------------------------------------------------- */

template <typename T>
Compute *Modify::compute_creator(MPM *mpm, vector<string> args)
{
  return new T(mpm, args);
}

/* ----------------------------------------------------------------------
   delete a Compute from list of Computes
------------------------------------------------------------------------- */

void Modify::delete_compute(string id)
{
  int icompute = find_compute(id);
  if (icompute < 0) {
    error->all(FLERR, "Could not find compute ID to delete.\n");
  }
  delete_compute(icompute);
}

void Modify::delete_compute(int icompute)
{
  if (compute[icompute]) delete compute[icompute];
  compute.erase(compute.begin()+icompute);
}

void Modify::list_init(int mask, vector<int> &list) {

  list.clear();
  for (int i = 0; i < fix.size(); i++) if (fix[i]->mask & mask) list.push_back(i);
}

/* ----------------------------------------------------------------------
   1st half of integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::initial_integrate()
{
  for (int i = 0; i < list_initial_integrate.size(); i++)
    fix[list_initial_integrate[i]]->initial_integrate();
}

/* ----------------------------------------------------------------------
   after post_particles_to_grid(), only for relevant fixes
------------------------------------------------------------------------- */

void Modify::post_particles_to_grid()
{
  for (int i = 0; i < list_post_particles_to_grid.size(); i++)
    fix[list_post_particles_to_grid[i]]->post_particles_to_grid();
}

/* ----------------------------------------------------------------------
   after update_grid_state(), only for relevant fixes
------------------------------------------------------------------------- */

void Modify::post_update_grid_state(){
  for (int i = 0; i < list_post_update_grid_state.size(); i++)
    fix[list_post_update_grid_state[i]]->post_update_grid_state();
}

/* ----------------------------------------------------------------------
   after grid_to_point(), only for relevant fixes
------------------------------------------------------------------------- */

void Modify::post_grid_to_point(){
  for (int i = 0; i < list_post_grid_to_point.size(); i++)
    fix[list_post_grid_to_point[i]]->post_grid_to_point();
}

/* ----------------------------------------------------------------------
   after advance_particles(), only for relevant fixes
------------------------------------------------------------------------- */

void Modify::post_advance_particles(){
  for (int i = 0; i < list_post_advance_particles.size(); i++)
    fix[list_post_advance_particles[i]]->post_advance_particles();
}

/* ----------------------------------------------------------------------
   after velocities_to_grid(), only for relevant fixes
------------------------------------------------------------------------- */

void Modify::post_velocities_to_grid(){
  for (int i = 0; i < list_post_velocities_to_grid.size(); i++)
    fix[list_post_velocities_to_grid[i]]->post_velocities_to_grid();
}

/* ----------------------------------------------------------------------
   final_integrate, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::final_integrate(){
  for (int i = 0; i < list_final_integrate.size(); i++)
    fix[list_final_integrate[i]]->final_integrate();
}


void Modify::run_computes(){
  for (int i = 0; i < compute.size(); i++)
    compute[i]->compute_value();
}
