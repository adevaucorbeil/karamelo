#include "modify.h"
#include "style_compute.h"
#include "style_fix.h"
#include "fix.h"
#include "compute.h"


Modify::Modify(MPM *mpm) : Pointers(mpm)
{
  n_initial_integrate = n_post_integrate = 0;
  n_pre_exchange = n_pre_neighbor = 0;
  n_pre_force = n_post_force = 0;
  n_final_integrate = n_end_of_step = n_thermo_energy = 0;
  n_initial_integrate_respa = n_post_integrate_respa = 0;
  n_pre_force_respa = n_post_force_respa = n_final_integrate_respa = 0;
  n_min_pre_exchange = n_min_pre_force = n_min_post_force = n_min_energy = 0;

  list_initial_integrate = list_post_integrate = NULL;
  list_pre_exchange = list_pre_neighbor = NULL;
  list_pre_force = list_post_force = NULL;
  list_final_integrate = list_end_of_step = NULL;
  list_thermo_energy = NULL;
  list_initial_integrate_respa = list_post_integrate_respa = NULL;
  list_pre_force_respa = list_post_force_respa = NULL;
  list_final_integrate_respa = NULL;
  list_min_pre_exchange = list_min_pre_neighbor = NULL;
  list_min_pre_force = list_min_post_force = NULL;
  list_min_energy = NULL;

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
  delete [] list_initial_integrate;
  delete [] list_post_integrate;
  delete [] list_pre_exchange;
  delete [] list_pre_neighbor;
  delete [] list_pre_force;
  delete [] list_post_force;
  delete [] list_final_integrate;
  delete [] list_end_of_step;
  delete [] list_thermo_energy;
  delete [] list_initial_integrate_respa;
  delete [] list_post_integrate_respa;
  delete [] list_pre_force_respa;
  delete [] list_post_force_respa;
  delete [] list_final_integrate_respa;
  delete [] list_min_pre_exchange;
  delete [] list_min_pre_neighbor;
  delete [] list_min_pre_force;
  delete [] list_min_post_force;
  delete [] list_min_energy;

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

  if (find_fix(args[0]) >= 0) {
    cout << "Error: reuse of fix ID" << endl;
    exit(1);
  }

    // create the Fix

  if (fix_map->find(args[1]) != fix_map->end()) {
    FixCreator fix_creator = (*fix_map)[args[1]];
    fix.push_back(fix_creator(mpm, args));
    fix.back()->init();
  }
  else {
    cout << "Unknown fix style " << args[1] << endl;
    exit(1);
  }
  
}

int Modify::find_fix(string name)
{
  for (int ifix = 0; ifix < fix.size(); ifix++) {
    cout << "fix["<< ifix <<"]->id=" << fix[ifix]->id << endl;
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
    cout << "Could not find fix ID to delete" << endl;
    exit(1);
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
    cout << "Error: reuse of compute ID" << endl;
    exit(1);
  }

    // create the Compute

  if (compute_map->find(args[1]) != compute_map->end()) {
    ComputeCreator compute_creator = (*compute_map)[args[1]];
    compute.push_back(compute_creator(mpm, args));
    compute.back()->init();
  }
  else {
    cout << "Unknown compute style " << args[1] << endl;
    exit(1);
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
    cout << "Could not find compute ID to delete" << endl;
    exit(1);
  }
  delete_compute(icompute);
}

void Modify::delete_compute(int icompute)
{
  if (compute[icompute]) delete compute[icompute];
  compute.erase(compute.begin()+icompute);
}
