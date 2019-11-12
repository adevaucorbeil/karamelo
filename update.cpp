#include <iostream>
#include "update.h"
#include "scheme.h"
#include "method.h"
#include "input.h"
#include "var.h"
#include "style_scheme.h"
#include "style_method.h"
#include <vector>

using namespace std;

Update::Update(MPM *mpm) : Pointers(mpm)
{
  cout << "In Update::Update()" << endl;

  atime = 0;
  atimestep = 0;
  ntimestep = 0; // set the current timestep to 0
  firststep = laststep = 0;
  beginstep = endstep = 0;
  first_update = 0;
  dt = 1e-16;
  dt_constant = false;
  dt_factor = 0.9;

  // Default scheme is MUSL:
  vector<string> scheme_args;
  scheme_args.push_back("musl");
  create_scheme(scheme_args);

  method = NULL;
  //vector<string> method_args;
  //method_args.push_back("tlmpm");
  //create_method(method_args);
}

Update::~Update()
{
  delete scheme;
  delete method;
}

void Update::set_dt_factor(vector<string> args){
  if (args.size()!=1) {
    cout << "Illegal dt_factor command: not enough arguments or too many arguments" << endl;
    exit(1);
  }
  dt_factor = input->parsev(args[0]);
}

void Update::set_dt(vector<string> args){
  if (args.size()!=1) {
    cout << "Illegal set_dt command: not enough arguments or too many arguments" << endl;
    exit(1);
  }
  dt = input->parsev(args[0]);
  dt_constant = true;
  (*input->vars)["dt"] = Var("dt", dt);
}

void Update::create_scheme(vector<string> args){
  if (args.size() < 1) {
    cout << "Illegal run_scheme command: not enough arguments" << endl;
    exit(1);
  }

  new_scheme(args);

  scheme_style = args[0];
}

void Update::new_scheme(vector<string> args){

  string style = args[0];

  if (0) return;

#define SCHEME_CLASS
#define SchemeStyle(key,Class) \
  else if (style.compare(#key) == 0) scheme = new Class(mpm,args);
#include "style_scheme.h"
#undef SchemeStyle
#undef SCHEME_CLASS

  else {
    cout << "Illegal scheme style" << endl;
    exit(1);
  }
}

void Update::create_method(vector<string> args){
  if (args.size() < 3) {
    cout << "Illegal method command: not enough arguments" << endl;
    exit(1);
  }

  new_method(args);

  method_type = args[0];
  method_shape_function = args[2];
}

void Update::new_method(vector<string> args){

  string style = args[0];

  if (0) return;

#define METHOD_CLASS
#define MethodStyle(key,Class) \
  else if (style.compare(#key) == 0) method = new Class(mpm,args);
#include "style_method.h"
#undef MethodStyle
#undef METHOD_CLASS

  else {
    cout << "Illegal method style" << endl;
    exit(1);
  }

  method->setup(args);
}

/* ----------------------------------------------------------------------
   update elapsed simulation time
   called at end of runs or when timestep size changes
------------------------------------------------------------------------- */

void Update::update_time()
{
  atime += (ntimestep-atimestep) * dt;
  atimestep = ntimestep;
  (*input->vars)["time"] = Var("time", atime);
}


int Update::update_timestep()
{
  update->ntimestep++;
  (*input->vars)["timestep"] = Var("timestep", ntimestep);
  return update->ntimestep;
}
