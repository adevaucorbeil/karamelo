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

#include <iostream>
#include "update.h"
#include "scheme.h"
#include "method.h"
#include "input.h"
#include "var.h"
#include "style_scheme.h"
#include "style_method.h"
#include <vector>
#include "error.h"

using namespace std;

Update::Update(MPM *mpm) : Pointers(mpm)
{
  // cout << "In Update::Update()" << endl;

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
  shape_function = ShapeFunctions::LINEAR;
  sub_method_type = SubMethodType::FLIP;
  alpha = 0.99;
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
    error->all(FLERR, "Illegal dt_factor command: not enough arguments or too many arguments.\n");
  }
  dt_factor = input->parsev(args[0]);
}

void Update::set_dt(vector<string> args){
  if (args.size()!=1) {
    error->all(FLERR, "Illegal dt_factor command: not enough arguments or too many arguments.\n");
  }
  dt = input->parsev(args[0]);
  dt_constant = true;
  (*input->vars)["dt"] = Var("dt", dt);
}


/*! This function is the C++ equivalent to the scheme() user function.\n
 * Syntax: scheme(type)\n
 * It points the pointer Update::scheme to the desired Scheme type selected from style_scheme.h
 */
void Update::create_scheme(vector<string> args){
  if (args.size() < 1) {
    error->all(FLERR, "Illegal scheme command: not enough arguments.\n");
  }

  scheme_style = args[0];

  if (0) return;

#define SCHEME_CLASS
#define SchemeStyle(key,Class) \
  else if (scheme_style.compare(#key) == 0) scheme = new Class(mpm,args);
#include "style_scheme.h"
#undef SchemeStyle
#undef SCHEME_CLASS

  else {
    error->all(FLERR, "Illegal scheme style.\n");
  }
}

/*! This function is the C++ equivalent to the method() user function.\n
 * Syntax: method(type, type specific arguments)\n
 * It points the pointer Update::method to the desired Method type selected from style_method.h
 */
void Update::create_method(vector<string> args){
  if (args.size() < 3) {
    error->all(FLERR, "Illegal method command: not enough arguments.\n");
  }

  int n = 0;
  method_type = args[0];
  

  if (0) return;

#define METHOD_CLASS
#define MethodStyle(key,Class) \
  else if (method_type.compare(#key) == 0) method = new Class(mpm,args);
#include "style_method.h"
#undef MethodStyle
#undef METHOD_CLASS

  else {
    error->all(FLERR, "Illegal method style.\n");
  }

  bool isFLIP = false;
  n++;
  // Method used: PIC, FLIP or APIC:
  if (map_sub_method_type.count(args[n]) > 0) {
    sub_method_type = map_sub_method_type.at(args[n]);
    if (sub_method_type == SubMethodType::PIC)
      alpha = 0;
    else if (sub_method_type == SubMethodType::FLIP) {
      isFLIP = true;
      if (args.size() < 4) {
	error->all(FLERR, "Illegal modify_method command: not enough arguments.\n");
      }
    }
  } else {
    error->all(FLERR, "Error: method type " + args[n] + " not understood. Expect: PIC, FLIP or APIC\n");
  }

  n++;
  
  // Check if the shape function given in the inputfile is recognized.
  if (args.size() > 1 + isFLIP) {
    if (map_shape_functions.count(args[n]) > 0) {
      shape_function = map_shape_functions.at(args[n]);
      cout << "Setting up " << args[n] << " basis functions\n";
      n++;
    } else {
      error->all(FLERR, "Illegal method_method argument: form function of type \033[1;31m" + args[2] + "\033[0m is unknown. Available options are:  \033[1;32mlinear\033[0m, \033[1;32mcubic-spline\033[0m, \033[1;32mquadratic-spline\033[0m, \033[1;32mBernstein-quadratic\033[0m.\n");
    }
  }

  if (args.size() > n + isFLIP) {
    error->all(FLERR, "Illegal modify_method command: too many arguments: " + to_string(n + isFLIP) + " expected, " + to_string(args.size()) + " received.\n");
  }

  if (isFLIP) alpha = input->parsev(args[n]);

  n++;

  method->setup(vector<string>(args.begin() + n, args.end()));
}

/*! Update elapsed simulation time.
 *  Called at end of runs or when timestep size changes.
 */
void Update::update_time()
{
  atime += (ntimestep-atimestep) * dt;
  atimestep = ntimestep;
  (*input->vars)["time"] = Var("time", atime);
}


/*! Update simulation timestep.
 */
int Update::update_timestep()
{
  update->ntimestep++;
  (*input->vars)["timestep"] = Var("timestep", ntimestep);
  return update->ntimestep;
}

/*! Write method, scheme, timestep, dt... to restart file.
 */
void Update::write_restart(ofstream *of) {
  // The informations to be stored in the restart file are:
  // - The method type
  // - The type of shape function
  // - FLIP and/or PIC, or APIC
  // - The scheme type
  // - The timestep
  // - dt
  // - dt_factor
  

  //of->write(reinterpret_cast<const char *>(&flag), sizeof(int));
}


/*! Write method, scheme, timestep, dt... to restart file.
 */
void Update::read_restart(ifstream *ifr) {
  // The informations to be stored in the restart file are:
  // - The method type
  // - The type of shape function
  // - FLIP and/or PIC, or APIC
  // - The scheme type
  // - The timestep
  // - dt
  // - dt_factor
}
