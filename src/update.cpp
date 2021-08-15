/* ----------------------------------------------------------------------
 *
 *                    ***       Karamelo       ***
 *               Parallel Material Point Method Simulator
 * 
 * Copyright (2019) Alban de Vaucorbeil, alban.devaucorbeil@monash.edu
 * Materials Science and Engineering, Monash University
 * Clayton VIC 3800, Australia
 * 
 * This software is distributed under the GNU General Public License.
 *
 * ----------------------------------------------------------------------- */

#include "update.h"
#include "error.h"
#include "input.h"
#include "method.h"
#include "scheme.h"
#include "style_method.h"
#include "style_scheme.h"
#include "universe.h"
#include "var.h"
#include <iostream>
#include <vector>

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

  method = nullptr;
  shape_function = ShapeFunctions::LINEAR;
  sub_method_type = SubMethodType::FLIP;
  PIC_FLIP = 0.99;
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
  else if (scheme_style.compare(#key) == 0) scheme = new Class(mpm);
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
  else if (method_type.compare(#key) == 0) method = new Class(mpm);
#include "style_method.h"
#undef MethodStyle
#undef METHOD_CLASS

  else {
    error->all(FLERR, "Illegal method style.\n");
  }

  bool isFLIP = false;
  n++;
  // Method used: PIC, FLIP, APIC, or ASFLIP:
  if (map_sub_method_type.count(args[n]) > 0) {
    sub_method_type = map_sub_method_type.at(args[n]);
    if (sub_method_type == SubMethodType::PIC)
      PIC_FLIP = 0;
    else if (sub_method_type == SubMethodType::FLIP ||
             sub_method_type == SubMethodType::AFLIP ||
             sub_method_type == SubMethodType::ASFLIP) {
      isFLIP = true;
      if (args.size() < 4) {
	error->all(FLERR, "Illegal modify_method command: not enough arguments.\n");
      }
    }
  } else {
    error->all(
        FLERR,
        "Error: method type " + args[n] +
            " not understood. Expect: PIC, FLIP, APIC, AFLIP, or ASFLIP\n");
  }

  n++;
  
  // Check if the shape function given in the inputfile is recognized.
  if (args.size() > n + isFLIP) {
    if (map_shape_functions.count(args[n]) > 0) {
      shape_function = map_shape_functions.at(args[n]);
      n++;
    } else {
      error->all(FLERR, "Illegal method_method argument: form function of type \033[1;31m" + args[2] + "\033[0m is unknown. Available options are:  \033[1;32mlinear\033[0m, \033[1;32mcubic-spline\033[0m, \033[1;32mquadratic-spline\033[0m, \033[1;32mBernstein-quadratic\033[0m.\n");
    }
  }

  if (isFLIP) {
    PIC_FLIP = input->parsev(args[n]);
    n++;
  }

  if (args.size() >= n + 1) {
    if (args[n].compare("thermo-mechanical") == 0) {
      method->temp = true;
    } else if (args[n].compare("mechanical") == 0) {
      method->temp = false;
    } else {
      error->all(
          FLERR,
          "Illegal modify_method command: keyword " + args[n] +
              " unknown. Expected \"thermo-mechanical\" or \"mechanical\".\n");
    }
  }
  n++;

  if (args.size() >= n + 1) {
    if (args[n].compare("gradient-enhanced") == 0) {
      method->ge = true;
    } else {
      error->all(FLERR,
                 "Illegal modify_method command: keyword " + args[n] +
                     " unknown. Expected \"gradient-enhanced\" or nothing.\n");
    }
  }
  n++;

  if (n < args.size()) {
    additional_args = vector<string>(args.begin() + n, args.end());
  } else {
    additional_args.clear();
  }
  method->setup(additional_args);
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
  // - If temperature is involved (method->temp)
  // - The scheme style
  // - FLIP and/or PIC, or APIC
  // - PIC_FLIP
  // - The type of shape function
  // - additional_args
  // - The timestep
  // - dt
  // - dt_factor
  // - dt_contant
  
  // Method type:
  size_t N = method_type.size();
  of->write(reinterpret_cast<const char *>(&N), sizeof(size_t));
  of->write(reinterpret_cast<const char *>(method_type.c_str()), N);

  // method->temp
  of->write(reinterpret_cast<const char *>(&method->temp), sizeof(bool));

  // Scheme style:
  N = scheme_style.size();
  of->write(reinterpret_cast<const char *>(&N), sizeof(size_t));
  of->write(reinterpret_cast<const char *>(scheme_style.c_str()), N);

  // Sub-method type (PIC and/or FLIP or APIC):
  of->write(reinterpret_cast<const char *>(&sub_method_type), sizeof(int));

  // PIC_FLIP
  of->write(reinterpret_cast<const char *>(&PIC_FLIP), sizeof(double));
  // cout << "PIC_FLIP=" << PIC_FLIP << endl;

  // The type of shape function:
  of->write(reinterpret_cast<const char *>(&shape_function), sizeof(int));

  // additional_args: 
  N = additional_args.size();
  of->write(reinterpret_cast<const char *>(&N), sizeof(size_t));
  for (size_t i = 0; i < N; i++) {
    size_t Ns = additional_args[i].size();
    of->write(reinterpret_cast<const char *>(&Ns), sizeof(size_t));
    of->write(reinterpret_cast<const char *>(additional_args[i].c_str()), Ns);    
  }

  // atime:
  of->write(reinterpret_cast<const char *>(&atime), sizeof(double));

  // Timestep:
  of->write(reinterpret_cast<const char *>(&ntimestep), sizeof(bigint));

  // dt:
  of->write(reinterpret_cast<const char *>(&dt), sizeof(double));

  // dt_factor:
  of->write(reinterpret_cast<const char *>(&dt_factor), sizeof(double));

  // dt_contant:
  of->write(reinterpret_cast<const char *>(&dt_constant), sizeof(bool));
}


/*! Write method, scheme, timestep, dt... to restart file.
 */
void Update::read_restart(ifstream *ifr) {
  // The informations to be stored in the restart file are:  // - The method type
  // - The scheme style
  // - FLIP and/or PIC, or APIC
  // - PIC_FLIP
  // - The type of shape function
  // - The timestep
  // - dt
  // - dt_factor
  // - dt_contant
  
  // Method type:
  size_t N = 0;
  ifr->read(reinterpret_cast<char *>(&N), sizeof(size_t));
  method_type.resize(N);
  ifr->read(reinterpret_cast<char *>(&method_type[0]), N);
  // cout << "method_type=" << method_type << endl;

  if (0) return;

#define METHOD_CLASS
#define MethodStyle(key,Class) \
  else if (method_type.compare(#key) == 0) method = new Class(mpm);
#include "style_method.h"
#undef MethodStyle
#undef METHOD_CLASS

  else {
    error->all(FLERR, "Illegal method style.\n");
  }


  // method->temp
  ifr->read(reinterpret_cast<char *>(&method->temp), sizeof(bool));

  // Scheme style:
  N = 0;
  ifr->read(reinterpret_cast<char *>(&N), sizeof(size_t));
  scheme_style.resize(N);
  ifr->read(reinterpret_cast<char *>(&scheme_style[0]), N);
  // cout << "scheme_style=" << scheme_style << endl;

  if (0) return;

#define SCHEME_CLASS
#define SchemeStyle(key,Class) \
  else if (scheme_style.compare(#key) == 0) scheme = new Class(mpm);
#include "style_scheme.h"
#undef SchemeStyle
#undef SCHEME_CLASS

  else {
    error->all(FLERR, "Illegal scheme style.\n");
  }
  // Sub-method type (PIC and/or FLIP or APIC):
  ifr->read(reinterpret_cast<char *>(&sub_method_type), sizeof(int));
  // cout << "sub_method_type=" << sub_method_type << endl;

  // PIC_FLIP
  ifr->read(reinterpret_cast<char *>(&PIC_FLIP), sizeof(double));
  // cout << "PIC_FLIP=" << PIC_FLIP << endl;

  // The type of shape function:
  ifr->read(reinterpret_cast<char *>(&shape_function), sizeof(int));
  // cout << "shape_functions=" << shape_function << endl;

  // additional_args: 
  ifr->read(reinterpret_cast<char *>(&N), sizeof(size_t));
  additional_args.resize(N);
  for (size_t i = 0; i < N; i++) {
    size_t Ns = 0;
    ifr->read(reinterpret_cast<char *>(&Ns), sizeof(size_t));
    additional_args[i].resize(Ns);
    ifr->read(reinterpret_cast<char *>(&additional_args[i][0]), Ns);
    // cout << "additional_args[" << i << "]=" << additional_args[i] << endl;
  }
  method->setup(additional_args);

  // atime:
  ifr->read(reinterpret_cast<char *>(&atime), sizeof(double));
  // cout << "atime=" << atime << endl;
  (*input->vars)["time"] = Var("time", atime);

  // Timestep:
  ifr->read(reinterpret_cast<char *>(&ntimestep), sizeof(bigint));
  // cout << "ntimestep=" << ntimestep << endl;
  atimestep = ntimestep;
  (*input->vars)["timestep"] = Var("timestep", ntimestep);

  // dt:
  ifr->read(reinterpret_cast<char *>(&dt), sizeof(double));
  // cout << "dt=" << dt << endl;

  // dt_factor:
  ifr->read(reinterpret_cast<char *>(&dt_factor), sizeof(double));
  // cout << "dt_factor=" << dt_factor << endl;

  // dt_contant:
  ifr->read(reinterpret_cast<char *>(&dt_constant), sizeof(bool));
  // cout << "dt_constant=" << dt_constant << endl;

  if (dt_constant)
    (*input->vars)["dt"] = Var("dt", dt);
  
}
