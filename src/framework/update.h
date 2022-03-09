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

#ifndef MPM_UPDATE_H
#define MPM_UPDATE_H

#include <pointers.h>
#include <vector>
#include <map>
/*! Stores everything related to time steps as well as the Scheme and Method classes.
 */
class Update : protected Pointers {
 public:
  enum class ShapeFunctions;
  enum class SubMethodType;

  double run_duration;                ///< Stop simulation if elapsed simulation time exceeds this.
  double elapsed_time_in_run;	      ///< Elapsed simulation time for a single run;
  double dt;                          ///< Timestep
  double dt_factor;                   ///< Timestep factor
  bool dt_constant;                   ///< is dt constant?
  bigint ntimestep;                   ///< current step
  int nsteps;                         ///< Number of steps to run
  double atime;                       ///< Simulation time at atime_step
  double maxtime;                     ///< Maximum simulation time (infinite if -1)
  bigint atimestep;                   ///< Last timestep atime was updated
  bigint firststep,laststep;          ///< 1st & last step of this run
  bigint beginstep,endstep;           ///< 1st and last step of multiple runs
  int first_update;                   ///< 0 before initial update, 1 after

  class Scheme *scheme;               ///< Pointer to the type of Scheme used
  string scheme_style;                ///< Name of the scheme style

  class Method *method;               ///< Pointer to the type of Method used
  string method_type;                 ///< Name of the method type
  SubMethodType sub_method_type;      ///< Name of the velocity updating method type
  ShapeFunctions shape_function;      ///< Type of shape function used
  double PIC_FLIP;                    ///< PIC/FLIP mixing factor
  bool temp;                          ///< True for thermo-mechanical simulations

  Update(class MPM *);
  ~Update();
  void set_dt_factor(vector<string>); ///< Sets the factor to be applied to the CFL timestep
  void set_dt(vector<string>);        ///< Sets the timestep
  void create_scheme(vector<string>); ///< Creates a scheme: USL, or MUSL.
  void create_method(vector<string>); ///< Creates a method: tlmpm, ulmpm, tlcpdi, ...
  void update_time();                 ///< Update elapsed time
  int update_timestep();              ///< Update timestep
  void write_restart(ofstream*);      ///< Write method, scheme, timestep, dt... to restart file
  void read_restart(ifstream*);       ///< Read method, scheme, timestep, dt... to restart file

  enum class SubMethodType {
    PIC,
    FLIP,
    APIC,
    MLS,
  };
  enum class ShapeFunctions {
    LINEAR,
    CUBIC_SPLINE,
    QUADRATIC_SPLINE,
    BERNSTEIN,
  };

private:
  const map<string, SubMethodType> map_sub_method_type{{"PIC", SubMethodType::PIC},
                                                       {"FLIP", SubMethodType::FLIP},
                                                       {"APIC", SubMethodType::APIC},
                                                       {"MLS", SubMethodType::MLS}};
  const map<string, ShapeFunctions> map_shape_functions{
      {"linear", ShapeFunctions::LINEAR},
      {"cubic-spline", ShapeFunctions::CUBIC_SPLINE},
      {"quadratic-spline", ShapeFunctions::QUADRATIC_SPLINE},
      {"Bernstein-quadratic", ShapeFunctions::BERNSTEIN}};

  vector<string> additional_args;     ///< Read method, scheme, timestep, dt... to restart file
};

#endif
