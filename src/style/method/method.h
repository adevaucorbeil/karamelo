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

#ifndef MPM_METHOD_H
#define MPM_METHOD_H

#include <pointers.h>
#include <matrix.h>
#include <vector>
#include <Kokkos_Core.hpp>

class Solid;
class Grid;

/*! Parent class of all the MPM methods supported: ULMPM, TLMPM, ULCPDI and TLCPDI, to date.
 *
 * Stores the type of method used, and virtual descriptions of the functions used outside of
 * the class such as grid_to_points() or update_stress().
 */
class Method : protected Pointers {
 public:
  int style = 0;

  Method(class MPM *);

  void setup(vector<string> args);
  void compute_grid_weight_functions_and_gradients(Solid &solid);

  bool apic();
  vector<Grid *> grids();
  
  void reset_mass_nodes(Grid &grid);
  void compute_mass_nodes(Solid &solid);

  void reset_velocity_nodes(Grid &grid);
  void compute_velocity_nodes(Solid &solid);

  void reset_force_nodes(Grid &grid);
  void compute_force_nodes(Solid &solid);

  void update_grid_velocities(Grid &grid);

  void compute_velocity_acceleration(Solid &solid);

  void update_position(Solid &solid);
  void advance_particles(Solid &solid);

  void compute_rate_deformation_gradient(bool doublemapping, Solid &solid);

  void update_deformation_gradient(Solid &solid);
  void Fbar_anti_vol_locking(Solid &solid);
  void update_stress(Solid &solid);
  
  void adjust_dt();
  void reset();
  void exchange_particles();

  bool is_TL;         ///< true: the method is total Lagrangian; false: it is updated Lagrangian
  bool is_CPDI;       ///< true if the method is a CPDI-like
  bool temp;          ///< true for thermo-mechanical simulations
  bool anti_volumetric_locking; ///< true when anti volumetric locking measures are used
  int rigid_solids;
  bool slip_contacts; ///< true when slip contacts are activated, false otherwise
};

#endif

/*! \defgroup method method
  
\section Syntax Syntax
\code
method(MPM_type, velocity_updating_method, shape_function, optional:method_specific_arguments)
\endcode

<ul>
<li>MPM_type: the type of MPM used: either tlmpm, ulmpm, tlcpdi, or ulcpdi.</li>
<li>velocity_updating_method: the type of method used to update/interpolate velocities: either PIC, FLIP, or APIC.</li>
<li>shape_function: type of shape functions to be used: either linear, cubic-spline or Bernstein-quadratic.</li>
<li>method_specific_arguments: optional list of arguments specific to the MPM or velocity updating method used. Typically, the PIC/FLIP mixing ratio (to be used with velocity_updating_method == FLIP). </li>
</ul>

\section Examples Examples
\code
alpha = 1
method(tlmpm, FLIP, Bernstein-quadratic, alpha)
\endcode
Defines a Total Lagrangian MPM simulation using a mixture of PIC and FLIP for the update of the particle velocities with a FLIP/PIC ratio of \f$alpha=1\f$ with Bernstein-quadratic shape functions. It is therefore TLMPM using pure FLIP.

\code
method(ulmpm, APIC, cubic-spline)
\endcode
Defines a Updated Lagrangian MPM (aka standard MPM) simulation using APIC for the velocity interpolation and update with cubic B-spline shape functions.

\section Description Description

This function sets:
<ol>
<li>the type of MPM algorithm to be used. To date total Lagrangian, updated MPM and CPDI are implemented. To use them, MPM_type has to be set to either: tlmpm, ulmpm, tlcpdi, ulcpdi, respectively.</li>
<li>the algorithm used to update/interpolate the velocities. APIC (APIC), PIC (PIC), FLIP (FLIP) or a mixture of FLIP and PIC (FLIP)can be used. The keyword FLIP is used for either pure FLIP, or PIC/FLIP mixture. Therefore, if specifying FLIP, an extra argument is necessary: the ratio between FLIP and PIC.</li>
<li>the type of shape function used. To date this includes: linear shape functions (linear), cubic B-splines (cubic-spline), and Bernstein quadratic functions (Bernstein-quadratic). Bernstein shape functions are only supported with TLMPM.</li>
</ol>

*/
