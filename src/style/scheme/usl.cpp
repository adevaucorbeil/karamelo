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

#include <usl.h>
#include <iostream>
#include <vector>
#include <scheme.h>
#include <method.h>
#include <update.h>
#include <output.h>
#include <modify.h>
#include <universe.h>
#include <domain.h>
#include <iomanip>

using namespace std;

USL::USL(MPM *mpm):
  Scheme(mpm)
{
  // cout << "In USL::USL()" << endl;
}

void USL::setup()
{
  output->setup();
}

void USL::run(Var condition)
{
  bigint ntimestep = update->ntimestep;
  
  // cout << "In USL::run" << endl;

  output->write(ntimestep);

  Method &method = *update->method;

  while (condition.result(mpm))
  {
    // prepare
    ntimestep = update->update_timestep();
    for (Solid *solid: domain->solids)
      method.compute_grid_weight_functions_and_gradients(*solid);
    modify->prepare();
    
    // grid reset
    method.reset();
    for (Grid *grid: method.grids())
    {
      method.reset_mass_nodes(*grid);
      method.reset_velocity_nodes(*grid);
      method.reset_force_nodes(*grid);
    }

    // p2g
    for (Solid *solid: domain->solids)
    {
      modify->initial_integrate(*solid);
      method.compute_mass_nodes(*solid);
    }

    for (Grid *grid: method.grids())
      grid->reduce_mass_ghost_nodes();

    for (Solid *solid: domain->solids)
    {
      method.compute_velocity_nodes(*solid);
      method.compute_force_nodes(*solid);
    }

    // grid update
    for (Grid *grid: method.grids())
    {
      grid->reduce_ghost_nodes(true, true, method.temp);

      modify->post_particles_to_grid(*grid);
      method.update_grid_velocities(*grid);
      modify->post_update_grid_state(*grid);
    }
    
    // g2p
    for (Solid *solid: domain->solids)
    {
      method.compute_rate_deformation_gradient(false, *solid);
      method.compute_velocity_acceleration(*solid);
      method.update_position(*solid);
      modify->post_grid_to_point(*solid);
      method.advance_particles(*solid);

      modify->post_advance_particles(*solid);
    }

    // grid update
    for (Grid *grid: method.grids())
    {
      grid->reduce_ghost_nodes(true, false, update->temp);
      
      method.update_grid_positions(*grid);
    }
    
    for (Solid *solid: domain->solids) {
      method.update_deformation_gradient(*solid);
      if (method.anti_volumetric_locking)
	method.Fbar_anti_vol_locking(*solid);
      method.update_stress(*solid);
    }

    method.exchange_particles();
    update->update_time();
    method.adjust_dt();

    for (Solid *solid: domain->solids)
      modify->final_integrate(*solid);

    modify->reduce();

    if (update->maxtime != -1 && update->atime > update->maxtime)
    {
      update->nsteps = ntimestep;
      output->write(ntimestep);

      break;
    }

    if (ntimestep == output->next || ntimestep == update->nsteps)
      output->write(ntimestep);
  }
}

