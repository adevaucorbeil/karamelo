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

#include <usf.h>
#include <iostream>
#include <vector>
#include <scheme.h>
#include <method.h>
#include <update.h>
#include <output.h>
#include <modify.h>
#include <universe.h>
#include <domain.h>

using namespace std;

USF::USF(MPM *mpm) : Scheme(mpm) {
  // cout << "In USF::USF()" << endl;
}

void USF::setup(){
  output->setup();
}

void USF::run(Var condition){

  bigint ntimestep = update->ntimestep;
  
  // cout << "In USF::run" << endl;

  output->write(ntimestep);

  Method &method = *update->method;

  while ((bool) condition.result(mpm))
  {
    // prepare
    ntimestep = update->update_timestep();
    for (Solid *solid: domain->solids)
      for (int ip = 0; ip < solid->np_local; ip++)
        method.compute_grid_weight_functions_and_gradients(*solid, ip);
    modify->prepare();
    
    // grid reset
    method.reset();
    for (Grid *grid: method.grids())
      for (int in = 0; in < grid->nnodes_local + grid->nnodes_ghost; in++)
      {
        method.reset_mass_nodes(*grid, in);
        method.reset_velocity_nodes(*grid, in);
      }

    // p2g 1
    for (Solid *solid: domain->solids)
      for (int ip = 0; ip < solid->np_local; ip++)
      {
        modify->initial_integrate(*solid, ip);
        method.compute_mass_nodes(*solid, ip);
      }

    for (Grid *grid: method.grids())
      grid->reduce_mass_ghost_nodes();

    for (Solid *solid: domain->solids)
      for (int ip = 0; ip < solid->np_local; ip++)
        method.compute_velocity_nodes(*solid, ip);
    
    // grid update 1
    for (Grid *grid: method.grids())
    {
      grid->reduce_ghost_nodes(true, false, method.temp);
      
      for (int in = 0; in < grid->nnodes_local + grid->nnodes_ghost; in++)
        modify->post_update_grid_state(*grid, in);
    }
    

    for (Solid *solid: domain->solids)
      for (int ip = 0; ip < solid->np_local; ip++)
      {
        method.compute_rate_deformation_gradient(true, *solid, ip);
        method.update_deformation_gradient_stress(true, *solid, ip);
      }

    for (Grid *grid: method.grids())
      for (int in = 0; in < grid->nnodes_local + grid->nnodes_ghost; in++)
        method.reset_force_nodes(*grid, in);

    // p2g 2
    for (Solid *solid: domain->solids)
      for (int ip = 0; ip < solid->np_local; ip++)
        method.compute_force_nodes(*solid, ip);
    
    // grid update 2
    for (Grid *grid: method.grids())
    {
      grid->reduce_ghost_nodes(false, true, method.temp);

      for (int in = 0; in < grid->nnodes_local + grid->nnodes_ghost; in++)
      {
        modify->post_particles_to_grid(*grid, in);
        method.update_grid_velocities(*grid, in);
        modify->post_update_grid_state(*grid, in);
      }
    }

    for (Solid *solid: domain->solids)
      for (int ip = 0; ip < solid->np_local; ip++)
      {
        method.compute_velocity_acceleration(*solid, ip);
        method.update_position(*solid, ip);
        modify->post_grid_to_point(*solid, ip);
        method.advance_particles(*solid, ip);
        modify->post_advance_particles(*solid, ip);
      }

    // grid update
    for (Grid *grid: method.grids())
    {
      grid->reduce_ghost_nodes(true, false, update->temp);
      
      for (int in = 0; in < grid->nnodes_local + grid->nnodes_ghost; in++)
      {
        modify->post_velocities_to_grid(*grid, in);
        method.update_grid_positions(*grid, in);
      }
    }

    method.exchange_particles();
    update->update_time();
    method.adjust_dt();
    
    for (Solid *solid: domain->solids)
      for (int ip = 0; ip < solid->np_local; ip++)
        modify->final_integrate(*solid, ip);

    if ((update->maxtime != -1) && (update->atime > update->maxtime)) {
      update->nsteps = ntimestep;
      output->write(ntimestep);
      break;
    }

    if (ntimestep == output->next || ntimestep == update->nsteps) {
      output->write(ntimestep);
    }
  }
}

