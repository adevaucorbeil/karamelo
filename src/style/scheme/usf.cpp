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
      method.compute_grid_weight_functions_and_gradients(*solid);
    modify->prepare();
    
    // grid reset
    method.reset();
    for (Grid *grid: method.grids())
    {
      method.reset_mass_nodes(*grid);
      method.reset_velocity_nodes(*grid);
    }

    // p2g 1
    for (Solid *solid: domain->solids)
    {
      modify->initial_integrate(*solid);
      method.compute_mass_nodes(*solid);
    }

    for (Grid *grid: method.grids())
      grid->reduce_mass_ghost_nodes();

    for (Solid *solid: domain->solids)
      method.compute_velocity_nodes(*solid);
    
    // grid update 1
    for (Grid *grid: method.grids())
    {
      grid->reduce_ghost_nodes(true, false, method.temp);
      
      modify->post_update_grid_state(*grid);
    }
    

    for (Solid *solid: domain->solids)
    {
      method.compute_rate_deformation_gradient(true, *solid);
      method.update_deformation_gradient_stress(true, *solid);
    }

    for (Grid *grid: method.grids())
      method.reset_force_nodes(*grid);

    // p2g 2
    for (Solid *solid: domain->solids)
      method.compute_force_nodes(*solid);
    
    // grid update 2
    for (Grid *grid: method.grids())
    {
      grid->reduce_ghost_nodes(false, true, method.temp);
      modify->post_particles_to_grid(*grid);
      method.update_grid_velocities(*grid);
      modify->post_update_grid_state(*grid);
    }

    for (Solid *solid: domain->solids)
    {
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
      modify->post_velocities_to_grid(*grid);
      method.update_grid_positions(*grid);
    }

    method.exchange_particles();
    update->update_time();
    method.adjust_dt();
    
    for (Solid *solid: domain->solids)
      modify->final_integrate(*solid);

    modify->reduce();

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

