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
    method.compute_grid_weight_functions_and_gradients();
    modify->prepare();
    
    // grid reset
    method.reset();
    for (Grid *grid: method.grids())
      for (int in = 0; in < grid->nnodes_local + grid->nnodes_ghost; in++)
      {
        grid->mass.at(in) = 0;
        grid->v.at(in) = Vector3d();
        grid->mb.at(in) = Vector3d();
        if (method.temp)
          grid->T.at(in) = 0;
      }

    // p2g 1
    for (Solid *solid: domain->solids)
      for (int ip = 0; ip < solid->np_local; ip++)
        modify->initial_integrate(*solid, ip);

    for (Solid *solid: domain->solids)
      for (int i = 0; i < solid->neigh_n.size(); i++)
        method.compute_mass_nodes(*solid,
                                   solid->neigh_n.at(i),
                                   solid->neigh_p.at(i),
                                   solid->wf.at(i));

    for (Grid *grid: method.grids())
      grid->reduce_mass_ghost_nodes();

    for (Solid *solid: domain->solids)
      for (int i = 0; i < solid->neigh_n.size(); i++)
      {
        int in = solid->neigh_n.at(i);
        int ip = solid->neigh_p.at(i);
        double wf = solid->wf.at(i);

        method.compute_velocity_nodes(*solid, in, ip, wf);
        if (update->temp)
          method.compute_temperature_nodes(*solid, in, ip, wf);
      }
    
    for (Grid *grid: method.grids())
      grid->reduce_ghost_nodes(true, false, method.temp);

    // grid update 1
    //modify->post_update_grid_state();

    method.compute_rate_deformation_gradient(true);
    method.update_deformation_gradient();
    method.update_stress(true);

    for (Grid *grid: method.grids())
      for (int in = 0; in < grid->nnodes_local + grid->nnodes_ghost; in++)
      {
        grid->mass.at(in) = 0;
        grid->f.at(in) = Vector3d();
        grid->mb.at(in) = Vector3d();
        if (method.temp)
        {
          grid->Qint.at(in) = 0;
          grid->Qext.at(in) = 0;
        }
      }

    // p2g 2
    for (Solid *solid: domain->solids)
      for (int i = 0; i < solid->neigh_n.size(); i++)
      {
        int in = solid->neigh_n.at(i);
        int ip = solid->neigh_p.at(i);
        double wf = solid->wf.at(i);
        const Vector3d &wfd = solid->wfd.at(i);

        method.compute_force_nodes(*solid, in, ip, wf, wfd);
        if (update->temp)
          method.compute_temperature_driving_force_nodes(*solid, in, ip, wf, wfd);
      }
    
    // grid update 2
    for (Grid *grid: method.grids())
    {
      grid->reduce_ghost_nodes(false, true, method.temp);

      for (int in = 0; in < grid->nnodes_local + grid->nnodes_ghost; in++)
      {
        modify->post_particles_to_grid(*grid, in);

        method.update_grid_velocities(*grid, in);
        if (method.temp)
          method.update_grid_temperature(*grid, in);

        modify->post_update_grid_state(*grid, in);
      }
    }

    method.grid_to_points();

    modify->post_grid_to_point();

    method.advance_particles();

    //modify->post_advance_particles();

    method.update_grid_positions();

    method.exchange_particles();

    update->update_time();
    method.adjust_dt();

    //modify->final_integrate();

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

