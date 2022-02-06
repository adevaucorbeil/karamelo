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

using namespace std;

USL::USL(MPM *mpm) : Scheme(mpm) {
  // cout << "In USL::USL()" << endl;
}

void USL::setup(){
  output->setup();
}

void USL::run(Var condition){

  bigint ntimestep = update->ntimestep;
  
  // cout << "In USL::run" << endl;

  output->write(ntimestep);

  Method &method = *update->method;

  while (condition.result(mpm))
  {
    // prepare
    ntimestep = update->update_timestep();
    method.compute_grid_weight_functions_and_gradients();
    modify->prepare();
    
    // grid reset
    method.reset();
    method.reset_mass_nodes();
    method.reset_nodes();

    // p2g
    for (Solid *solid: domain->solids)
      for (int ip = 0; ip < solid->np_local; ip++)
        modify->initial_integrate(*solid, ip);

    for (Solid *solid: domain->solids)
      for (int i = 0; i < solid->neigh_n.size(); i++)
        method.compute_mass_nodes(*solid,
                                   solid->neigh_n.at(i),
                                   solid->neigh_p.at(i),
                                   solid->wf.at(i));

    method.reduce_mass_ghost_nodes();

    for (Solid *solid: domain->solids)
      for (int i = 0; i < solid->neigh_n.size(); i++)
      {
        int in = solid->neigh_n.at(i);
        int ip = solid->neigh_p.at(i);
        double wf = solid->wf.at(i);
        const Vector3d &wfd = solid->wfd.at(i);

        method.compute_velocity_nodes(*solid, in, ip, wf);
        method.compute_force_nodes(*solid, in, ip, wf, wfd);
        if (update->temp)
        {
          method.compute_temperature_nodes(*solid, in, ip, wf);
          method.compute_temperature_driving_force_nodes(*solid, in, ip, wf, wfd);
        }
      }

    method.reduce_ghost_nodes();

    // grid update
    //modify->post_particles_to_grid();

    method.update_grid_state();

    //modify->post_update_grid_state();

    method.compute_rate_deformation_gradient(false);

    // g2p
    method.grid_to_points();

    modify->post_grid_to_point();

    method.advance_particles();

    //modify->post_advance_particles();
    
    //method.velocities_to_grid();

    //modify->post_velocities_to_grid();

    method.update_grid_positions();

    method.update_deformation_gradient();
    method.update_stress(false);

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

