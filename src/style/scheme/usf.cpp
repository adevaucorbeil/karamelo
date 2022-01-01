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

  //for (int i=0; i<nsteps; i++){
  while ((bool) condition.result(mpm)) {
    ntimestep = update->update_timestep();

    update->method->compute_grid_weight_functions_and_gradients();

    update->method->reset();
    modify->initial_integrate();

    update->method->particles_to_grid_USF_1();

    modify->post_update_grid_state();

    update->method->compute_rate_deformation_gradient(true);
    update->method->update_deformation_gradient();
    update->method->update_stress(true);

    update->method->particles_to_grid_USF_2();

    modify->post_particles_to_grid();

    update->method->update_grid_state();

    modify->post_update_grid_state();

    update->method->grid_to_points();

    modify->post_grid_to_point();

    update->method->advance_particles();

    modify->post_advance_particles();

    update->method->update_grid_positions();

    update->method->exchange_particles();

    update->update_time();
    update->method->adjust_dt();

    modify->final_integrate();

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

