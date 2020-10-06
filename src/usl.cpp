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

#include "usl.h"
#include <iostream>
#include <vector>
#include "scheme.h"
#include "method.h"
#include "update.h"
#include "output.h"
#include "modify.h"
#include "universe.h"

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

  //for (int i=0; i<nsteps; i++){
  while ((bool) condition.result(mpm)) {
    ntimestep = update->update_timestep();

    update->method->compute_grid_weight_functions_and_gradients();

    update->method->reset();
    modify->initial_integrate();

    update->method->particles_to_grid();

    modify->post_particles_to_grid();

    update->method->update_grid_state();

    modify->post_update_grid_state();

    update->method->grid_to_points();

    modify->post_grid_to_point();

    update->method->advance_particles();

    modify->post_advance_particles();
    
    //update->method->velocities_to_grid();

    //modify->post_velocities_to_grid();

    update->method->update_grid_positions();

    update->method->compute_rate_deformation_gradient();
    update->method->update_deformation_gradient();
    update->method->update_stress();

    update->method->exchange_particles();

    update->method->adjust_dt();
    update->update_time();

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
  output->show_plot();
}

