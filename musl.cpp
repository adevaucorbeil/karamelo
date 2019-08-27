#include "musl.h"
#include <iostream>
#include <vector>
#include "scheme.h"
#include "method.h"
#include "update.h"
#include "output.h"
#include "modify.h"

using namespace std;

MUSL::MUSL(MPM *mpm, vector<string> args) : Scheme(mpm) {
  cout << "In MUSL::MUSL()" << endl;
}

void MUSL::setup(){
  output->setup();
}

void MUSL::run(int n){

  bigint ntimestep = update->ntimestep;
  int nsteps = update->nsteps;
  nsteps = n;

  // cout << "In MUSL::run" << endl;

  output->write(ntimestep);

  for (int i=0; i<nsteps; i++){
    ntimestep = ++update->ntimestep;

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
    
    update->method->velocities_to_grid();

    modify->post_velocities_to_grid();

    update->method->compute_rate_deformation_gradient();
    update->method->update_deformation_gradient();
    update->method->update_stress();

    update->method->adjust_dt();

    modify->final_integrate();

    if ((update->maxtime != -1) && (update->atime > update->maxtime)) {
      update->nsteps = ntimestep;
      output->write(ntimestep);
      break;
    }

    if (ntimestep == output->next || ntimestep == nsteps) {
      output->write(ntimestep);
    }
    
  }
}

