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
  update->method->setup();
}

void MUSL::run(int n){

  bigint ntimestep = update->ntimestep;

  // cout << "In MUSL::run" << endl;

  update->method->compute_grid_weight_functions_and_gradients();
  output->write(ntimestep);

  for (int i=0; i<n; i++){
    ntimestep = ++update->ntimestep;

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


    modify->final_integrate();

    if (ntimestep == output->next) {
      output->write(ntimestep);
    }
    
    update->method->adjust_dt();
  }
}

