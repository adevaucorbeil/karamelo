#include "musl.h"
#include <iostream>
#include <vector>
#include "scheme.h"
#include "method.h"
#include "update.h"
#include "output.h"

using namespace std;

MUSL::MUSL(MPM *mpm, vector<string> args) : Scheme(mpm) {
  cout << "In MUSL::MUSL()" << endl;
}

void MUSL::setup(){
  output->setup();
  update->method->setup();
}

void MUSL::run(int n){

  bigint ntimestep;

  cout << "In MUSL::run" << endl;

  for (int i=0; i<n; i++){
    ntimestep = ++update->ntimestep;

    update->method->compute_grid_weight_functions_and_gradients();
    update->method->particles_to_grid();
    update->method->update_grid_state();
    update->method->grid_to_points();
    update->method->advance_particles();
    update->method->velocities_to_grid();
    update->method->compute_rate_deformation_gradient();
    update->method->update_deformation_gradient();
    update->method->update_stress();

    if (ntimestep == output->next) {
      output->write(ntimestep);
    }
    
  }
}

