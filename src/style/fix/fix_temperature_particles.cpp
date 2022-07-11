/* ----------------------------------------------------------------------
 *
 *                    ***       Karamelo       ***
 *               Parallel Material Point Method Simulator
 * 
 * Copyright (2020) Alban de Vaucorbeil, alban.devaucorbeil@deakin.edu.au
 * Institute for Frontier Materials, Deakin University
 * Geelong VIC 3216, Australia

 * This software is distributed under the GNU General Public License.
 *
 * ----------------------------------------------------------------------- */

#include <fix_temperature_particles.h>
#include <domain.h>
#include <error.h>
#include <group.h>
#include <input.h>
#include <special_functions.h>
#include <universe.h>
#include <update.h>
#include <expression_operation.h>


using namespace std;
using namespace FixConst;



FixTemperatureParticles::FixTemperatureParticles(MPM *mpm, vector<string> args):
  Fix(mpm, args, INITIAL_INTEGRATE | POST_ADVANCE_PARTICLES)
{
  if (args.size() < 3) {
    error->all(FLERR, "Error: not enough arguments.\n");
  }

  if (args[2].compare("restart") ==
      0) { // If the keyword restart, we are expecting to have read_restart()
           // launched right after.
    igroup = stoi(args[3]);
    if (igroup == -1 && universe->me == 0) {
      cout << "Could not find group number " << args[3] << endl;
    }
    groupbit = group->bitmask[igroup];

    return;
  }

  if (args.size() < Nargs) {
    error->all(FLERR, "Error: too few arguments for fix_temperature_nodes.\n" +
                          usage);
  }

  if (group->pon[igroup].compare("particles") !=0 ) {
    error->one(FLERR, "fix_temperature_nodes needs to be given a group of nodes" +
                          group->pon[igroup] + ", " + args[2] +
                          " is a group of " + group->pon[igroup] + ".\n");
  }
  if (universe->me == 0) {
    cout << "Creating new fix FixTemperatureParticles with ID: " << args[0] << endl;
  }
  id = args[0];

  string time = "time";
  string previous = args[3];

  input->parsev(args[3]);
  Tvalue = &input->expressions[args[3]];


  // Replace "time" by "time - dt" in the x argument:
  previous = SpecialFunc::replace_all(input->parsev(previous).str(), "time", "(time - dt)");
  input->parsev(previous);
  Tprevvalue = &input->expressions[previous];
}

void FixTemperatureParticles::initial_integrate(Solid &solid)
{
  // Update the temperatures:
  Tprevvalue->evaluate(solid);

  int groupbit = this->groupbit;
  Kokkos::View<int*> mask = solid.mask;
  Kokkos::View<float*> T = solid.T;

  Kokkos::View<float **> Tpv = Tprevvalue->registers;

  Kokkos::parallel_for("FixTemperatureParticles::initial_integrate", solid.np_local,
		       KOKKOS_LAMBDA(const int &ip)
      {
        if (!(mask[ip] & groupbit))
          return;

        T[ip]        = Tpv(0, ip);
      });
}

void FixTemperatureParticles::post_advance_particles(Solid &solid)
{
  // Update the temperatures:
  Tvalue->evaluate(solid);

  int groupbit = this->groupbit;
  Kokkos::View<int*> mask = solid.mask;
  Kokkos::View<float*> T = solid.T;

  Kokkos::View<float **> Tv = Tvalue->registers;

  Kokkos::parallel_for("FixTemperatureParticles::post_advance_particles", solid.np_local,
		       KOKKOS_LAMBDA(const int &ip)
      {
        if (!(mask[ip] & groupbit))
          return;

        T[ip]        = Tv(0, ip);
      });
}

void FixTemperatureParticles::write_restart(ofstream *of) {
  // Tvalue.write_to_restart(of);
  // Tprevvalue.write_to_restart(of);
}

void FixTemperatureParticles::read_restart(ifstream *ifr) {
  // Tvalue.read_from_restart(ifr);
  // Tprevvalue.read_from_restart(ifr);
}
