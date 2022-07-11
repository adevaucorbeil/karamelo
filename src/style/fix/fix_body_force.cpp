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

#include <fix_body_force.h>
#include <domain.h>
#include <error.h>
#include <expression_operation.h>
#include <grid.h>
#include <group.h>
#include <input.h>
#include <universe.h>

using namespace std;
using namespace FixConst;


FixBodyForce::FixBodyForce(MPM *mpm, vector<string> args):
  Fix(mpm, args, POST_PARTICLES_TO_GRID)
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

    xset = yset = zset = false;
    return;
  }

  if (args.size() < Nargs.find(domain->dimension)->second) {
    error->all(FLERR, "Error: too few arguments for fix_body_force.\n" +
                          usage.find(domain->dimension)->second);
  }

  if (group->pon[igroup].compare("nodes") !=0  && group->pon[igroup].compare("all") !=0) {
    error->one(FLERR, "fix_body_force needs to be given a group of nodes" + group->pon[igroup] + ", " + args[2] + " is a group of "+ group->pon[igroup] + ".\n");
  }
  if (universe->me == 0) {
    cout << "Creating new fix FixBodyForce with ID: " << args[0] << endl;
  }
  id = args[0];

  xset = yset = zset = false;

  if (args[3] != "NULL") {
    xset = true;
    input->parsev(args[3]);

    fb[0] = &input->expressions[args[3]];
  }

  if (domain->dimension >= 2 && args[4] != "NULL")
  {
      yset = true;
      input->parsev(args[4]);

      fb[1] = &input->expressions[args[4]];
  } else {
    fb[1] = nullptr;
  }

  if (domain->dimension == 3 && args[5] != "NULL")
    {
      input->parsev(args[5]);
      zset = true;

      fb[2] = &input->expressions[args[5]];
    }
  else {
    fb[2] = nullptr;
  }
}

void FixBodyForce::prepare()
{
  ftot = Vector3d();
}

void FixBodyForce::reduce()
{
  Vector3d ftot_reduced;
  MPI_Allreduce(ftot.elements, ftot_reduced.elements, 3, MPI_FLOAT, MPI_SUM, universe->uworld);
}

void FixBodyForce::post_particles_to_grid(Grid &grid)
{
  
  // Go through all the nodes in the group and set fb to the right value:

  for (int i = 0; i < 3; i++)
    if (fb[i])
      fb[i]->evaluate(grid);

  int groupbit = this->groupbit;
  Kokkos::View<float*> mass = grid.mass;
  Kokkos::View<int*> mask = grid.mask;
  Kokkos::View<Vector3d*> mb = grid.mb;

  for (int i = 0; i < 3; i++)
    if (fb[i])
    {
      Kokkos::View<float **> fb_i = fb[i]->registers;

      Kokkos::parallel_reduce("FixBodyForce::post_particles_to_grid", grid.nnodes_local + grid.nnodes_ghost,
                              KOKKOS_LAMBDA(const int &in, float &ftot_i)
      {
	if (!mass[in] || !(mask[in] & groupbit))
          return;

	mb[in][i] += mass[in]*fb_i(0, in);

        ftot_i += fb_i(0, in);
      }, ftot[i]);
    }
}

void FixBodyForce::write_restart(ofstream *of) {
  of->write(reinterpret_cast<const char *>(&xset), sizeof(bool));
  of->write(reinterpret_cast<const char *>(&yset), sizeof(bool));
  of->write(reinterpret_cast<const char *>(&zset), sizeof(bool));

  // if (xset)
  //   xvalue.write_to_restart(of);
  // if (yset)
  //   yvalue.write_to_restart(of);
  // if (zset)
  //   zvalue.write_to_restart(of);
}

void FixBodyForce::read_restart(ifstream *ifr) {
  ifr->read(reinterpret_cast<char *>(&xset), sizeof(bool));
  ifr->read(reinterpret_cast<char *>(&yset), sizeof(bool));
  ifr->read(reinterpret_cast<char *>(&zset), sizeof(bool));

  // if (xset)
  //   xvalue.read_from_restart(ifr);
  // if (yset)
  //   yvalue.read_from_restart(ifr);
  // if (zset)
  //   zvalue.read_from_restart(ifr);
}
