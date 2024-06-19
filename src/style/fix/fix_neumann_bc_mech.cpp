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

#include <fix_neumann_bc_mech.h>
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


FixNeumannBCMech::FixNeumannBCMech(MPM *mpm, vector<string> args):
  Fix(mpm, args, INITIAL_INTEGRATE)
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

  if (args.size() < Nargs.find(domain->dimension)->second) {
    error->all(FLERR, "Error: too few arguments for fix_heat_flux_particles.\n" +
                          usage.find(domain->dimension)->second);
  }

  if (group->pon[igroup].compare("particles") !=0 ) {
    error->one(FLERR, "fix_convection_nodes needs to be given a group of particles" +
                          group->pon[igroup] + ", " + args[2] +
                          " is a group of " + group->pon[igroup] + ".\n");
  }
  if (universe->me == 0) {
    cout << "Creating new fix FixNeumannBCMech with ID: " << args[0] << endl;
  }

  id = args[0];

  input->parsev(args[3]);
  t[0] = &input->expressions[args[3]];

  if (domain->dimension >= 2) {
    input->parsev(args[4]);
    t[1] = &input->expressions[args[4]];
  }
  if (domain->dimension >= 3) {
    input->parsev(args[5]);
    t[2] = &input->expressions[args[5]];
  }
}

void FixNeumannBCMech::prepare()
{
  ftot = Vector3d();
}

void FixNeumannBCMech::reduce()
{
  Vector3d ftot_reduced;

  // Reduce ftot:
  MPI_Allreduce(&ftot, &ftot_reduced, 3, MPI_FLOAT, MPI_SUM,
                universe->uworld);

  (*input->vars)[id + "_x"] = Var(id + "_x", ftot_reduced[0]);
  (*input->vars)[id + "_y"] = Var(id + "_y", ftot_reduced[1]);
  (*input->vars)[id + "_z"] = Var(id + "_z", ftot_reduced[2]);
}

void FixNeumannBCMech::initial_integrate(Solid &solid)
{
  // Go through all the particles in the group and set v_update to the right value:
  for (int i = 0; i < domain->dimension; i++)
    t[i]->evaluate(solid);

  int groupbit = this->groupbit;
  int dimension = domain->dimension;
  Kokkos::View<int*> mask = solid.mask;
  Kokkos::View<Vector3d*> mbp = solid.mbp;
  Kokkos::View<float*> vol = solid.vol;

  for (int i = 0; i < domain->dimension; i++)
    if (t[i])
      {
	Kokkos::View<float **> t_i = t[i]->registers;
	Kokkos::parallel_reduce("FixNeumannBCMech::initial_integrate", solid.np_local,
			     KOKKOS_LAMBDA(const int &ip, float &ftot_i)
			     {
			       if (!(mask[ip] & groupbit))
				 return;

			       float Ap;
			       if (dimension == 1)
				 Ap = 1;
			       else if (dimension == 2)
				 Ap = Kokkos::sqrt(vol[ip]);
			       else         
				 Ap = Kokkos::pow(vol[ip], 2/3);

			       const float &f = Ap*t_i(0, ip);

			       mbp[ip][i] += f;
			       ftot_i += f;

			     },ftot[i]);
      }
}

void FixNeumannBCMech::write_restart(ofstream *of) {
  // for (int i = 0; i < domain->dimension; i++)
  //   t[i].write_to_restart(of);
}

void FixNeumannBCMech::read_restart(ifstream *ifr) {
  // for (int i = 0; i < domain->dimension; i++)
  //   t[i].read_from_restart(ifr);
}
