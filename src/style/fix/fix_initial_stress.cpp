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

#include <fix_initial_stress.h>
#include <domain.h>
#include <error.h>
#include <group.h>
#include <input.h>
#include <solid.h>
#include <universe.h>
#include <update.h>
#include <expression_operation.h>


using namespace std;
using namespace FixConst;


FixInitialStress::FixInitialStress(MPM *mpm, vector<string> args):
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

  if (args.size() < Nargs) {
    error->all(FLERR, "Error: not enough arguments.\n" + usage);
  }

  if (args.size() > Nargs) {
    error->all(FLERR, "Error: too many arguments.\n" + usage);
  }

  if (group->pon[igroup].compare("particles") !=0 && group->pon[igroup].compare("all") !=0) {
    error->all(FLERR, "fix_initial_stress needs to be given a group of particles" + group->pon[igroup] + ", " + args[2] + " is a group of " + group->pon[igroup] + ".\n");
  }
  if (universe->me == 0) {
    cout << "Creating new fix FixInitialStress with ID: " << args[0] << endl;
  }
  id = args[0];

  for (int i = 0; i < 6; i++) {
    if (args[i + 3] != "NULL") {
      input->parsev(args[i + 3]);
      s_value[i] = &input->expressions[args[i + 3]];
    } else {
      s_value[i] = nullptr;
    }
  }
}

void FixInitialStress::prepare()
{
  // for (Var &s_value: s_value)
  //   s_value.result(mpm);
}

void FixInitialStress::initial_integrate(Solid &solid)
{
  // cout << "In FixInitialStress::initial_integrate()" << endl;
  // Go through all the nodes in the group and set v to the right value:

  for (int i = 0; i < 6; i++)
    if (s_value[i])
      s_value[i]->evaluate(solid);


  int groupbit = this->groupbit;
  double ntimestep = update->ntimestep;
  Kokkos::View<int*> mask = solid.mask;
  Kokkos::View<Matrix3d*> sigma = solid.sigma;

  for (int i = 0; i < 6; i++)
    if (s_value[i])
      {
	Kokkos::View<double **> s_value_i = s_value[i]->registers;
	Kokkos::parallel_for("FixInitialStress::initial_integrate", solid.np_local,
			     KOKKOS_LAMBDA(const int &ip)
			     {
			       if (ntimestep != 1 || !(mask[ip] & groupbit))
				 return;

			       if      (i<=2) sigma[ip](i,i) = s_value_i(0, ip);
			       else if (i==3) sigma[ip](1,2) = sigma[ip](2,1) = s_value_i(0, ip);
			       else if (i==4) sigma[ip](0,2) = sigma[ip](2,0) = s_value_i(0, ip);
			       else if (i==5) sigma[ip](0,1) = sigma[ip](1,0) = s_value_i(0, ip);
			     });
      }

  
  if (update->method_type == "tlmpm" || update->method_type == "tlcpdi") {
    Kokkos::View<double*> vol0 = solid.vol0;
    Kokkos::View<Matrix3d*> vol0PK1 = solid.vol0PK1;

    Kokkos::parallel_for("FixInitialStress::initial_integrate_tlmpm", solid.np_local,
			 KOKKOS_LAMBDA(const int &ip)
			 {
			   solid.vol0PK1[ip] = vol0[ip]*sigma[ip];
			 });
  }
}
