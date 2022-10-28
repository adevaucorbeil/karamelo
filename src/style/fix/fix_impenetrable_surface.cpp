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

#include <fix_impenetrable_surface.h>
#include <domain.h>
#include <error.h>
#include <expression_operation.h>
#include <group.h>
#include <input.h>
#include <solid.h>
#include <universe.h>
#include <update.h>


using namespace std;
using namespace FixConst;


FixImpenetrableSurface::FixImpenetrableSurface(MPM *mpm, vector<string> args)
    : Fix(mpm, args, INITIAL_INTEGRATE) {

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

  if (group->pon[igroup] != "particles" &&
      group->pon[igroup] != "all") {
    error->all(FLERR,
               "fix_impenetrablesurface needs to be given a group of nodes" +
                   group->pon[igroup] + ", " + args[2] + " is a group of " +
                   group->pon[igroup] + ".\n");
  }
  if (universe->me == 0) {
    cout << "Creating new fix FixImpenetrableSurface with ID: " << args[0]
     << endl;
  }
  id = args[0];

  int k = 2;

  K = input->parsev(args[++k]).result(mpm);

  for (int dim = 0; dim < 3; dim++) {
    k++;
    input->parsev(args[k]);
    xs[dim] = &input->expressions[args[k]];
  }
  for (int dim = 0; dim < 3; dim++) {
    k++;
    input->parsev(args[k]);
    normal[dim] = &input->expressions[args[k]];
    cout << "normal[" << dim << "]=" << args[k] << endl;
  }
}

void FixImpenetrableSurface::prepare()
{
  ftot = Vector3d();
}

void FixImpenetrableSurface::reduce()
{
  Vector3d ftot_reduced;

  // Reduce ftot:
  MPI_Allreduce(ftot.elements, ftot_reduced.elements, 3, MPI_FLOAT, MPI_SUM,
                universe->uworld);

  // (*input->vars)[id + "_x"] = Var(id + "_x", ftot_reduced[0]);
  // (*input->vars)[id + "_y"] = Var(id + "_y", ftot_reduced[1]);
  // (*input->vars)[id + "_z"] = Var(id + "_z", ftot_reduced[2]);
  input->parsev(id + "_x", ftot_reduced[0]);
  input->parsev(id + "_y", ftot_reduced[1]);
  input->parsev(id + "_z", ftot_reduced[2]);
}

void FixImpenetrableSurface::initial_integrate(Solid &solid) {
  // cout << "In FixImpenetrableSurface::initial_integrate()\n";

  // Go through all the particles in the group and set b to the right value:
  for (int i = 0; i < 3; i++) {
    if (xs[i])
      xs[i]->evaluate(solid);
    if (normal[i])
      normal[i]->evaluate(solid);
  }

  int K = this->K;
  int groupbit = this->groupbit;
  Kokkos::View<int*> mask = solid.mask;
  Kokkos::View<float*> smass = solid.mass;
  Kokkos::View<Vector3d*> sx = solid.x;
  Kokkos::View<Vector3d*> smbp = solid.mbp;
  Kokkos::View<float*> sdamage = solid.damage;

  float G = solid.mat->G;

  Kokkos::View<float **> xs0_i = xs[0]->registers;
  Kokkos::View<float **> xs1_i = xs[1]->registers;
  Kokkos::View<float **> xs2_i = xs[2]->registers;

  Kokkos::View<float **> normal0_i = normal[0]->registers;
  Kokkos::View<float **> normal1_i = normal[1]->registers;
  Kokkos::View<float **> normal2_i = normal[2]->registers;

  float ftot_0, ftot_1, ftot_2;

  Kokkos::parallel_reduce("FixImpenetrableSurface::initial_integrate", solid.np_local,
  KOKKOS_LAMBDA(const int &ip, float &lftot_0, float &lftot_1, float &lftot_2)
  {
    if (!smass[ip] || !(mask[ip] & groupbit))
      return;

    Vector3d xs_ = Vector3d(xs0_i(0, ip), xs1_i(0, ip), xs2_i(0, ip));
    Vector3d normal_ = Vector3d(normal0_i(0, ip), normal1_i(0, ip), normal2_i(0, ip));
    normal_.normalize();

    float p = normal_.dot(xs_ - sx[ip]);

    if (p < 0)
      return;

    const Vector3d &f = K*G*p*(1 - sdamage[ip])*normal_;

    smbp[ip] += f;
    lftot_0 += f[0];
    lftot_1 += f[1];
    lftot_2 += f[2];
  }, ftot_0, ftot_1, ftot_2);

  ftot = Vector3d(ftot_0, ftot_1, ftot_2);
}

void FixImpenetrableSurface::write_restart(ofstream *of) {
  // of->write(reinterpret_cast<const char *>(&K), sizeof(float));

  // xs_x.write_to_restart(of);
  // xs_y.write_to_restart(of);
  // xs_z.write_to_restart(of);

  // nx.write_to_restart(of);
  // ny.write_to_restart(of);
  // nz.write_to_restart(of);
}

void FixImpenetrableSurface::read_restart(ifstream *ifr) {
  // ifr->read(reinterpret_cast<char *>(&K), sizeof(float));

  // xs_x.read_from_restart(ifr);
  // xs_y.read_from_restart(ifr);
  // xs_z.read_from_restart(ifr);

  // nx.read_from_restart(ifr);
  // ny.read_from_restart(ifr);
  // nz.read_from_restart(ifr);
}
