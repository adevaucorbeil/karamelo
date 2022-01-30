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

  int k = 2;

  K = input->parsev(args[++k]).result(mpm);

  xs_x = input->parsev(args[++k]);
  xs_y = input->parsev(args[++k]);
  xs_z = input->parsev(args[++k]);
  nx = input->parsev(args[++k]);
  ny = input->parsev(args[++k]);
  nz = input->parsev(args[++k]);

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
}

void FixImpenetrableSurface::prepare()
{
  xs_x.result(mpm);
  xs_y.result(mpm);
  xs_z.result(mpm);
  nx  .result(mpm);
  ny  .result(mpm);
  nz  .result(mpm);

  ftot = Vector3d();
}

void FixImpenetrableSurface::reduce()
{
  Vector3d ftot_reduced;

  // Reduce ftot:
  MPI_Allreduce(ftot.elements, ftot_reduced.elements, 3, MPI_DOUBLE, MPI_SUM,
                universe->uworld);

  (*input->vars)[id + "_x"] = Var(id + "_x", ftot_reduced[0]);
  (*input->vars)[id + "_y"] = Var(id + "_y", ftot_reduced[1]);
  (*input->vars)[id + "_z"] = Var(id + "_z", ftot_reduced[2]);
}

void FixImpenetrableSurface::initial_integrate(Solid &solid, int ip) {
  // cout << "In FixImpenetrableSurface::initial_integrate()\n";

  // Go through all the particles in the group and set b to the right value:
  Vector3d xs(xs_x.result(mpm, true),
              xs_y.result(mpm, true),
              xs_z.result(mpm, true));
  Vector3d n(nx.result(mpm, true),
             ny.result(mpm, true),
             nz.result(mpm, true)); // Outgoing normal
  n.normalize();

  // cout << "line 1: " << line1[0] << "x + " << line1[1] << "y + " << line1[2]
  // << endl; cout << "line 2: " << line2[0] << "x + " << line2[1] << "y + " <<
  // line2[2] << endl;

  if (!solid.mass.at(ip) || !(solid.mask.at(ip) & groupbit))
    return;

  double p = n.dot(solid.x.at(ip) - xs);
  // if (s->ptag[ip] == 1) {
  //   cout << id << "- Particle " << s->ptag[ip] << "\t";
  //   cout << "p = " << p << "\t";
  //   cout << "xp = [" << s->x[ip][0] << "," << s->x[ip][1] << ","
  //        << s->x[ip][2] << "]\t" << endl;
  //   cout << "xs = [" << xs[0] << "," << xs[1] << "," << xs[2] << "]\t"
  //        << endl;
  //   cout << "n = [" << n[0] << "," << n[1] << "," << n[2] << "]"
  //        << endl;
  // }

  if (p < 0)
    return;

    // cout << "Particle " << s->ptag[ip] << " is inside\n";

  const Vector3d &f = K*solid.mat->G*p*(1 - solid.damage.at(ip))*n;
  solid.mbp.at(ip) += f;
  ftot += f;
}

void FixImpenetrableSurface::write_restart(ofstream *of) {
  of->write(reinterpret_cast<const char *>(&K), sizeof(double));

  xs_x.write_to_restart(of);
  xs_y.write_to_restart(of);
  xs_z.write_to_restart(of);

  nx.write_to_restart(of);
  ny.write_to_restart(of);
  nz.write_to_restart(of);
}

void FixImpenetrableSurface::read_restart(ifstream *ifr) {
  ifr->read(reinterpret_cast<char *>(&K), sizeof(double));

  xs_x.read_from_restart(ifr);
  xs_y.read_from_restart(ifr);
  xs_z.read_from_restart(ifr);

  nx.read_from_restart(ifr);
  ny.read_from_restart(ifr);
  nz.read_from_restart(ifr);
}
