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

#include <fix_indent_minimize_penetration.h>
#include <domain.h>
#include <error.h>
#include <group.h>
#include <input.h>
#include <solid.h>
#include <universe.h>
#include <update.h>


using namespace std;
using namespace FixConst;


#define four_thirds 1.333333333

FixIndentMinimizePenetration::FixIndentMinimizePenetration(MPM *mpm, vector<string> args)
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

    R = mu = 0;
    return;
  }

  if (args.size() < Nargs) {
    error->all(FLERR, "Error: not enough arguments.\n" + usage);
  }
  if (args.size() > Nargs) {
    error->all(FLERR, "Error: too many arguments.\n" + usage);
  }

  if (group->pon[igroup] != "particles" &&
      group->pon[igroup] != "all") {
    error->all(FLERR, "fix_indent_hertz needs to be given a group of particles" +
                          group->pon[igroup] + ", " + args[2] +
                          " is a group of " + group->pon[igroup] + ".\n");
  }
  if (universe->me == 0) {
    cout << "Creating new fix FixIndentMinimizePenetration with ID: " << args[0] << endl;
  }
  id = args[0];

  type = args[3];
  if (args[3] == "sphere") {
    type = "sphere";
  } else {
    error->all(FLERR, "Error indent type " + args[3] +
                          " unknown. Only type sphere is supported.\n");
  }

  R = input->parsev(args[4]);
  xvalue = input->parsev(args[5]);
  yvalue = input->parsev(args[6]);
  zvalue = input->parsev(args[7]);
  vxvalue = input->parsev(args[8]);
  vyvalue = input->parsev(args[9]);
  vzvalue = input->parsev(args[10]);
  mu = input->parsev(args[11]);
}

void FixIndentMinimizePenetration::prepare()
{
  xvalue.result(mpm);
  yvalue.result(mpm);
  zvalue.result(mpm);
  vxvalue.result(mpm);
  vyvalue.result(mpm);
  vzvalue.result(mpm);
  
  int solid = group->solid[igroup];

  cellsizeSq = 0;
  if (solid == -1)
    for (const Solid *solid: domain->solids)
      cellsizeSq += solid->grid->cellsize*solid->grid->cellsize;
  else
    cellsizeSq += domain->solids.at(solid)->grid->cellsize*domain->solids.at(solid)->grid->cellsize;

  A = 0;
  ftot = Vector3d();
}

void FixIndentMinimizePenetration::reduce()
{
  double A_reduced;
  Vector3d ftot_reduced;

  // Reduce ftot:
  MPI_Allreduce(&A, &A_reduced, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  MPI_Allreduce(ftot.elements, ftot_reduced.elements, 3, MPI_DOUBLE, MPI_SUM,
                universe->uworld);

  (*input->vars)[id + "_s"] = Var(id + "_s", A_reduced);
  (*input->vars)[id + "_x"] = Var(id + "_x", ftot_reduced[0]);
  (*input->vars)[id + "_y"] = Var(id + "_y", ftot_reduced[1]);
  (*input->vars)[id + "_z"] = Var(id + "_z", ftot_reduced[2]);
}

void FixIndentMinimizePenetration::initial_integrate(Solid &solid, int ip)
{
  // Go through all the particles in the group and set b to the right value:
  if (!solid.mass.at(ip) || !(solid.mask.at(ip) & groupbit))
    return;

  Vector3d xs(xvalue .result(mpm, true),
              yvalue .result(mpm, true),
              zvalue .result(mpm, true));
  Vector3d vs(vxvalue.result(mpm, true),
              vyvalue.result(mpm, true),
              vzvalue.result(mpm, true));
  
  // Gross screening:
  Vector3d xsp = solid.x.at(ip) - xs;

  double Rs = 0;
  if (domain->dimension == 2)
  {
	if (domain->axisymmetric)
	  Rs = 0.5*sqrt(solid.vol0.at(ip)/solid.x0.at(ip)[0]);
	else
	  Rs = 0.5*sqrt(solid.vol0.at(ip));
	Rs += R;

	if (xsp[0] >=  Rs || xsp[1] >=  Rs ||
        xsp[0] <= -Rs || xsp[1] <= -Rs)
      return;
  }
  else if (domain->dimension == 3)
  {
	Rs = R + 0.5*cbrt(solid.vol0.at(ip));

	if (xsp[0] >=  Rs || xsp[1] >=  Rs || xsp[2] >=  Rs ||
        xsp[0] <= -Rs || xsp[1] <= -Rs || xsp[2] <= -Rs)
      return;
  }

  // Finer screening:
  double r = xsp.norm();

  if (r >= Rs)
    return;

  // penetration
  double p = Rs - r;

  if (p <= 0)
    return;

  xsp /= r;
  double fmag = solid.mass.at(ip)*p/update->dt/update->dt;
  Vector3d f = fmag*xsp;

  if (mu)
  {
	const Vector3d &vps = vs - solid.v.at(ip);
	Vector3d vt = vps - vps.dot(xsp)*xsp;
	double vtnorm = vt.norm();

	if (vtnorm)
    {
	  vt /= vtnorm;
	  f += MIN(solid.mass.at(ip)*vtnorm/update->dt, mu*fmag)*vt;
	}
  }

  A += cellsizeSq*solid.F.at(ip)(0, 0)*solid.F.at(ip)(2, 2);
  solid.mbp.at(ip) += f;
  ftot += f;
}

void FixIndentMinimizePenetration::write_restart(ofstream *of) {
  of->write(reinterpret_cast<const char *>(&R), sizeof(double));
  of->write(reinterpret_cast<const char *>(&mu), sizeof(double));
  xvalue.write_to_restart(of);
  yvalue.write_to_restart(of);
  zvalue.write_to_restart(of);

  vxvalue.write_to_restart(of);
  vyvalue.write_to_restart(of);
  vzvalue.write_to_restart(of);
}

void FixIndentMinimizePenetration::read_restart(ifstream *ifr) {
  ifr->read(reinterpret_cast<char *>(&R), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&mu), sizeof(double));
  xvalue.read_from_restart(ifr);
  yvalue.read_from_restart(ifr);
  zvalue.read_from_restart(ifr);

  vxvalue.read_from_restart(ifr);
  vyvalue.read_from_restart(ifr);
  vzvalue.read_from_restart(ifr);
  type = "sphere";
}
