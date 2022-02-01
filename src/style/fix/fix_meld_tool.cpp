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

#include <fix_meld_tool.h>
#include <domain.h>
#include <error.h>
#include <group.h>
#include <input.h>
#include <solid.h>
#include <universe.h>
#include <update.h>


using namespace std;
using namespace FixConst;


FixMeldTool::FixMeldTool(MPM *mpm, vector<string> args)
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
    error->all(FLERR, "fix_meldtool needs to be given a group of nodes" +
                          group->pon[igroup] + ", " + args[2] +
                          " is a group of " + group->pon[igroup] + ".\n");
  }
  if (universe->me == 0) {
    cout << "Creating new fix FixMeldTool with ID: " << args[0] << endl;
  }
  id = args[0];

  int k = 2;
  
  K = input->parsev(args[++k]).result(mpm);
  k++;
  if (args[k] == "x") {
    dim = X;
    axis0 = Y;
    axis1 = Z;
  } else if (args[k] == "y") {
    dim = Y;
    axis0 = X;
    axis1 = Z;
  } else if (args[k] == "z") {
    dim = Z;
    axis0 = X;
    axis1 = Y;
  } else {
    error->all(FLERR, "Unknown dim: " + args[k] + "! Options are : x, y, and z.\n");
  }

  w = input->parsev(args[++k]).result(mpm);
  c1 = input->parsev(args[++k]);
  c2 = input->parsev(args[++k]);
  theta = input->parsev(args[++k]);
  lo = input->parsev(args[++k]).result(mpm);
  hi = input->parsev(args[++k]).result(mpm);
  Rmax = input->parsev(args[++k]).result(mpm);
  RmaxSq = Rmax * Rmax;
}

void FixMeldTool::prepare()
{
  theta.result(mpm);
  c1.result(mpm);
  c2.result(mpm);

  ftot = Vector3d();
}

void FixMeldTool::reduce()
{
  Vector3d ftot_reduced;

  // Reduce ftot:
  MPI_Allreduce(ftot.elements, ftot_reduced.elements, 3, MPI_DOUBLE, MPI_SUM,
                universe->uworld);

  (*input->vars)[id + "_x"] = Var(id + "_x", ftot_reduced[0]);
  (*input->vars)[id + "_y"] = Var(id + "_y", ftot_reduced[1]);
  (*input->vars)[id + "_z"] = Var(id + "_z", ftot_reduced[2]);
}

void FixMeldTool::initial_integrate(Solid &solid, int ip)
{
  // cout << "In FixMeldTool::initial_integrate()\n";

  // Go through all the particles in the group and set b to the right value:
  if (!solid.mass.at(ip) || !(solid.mask.at(ip) & groupbit))
    return;

  double theta_ = theta.result(mpm, true);
  double c = cos(theta_);
  double s = sin(theta_);

  double c1_ = c1.result(mpm, true);
  double c2_ = c2.result(mpm, true);

  Vector3d xprime;
  Matrix3d R;

  if (dim == X)
  {
    R = Matrix3d(1, 0, 0,
                 0, c, s,
                 0, -s, c);

    xprime = Vector3d(lo, c1_, c2_);
  }
  else if (dim == Y)
  {
    R = Matrix3d(c, 0, s,
                 0, 1, 0,
                 -s, 0, c);

    xprime = Vector3d(c1_, lo, c2_);
  }
  else if (dim == Z)
  {
    R = Matrix3d(1, 0, 0,
                 0, c, s,
                 0, -s, c);

    xprime = Vector3d(c1_, c2_, lo);
  }

  xprime = solid.x.at(ip) - xprime;
  // if (update->ntimestep > 89835 && (solid.ptag.at(ip)==12 || solid.ptag.at(ip)==21)) {
  //   cout << "Check Particle " << solid.ptag.at(ip) << "\txprime=[" << xprime[0] << "," << xprime[1] << "," << xprime[2] << "]\n";
  //   cout << "R=\n" << R << endl;
  // }
  if (xprime(dim)   < 0    || xprime(dim)   > hi - lo ||
      xprime(axis0) > Rmax || xprime(axis0) < -Rmax   ||
      xprime(axis1) > Rmax || xprime(axis1) < -Rmax)
    return;

  // if (solid.ptag.at(ip)==12 || solid.ptag.at(ip)==21) {
  //   cout << "Particle " << solid.ptag.at(ip) << " in 1\n";
  // }

  double rSq = xprime(axis0) * xprime(axis0) + xprime(axis1) * xprime(axis1);

  if (rSq > RmaxSq)
    return;

  // if (solid.ptag.at(ip)==12 || solid.ptag.at(ip)==21) {
  //     cout << "Particle " << solid.ptag.at(ip) << " in 2\n";
  // }
  xprime = R*xprime;
  double p0 = xprime[axis0];
  double p1 = xprime[axis1];
  double p2 = xprime[dim];
  double pext = Rmax - sqrt(p0*p0 + p1*p1);
  double p;

  Vector3d f;

  if (p0 > w)
  {
    p = p0 - w;
    f[axis0] = -p;
  }
  else if (p0 < -w)
  {
    p = -w - p0;
    f[axis0] = p;
  }

  if (p1 > w)
  {
    p = p1 - w;
    f[axis1] = -p;
  }
  else if (p1 < -w)
  {
    p = -p1 - w;
    f[axis1] = p;
  }

  if (pext > 0 && pext < f.norm())
  {
    f = Vector3d();
    double r = sqrt(rSq);
    f[axis0] = pext*xprime[axis0]/r;
    f[axis1] = pext*xprime[axis1]/r;
  }

  if (p2 > 0 || p2 < f.norm())
  {
    f = Vector3d();
    f[dim] = -p2;
  }

  f = K*solid.mat->G*(1 - solid.damage.at(ip))*R.transpose()*f;
  solid.mbp.at(ip) += f;
  // if (solid.ptag.at(ip)==12 || solid.ptag.at(ip)==21) {
  //     Vector3d dx = solid.x.at(ip) - c;
  //     cout << "Particle " << solid.ptag.at(ip) << " f=[" << f[0] << "," << f[1] << "," << f[2] << "]\tw=" << w << " p0=" << p0 << " p1=" << p1 << " p2=" << p2 << "\txprime=[" << xprime[0] << "," << xprime[1] << "," << xprime[2] << "]\tdx=[" << dx(0) << "," << dx(1) << "," << dx(2) << "]\n";
  //     cout << "R=\n" << R << endl;
  // }
  ftot += f;
  // if (f[dim] != 0) {
  //     cout << "particle " << solid.ptag.at(ip) << " force:" << f[0] << ", " << f[1] << ", " << f[2] << endl;
  // }
}


void FixMeldTool::write_restart(ofstream *of) {
  of->write(reinterpret_cast<const char *>(&dim), sizeof(int));
  of->write(reinterpret_cast<const char *>(&axis0), sizeof(int));
  of->write(reinterpret_cast<const char *>(&axis1), sizeof(int));
  of->write(reinterpret_cast<const char *>(&K), sizeof(double));
  of->write(reinterpret_cast<const char *>(&w), sizeof(double));
  of->write(reinterpret_cast<const char *>(&lo), sizeof(double));
  of->write(reinterpret_cast<const char *>(&hi), sizeof(double));
  of->write(reinterpret_cast<const char *>(&Rmax), sizeof(double));
  
  c1.write_to_restart(of);
  c2.write_to_restart(of);
  theta.write_to_restart(of);
}

void FixMeldTool::read_restart(ifstream *ifr) {
  ifr->read(reinterpret_cast<char *>(&dim), sizeof(int));
  ifr->read(reinterpret_cast<char *>(&axis0), sizeof(int));
  ifr->read(reinterpret_cast<char *>(&axis1), sizeof(int));
  ifr->read(reinterpret_cast<char *>(&K), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&w), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&lo), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&hi), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&Rmax), sizeof(double));
  RmaxSq = Rmax * Rmax;

  c1.read_from_restart(ifr);
  c2.read_from_restart(ifr);
  theta.read_from_restart(ifr);
}
