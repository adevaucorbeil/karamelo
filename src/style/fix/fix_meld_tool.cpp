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
#include <expression_operation.h>
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

  input->parsev(args[++k]);
  c1 = &input->expressions[args[k]];
  input->parsev(args[++k]);
  c2 = &input->expressions[args[k]];
  input->parsev(args[++k]);
  theta = &input->expressions[args[k]];

  lo = input->parsev(args[++k]).result(mpm);
  hi = input->parsev(args[++k]).result(mpm);
  Rmax = input->parsev(args[++k]).result(mpm);
  RmaxSq = Rmax * Rmax;
}

void FixMeldTool::prepare()
{
  ftot = Vector3d();
}

void FixMeldTool::reduce()
{
  Vector3d ftot_reduced;

  // Reduce ftot:
  MPI_Allreduce(ftot.elements, ftot_reduced.elements, 3, MPI_FLOAT, MPI_SUM,
                universe->uworld);

  input->parsev(id + "_x", ftot_reduced[0]);
  input->parsev(id + "_y", ftot_reduced[1]);
  input->parsev(id + "_z", ftot_reduced[2]);
}

void FixMeldTool::initial_integrate(Solid &solid)
{
  // cout << "In FixMeldTool::initial_integrate()\n";

  // Go through all the particles in the group and set b to the right value:

  c1->evaluate(solid);
  c2->evaluate(solid);
  theta->evaluate(solid);

  Kokkos::View<float **> c1_ = c1->registers;
  Kokkos::View<float **> c2_ = c2->registers;
  Kokkos::View<float **> theta_ = theta->registers;

  int groupbit = this->groupbit;
  int dim = this->dim;
  int axis0 = this->axis0;
  int axis1 = this->axis1;
  float K = this->K;
  float w = this->w;
  float lo = this->lo;
  float hi = this->hi;
  float Rmax = this->Rmax;
  float RmaxSq = this->RmaxSq;
  Kokkos::View<float*> mass = solid.mass;
  Kokkos::View<int*> mask = solid.mask;
  Kokkos::View<Vector3d*> sx = solid.x;
  Kokkos::View<Vector3d*> smbp = solid.mbp;
  Kokkos::View<float*> sdamage = solid.damage;
  float G = solid.mat->G;

  float f0, f1, f2;
  int n;

  Kokkos::parallel_reduce("FixMeldTool::initial_integrate", solid.np_local,
			  KOKKOS_LAMBDA(const int &ip, float &ftot0, float &ftot1, float &ftot2, int &n_)
  {
    if (!mass[ip] || !(mask[ip] & groupbit))
      return;

    const float &c = Kokkos::Experimental::cos(theta_(0, ip));
    const float &s = Kokkos::Experimental::sin(theta_(0, ip));

    Vector3d xprime, xtool;
    Matrix3d R;

    if (dim == 0)
    {
      R = Matrix3d(1, 0, 0,
                   0, c, s,
                   0, -s, c);

      xtool = Vector3d(lo, c1_(0, ip), c2_(0, ip));
    }
    else if (dim == 1)
    {
      R = Matrix3d(c, 0, s,
                   0, 1, 0,
                   -s, 0, c);

      xtool = Vector3d(c1_(0, ip), lo, c2_(0, ip));
    }
    else if (dim == 2)
    {
      R = Matrix3d(1, 0, 0,
                   0, c, s,
                   0, -s, c);

      xtool = Vector3d(c1_(0, ip), c2_(0, ip), lo);
    }

    xprime = sx[ip] - xtool;
    if (xprime[dim]   < 0    || xprime[dim]   > hi - lo ||
        xprime[axis0] > Rmax || xprime[axis0] < -Rmax   ||
        xprime[axis1] > Rmax || xprime[axis1] < -Rmax)
      return;

    const float &p0 = xprime[axis0];
    const float &p1 = xprime[axis1];
    const float &p2 = xprime[dim];

    const float &rSq = p0 * p0 + p1 * p1;

    if (rSq > RmaxSq)
      return;

    xprime = R*xprime;
    const float &pext = Rmax - Kokkos::Experimental::sqrt(rSq);

    Vector3d f;

    if (p0 > w) {
      f[axis0] = w - p0;
    } else if (p0 < -w) {
      f[axis0] = -w - p0;
    }

    if (p1 > w) {
      f[axis1] = w - p1;
    } else if (p1 < -w) {
      f[axis1] = -w - p1;
    }

    if (pext > 0 && pext < f.norm()) {
      f = Vector3d();
      const float &r = Kokkos::Experimental::sqrt(rSq);
      f[axis0] = pext * xprime[axis0] / r;
      f[axis1] = pext * xprime[axis1] / r;
    }

    if (p2 > 0 && p2 < f.norm()) {
      f = Vector3d();
      f[dim] = -p2;
    }

    f = K * G * (1 - sdamage[ip]) * R.transpose() * f;
    smbp[ip] += f;

    ftot0 += f[0];
    ftot1 += f[1];
    ftot2 += f[2];
    n_ += 1;
  }, f0, f1, f2, n);

  ftot += Vector3d(f0, f1, f2);
}


void FixMeldTool::write_restart(ofstream *of) {
  // of->write(reinterpret_cast<const char *>(&dim), sizeof(int));
  // of->write(reinterpret_cast<const char *>(&axis0), sizeof(int));
  // of->write(reinterpret_cast<const char *>(&axis1), sizeof(int));
  // of->write(reinterpret_cast<const char *>(&K), sizeof(float));
  // of->write(reinterpret_cast<const char *>(&w), sizeof(float));
  // of->write(reinterpret_cast<const char *>(&lo), sizeof(float));
  // of->write(reinterpret_cast<const char *>(&hi), sizeof(float));
  // of->write(reinterpret_cast<const char *>(&Rmax), sizeof(float));
  
  // c1.write_to_restart(of);
  // c2.write_to_restart(of);
  // theta.write_to_restart(of);
}

void FixMeldTool::read_restart(ifstream *ifr) {
  // ifr->read(reinterpret_cast<char *>(&dim), sizeof(int));
  // ifr->read(reinterpret_cast<char *>(&axis0), sizeof(int));
  // ifr->read(reinterpret_cast<char *>(&axis1), sizeof(int));
  // ifr->read(reinterpret_cast<char *>(&K), sizeof(float));
  // ifr->read(reinterpret_cast<char *>(&w), sizeof(float));
  // ifr->read(reinterpret_cast<char *>(&lo), sizeof(float));
  // ifr->read(reinterpret_cast<char *>(&hi), sizeof(float));
  // ifr->read(reinterpret_cast<char *>(&Rmax), sizeof(float));
  // RmaxSq = Rmax * Rmax;

  // c1.read_from_restart(ifr);
  // c2.read_from_restart(ifr);
  // theta.read_from_restart(ifr);
}
