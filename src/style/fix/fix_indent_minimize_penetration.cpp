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
#include <expression_operation.h>
#include <group.h>
#include <input.h>
#include <method.h>
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

  R = input->parsev(args[4]).result(mpm);

  input->parsev(args[5]);
  xvalue = &input->expressions[args[5]];
  input->parsev(args[6]);
  yvalue = &input->expressions[args[6]];
  input->parsev(args[7]);
  zvalue = &input->expressions[args[7]];

  input->parsev(args[8]);
  vxvalue = &input->expressions[args[8]];
  input->parsev(args[9]);
  vyvalue = &input->expressions[args[9]];
  input->parsev(args[10]);
  vzvalue = &input->expressions[args[10]];

  mu = input->parsev(args[11]).result(mpm);
}

void FixIndentMinimizePenetration::prepare()
{
  A = 0;
  ftot = Vector3d();
}

void FixIndentMinimizePenetration::reduce()
{
  float A_reduced;
  Vector3d ftot_reduced;

  // Reduce ftot:
  MPI_Allreduce(&A, &A_reduced, 1, MPI_FLOAT, MPI_SUM, universe->uworld);
  MPI_Allreduce(ftot.elements, ftot_reduced.elements, 3, MPI_FLOAT, MPI_SUM,
                universe->uworld);

  input->parsev(id + "_s", A_reduced);
  input->parsev(id + "_x", ftot_reduced[0]);
  input->parsev(id + "_y", ftot_reduced[1]);
  input->parsev(id + "_z", ftot_reduced[2]);
}

void FixIndentMinimizePenetration::initial_integrate(Solid &solid)
{
  xvalue->evaluate(solid);
  yvalue->evaluate(solid);
  zvalue->evaluate(solid);

  vxvalue->evaluate(solid);
  vyvalue->evaluate(solid);
  vzvalue->evaluate(solid);

  Kokkos::View<float **> xvalue_ = xvalue->registers;
  Kokkos::View<float **> yvalue_ = yvalue->registers;
  Kokkos::View<float **> zvalue_ = zvalue->registers;

  Kokkos::View<float **> vxvalue_ = vxvalue->registers;
  Kokkos::View<float **> vyvalue_ = vyvalue->registers;
  Kokkos::View<float **> vzvalue_ = vzvalue->registers;

  // Go through all the particles in the group and set b to the right value:

  int groupbit = this->groupbit;
  float R      = this->R;
  float mu     = this->mu;
  int dimension = domain->dimension;
  bool axisymmetric = domain->axisymmetric;

  Kokkos::View<float*>    mass = solid.mass;
  Kokkos::View<int*>      mask = solid.mask;
  Kokkos::View<Vector3d*>  sx  = solid.x;
  Kokkos::View<Vector3d*>  sx0 = solid.x0;
  Kokkos::View<Vector3d*>   sv = solid.v;
  Kokkos::View<float*>    svol = update->method->is_TL ? solid.vol0 : solid.vol;
  Kokkos::View<Vector3d*> smbp = solid.mbp;
  Kokkos::View<Matrix3d*> sF = solid.F;

  const float &cellsizeSq = solid.grid->cellsize * solid.grid->cellsize;
  float dt = update->dt;
  float dtSq = dt*dt;
  
  float f0 = 0, f1 = 0, f2 = 0;

  Kokkos::parallel_reduce("FixIndentMinimizePenetration::initial_integrate", solid.np_local,
  KOKKOS_LAMBDA(const int &ip, float &ftot0, float &ftot1, float &ftot2, float &A_)
  {
    if (!mass[ip] || !(mask[ip] & groupbit))
      return;

    const Vector3d xs( xvalue_(0, ip),  yvalue_(0, ip),  zvalue_(0, ip));
    const Vector3d vs(vxvalue_(0, ip), vyvalue_(0, ip), vzvalue_(0, ip));

    // Gross screening:
    Vector3d xsp = sx[ip] - xs;

    float Rs = 0;
    if (dimension == 2) {
      if (axisymmetric)
        Rs = 0.5 * Kokkos::Experimental::sqrt(svol[ip] / sx0[ip][0]);
      else
        Rs = 0.5 * Kokkos::Experimental::sqrt(svol[ip]);
      Rs += R;

      if (xsp[0] >= Rs || xsp[1] >= Rs || xsp[0] <= -Rs || xsp[1] <= -Rs)
        return;
    } else if (dimension == 3) {
      Rs = R + 0.5 * Kokkos::Experimental::cbrt(svol[ip]);

      if (xsp[0] >= Rs || xsp[1] >= Rs || xsp[2] >= Rs || xsp[0] <= -Rs ||
          xsp[1] <= -Rs || xsp[2] <= -Rs)
        return;
    }

    // Finer screening:
    const float &r = xsp.norm();

    if (r >= Rs)
      return;

    // penetration
    const float &p = Rs - r;

    if (p <= 0)
      return;

    xsp /= r;
    const float &fmag = mass[ip]*p/dtSq;
    Vector3d f = fmag*xsp;

    if (mu) {
      const Vector3d &vps = vs - sv[ip];
      Vector3d vt = vps - vps.dot(xsp) * xsp;
      const float &vtnorm = vt.norm();

      if (vtnorm) {
        vt /= vtnorm;
        f += MIN(mass[ip] * vtnorm / dt, mu * fmag) * vt;
      }
    }

    A_ += cellsizeSq*sF[ip](0, 0)*sF[ip](2, 2);
    smbp[ip] += f;

    ftot0 += f[0];
    ftot1 += f[1];
    ftot2 += f[2];
  }, f0, f1, f2, A);

  ftot += Vector3d(f0, f1, f2);
}

void FixIndentMinimizePenetration::write_restart(ofstream *of) {
  // of->write(reinterpret_cast<const char *>(&R), sizeof(float));
  // of->write(reinterpret_cast<const char *>(&mu), sizeof(float));
  // xvalue.write_to_restart(of);
  // yvalue.write_to_restart(of);
  // zvalue.write_to_restart(of);

  // vxvalue.write_to_restart(of);
  // vyvalue.write_to_restart(of);
  // vzvalue.write_to_restart(of);
}

void FixIndentMinimizePenetration::read_restart(ifstream *ifr) {
  // ifr->read(reinterpret_cast<char *>(&R), sizeof(float));
  // ifr->read(reinterpret_cast<char *>(&mu), sizeof(float));
  // xvalue.read_from_restart(ifr);
  // yvalue.read_from_restart(ifr);
  // zvalue.read_from_restart(ifr);

  // vxvalue.read_from_restart(ifr);
  // vyvalue.read_from_restart(ifr);
  // vzvalue.read_from_restart(ifr);
  // type = "sphere";
}
