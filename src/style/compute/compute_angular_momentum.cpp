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

#include <compute_angular_momentum.h>
#include <domain.h>
#include <error.h>
#include <group.h>
#include <input.h>
#include <math_special.h>
#include <output.h>
#include <universe.h>
#include <update.h>
#include <matrix.h>
#include <iostream>
#include <string>
#include <vector>
#include <expression_operation.h>

using namespace std;
using namespace MathSpecial;


ComputeAngularMomentum::ComputeAngularMomentum(MPM *mpm, vector<string> args)
    : Compute(mpm, args) {

  
  if (args.size() < Nargs) {
    error->all(FLERR, "Error: not enough arguments.\n" + usage);
  }
  
  if (group->pon[igroup] != "particles" &&
      group->pon[igroup] != "all") {
    error->all(FLERR,
               "compute_angular_momentum needs to be given a group of particles" +
                   group->pon[igroup] + ", " + args[2] + " is a group of " +
                   group->pon[igroup] + ".\n");
  }

  
  if (universe->me == 0) 
    cout << "Creating new compute ComputeAngularMomentum with ID: " << args[0] << endl;
  id = args[0];

  int k = 2;
  
  x0[0] = input->parsev(args[++k]).result(mpm);
  x0[1] = input->parsev(args[++k]).result(mpm);
  x0[2] = input->parsev(args[++k]).result(mpm);

  input->parsev(id, 0);

  t = update->ntimestep;
  J = Vector3d();
  input->parsev(id + "_x", J(0));
  input->parsev(id + "_y", J(1));
  input->parsev(id + "_z", J(2));
}

ComputeAngularMomentum::~ComputeAngularMomentum() {}

void ComputeAngularMomentum::compute_value(Solid &solid) {
  // double Ek_tmp, Ek_reduced = 0;
  double Jx, Jy, Jz;
  Jx = Jy = Jz = 0;

  if (t != update->ntimestep) {
    t = update->ntimestep;
    J = Vector3d();
  }

  Kokkos::View<Vector3d*> x = solid.x;
  Kokkos::View<Vector3d*> v = solid.v;
  Kokkos::View<double*> mass = solid.mass;
  Kokkos::View<int*> mask = solid.mask;

  int groupbit = this->groupbit;
  Vector3d x0 = this->x0;

  if (update->ntimestep == output->next ||
      update->ntimestep == update->nsteps){
    Kokkos::parallel_reduce("ComputeAngularMomentum::compute_value", solid.np_local,
			    KOKKOS_LAMBDA(const int &ip, double &lJx, double &lJy, double &lJz)
			    {
			      if (mask[ip] & groupbit) {
				const Vector3d &Jp = -mass[ip] * v[ip].cross(x[ip] - x0);
				lJx += Jp(0);
				lJy += Jp(1);
				lJz += Jp(2);
			      }
			    }, Jx, Jy, Jz);
    J += Vector3d(Jx, Jy, Jz);
  }

  // Reduce J:

  Vector3d J_reduced;
  MPI_Allreduce(J.elements, J_reduced.elements, 3, MPI_DOUBLE, MPI_SUM, universe->uworld);
  input->parsev(id + "_x", J_reduced(0));
  input->parsev(id + "_y", J_reduced(1));
  input->parsev(id + "_z", J_reduced(2));
}
