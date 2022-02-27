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

#include <populate_cylindrical_coordinates.h>
#include <domain.h>
#include <error.h>
#include <group.h>
#include <input.h>
#include <method.h>
#include <solid.h>
#include <universe.h>
#include <update.h>
#include <var.h>
#include <iostream>

using namespace std;

PopulateCylindricalCoordinates::PopulateCylindricalCoordinates(MPM *mpm) : Pointers(mpm) {}

Var PopulateCylindricalCoordinates::command(vector<string> args) {
  // cout << "In PopulateCylindricalCoordinates::command()" << endl;

  if (args.size() < Nargs.find(domain->dimension)->second) {
      error->all(FLERR,
		 "Error: not enough arguments.\n" + usage.find(domain->dimension)->second);
    }

  int isolid = domain->find_solid(args[0]);

  if (isolid < 0) {
    if (args[0].compare("all")!=0) {
      error->all(FLERR, "Error: solid " + args[0] + " unknown.\n");
    }
  }

  int iregion = domain->find_region(args[1]);

  if (iregion < 0)
    {
      error->all(FLERR, "Error: region " + args[1] + " unknown.\n");
    }

  double c1, c2, R, lo, hi, T0;
  int Ndx, Ndr, Ndtheta;
  char axis;
  int ii, jj, kk;

  if (domain->dimension == 3) {
    if (args[2] == "x") {
      axis = 'x';
      ii = 1;
      jj = 2;
      kk = 0;      
    } else if (args[2] == "y") {
      axis = 'y';
      ii = 2;
      jj = 0;
      kk = 1;
    } else if (args[2] == "z") {
      axis = 'z';
      ii = 0;
      jj = 1;
      kk = 2;
    } else {
      error->all(FLERR, "Error: populate_cylindrical_coordinates, cylinder axis not understood, expect x, y, or z, received " + args[2]+".\n");
    }

    c1 = input->parsev(args[3]);
    c2 = input->parsev(args[4]);
    R = input->parsev(args[5]);
    lo = input->parsev(args[6]);
    hi = input->parsev(args[7]);
    Ndx = (int) input->parsev(args[8]);
    Ndr = (int) input->parsev(args[9]);
    Ndtheta = (int) input->parsev(args[10]);
    T0 = input->parsev(args[11]);
  } else {
    axis = 'z';

    c1 = input->parsev(args[2]);
    c2 = input->parsev(args[3]);
    R = input->parsev(args[4]);
    lo = hi = 0;
    Ndx = 1;
    Ndr = (int) input->parsev(args[5]);
    Ndtheta = (int) input->parsev(args[6]);
    T0 = input->parsev(args[7]);
    ii = 0;
    jj = 1;
    kk = 2;
  }

  if (universe->me == 0) {
    cout << "axis, c1, c2, R = " << axis << "\t" << c1 << "\t" << c2 << "\t"
         << R << endl;
  }

  int N = Ndx * Ndr * Ndtheta;

  Solid *s = domain->solids[isolid];

  int np_local_reduced;
  domain->np_total -= s->np;
  s->grow(N);

  int l = 0;

  double Ndr_ = R/Ndr;
  double Ndtheta_ = (2*M_PI)/Ndtheta;
  double Ndx_ = (hi - lo)/Ndx;
  if (Ndx_ == 0)
    Ndx_++;
  double vol_ = Ndr_ * Ndtheta_ * Ndx_;

  for (int i = 0; i < Ndx; i++) {
    for (int j = 0; j < Ndr; j++) {
      for (int k = 0; k < Ndtheta; k++) {
        if (l >= N) {
          error->all(FLERR, "Error in Solid::populate(), exceeding the "
                            "allocated number of particles.\n");
        }

        s->x0[l][ii] = s->x[l][ii] = c1 + j * Ndr_ * cos(k * Ndtheta_);
        s->x0[l][jj] = s->x[l][jj] = c2 + j * Ndr_ * sin(k * Ndtheta_);
        s->x0[l][kk] = s->x[l][kk] = lo + i * Ndx_;

        // Check if the particle is inside the region:
        if (domain->inside_subdomain(s->x0[l][0], s->x0[l][1], s->x0[l][2]) &&
            domain->regions[iregion]->inside(s->x0[l][0], s->x0[l][1], s->x0[l][2]) ==
                1) {
          s->a[l] = Vector3d();
          s->v[l] = Vector3d();
          s->f[l] = Vector3d();
          s->mbp[l] = Vector3d();
          s->v_update[l] = Vector3d();
          s->rho0[l] = s->rho[l] = s->mat->rho0;
	  if (domain->axisymmetric == true) {
	    if (j==0) {
	      s->vol0[l] = s->vol[l] = s->x0[l][0] * 2 * M_PI * Ndr_ * Ndx_;
	    } else {
	      s->vol0[l] = s->vol[l] = s->x0[l][0] * j * Ndr_ * Ndtheta_ * Ndr_ * Ndx_;
	    }
	  } else {
	    if (j==0) {
	      s->vol0[l] = s->vol[l] = 2 * M_PI * Ndr_ * Ndx_;
	    } else {
	      s->vol0[l] = s->vol[l] = j * Ndr_ * Ndtheta_ * Ndr_ * Ndx_;
	    }
	  }
          s->mass[l] = s->rho0[l] * s->vol0[l];

          s->eff_plastic_strain[l] = 0;
          s->eff_plastic_strain_rate[l] = 0;
          s->damage[l] = 0;
          s->damage_init[l] = 0;
          s->T[l] = T0;
          s->ienergy[l] = 0;
          s->strain_el[l] = Matrix3d();
          s->sigma[l] = Matrix3d();
          s->vol0PK1[l] = Matrix3d();
          s->L[l] = Matrix3d();
          s->F[l] = Matrix3d::identity();
          s->R[l] = Matrix3d::identity();
          s->D[l] = Matrix3d();
          s->Finv[l] = Matrix3d();
          s->Fdot[l] = Matrix3d();
          s->J[l] = 1;
          s->mask[l] = 1;
          s->dtCFL[l] = 0;
          l++;
        }
        if (j == 0)
          break;
      }
    }
  }
  if (s->np_local > l)
    {
      s->grow(l);
    }
  s->np_local = l; // Adjust np_local to account for the particles outside the domain

  tagint ptag0 = 0;

  for (int proc=0; proc<universe->nprocs; proc++){
    int np_local_bcast;
    if (proc == universe->me) {
      // Send np_local
      np_local_bcast = s->np_local;
    } else {
      // Receive np_local
      np_local_bcast = 0;
    }
    MPI_Bcast(&np_local_bcast,1,MPI_INT,proc,universe->uworld);
    if (universe->me > proc) ptag0 += np_local_bcast;
  }

  s->np_local = l; // Adjust np to account for the particles outside the domain

  // cout << "ptag0=" << ptag0 << endl;
  // cout << "s->np_local=" << s->np_local << endl;
  for (int i = 0; i < s->np_local; i++)
  {
    s->ptag[i] = ptag0 + i + 1 + domain->np_total;
  }

  if (l != s->np_local)
  {
    cout << "Error l=" << l << " != np_local=" << s->np_local << endl;
    exit(1);
  }

  MPI_Allreduce(&s->np_local, &np_local_reduced, 1, MPI_INT, MPI_SUM, universe->uworld);
  s->np += np_local_reduced;
  domain->np_total += s->np;

  return Var(0);
}
