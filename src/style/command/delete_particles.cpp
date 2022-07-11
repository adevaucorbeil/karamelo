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

#include <delete_particles.h>
#include <domain.h>
#include <error.h>
#include <group.h>
#include <iostream>
#include <method.h>
#include <universe.h>
#include <update.h>
#include <var.h>

using namespace std;

DeleteParticles::DeleteParticles(MPM *mpm) : Pointers(mpm) {}

Var DeleteParticles::command(vector<string> args) {
  // cout << "In DeleteParticles::command()" << endl;

  if (args.size() < 3) error->all(FLERR, "Illegal delete_particles command\n");

  int ns = domain->solids.size();

  int isolid = domain->find_solid(args[0]);

  if (isolid < 0) {
    if (args[0].compare("all")!=0) error->all(FLERR, "Error: solid " + args[0] + " unknown.\n");
  }

  if (args[1].compare("region")==0) delete_region(args, isolid);
  else error->all(FLERR, "Error: use of illegal keyword for delete_particles command: " + args[1] + "\n");

  return Var(0);
}

void DeleteParticles::delete_region(vector<string> args, int isolid) {
  int iregion = domain->find_region(args[2]);

  if (iregion < 0)
    error->all(FLERR, "Error: region " + args[2] + " unknown.\n");

  const int &ns = domain->solids.size();

  for (int is = 0; is < ns; is++) {
    if ((isolid < 0) || (is == isolid)) {

      Solid &s = *domain->solids[is];

      int &np_local = s.np_local;

      /* Create the list of particles to remove*/
      Kokkos::View<bool*> dlist("dlist", np_local);

      Kokkos::View<bool*>::HostMirror     dlist_mirror = create_mirror(dlist);
      Kokkos::View<Vector3d*>::HostMirror x0_mirror    = create_mirror(s.x0);
      deep_copy(x0_mirror, s.x0);

      for (int ip = 0; ip < np_local; ip++) {
        if (domain->regions[iregion]->inside(
                x0_mirror[ip][0], x0_mirror[ip][1],
                x0_mirror[ip][2]) == 1)
          dlist_mirror[ip] = true;
        else
          dlist_mirror[ip] = false;
      }

      deep_copy(dlist, dlist_mirror);
      int ndelete;


      /* Create a View containing the indexes of the particles to keep.*/
      Kokkos::View<int*> sum_post("sum_post", np_local);
      Kokkos::View<int*> inew    ("inew"    , np_local);

      Kokkos::parallel_scan("scan", np_local,
			    KOKKOS_LAMBDA(int ip, int &partial_sum, bool is_final)
      {
	partial_sum += dlist[ip];
	if(is_final) {
	  sum_post[ip] = partial_sum;
	  if (!dlist[ip])
	    inew[ip - partial_sum] = ip;
	}
      }, ndelete);

      bool temp = update->method->temp;

      Kokkos::View<tagint*> ptag = s.ptag;

      Kokkos::View<Vector3d*> x = s.x;
      Kokkos::View<Vector3d*> x0 = s.x0;
      Kokkos::View<Vector3d*> v = s.v;
      Kokkos::View<Vector3d*> v_update = s.v_update;
      Kokkos::View<Vector3d*> a = s.a;

      Kokkos::View<Vector3d*> mbp = s.mbp;
      Kokkos::View<Vector3d*> f = s.f;

      Kokkos::View<Matrix3d*> sigma = s.sigma;
      Kokkos::View<Matrix3d*> strain_el = s.strain_el;
      Kokkos::View<Matrix3d*> vol0PK1 = s.vol0PK1;
      Kokkos::View<Matrix3d*> L = s.L;
      Kokkos::View<Matrix3d*> F = s.F;
      Kokkos::View<Matrix3d*> R = s.R;
      Kokkos::View<Matrix3d*> D = s.D;
      Kokkos::View<Matrix3d*> Finv = s.Finv;
      Kokkos::View<Matrix3d*> Fdot = s.Fdot;

      Kokkos::View<double*> J = s.J;
      Kokkos::View<double*> vol = s.vol;
      Kokkos::View<double*> vol0 = s.vol0;
      Kokkos::View<double*> rho = s.rho;
      Kokkos::View<double*> rho0 = s.rho0;
      Kokkos::View<double*> mass = s.mass;
      Kokkos::View<double*> eff_plastic_strain = s.eff_plastic_strain;
      Kokkos::View<double*> eff_plastic_strain_rate = s.eff_plastic_strain_rate;
      Kokkos::View<double*> damage = s.damage;
      Kokkos::View<double*> damage_init = s.damage_init;
      Kokkos::View<double*> ienergy = s.ienergy;
      Kokkos::View<int*> mask = s.mask;

      Kokkos::View<double*> T = s.T;
      Kokkos::View<double*> gamma = s.gamma;
      Kokkos::View<Vector3d*> q = s.q;

      Kokkos::View<int*> error_flag = s.error_flag;
      Kokkos::View<double*> dtCFL = s.dtCFL;

      double vtot_local;
      double mtot_local;

      np_local -= ndelete;

      Kokkos::parallel_reduce("copy_particles", np_local,
			      KOKKOS_LAMBDA(int ip, double &vtot_, double &mtot_)
      {
	if (dlist[ip]) {
	  const int &j = inew[np_local - sum_post[ip]];
	  ptag[ip]     = ptag[j];
	  x0[ip]       = x0[j];
	  x[ip]        = x[j];
	  v[ip]        = v[j];
	  v_update[ip] = v_update[j];
	  a[ip]        = a[j];
	  mbp[ip]      = mbp[j];
	  f[ip]        = f[j];
	  vol0[ip]     = vol0[j];
	  vol[ip]      = vol[j];
	  rho0[ip]     = rho0[j];
	  rho[ip]      = rho[j];
	  mass[ip]     = mass[j];
	  eff_plastic_strain[ip]      = eff_plastic_strain[j];
	  eff_plastic_strain_rate[ip] = eff_plastic_strain_rate[j];
	  damage[ip]      = damage[j];
	  damage_init[ip] = damage_init[j];
	  if (temp) {
	    T[ip]        = T[j];
	    gamma[ip]    = gamma[j];
	    q[ip]        = q[j];
	  }
	  ienergy[ip]    = ienergy[j];
	  mask[ip]       = mask[j];
	  sigma[ip]      = sigma[j];
	  strain_el[ip]  = strain_el[j];
	  vol0PK1[ip]    = vol0PK1[j];
	  L[ip]          = L[j];
	  F[ip]          = F[j];
	  R[ip]          = R[j];
	  D[ip]          = D[j];
	  Finv[ip]       = Finv[j];
	  dtCFL[ip]      = dtCFL[j];
	  Fdot[ip]       = Fdot[j];
	  J[ip]          = J[j];
	  error_flag[ip] = error_flag[j];
	}
	vtot_ += vol[ip];
	mtot_ += mass[ip];
      }, vtot_local, mtot_local);

        resize(ptag, np_local);
	resize(x0, np_local);
	resize(x, np_local);
	resize(v, np_local);
	resize(v_update, np_local);
	resize(a, np_local);
	resize(mbp, np_local);
	resize(f, np_local);
	resize(vol0, np_local);
	resize(vol, np_local);
	resize(rho0, np_local);
	resize(rho, np_local);
	resize(mass, np_local);
	resize(eff_plastic_strain, np_local);
	resize(eff_plastic_strain_rate, np_local);
	resize(damage, np_local);
	resize(damage_init, np_local);
	if (temp) {
	  resize(T, np_local);
	  resize(gamma, np_local);
	  resize(q, np_local);
	}
	resize(ienergy, np_local);
	resize(mask, np_local);
	resize(sigma, np_local);
	resize(strain_el, np_local);
	resize(vol0PK1, np_local);
	resize(L, np_local);
	resize(F, np_local);
	resize(R, np_local);
	resize(D, np_local);
	resize(Finv, np_local);
	resize(dtCFL, np_local);
	resize(Fdot, np_local);
	resize(J, np_local);
	resize(error_flag, np_local);

	double values[3] = {(double) ndelete, vtot_local, mtot_local};
	double reduced_values[3] = {0, 0, 0};

	MPI_Allreduce(values, reduced_values, 3, MPI_DOUBLE, MPI_SUM, universe->uworld);

	int ndelete_reduced = (int) reduced_values[0];
	s.vtot = reduced_values[1];
	s.mtot = reduced_values[2];
	if (universe->me == 0) {
	  cout << "Deleting " << ndelete_reduced
	       << " particles from solid " << s.id << endl;
	  cout << "Solid " << s.id
	       << " new total volume = " << s.vtot << endl;
	}

	domain->np_total -= ndelete_reduced;
    }
  }
}
