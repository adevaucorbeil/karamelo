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

#include "ulmpm.h"
#include "basis_functions.h"
#include "domain.h"
#include "error.h"
#include "grid.h"
#include "input.h"
#include "solid.h"
#include "universe.h"
#include "update.h"
#include "var.h"
#include <Eigen/Eigen>
#include <algorithm>
#include <iostream>
#include <math.h>
#include <vector>

using namespace std;

ULMPM::ULMPM(MPM *mpm) : Method(mpm)
{
  cout << "In ULMPM::ULMPM()" << endl;

  update_wf   = 1;
  update->PIC_FLIP = 0.99;
  apic        = false;

  // Default base function (linear):
  basis_function            = &BasisFunction::linear;
  derivative_basis_function = &BasisFunction::derivative_linear;

  rigid_solids = 0;
}

ULMPM::~ULMPM() {}

void ULMPM::setup(vector<string> args)
{
  if (args.size() > 0) {
    error->all(FLERR, "Illegal modify_method command: too many arguments.\n");
  }

  if (update->shape_function == update->ShapeFunctions::LINEAR) {
    cout << "Setting up linear basis functions\n";
    basis_function = &BasisFunction::linear;
    derivative_basis_function = &BasisFunction::derivative_linear;
  } else if (update->shape_function == update->ShapeFunctions::CUBIC_SPLINE) {
    cout << "Setting up cubic-spline basis functions\n";
    basis_function = &BasisFunction::cubic_spline;
    derivative_basis_function = &BasisFunction::derivative_cubic_spline;
  } else if (update->shape_function == update->ShapeFunctions::QUADRATIC_SPLINE) {
    cout << "Setting up quadratic-spline basis functions\n";
    basis_function = &BasisFunction::quadratic_spline;
    derivative_basis_function = &BasisFunction::derivative_quadratic_spline;
  } else if (update->shape_function == update->ShapeFunctions::BERNSTEIN) {
    cout << "Setting up Bernstein-quadratic basis functions\n";
    basis_function = &BasisFunction::bernstein_quadratic;
    derivative_basis_function = &BasisFunction::derivative_bernstein_quadratic;
  } else {
    error->all(FLERR, "Error: shape function not supported! Supported functions are:  \033[1;32mlinear\033[0m, \033[1;32mcubic-spline\033[0m, \033[1;32mquadratic-spline\033[0m, \033[1;32mBernstein-quadratic\033[0m.\n");
  }

  if (update->sub_method_type == update->SubMethodType::APIC) {
    apic = true;
    update->PIC_FLIP = 0;
  }
}

void ULMPM::compute_grid_weight_functions_and_gradients()
{
  if (!update_wf)
    return;

  bigint nsolids, np_local, nnodes_local, nnodes_ghost;

  nsolids = domain->solids.size();

  if (nsolids)
  {
    for (int isolid = 0; isolid < nsolids; isolid++)
    {
      if (update->ntimestep == 0 && domain->solids[isolid]->mat->rigid)
        rigid_solids = 1;

      np_local = domain->solids[isolid]->np_local;
      nnodes_local = domain->solids[isolid]->grid->nnodes_local;
      nnodes_ghost = domain->solids[isolid]->grid->nnodes_ghost;

      vector<int> *numneigh_pn = &domain->solids[isolid]->numneigh_pn;
      vector<int> *numneigh_np = &domain->solids[isolid]->numneigh_np;

      vector<vector<int>> *neigh_pn = &domain->solids[isolid]->neigh_pn;
      vector<vector<int>> *neigh_np = &domain->solids[isolid]->neigh_np;

      vector<vector<double>> *wf_pn = &domain->solids[isolid]->wf_pn;
      vector<vector<double>> *wf_np = &domain->solids[isolid]->wf_np;

      vector<vector<Eigen::Vector3d>> *wfd_pn = &domain->solids[isolid]->wfd_pn;
      vector<vector<Eigen::Vector3d>> *wfd_np = &domain->solids[isolid]->wfd_np;

      Eigen::Vector3d r;
      double s[3], sd[3];
      vector<Eigen::Vector3d> *xp = &domain->solids[isolid]->x;
      vector<Eigen::Vector3d> *xn = &domain->solids[isolid]->grid->x0;
      double inv_cellsize = 1.0 / domain->solids[isolid]->grid->cellsize;
      double wf;
      Eigen::Vector3d wfd;

      vector<array<int, 3>> *ntype = &domain->solids[isolid]->grid->ntype;
      vector<bool> *nrigid = &domain->solids[isolid]->grid->rigid;

      unordered_map<int, int> *map_ntag = &domain->solids[isolid]->grid->map_ntag;
      unordered_map<int, int>::iterator it;

      r.setZero();

      for (int in = 0; in < nnodes_local + nnodes_ghost; in++)
      {
        (*neigh_np)[in].clear();
        (*numneigh_np)[in] = 0;
        (*wf_np)[in].clear();
        (*wfd_np)[in].clear();
      }

      if (np_local && (nnodes_local + nnodes_ghost))
      {
        for (int ip = 0; ip < np_local; ip++)
        {
          (*neigh_pn)[ip].clear();
          (*numneigh_pn)[ip] = 0;
          (*wf_pn)[ip].clear();
          (*wfd_pn)[ip].clear();

          // Calculate what nodes particle ip will interact with:
	  int nx = domain->solids[isolid]->grid->nx_global;
	  int ny = domain->solids[isolid]->grid->ny_global;
	  int nz = domain->solids[isolid]->grid->nz_global;

          vector<int> n_neigh;

          if (update->shape_function == update->ShapeFunctions::LINEAR)
          {
	    int i0 = (int) (((*xp)[ip][0] - domain->boxlo[0])*inv_cellsize);
	    int j0 = (int) (((*xp)[ip][1] - domain->boxlo[1])*inv_cellsize);
	    int k0 = (int) (((*xp)[ip][2] - domain->boxlo[2])*inv_cellsize);

            for (int i = i0; i < i0 + 2; i++)
            {
              if (ny > 1)
              {
                for (int j = j0; j < j0 + 2; j++)
                {
                  if (nz > 1)
                  {
                    for (int k = k0; k < k0 + 2; k++)
                    {
		      it = (*map_ntag).find(nz * ny * i + nz * j + k);
		      if (it != (*map_ntag).end())
			{
			  n_neigh.push_back(it->second);
			}
                    }
                  }
                  else
                  {
		    it = (*map_ntag).find(ny * i + j);
		    if (it != (*map_ntag).end())
		      {
			n_neigh.push_back(it->second);
		      }
                  }
                }
              }
              else
              {
		if (i < nnodes_local + nnodes_ghost)
                  n_neigh.push_back(i);
              }
            }
	  } else {
	    // cubic and quadratic B-splines
            int i0 =
                (int)(((*xp)[ip][0] - domain->boxlo[0]) * inv_cellsize - 1);
            int j0 =
                (int)(((*xp)[ip][1] - domain->boxlo[1]) * inv_cellsize - 1);
            int k0 =
                (int)(((*xp)[ip][2] - domain->boxlo[2]) * inv_cellsize - 1);

            for (int i = i0; i < i0 + 4; i++) {
              if (ny > 1) {
                for (int j = j0; j < j0 + 4; j++) {
                  if (nz > 1) {
                    for (int k = k0; k < k0 + 4; k++) {
                      it = (*map_ntag).find(nz * ny * i + nz * j + k);
                      if (it != (*map_ntag).end()) {
                        n_neigh.push_back(it->second);
                      }
                    }
                  }
                  else
                  {
		    it = (*map_ntag).find(ny * i + j);
		    if (it != (*map_ntag).end())
		      {
			n_neigh.push_back(it->second);
		      }
                  }
                }
              }
              else
              {
		if (i < nnodes_local + nnodes_ghost)
                  n_neigh.push_back(i);
              }
            }
          }

          // cout << "[";
          // for (auto ii: n_neigh)
          //   cout << domain->solids[isolid]->grid->ntag[ii] << ' ';
          // cout << "]\n";

          // for (int in=0; in<nnodes; in++) {
          for (auto in : n_neigh)
          {

            // Calculate the distance between each pair of particle/node:
	    r = ((*xp)[ip] - (*xn)[in]) * inv_cellsize;

	    s[0] = basis_function(r[0], (*ntype)[in][0]);
	    wf = s[0];
	    if (wf != 0) {
	      if (domain->dimension >= 2) {
		s[1] = basis_function(r[1], (*ntype)[in][1]);
		wf *= s[1];
	      }
	      else s[1] = 1;
	      if (domain->dimension == 3 && wf != 0) {
		s[2] = basis_function(r[2], (*ntype)[in][2]);
		wf *= s[2];
	      }
	      else s[2] = 1;
	    }

	    if (wf != 0)
            {
              if (domain->solids[isolid]->mat->rigid)
                (*nrigid)[in] = true;
              // // cout << in << "\t";
              // // Check if this node is in n_neigh:
              // if (find(n_neigh.begin(), n_neigh.end(), in) == n_neigh.end())
              // {
              // 	// in is not in n_neigh
              //  	cout << "in=" << in << " not found in n_neigh for ip=" << ip
              //  << " which is :["; 	for (auto ii: n_neigh) 	  cout << ii << ' ';
              //  	cout << "]\n";
              // }

	      sd[0] = derivative_basis_function(r[0], (*ntype)[in][0], inv_cellsize);
	      if (domain->dimension >= 2) sd[1] = derivative_basis_function(r[1], (*ntype)[in][1], inv_cellsize);
	      if (domain->dimension == 3) sd[2] = derivative_basis_function(r[2], (*ntype)[in][2], inv_cellsize);

	      (*neigh_pn)[ip].push_back(in);
              (*neigh_np)[in].push_back(ip);
              (*numneigh_pn)[ip]++;
              (*numneigh_np)[in]++;

              (*wf_pn)[ip].push_back(wf);
              (*wf_np)[in].push_back(wf);

              
              if (domain->dimension == 3)
              {
                wfd[0] = sd[0] * s[1] * s[2];
                wfd[1] = s[0] * sd[1] * s[2];
                wfd[2] = s[0] * s[1] * sd[2];
              }
              else if (domain->dimension == 2)
              {
                wfd[0] = sd[0] * s[1];
                wfd[1] = s[0] * sd[1];
                wfd[2] = 0;
              }
	      else
              {
                wfd[0] = sd[0];
                wfd[1] = 0;
                wfd[2] = 0;
              }
              (*wfd_pn)[ip].push_back(wfd);
              (*wfd_np)[in].push_back(wfd);
              // cout << "ip=" << ip << ", in=" << in << ", wf=" << wf << ",
              // wfd=[" << wfd[0] << "," << wfd[1] << "," << wfd[2] << "]" <<
              // endl;
            }
          }
          // cout << endl;
        }
      }
      if (apic)
        domain->solids[isolid]->compute_inertia_tensor();
    }
  } // end if (nsolids)

  if (update->ntimestep == 0) {
    // Reduce rigid_solids
    int rigid_solids_reduced = 0;

    MPI_Allreduce(&rigid_solids, &rigid_solids_reduced, 1, MPI_INT, MPI_LOR,
		  universe->uworld);
    rigid_solids = rigid_solids_reduced;
  }
  if (rigid_solids) {
    domain->grid->reduce_rigid_ghost_nodes();
  }
}

void ULMPM::particles_to_grid() {
  bool grid_reset = false; // Indicate if the grid quantities have to be reset
  for (int isolid = 0; isolid < domain->solids.size(); isolid++) {

    if (isolid == 0)
      grid_reset = true;
    else
      grid_reset = false;

    domain->solids[isolid]->compute_mass_nodes(grid_reset);
  }

  domain->grid->reduce_mass_ghost_nodes();

  for (int isolid = 0; isolid < domain->solids.size(); isolid++) {

    if (isolid == 0)
      grid_reset = true;
    else
      grid_reset = false;

    if (apic)
      domain->solids[isolid]->compute_velocity_nodes_APIC(grid_reset);
    else
      domain->solids[isolid]->compute_velocity_nodes(grid_reset);
    domain->solids[isolid]->compute_external_forces_nodes(grid_reset);
    domain->solids[isolid]->compute_internal_forces_nodes_UL(grid_reset);

    if (temp) {
      domain->solids[isolid]->compute_temperature_nodes(grid_reset);
      domain->solids[isolid]->compute_external_temperature_driving_forces_nodes(grid_reset);
      domain->solids[isolid]->compute_internal_temperature_driving_forces_nodes();
    }
  }
  domain->grid->reduce_ghost_nodes(false, temp);
}

void ULMPM::update_grid_state() {
  domain->grid->update_grid_velocities();
  if (temp)
    domain->grid->update_grid_temperature();
}

void ULMPM::grid_to_points()
{
  for (int isolid = 0; isolid < domain->solids.size(); isolid++)
  {
    if (apic)
      domain->solids[isolid]->compute_rate_deformation_gradient_UL_APIC();
    domain->solids[isolid]->compute_particle_velocities_and_positions();
    domain->solids[isolid]->compute_particle_acceleration();
    if (temp) {
      domain->solids[isolid]->update_particle_temperature();
    }
  }
}

void ULMPM::advance_particles()
{
  for (int isolid = 0; isolid < domain->solids.size(); isolid++)
  {
    domain->solids[isolid]->update_particle_velocities(update->PIC_FLIP);
  }
}

void ULMPM::velocities_to_grid()
{
  bool grid_reset = false; // Indicate if the grid quantities have to be reset
  for (int isolid = 0; isolid < domain->solids.size(); isolid++)
  {

    if (isolid == 0)
      grid_reset = true;
    else
      grid_reset = false;

    if (!apic)
    {
      // domain->solids[isolid]->compute_mass_nodes(grid_reset);
      domain->solids[isolid]->compute_velocity_nodes(grid_reset);
      if (temp) {
	domain->solids[isolid]->compute_temperature_nodes(grid_reset);
      }
    }
  }
  domain->grid->reduce_ghost_nodes(true, temp);
}

void ULMPM::compute_rate_deformation_gradient() {
  if (!apic) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      domain->solids[isolid]->compute_rate_deformation_gradient_UL_MUSL();
      // domain->solids[isolid]->compute_deformation_gradient();
    }
  }
}

void ULMPM::update_deformation_gradient()
{
  for (int isolid = 0; isolid < domain->solids.size(); isolid++)
  {
    domain->solids[isolid]->update_deformation_gradient();
  }
}

void ULMPM::update_stress()
{
  for (int isolid = 0; isolid < domain->solids.size(); isolid++)
  {
    domain->solids[isolid]->update_stress();
    if (temp) {
      domain->solids[isolid]->update_heat_flux();
    }
  }
}

void ULMPM::adjust_dt()
{
  if (update->dt_constant) return; // dt is set as a constant, do not update

  double dtCFL = 1.0e22;
  double dtCFL_reduced = 1.0e22;

  for (int isolid = 0; isolid < domain->solids.size(); isolid++)
  {
    dtCFL = MIN(dtCFL, domain->solids[isolid]->dtCFL);
    if (dtCFL == 0)
    {
      cout << "Error: dtCFL == 0\n";
      cout << "domain->solids[" << isolid << "]->dtCFL == 0\n";
      error->one(FLERR, "");
    } else if (std::isnan(dtCFL)) {
      cout << "Error: dtCFL = " << dtCFL << "\n";
      cout << "domain->solids[" << isolid << "]->dtCFL == " << domain->solids[isolid]->dtCFL << "\n";
      error->one(FLERR, "");
    }
  }

  MPI_Allreduce(&dtCFL, &dtCFL_reduced, 1, MPI_DOUBLE, MPI_MIN, universe->uworld);

  update->dt = dtCFL_reduced * update->dt_factor;
  (*input->vars)["dt"] = Var("dt", update->dt);
}

void ULMPM::reset()
{
  int np_local;

  for (int isolid = 0; isolid < domain->solids.size(); isolid++)
  {
    domain->solids[isolid]->dtCFL = 1.0e22;
    np_local = domain->solids[isolid]->np_local;
    for (int ip = 0; ip < np_local; ip++) domain->solids[isolid]->mbp[ip].setZero();
  }
}

void ULMPM::exchange_particles() {
  int ip, np_local_old, size_buf;
  vector<Eigen::Vector3d> *xp;
  vector<double> buf_send;
  vector<double> buf_send_vect[universe->nprocs];
  vector<double> buf_recv_vect[universe->nprocs];
  vector<int> unpack_list;

  // Identify the particles that are not in the subdomain
  // and transfer their variables to the buffer:

  for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
    buf_send_vect[universe->me].clear();
    np_local_old = domain->solids[isolid]->np_local;
    xp = &domain->solids[isolid]->x;

    ip = 0;
    while (ip < domain->solids[isolid]->np_local) {
      if (!domain->inside_subdomain((*xp)[ip][0], (*xp)[ip][1], (*xp)[ip][2])) {
        // The particle is not located in the subdomain anymore:
        // transfer it to the buffer
        domain->solids[isolid]->pack_particle(ip, buf_send_vect[universe->me]);
        domain->solids[isolid]->copy_particle(
            domain->solids[isolid]->np_local - 1, ip);
        domain->solids[isolid]->np_local--;
      } else {
        ip++;
      }
    }

    // Resize particle variables:
    if (np_local_old - domain->solids[isolid]->np_local !=
        buf_send_vect[universe->me].size() / domain->solids[isolid]->comm_n) {
      error->one(
          FLERR,
          "Size of buffer does not match the number of particles that left the "
          "domain: " +
              to_string(np_local_old - domain->solids[isolid]->np_local) +
              "!=" + to_string(buf_send_vect[universe->me].size()) + "\n");
    }
    if (buf_send_vect[universe->me].size()) {
      domain->solids[isolid]->grow(domain->solids[isolid]->np_local);
    }

    // Exchange buffers:
    for (int sproc = 0; sproc < universe->nprocs; sproc++) {
      if (sproc == universe->me) {
        size_buf = buf_send_vect[sproc].size();
      }

      MPI_Bcast(&size_buf, 1, MPI_INT, sproc, universe->uworld);

      if (sproc != universe->me) {
        buf_send_vect[sproc].resize(size_buf);
      }

      MPI_Bcast(&buf_send_vect[sproc][0], size_buf, MPI_DOUBLE, sproc,
                universe->uworld);
    }

    // Check what particles are within the subdomain:
    for (int sproc = 0; sproc < universe->nprocs; sproc++) {
      if (sproc != universe->me) {
        unpack_list.clear();
        ip = 0;
        while (ip < buf_send_vect[sproc].size()) {
          if (domain->inside_subdomain(buf_send_vect[sproc][ip + 1],
                                       buf_send_vect[sproc][ip + 2],
                                       buf_send_vect[sproc][ip + 3])) {
            unpack_list.push_back(ip);
          }
          ip += domain->solids[isolid]->comm_n;
        }

        domain->solids[isolid]->grow(domain->solids[isolid]->np_local +
                                     unpack_list.size());

        // Unpack buffer:
        domain->solids[isolid]->unpack_particle(
            domain->solids[isolid]->np_local, unpack_list,
            buf_send_vect[sproc]);
      }
    }
  }
}
