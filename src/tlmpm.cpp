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

#include "tlmpm.h"
#include "basis_functions.h"
#include "domain.h"
#include "error.h"
#include "grid.h"
#include "input.h"
#include "method.h"
#include "solid.h"
#include "universe.h"
#include "update.h"
#include "var.h"
#include <Eigen/Eigen>
#include <algorithm>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>

using namespace std;

TLMPM::TLMPM(MPM *mpm) : Method(mpm) {
  // cout << "In TLMPM::TLMPM()" << endl;

  update_wf = true;
  update_mass_nodes = true;
  update->PIC_FLIP = 0.99;
  is_TL = true;

  // Default base function (linear):
  basis_function = &BasisFunction::linear;
  derivative_basis_function = &BasisFunction::derivative_linear;
}

TLMPM::~TLMPM()
{
}

void TLMPM::setup(vector<string> args)
{

  if (args.size() > 0) {
    error->all(FLERR, "Illegal modify_method command: too many arguments.\n");
  }
  
  if (update->shape_function == Update::ShapeFunctions::LINEAR) {
    if (universe->me == 0)
      cout << "Setting up linear basis functions\n";
    basis_function = &BasisFunction::linear;
    derivative_basis_function = &BasisFunction::derivative_linear;
  } else if (update->shape_function == Update::ShapeFunctions::CUBIC_SPLINE) {
    if (universe->me == 0)
      cout << "Setting up cubic-spline basis functions\n";
    basis_function = &BasisFunction::cubic_spline;
    derivative_basis_function = &BasisFunction::derivative_cubic_spline;
  } else if (update->shape_function == Update::ShapeFunctions::QUADRATIC_SPLINE) {
    if (universe->me == 0)
      cout << "Setting up quadratic-spline basis functions\n";
    basis_function = &BasisFunction::quadratic_spline;
    derivative_basis_function = &BasisFunction::derivative_quadratic_spline;
  } else if (update->shape_function == Update::ShapeFunctions::BERNSTEIN) {
    if (universe->me == 0)
      cout << "Setting up Bernstein-quadratic basis functions\n";
    basis_function = &BasisFunction::bernstein_quadratic;
    derivative_basis_function = &BasisFunction::derivative_bernstein_quadratic;
  } else {
    error->all(FLERR, "Error: shape function not supported! Supported functions are:  \033[1;32mlinear\033[0m, \033[1;32mcubic-spline\033[0m, \033[1;32mquadratic-spline\033[0m, \033[1;32mBernstein-quadratic\033[0m.\n");
  }

  if (update->sub_method_type == Update::SubMethodType::APIC) {
    update->PIC_FLIP = 0;
  }
}

void TLMPM::compute_grid_weight_functions_and_gradients()
{
  if (!update_wf) return;

  if (domain->np_local == 0) {
    error->one(
        FLERR,
        "Bad domain decomposition, some CPUs (at least CPU #" +
            to_string(universe->me) +
            ") do not have any particles "
            "attached to.\nTry to increase or decrease the number of CPUs "
            "(using prime numbers might help).\n");
  }

  // cout << "In TLMPM::compute_grid_weight_functions_and_gradients()\n";
  bigint nsolids, np_local, nnodes_local, nnodes_ghost, nnodes;

  nsolids = domain->solids.size();

  if (nsolids) {
    for (int isolid=0; isolid<nsolids; isolid++){

      np_local = domain->solids[isolid]->np_local;
      nnodes = domain->solids[isolid]->grid->nnodes;
      nnodes_local = domain->solids[isolid]->grid->nnodes_local;
      nnodes_ghost = domain->solids[isolid]->grid->nnodes_ghost;

      vector<int> *numneigh_pn = &domain->solids[isolid]->numneigh_pn;
      vector<int> *numneigh_np = &domain->solids[isolid]->numneigh_np;

      vector<vector<int>> *neigh_pn = &domain->solids[isolid]->neigh_pn;
      vector<vector<int>> *neigh_np = &domain->solids[isolid]->neigh_np;

      vector<vector< double >> *wf_pn = &domain->solids[isolid]->wf_pn;
      vector<vector< double >> *wf_np = &domain->solids[isolid]->wf_np;

      vector<vector< Eigen::Vector3d >> *wfd_pn = &domain->solids[isolid]->wfd_pn;
      vector<vector< Eigen::Vector3d >> *wfd_np = &domain->solids[isolid]->wfd_np;

      Eigen::Vector3d r;
      double s[3], sd[3];
      vector<Eigen::Vector3d> *xp = &domain->solids[isolid]->x0;
      vector<Eigen::Vector3d> *xn = &domain->solids[isolid]->grid->x0;
      double inv_cellsize = 1.0 / domain->solids[isolid]->grid->cellsize;
      double wf;
      Eigen::Vector3d wfd;

      vector<array<int, 3>> *ntype = &domain->solids[isolid]->grid->ntype;
      vector<bool> *nrigid = &domain->solids[isolid]->grid->rigid;

      vector<tagint> *map_ntag = &domain->solids[isolid]->grid->map_ntag;
      int inn;
      tagint tag = 0;
      
      r.setZero();
      if (np_local && (nnodes_local + nnodes_ghost)) {

	int nx = domain->solids[isolid]->grid->nx_global;
	int ny = domain->solids[isolid]->grid->ny_global;
	int nz = domain->solids[isolid]->grid->nz_global;

	for (int ip=0; ip<np_local; ip++) {
	  // Calculate what nodes particle ip will interact with:

	  vector<int> n_neigh;

	  if (update->shape_function == Update::ShapeFunctions::LINEAR) {
	    int i0 = (int) (((*xp)[ip][0] - domain->solids[isolid]->solidlo[0])*inv_cellsize);
	    int j0 = (int) (((*xp)[ip][1] - domain->solids[isolid]->solidlo[1])*inv_cellsize);
	    int k0 = (int) (((*xp)[ip][2] - domain->solids[isolid]->solidlo[2])*inv_cellsize);

	    // cout << "(" << i0 << "," << j0 << "," << k0 << ")\t";

	    for(int i=i0; i<i0+2;i++){
	      if (ny>1){
		for(int j=j0; j<j0+2;j++){
		  if (nz>1){
		    for(int k=k0; k<k0+2;k++){
		      tag = nz * ny * i + nz * j + k;
		      if (tag < nnodes) {
			inn = (*map_ntag)[tag];
			if (inn != -1) {
			  n_neigh.push_back(inn);
			}
		      }
		    }
		  } else {
		    tag = ny * i + j;
		    if (tag < nnodes) {
		      inn = (*map_ntag)[tag];
		      if (inn != -1) {
			n_neigh.push_back(inn);
		      }
		    }
		  }
		}
	      } else {
		if (i < nnodes_local + nnodes_ghost)
		  n_neigh.push_back(i);
	      }
	    }
	  } else if (update->shape_function == Update::ShapeFunctions::BERNSTEIN){
	    int i0 = 2*(int) (((*xp)[ip][0] - domain->solids[isolid]->solidlo[0])*inv_cellsize);
	    int j0 = 2*(int) (((*xp)[ip][1] - domain->solids[isolid]->solidlo[1])*inv_cellsize);
	    int k0 = 2*(int) (((*xp)[ip][2] - domain->solids[isolid]->solidlo[2])*inv_cellsize);

	    if ((i0 >= 1) && (i0 % 2 != 0)) i0--;
	    if ((j0 >= 1) && (j0 % 2 != 0)) j0--;
	    if (nz>1) if ((k0 >= 1) && (k0 % 2 != 0)) k0--;

	    // cout << "(" << i0 << "," << j0 << "," << k0 << ")\t";

	    for(int i=i0; i<i0+3;i++){
	      if (ny>1){
		for(int j=j0; j<j0+3;j++){
		  if (nz>1){
		    for(int k=k0; k<k0+3;k++){
		      tag = nz * ny * i + nz * j + k;
		      if (tag < nnodes) {
			inn = (*map_ntag)[tag];
			if (inn != -1) {
			  n_neigh.push_back(inn);
			}
		      }
		    }
		  } else {
		    tag = ny * i + j;
		    if (tag < nnodes) {
		      inn = (*map_ntag)[tag];
		      if (inn != -1) {
			n_neigh.push_back(inn);
		      }
		    }
		  }
		}
	      } else {
		if (i < nnodes_local + nnodes_ghost)
		  n_neigh.push_back(i);
	      }
	    }
          } else {
	    //(update->shape_function == Update::ShapeFunctions::CUBIC_SPLINE || update->shape_function == Update::ShapeFunctions::QUADRATIC_SPLINE) {
            int i0 = (int) (((*xp)[ip][0] - domain->solids[isolid]->solidlo[0])*inv_cellsize - 1);
	    int j0 = (int) (((*xp)[ip][1] - domain->solids[isolid]->solidlo[1])*inv_cellsize - 1);
	    int k0 = (int) (((*xp)[ip][2] - domain->solids[isolid]->solidlo[2])*inv_cellsize - 1);

	    // cout << "(" << i0 << "," << j0 << "," << k0 << ")\t";

	    for(int i=i0; i<i0+4;i++){
	      if (ny>1) {
		for(int j=j0; j<j0+4;j++){
		  if (nz>1){
		    for(int k=k0; k<k0+4;k++){
		      tag = nz * ny * i + nz * j + k;
			if (tag < nnodes) {
			  inn = (*map_ntag)[nz*ny*i+nz*j+k];
			  if (inn != -1) {
			    n_neigh.push_back(inn);
			  }
			}
		    }
		  } else {
		    tag = ny * i + j;
		    if (tag < nnodes) {
		      inn = (*map_ntag)[ny*i+j];
		      if (inn != -1) {
			n_neigh.push_back(inn);
		      }
		    }
		  }
		}
	      } else {
		if (i < nnodes_local + nnodes_ghost)
		  n_neigh.push_back(i);
	      }
	    }
          }

          // cout << "ip: "<< ip << "\t";
	  // cout << "[";
	  // for (auto ii: n_neigh)
	  //   cout << domain->solids[isolid]->grid->ntag[ii] << ' ';
	  // cout << "]\n";

	  //for (int in=0; in<nnodes; in++) {
	  for (auto in: n_neigh) {
	    // Calculate the distance between each pair of particle/node:
	    r = ((*xp)[ip] - (*xn)[in]) * inv_cellsize;

	    s[0] = basis_function(r[0], (*ntype)[in][0]);
	    if (domain->dimension >= 2) s[1] = basis_function(r[1], (*ntype)[in][1]);
	    else s[1] = 1;
	    if (domain->dimension == 3) s[2] = basis_function(r[2], (*ntype)[in][2]);
	    else s[2] = 1;

	    if (s[0] != 0 && s[1] != 0 && s[2] != 0) {
	      if (domain->solids[isolid]->mat->rigid) (*nrigid)[in] = true;
	      // // cout << in << "\t";
	      // // Check if this node is in n_neigh:
	      // if (find(n_neigh.begin(), n_neigh.end(), in) == n_neigh.end()) {
	      // 	// in is not in n_neigh
	      //  	cout << "in=" << in << " not found in n_neigh for ip=" << ip << " which is :[";
	      //  	for (auto ii: n_neigh)
	      //  	  cout << ii << ' ';
	      //  	cout << "]\n";
	      // }

	      sd[0] = derivative_basis_function(r[0], (*ntype)[in][0], inv_cellsize);
	      if (domain->dimension >= 2) sd[1] = derivative_basis_function(r[1], (*ntype)[in][1], inv_cellsize);
	      if (domain->dimension == 3) sd[2] = derivative_basis_function(r[2], (*ntype)[in][2], inv_cellsize);

	      (*neigh_pn)[ip].push_back(in);
	      (*neigh_np)[in].push_back(ip);
	      (*numneigh_pn)[ip]++;
	      (*numneigh_np)[in]++;
	      if (domain->dimension == 1) wf = s[0];
	      if (domain->dimension == 2) wf = s[0]*s[1];
	      if (domain->dimension == 3) wf = s[0]*s[1]*s[2];

	      (*wf_pn)[ip].push_back(wf);
	      (*wf_np)[in].push_back(wf);

	      if (domain->dimension == 1)
		{
		  wfd[0] = sd[0];
		  wfd[1] = 0;
		  wfd[2] = 0;
		}
	      else if (domain->dimension == 2)
		{
		  wfd[0] = sd[0]*s[1];
		  wfd[1] = s[0]*sd[1];
		  wfd[2] = 0;
		}
	      else if (domain->dimension == 3)
		{
		  wfd[0] = sd[0]*s[1]*s[2];
		  wfd[1] = s[0]*sd[1]*s[2];
		  wfd[2] = s[0]*s[1]*sd[2];
		}
	      (*wfd_pn)[ip].push_back(wfd);
	      (*wfd_np)[in].push_back(wfd);
	      // cout << "ip=" << ip << ", in=" << in << ", wf=" << wf << ", wfd=[" << wfd[0] << "," << wfd[1] << "," << wfd[2] << "]" << endl;
	    }
	  }
	  // cout << endl;
	}
      }
      if (update->sub_method_type == Update::SubMethodType::APIC) domain->solids[isolid]->compute_inertia_tensor();
    }
  }

  update_wf = false;
}

void TLMPM::particles_to_grid()
{
  bool grid_reset = true; // Indicate if the grid quantities have to be reset
  if (update_mass_nodes) {
    for (int isolid=0; isolid<domain->solids.size(); isolid++){
      domain->solids[isolid]->compute_mass_nodes(grid_reset);
      domain->solids[isolid]->grid->reduce_mass_ghost_nodes();
    }
    update_mass_nodes = false;
  }

  for (int isolid=0; isolid<domain->solids.size(); isolid++){
    if (update->sub_method_type == Update::SubMethodType::APIC) domain->solids[isolid]->compute_velocity_nodes_APIC(grid_reset);
    else domain->solids[isolid]->compute_velocity_nodes(grid_reset);
    domain->solids[isolid]->compute_external_forces_nodes(grid_reset);
    domain->solids[isolid]->compute_internal_forces_nodes_TL();

    if (temp) {
      domain->solids[isolid]->compute_temperature_nodes(grid_reset);
      domain->solids[isolid]->compute_external_temperature_driving_forces_nodes(grid_reset);
      domain->solids[isolid]->compute_internal_temperature_driving_forces_nodes();
    }
    domain->solids[isolid]->grid->reduce_ghost_nodes(true, true, temp);
  }
}

void TLMPM::particles_to_grid_USF_1()
{
  bool grid_reset = true; // Indicate if the grid quantities have to be reset
  if (update_mass_nodes) {
    for (int isolid=0; isolid<domain->solids.size(); isolid++){
      domain->solids[isolid]->compute_mass_nodes(grid_reset);
      domain->solids[isolid]->grid->reduce_mass_ghost_nodes();
    }
    update_mass_nodes = false;
  }

  for (int isolid=0; isolid<domain->solids.size(); isolid++){
    if (update->sub_method_type == Update::SubMethodType::APIC)
      domain->solids[isolid]->compute_velocity_nodes_APIC(grid_reset);
    else
      domain->solids[isolid]->compute_velocity_nodes(grid_reset);

    if (temp)
      domain->solids[isolid]->compute_temperature_nodes(grid_reset);

    domain->solids[isolid]->grid->reduce_ghost_nodes(true, false, temp);
  }
}

void TLMPM::particles_to_grid_USF_2()
{
  bool grid_reset = true; // Indicate if the grid quantities have to be reset

  for (int isolid=0; isolid<domain->solids.size(); isolid++){
    domain->solids[isolid]->compute_external_forces_nodes(grid_reset);
    domain->solids[isolid]->compute_internal_forces_nodes_TL();

    if (temp) {
      domain->solids[isolid]->compute_external_temperature_driving_forces_nodes(grid_reset);
      domain->solids[isolid]->compute_internal_temperature_driving_forces_nodes();
    }
    domain->solids[isolid]->grid->reduce_ghost_nodes(false, true, temp);
  }
}

void TLMPM::update_grid_state()
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    domain->solids[isolid]->grid->update_grid_velocities();
    if (temp) {
      domain->solids[isolid]->grid->update_grid_temperature();
    }
  }
}

void TLMPM::grid_to_points()
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    if (domain->solids[isolid]->mat->rigid) {
      domain->solids[isolid]->compute_particle_velocities_and_positions();
      domain->solids[isolid]->compute_particle_acceleration();
    } else {
      domain->solids[isolid]->compute_particle_accelerations_velocities_and_positions();
    }
    if (temp) {
      domain->solids[isolid]->update_particle_temperature();
    }
  }
}

void TLMPM::advance_particles()
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    domain->solids[isolid]->update_particle_velocities(update->PIC_FLIP);
  }
}

void TLMPM::velocities_to_grid()
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    if (update->sub_method_type == Update::SubMethodType::APIC)
      domain->solids[isolid]->compute_velocity_nodes_APIC(true);
    else
      domain->solids[isolid]->compute_velocity_nodes(true);
    if (temp) {
      domain->solids[isolid]->compute_temperature_nodes(true);
    }
    domain->solids[isolid]->grid->reduce_ghost_nodes(true, false, temp);
  }
}

void TLMPM::update_grid_positions()
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    domain->solids[isolid]->grid->update_grid_positions();
  }
}

void TLMPM::compute_rate_deformation_gradient(bool doublemapping) {
  for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
    if (update->sub_method_type == Update::SubMethodType::APIC)
      domain->solids[isolid]->compute_rate_deformation_gradient_TL_APIC(
          doublemapping);
    else
      domain->solids[isolid]->compute_rate_deformation_gradient_TL(
          doublemapping);
  }
}

void TLMPM::update_deformation_gradient()
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    domain->solids[isolid]->update_deformation_gradient();
  }
}

void TLMPM::update_stress(bool doublemapping)
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    domain->solids[isolid]->update_stress();
    if (temp) {
      domain->solids[isolid]->update_heat_flux(doublemapping);
    }
  }
}

void TLMPM::adjust_dt()
{
  if (update->dt_constant) return; // dt is set as a constant, do not update


  double dtCFL = 1.0e22;
  double dtCFL_reduced = 1.0e22;

  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    dtCFL = MIN(dtCFL, domain->solids[isolid]->dtCFL);
    if (dtCFL == 0) {
      cout << "Error: dtCFL == 0\n";
      cout << "domain->solids[" << isolid << "]->dtCFL == 0\n";
      error->all(FLERR, "");
    } else if (std::isnan(dtCFL)) {
      cout << "Error: dtCFL = " << dtCFL << "\n";
      cout << "domain->solids[" << isolid << "]->dtCFL == " << domain->solids[isolid]->dtCFL << "\n";
      error->all(FLERR, "");
    }
  }

  MPI_Allreduce(&dtCFL, &dtCFL_reduced, 1, MPI_DOUBLE, MPI_MIN, universe->uworld);

  update->dt = dtCFL_reduced * update->dt_factor;
  (*input->vars)["dt"] = Var("dt", update->dt);
}

void TLMPM::reset()
{
  int np_local;

  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    domain->solids[isolid]->dtCFL = 1.0e22;
    np_local = domain->solids[isolid]->np_local;
    for (int ip = 0; ip < np_local; ip++) domain->solids[isolid]->mbp[ip].setZero();
  }
}
