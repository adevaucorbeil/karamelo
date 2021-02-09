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
#include <map>
#include <math.h>
#include <string>
#include <vector>

using namespace std;

TLMPM::TLMPM(MPM *mpm, vector<string> args) : Method(mpm) {
  cout << "In TLMPM::TLMPM()" << endl;

  update_wf = 1;
  method_type = "FLIP";
  FLIP = 0.99;
  is_TL = true;

  // Default base function (linear):
  shape_function = "linear";
  basis_function = &BasisFunction::linear;
  derivative_basis_function = &BasisFunction::derivative_linear;
}

TLMPM::~TLMPM()
{
}

void TLMPM::setup(vector<string> args)
{
  int n = 1;
  bool isFLIP = false;
  // Method used: PIC, FLIP or APIC:
  if (args[n].compare("PIC") == 0) {
    method_type = "PIC";
    FLIP = 0;
  } else if (args[n].compare("FLIP") == 0) {
    method_type = "FLIP";
    isFLIP = true;

    if (args.size() < 2) {
      error->all(FLERR, "Illegal modify_method command: not enough arguments.\n");
    }

  } else if (args[n].compare("APIC") == 0) {
    method_type = "APIC";
    FLIP = 0;
  } else {
    error->all(FLERR, "Error: method type " + args[n] + " not understood. Expect: PIC, FLIP or APIC\n");
  }

  n++;
  
  if (args.size() > 1 + isFLIP) {
    if (args[n].compare("linear") == 0) {
      shape_function = "linear";
      cout << "Setting up linear basis functions\n";
      basis_function = &BasisFunction::linear;
      derivative_basis_function = &BasisFunction::derivative_linear;
      n++;
    } else if (args[n].compare("cubic-spline") == 0) {
      shape_function = "cubic-spline";
      cout << "Setting up cubic-spline basis functions\n";
      basis_function = &BasisFunction::cubic_spline;
      derivative_basis_function = &BasisFunction::derivative_cubic_spline;
      n++;
    } else if (args[n].compare("quadratic-spline") == 0) {
      shape_function = "quadratic-spline";
      cout << "Setting up quadratic-spline basis functions\n";
      basis_function = &BasisFunction::quadratic_spline;
      derivative_basis_function = &BasisFunction::derivative_quadratic_spline;
      n++;
    } else if (args[n].compare("Bernstein-quadratic") == 0) {
      shape_function = "Bernstein-quadratic";
      cout << "Setting up Bernstein-quadratic basis functions\n";
      basis_function = &BasisFunction::bernstein_quadratic;
      derivative_basis_function = &BasisFunction::derivative_bernstein_quadratic;
      n++;
    } else {
      error->all(FLERR, "Illegal method_method argument: form function of type \033[1;31m" + args[n] + "\033[0m is unknown. Available options are:  \033[1;32mlinear\033[0m, \033[1;32mcubic-spline\033[0m, \033[1;32mquadratic-spline\033[0m, \033[1;32mBernstein-quadratic\033[0m.\n");
    }
  }

  if (args.size() > n + isFLIP + 1) {
    error->all(FLERR,
               "Illegal modify_method command: too many arguments: maximum " +
                   to_string(n + isFLIP + 1) + " expected, " +
                   to_string(args.size()) + " received.\n");
  }

  if (isFLIP)
    FLIP = input->parsev(args[n]);

  if (args.size() == n + isFLIP + 1) {
    if (args[n + 1].compare("thermo-mechanical") == 0) {
      temp = true;
    } else if (args[n + 1].compare("mechanical") == 0) {
      temp = false;
    } else {
      error->all(
          FLERR,
          "Illegal modify_method command: keyword " + args[n + 1] +
              " unknown. Expected \"thermo-mechanical\" or \"mechanical\".\n");
    }
  }
}

void TLMPM::compute_grid_weight_functions_and_gradients()
{
  if (!update_wf) return;

  // cout << "In TLMPM::compute_grid_weight_functions_and_gradients()\n";
  bigint nsolids, np_local, nnodes_local, nnodes_ghost;

  nsolids = domain->solids.size();

  if (nsolids) {
    for (int isolid=0; isolid<nsolids; isolid++){

      np_local = domain->solids[isolid]->np_local;
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
      vector<Eigen::Vector3d> *xp = &domain->solids[isolid]->x;
      vector<Eigen::Vector3d> *xn = &domain->solids[isolid]->grid->x0;
      double inv_cellsize = 1.0 / domain->solids[isolid]->grid->cellsize;
      double wf;
      Eigen::Vector3d wfd;

      vector<array<int, 3>> *ntype = &domain->solids[isolid]->grid->ntype;
      vector<bool> *nrigid = &domain->solids[isolid]->grid->rigid;

      map<int, int> *map_ntag = &domain->solids[isolid]->grid->map_ntag;
      map<int, int>::iterator it;

      bool linear, cubic, quadratic, bernstein;
      linear = cubic = quadratic = bernstein = false;

      if (update->method_shape_function.compare("linear")==0) linear = true;
      if (update->method_shape_function.compare("cubic-spline")==0) cubic = true;
      if (update->method_shape_function.compare("quadratic-spline")==0) quadratic = true;
      if (update->method_shape_function.compare("Bernstein-quadratic")==0) bernstein = true;
      
      r.setZero();
      if (np_local && (nnodes_local + nnodes_ghost)) {

	int nx = domain->solids[isolid]->grid->nx_global;
	int ny = domain->solids[isolid]->grid->ny_global;
	int nz = domain->solids[isolid]->grid->nz_global;

	for (int ip=0; ip<np_local; ip++) {
	  // Calculate what nodes particle ip will interact with:

	  vector<int> n_neigh;

	  if (linear) {
	    int i0 = (int) (((*xp)[ip][0] - domain->solids[isolid]->solidlo[0])*inv_cellsize);
	    int j0 = (int) (((*xp)[ip][1] - domain->solids[isolid]->solidlo[1])*inv_cellsize);
	    int k0 = (int) (((*xp)[ip][2] - domain->solids[isolid]->solidlo[2])*inv_cellsize);

	    // cout << "(" << i0 << "," << j0 << "," << k0 << ")\t";

	    for(int i=i0; i<i0+2;i++){
	      if (ny>1){
		for(int j=j0; j<j0+2;j++){
		  if (nz>1){
		    for(int k=k0; k<k0+2;k++){
		      it = (*map_ntag).find(nz*ny*i+nz*j+k);
		      if (it != (*map_ntag).end()) {
			n_neigh.push_back(it->second);
		      }
		    }
		  } else {
		    it = (*map_ntag).find(ny*i+j);
		    if (it != (*map_ntag).end()) {
		      n_neigh.push_back(it->second);
		    }
		  }
		}
	      } else {
		if (i < nnodes_local + nnodes_ghost)
		  n_neigh.push_back(i);
	      }
	    }
	  } else if (bernstein){
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
		      it = (*map_ntag).find(nz*ny*i+nz*j+k);
		      if (it != (*map_ntag).end()) {
			n_neigh.push_back(it->second);
		      }
		    }
		  } else {
		    it = (*map_ntag).find(ny*i+j);
		    if (it != (*map_ntag).end()) {
		      n_neigh.push_back(it->second);
		    }
		  }
		}
	      } else {
		if (i < nnodes_local + nnodes_ghost)
		  n_neigh.push_back(i);
	      }
	    }
	  } else if (cubic || quadratic) {
	    int i0 = (int) (((*xp)[ip][0] - domain->solids[isolid]->solidlo[0])*inv_cellsize - 1);
	    int j0 = (int) (((*xp)[ip][1] - domain->solids[isolid]->solidlo[1])*inv_cellsize - 1);
	    int k0 = (int) (((*xp)[ip][2] - domain->solids[isolid]->solidlo[2])*inv_cellsize - 1);

	    // cout << "(" << i0 << "," << j0 << "," << k0 << ")\t";

	    for(int i=i0; i<i0+4;i++){
	      if (ny>1) {
		for(int j=j0; j<j0+4;j++){
		  if (nz>1){
		    for(int k=k0; k<k0+4;k++){
		      it = (*map_ntag).find(nz*ny*i+nz*j+k);
		      if (it != (*map_ntag).end()) {
			n_neigh.push_back(it->second);
		      }
		    }
		  } else {
		    it = (*map_ntag).find(ny*i+j);
		    if (it != (*map_ntag).end()) {
		      n_neigh.push_back(it->second);
		    }
		  }
		}
	      } else {
		if (i < nnodes_local + nnodes_ghost)
		  n_neigh.push_back(i);
	      }
	    }
	  } else {
	    error->all(FLERR, "Shape function type not supported by TLMPM::compute_grid_weight_functions_and_gradients(): " + update->method_shape_function + ".\n");
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
      if (method_type.compare("APIC") == 0) domain->solids[isolid]->compute_inertia_tensor();
    }
  }

  update_wf = 0;
}

void TLMPM::particles_to_grid()
{
  bool grid_reset = true; // Indicate if the grid quantities have to be reset
  for (int isolid=0; isolid<domain->solids.size(); isolid++){
    if (update->ntimestep==1)
      {
	domain->solids[isolid]->compute_mass_nodes(grid_reset);
	domain->solids[isolid]->grid->reduce_mass_ghost_nodes();
      }
  }

  for (int isolid=0; isolid<domain->solids.size(); isolid++){
    if (method_type.compare("APIC") == 0) domain->solids[isolid]->compute_velocity_nodes_APIC(grid_reset);
    else domain->solids[isolid]->compute_velocity_nodes(grid_reset);
    domain->solids[isolid]->compute_external_forces_nodes(grid_reset);
    domain->solids[isolid]->compute_internal_forces_nodes_TL();

    if (temp) {
      domain->solids[isolid]->compute_temperature_nodes(grid_reset);
      domain->solids[isolid]->compute_external_temperature_driving_forces_nodes(grid_reset);
      domain->solids[isolid]->compute_internal_temperature_driving_forces_nodes();
    }
    domain->solids[isolid]->grid->reduce_ghost_nodes(false, temp);
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
    domain->solids[isolid]->compute_particle_velocities_and_positions();
    domain->solids[isolid]->compute_particle_acceleration();
    if (temp) {
      domain->solids[isolid]->update_particle_temperature();
    }
  }
}

void TLMPM::advance_particles()
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    domain->solids[isolid]->update_particle_velocities(FLIP);
  }
}

void TLMPM::velocities_to_grid()
{
  if (method_type.compare("APIC") != 0) {
    for (int isolid=0; isolid<domain->solids.size(); isolid++) {
      //domain->solids[isolid]->compute_mass_nodes();
      domain->solids[isolid]->compute_velocity_nodes(true);
      if (temp) {
	domain->solids[isolid]->compute_temperature_nodes(true);
      }
      domain->solids[isolid]->grid->reduce_ghost_nodes(true, temp);
    }
  }
}

void TLMPM::update_grid_positions()
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    domain->solids[isolid]->grid->update_grid_positions();
  }
}

void TLMPM::compute_rate_deformation_gradient()
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    if (method_type.compare("APIC") == 0) domain->solids[isolid]->compute_rate_deformation_gradient_TL_APIC();
    else domain->solids[isolid]->compute_rate_deformation_gradient_TL();
    //domain->solids[isolid]->compute_deformation_gradient();
  }
}

void TLMPM::update_deformation_gradient()
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    domain->solids[isolid]->update_deformation_gradient();
  }
}

void TLMPM::update_stress()
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    domain->solids[isolid]->update_stress();
    if (temp) {
      domain->solids[isolid]->update_heat_flux();
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
