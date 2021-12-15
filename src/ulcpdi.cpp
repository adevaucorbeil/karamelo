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

#include <ulcpdi.h>
#include <basis_functions.h>
#include <domain.h>
#include <error.h>
#include <grid.h>
#include <input.h>
#include <solid.h>
#include <universe.h>
#include <update.h>
#include <var.h>
#include <matrix.h>
#include <algorithm>
#include <iostream>
#include <math.h>
#include <vector>

using namespace std;

ULCPDI::ULCPDI(MPM *mpm) : Method(mpm) {

  // cout << "In ULCPDI::ULCPDI()" << endl;

  update_wf = 1;
  is_CPDI = true;
  style = 0;   //Default CPDI style is known_styles[style]="R4";

  // Default base function (linear):
  basis_function = &BasisFunction::linear;
  derivative_basis_function = &BasisFunction::derivative_linear;
}

ULCPDI::~ULCPDI()
{
}

void ULCPDI::setup(vector<string> args)
{
  if (args.size() > 1) {
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

  if (args.size() > 0) {
    bool found_style = false;
    for (int i=0; i<sizeof(known_styles)/sizeof(string); i++) {
      if (known_styles[i].compare(args[0])==0) {
        style = i;
	found_style = true;
        break;
      }
    }
    if (!found_style) {
      string error_str = "CPDI style \033[1;31m" + args[0] + "\033[0m unknown. Available options are:";
      for (int i=0; i<sizeof(known_styles)/sizeof(string); i++) {
	if (i) error_str += ",";
	error_str += " \033[1;32m" + known_styles[i] + "\033[0m";
      }
      error_str += "\n";
      error->all(FLERR, error_str);
    }
  } else {
    string error_str = "Error: CPDI style unspecified in the method() command.\n";
    for (int i=0; i<sizeof(known_styles)/sizeof(string); i++) {
	if (i) error_str += ",";
	error_str += " \033[1;32m" + known_styles[i] + "\033[0m";
      }
    error_str += "\n";
    error->all(FLERR, error_str);
  }

  if (universe->me == 0)
    cout << "Using CPDI-" << known_styles[style] << endl;
}

void ULCPDI::compute_grid_weight_functions_and_gradients()
{
  if (!update_wf) return;

  if (domain->dimension !=2)
    {
      error->all(FLERR, "Error: ULCPDI is only 2D....\n");
    }

  bigint nsolids, np_local, nnodes, nc;

  nsolids = domain->solids.size();

  if (nsolids)
    {
    for (int isolid=0; isolid<nsolids; isolid++)
      {
	Solid *s = domain->solids[isolid];

	np_local = s->np_local;
	nc = s->nc;
	nnodes = s->grid->nnodes_local + s->grid->nnodes_ghost;

	vector<int> *numneigh_pn = &s->numneigh_pn;
	vector<int> *numneigh_np = &s->numneigh_np;

	vector<vector<int>> *neigh_pn = &s->neigh_pn;
	vector<vector<int>> *neigh_np = &s->neigh_np;

	vector<vector< double >> *wf_pn = &s->wf_pn;
	vector<vector< double >> *wf_pn_corners = &s->wf_pn_corners;
	vector<vector< double >> *wf_np = &s->wf_np;

	vector<vector< Vector3d >> *wfd_pn = &s->wfd_pn;
	vector<vector< Vector3d >> *wfd_np = &s->wfd_np;

	vector<Vector3d> *xp  = &s->x;
	vector<Vector3d> *xpc = &s->xpc;
	vector<Vector3d> *xn  = &s->grid->x0;
	vector<Vector3d> *rp  = &s->rp;

	double inv_cellsize          = 1.0 / s->grid->cellsize;
	vector<array<int, 3>> *ntype = &s->grid->ntype;

	double wf;
	double phi[3];
	Vector3d r, wfd;
	vector<Vector3d> xcorner(nc);
	vector<double> wfc(nc, 0);

	double a, b, inv_Vp, alpha_over_Vp, sixVp;

	for (int in = 0; in < nnodes; in++)
	  {
	    (*neigh_np)[in].clear();
	    (*numneigh_np)[in] = 0;
	    (*wf_np)[in].clear();
	    (*wfd_np)[in].clear();
	  }

	if (np_local && nnodes)
	  {
	  for (int ip = 0; ip < np_local; ip++)
	    {
	      (*neigh_pn)[ip].clear();
	      (*numneigh_pn)[ip] = 0;
	      (*wf_pn)[ip].clear();
	      for(int ic=0; ic<nc; ic++) (*wf_pn_corners)[nc*ip+ic].clear();
	      (*wfd_pn)[ip].clear();

	      // Calculate what nodes the corner of Omega_p will interact with:
	      int nx = s->grid->nx;
	      int ny = s->grid->ny;
	      int nz = s->grid->nz;

	      vector<int> n_neigh;
	      int m;


	      if (style==0)
		{ //CPDI-R4
		  // Calculate the coordinates of the particle domain's corners:
		  if (domain->dimension == 1)
		    {
		      xcorner[0] = (*xp)[ip] - (*rp)[ip];
		      xcorner[1] = (*xp)[ip] + (*rp)[ip];
		    }

		if (domain->dimension == 2)
		  {
		    xcorner[0] = (*xp)[ip] - (*rp)[2*ip] - (*rp)[2*ip+1];
		    xcorner[1] = (*xp)[ip] + (*rp)[2*ip] - (*rp)[2*ip+1];
		    xcorner[2] = (*xp)[ip] + (*rp)[2*ip] + (*rp)[2*ip+1];
		    xcorner[3] = (*xp)[ip] - (*rp)[2*ip] + (*rp)[2*ip+1];
		  }

		if (domain->dimension == 3)
		  {
		    error->all(FLERR, "Unsupported!\n");
		  }
		}

	      for (int ic=0; ic<nc; ic++)
		{ // Do this for all corners

		  if (style==1)
		    { // CPDI-Q4
		      xcorner[ic][0] = (*xpc)[nc*ip+ic][0];
		      xcorner[ic][1] = (*xpc)[nc*ip+ic][1];
		      xcorner[ic][2] = (*xpc)[nc*ip+ic][2];
		    }

		  int i0, j0, k0;

		  if (update->shape_function == Update::ShapeFunctions::LINEAR)
		    {
		      i0 = (int) ((xcorner[ic][0] - domain->boxlo[0])*inv_cellsize);
		      j0 = (int) ((xcorner[ic][1] - domain->boxlo[1])*inv_cellsize);
		      k0 = (int) ((xcorner[ic][2] - domain->boxlo[2])*inv_cellsize);

		      m = 2;

		    }
		  else if (update->shape_function == Update::ShapeFunctions::BERNSTEIN)
		    {
		      i0 = 2*(int) ((xcorner[ic][0] - domain->boxlo[0])*inv_cellsize);
		      j0 = 2*(int) ((xcorner[ic][1] - domain->boxlo[1])*inv_cellsize);
		      k0 = 2*(int) ((xcorner[ic][2] - domain->boxlo[2])*inv_cellsize);

		      if ((i0 >= 1) && (i0 % 2 != 0)) i0--;
		      if ((j0 >= 1) && (j0 % 2 != 0)) j0--;
		      if (nz>1) if ((k0 >= 1) && (k0 % 2 != 0)) k0--;

		      m = 3;

		    }
		  else
		    { // cubic and quadratic B-splines
		      i0 = (int) ((xcorner[ic][0] - domain->boxlo[0])*inv_cellsize - 1);
		      j0 = (int) ((xcorner[ic][1] - domain->boxlo[1])*inv_cellsize - 1);
		      k0 = (int) ((xcorner[ic][2] - domain->boxlo[2])*inv_cellsize - 1);

		      m = 4;

		    }

		  for(int i = i0; i < i0 + m; i++)
		    {
		      if (ny > 1)
			{
			  for(int j = j0; j < j0 + m; j++)
			    {
			      if (nz > 1)
				{
				  for(int k = k0; k < k0 + m; k++)
				    {
				      int n = nz * ny * i + nz * j + k;
				      if (n < nnodes)
					n_neigh.push_back(n);
				    }
				}
			      else
				{
				  int n = ny * i + j;
				  if (n < nnodes)
				    {
				      n_neigh.push_back(n);
				    }
				}
			    }
			}
		      else
			{
			  if (i < nnodes)
			    n_neigh.push_back(i);
			}
		    }
		}

	      // Keep only unique values of in in n_neigh:
	      vector<int>::iterator it;
	      sort(n_neigh.begin(), n_neigh.end());           // Sort all values
	      it = unique (n_neigh.begin(), n_neigh.end());   // Remove indentical consecutive values
	      n_neigh.resize( distance(n_neigh.begin(),it) ); // Resize vector.

	      // cout << "[";
	      // for (auto ii: n_neigh)
	      //   cout << ii << ' ';
	      // cout << "]\n";

	      inv_Vp = 1.0 / s->vol[ip];

	      if (style==1)
		{
		  a = (xcorner[3][0]-xcorner[0][0])*(xcorner[1][1]-xcorner[2][1])
		    - (xcorner[1][0]-xcorner[2][0])*(xcorner[3][1]-xcorner[0][1]);

		  b = (xcorner[2][0]-xcorner[3][0])*(xcorner[0][1]-xcorner[1][1])
		    - (xcorner[0][0]-xcorner[1][0])*(xcorner[2][1]-xcorner[3][1]);

		  alpha_over_Vp = 0.0417 * inv_Vp;
		  sixVp = 6 * s->vol[ip];
		}

	      //for (int in=0; in<nnodes; in++) {
	      for (auto in: n_neigh)
		{
		  wf = 0;

		  for(int ic=0; ic<nc; ic++)
		    {
		      // Calculate the distance between each pair of particle/node:
		      r = (xcorner[ic] - (*xn)[in]) * inv_cellsize;

		      phi[0] = basis_function(r[0], (*ntype)[in][0]);
		      phi[1] = basis_function(r[1], (*ntype)[in][1]);
		      if (domain->dimension == 3) phi[2] = basis_function(r[2], (*ntype)[in][2]);
		      else phi[2] = 1;

		      wfc[ic] = phi[0]*phi[1]*phi[2]; // Shape function of the corner node

		      if (style==0 && wfc[ic] > 1.0e-12) wf += wfc[ic];
		    }

		  if (style==0) wf *= 0.25;
		  if (style==1)
		    {
		      wf = alpha_over_Vp*((sixVp - a - b)*wfc[0]
					  + (sixVp - a + b)*wfc[1]
					  + (sixVp + a + b)*wfc[2]
					  + (sixVp + a - b)*wfc[3]);
		    }

		  if (wf > 1.0e-12)
		    {
		      if (style==0)
			{
			  wfd[0] = (wfc[0] - wfc[2]) * ((*rp)[domain->dimension*ip][1] - (*rp)[domain->dimension*ip+1][1])
			    + (wfc[1] - wfc[3]) * ((*rp)[domain->dimension*ip][1] + (*rp)[domain->dimension*ip+1][1]);

			  wfd[1] = (wfc[0] - wfc[2]) * ((*rp)[domain->dimension*ip+1][0] - (*rp)[domain->dimension*ip][0])
			    - (wfc[1] - wfc[3]) * ((*rp)[domain->dimension*ip][0] + (*rp)[domain->dimension*ip+1][0]);

			  wfd[2] = 0;

			  wfd *= inv_Vp;
			}

		      if (style==1)
			{
			  wfd[0] = wfc[0] * (xcorner[1][1]-xcorner[3][1])
			    + wfc[1] * (xcorner[2][1]-xcorner[0][1])
			    + wfc[2] * (xcorner[3][1]-xcorner[1][1])
			    + wfc[3] * (xcorner[0][1]-xcorner[2][1]);

			  wfd[1] = wfc[0] * (xcorner[3][0]-xcorner[1][0])
			    + wfc[1] * (xcorner[0][0]-xcorner[2][0])
			    + wfc[2] * (xcorner[1][0]-xcorner[3][0])
			    + wfc[3] * (xcorner[2][0]-xcorner[0][0]);

			  wfd[2] = 0;

			  wfd *= 0.5*inv_Vp;
			  for(int ic=0; ic<nc; ic++) (*wf_pn_corners)[nc*ip+ic].push_back(wfc[ic]);
			}

		      (*neigh_pn)[ip].push_back(in);
		      (*neigh_np)[in].push_back(ip);
		      (*numneigh_pn)[ip]++;
		      (*numneigh_np)[in]++;

		      (*wf_pn)[ip].push_back(wf);
		      (*wf_np)[in].push_back(wf);
		      (*wfd_pn)[ip].push_back(wfd);
		      (*wfd_np)[in].push_back(wfd);
		    }
		}
	    }
	  }
	if (method_type.compare("APIC") == 0) s->compute_inertia_tensor();
      }
    }
}

void ULCPDI::particles_to_grid()
{
  bool grid_reset = false; // Indicate if the grid quantities have to be reset
  for (int isolid=0; isolid<domain->solids.size(); isolid++){

    if (isolid == 0) grid_reset = true;
    else grid_reset = false;

    domain->solids[isolid]->compute_mass_nodes(grid_reset);
  }
  for (int isolid=0; isolid<domain->solids.size(); isolid++){

    if (isolid == 0) grid_reset = true;
    else grid_reset = false;

    if (method_type.compare("APIC") == 0) domain->solids[isolid]->compute_velocity_nodes_APIC(grid_reset);
    else domain->solids[isolid]->compute_velocity_nodes(grid_reset);
    domain->solids[isolid]->compute_external_and_internal_forces_nodes_UL(grid_reset);
    /*compute_thermal_energy_nodes();*/
  }
}

void ULCPDI::particles_to_grid_USF_1()
{
  bool grid_reset = false; // Indicate if the grid quantities have to be reset
  for (int isolid=0; isolid<domain->solids.size(); isolid++){

    if (isolid == 0) grid_reset = true;
    else grid_reset = false;

    domain->solids[isolid]->compute_mass_nodes(grid_reset);
  }
  for (int isolid=0; isolid<domain->solids.size(); isolid++){

    if (isolid == 0) grid_reset = true;
    else grid_reset = false;

    if (method_type.compare("APIC") == 0)
      domain->solids[isolid]->compute_velocity_nodes_APIC(grid_reset);
    else
      domain->solids[isolid]->compute_velocity_nodes(grid_reset);
  }
}

void ULCPDI::particles_to_grid_USF_2()
{
  bool grid_reset = false; // Indicate if the grid quantities have to be reset
  for (int isolid=0; isolid<domain->solids.size(); isolid++){

    if (isolid == 0) grid_reset = true;
    else grid_reset = false;

    domain->solids[isolid]->compute_external_and_internal_forces_nodes_UL(grid_reset);
    /*compute_thermal_energy_nodes();*/
  }
}

void ULCPDI::update_grid_state()
{
  domain->grid->update_grid_velocities();
}

void ULCPDI::grid_to_points()
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    domain->solids[isolid]->compute_particle_velocities_and_positions();
    domain->solids[isolid]->compute_particle_acceleration();
  }
}

void ULCPDI::advance_particles()
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    domain->solids[isolid]->update_particle_velocities(FLIP);
  }
}

void ULCPDI::velocities_to_grid()
{
  bool grid_reset = false; // Indicate if the grid quantities have to be reset
  for (int isolid=0; isolid<domain->solids.size(); isolid++){

    if (isolid == 0) grid_reset = true;
    else grid_reset = false;

    if (method_type.compare("APIC") != 0) { 
      //domain->solids[isolid]->compute_mass_nodes(grid_reset);
      domain->solids[isolid]->compute_velocity_nodes(grid_reset);
    }
    // domain->solids[isolid]->grid->update_grid_positions();
  }
}

void ULCPDI::compute_rate_deformation_gradient(bool doublemapping) {
  for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
    if (update->sub_method_type == Update::SubMethodType::APIC)
      domain->solids[isolid]->compute_rate_deformation_gradient_UL_APIC(doublemapping);
    else
      domain->solids[isolid]->compute_rate_deformation_gradient_UL(doublemapping);
  }
}

void ULCPDI::update_deformation_gradient()
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    domain->solids[isolid]->update_deformation_gradient();
    domain->solids[isolid]->update_particle_domain();
  }
}

void ULCPDI::update_stress(bool doublemapping)
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    domain->solids[isolid]->update_stress();
  }
}

void ULCPDI::adjust_dt()
{
  update->update_time();
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

void ULCPDI::reset()
{
  int np_local;

  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    domain->solids[isolid]->dtCFL = 1.0e22;
    np_local = domain->solids[isolid]->np_local;
    for (int ip = 0; ip < np_local; ip++) domain->solids[isolid]->mbp[ip] = Vector3d();
  }
}


void ULCPDI::exchange_particles()
{
  cout << "Error: ULCPDI::exchange_particles() not coded\n";
  error->all(FLERR, "");
}
