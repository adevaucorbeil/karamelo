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

#include <tlmpm.h>
#include <basis_functions.h>
#include <domain.h>
#include <error.h>
#include <grid.h>
#include <input.h>
#include <method.h>
#include <solid.h>
#include <universe.h>
#include <update.h>
#include <var.h>
#include <matrix.h>
#include <algorithm>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include <mpm_math.h>

using namespace std;

TLMPM::TLMPM(MPM *mpm):
  Method(mpm)
{
  // cout << "In TLMPM::TLMPM()" << endl;

  update->PIC_FLIP = 0.99;
  is_TL = true;

  // Default base function (linear):
  basis_function = &BasisFunction::linear;
  derivative_basis_function = &BasisFunction::derivative_linear;
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

void TLMPM::compute_grid_weight_functions_and_gradients(Solid &solid, int ip)
{
  if (!update->atimestep)
    Method::compute_grid_weight_functions_and_gradients(solid, ip);
}

vector<Grid *> TLMPM::grids()
{
  vector<Grid *> grids;
  grids.reserve(domain->solids.size());

  for (Solid *solid: domain->solids)
    grids.push_back(solid->grid);

  return grids;
}

void TLMPM::reset_mass_nodes(Grid &grid, int in)
{
  if (!update->atimestep)
    Method::reset_mass_nodes(grid, in);
}

void TLMPM::compute_mass_nodes(Solid &solid, int ip)
{
  if (!update->atimestep)
    Method::compute_mass_nodes(solid, ip);
}

void TLMPM::compute_internal_force_nodes(Solid &solid, int ip)
{
  for (int i = 0; i < solid.neigh_n.at(ip).size(); i++)
  {
    int in = solid.neigh_n.at(ip).at(i);
    double wf = solid.wf.at(ip).at(i);
    const Vector3d &wfd = solid.wfd.at(ip).at(i);

    Vector3d &f = solid.grid->f.at(in);
    const Matrix3d &vol0PK1 = solid.vol0PK1.at(ip);
    const Vector3d &x0 = solid.x0.at(ip);

    if (update->sub_method_type == Update::SubMethodType::MLS)
      f -= vol0PK1*wf*solid.Di*(solid.grid->x0.at(in) - x0);
    else
      f -= vol0PK1*wfd;

    if (domain->axisymmetric)
      f[0] -= vol0PK1(2, 2)*wf/x0[0];
  }
}

void TLMPM::update_grid_positions(Grid &grid, int in)
{
  grid.x.at(in) += update->dt*grid.v.at(in);
}

vector<Matrix3d> &TLMPM::get_gradients(Solid &solid)
{
  return solid.Fdot;
}

void TLMPM::update_deformation_gradient_matrix(Solid &solid, int ip)
{
  solid.F.at(ip) += update->dt*solid.Fdot.at(ip);
}

void TLMPM::update_velocity_gradient_matrix(Solid &solid, int ip)
{
  // polar decomposition of the deformation gradient, F = R*U
  if (!MPM_Math::PolDec(solid.F.at(ip), solid.R.at(ip)))
  {
    cout << "Polar decomposition of deformation gradient failed for particle " << ip << ".\n";
    cout << "F:" << endl << solid.F.at(ip) << endl;
    cout << "timestep:" << endl << update->ntimestep << endl;
    error->one(FLERR, "");
  }

  // In TLMPM. L is computed from Fdot:
  solid.L.at(ip) = solid.Fdot.at(ip)*solid.Finv.at(ip);
  solid.D.at(ip) = 0.5*(solid.R.at(ip).transpose()*(solid.L.at(ip) + solid.L.at(ip).transpose())*solid.R.at(ip));
}
