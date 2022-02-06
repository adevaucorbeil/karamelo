/* -*- c++ -*- ----------------------------------------------------------
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

#include <method.h>

#include <solid.h>
#include <update.h>
#include <domain.h>
#include <error.h>

using namespace std;

Method::Method(MPM *mpm) : Pointers(mpm)
{
  is_TL = false;
  is_CPDI = false;
  ge = false;
  temp = false;
}

bool Method::apic()
{
  switch (update->sub_method_type)
  {
  case Update::SubMethodType::APIC:
  case Update::SubMethodType::AFLIP:
  case Update::SubMethodType::ASFLIP:
  case Update::SubMethodType::MLS:
    return true;
  }

  return false;
}

void Method::compute_mass_nodes(Solid &solid, int in, int ip, double wf)
{
  if (should_compute_mass_nodes() && !solid.grid->rigid.at(in) || solid.mat->rigid)
    solid.grid->mass.at(in) += wf*solid.mass.at(ip);
}

void Method::compute_velocity_nodes(Solid &solid, int in, int ip, double wf)
{
  if (solid.grid->rigid.at(in) && !solid.mat->rigid)
    return;

  if (double node_mass = solid.grid->mass.at(in))
  {
    double wf_mass = wf*solid.mass.at(ip);
    const Vector3d &dx = solid.grid->x0.at(in) - solid.x.at(ip);
    
    Vector3d vtemp = solid.v.at(ip);

    if (apic())
    {
      if (is_TL)
        vtemp += solid.Fdot.at(ip)*(solid.grid->x0.at(in) - solid.x0.at(ip));
      else
        vtemp += solid.L.at(ip)*dx;

      solid.grid->v.at(in) += vtemp/node_mass;
    }
    else
    {
      if (update->method->ge)
        vtemp += solid.L.at(ip)*dx;

      solid.grid->v.at(in) += wf_mass*vtemp/node_mass;
    }

    if (solid.grid->rigid.at(in))
      solid.grid->mb.at(in) += wf_mass*solid.v_update.at(ip)/node_mass;
  }
}

void Method::compute_force_nodes(Solid &solid, int in, int ip, double wf, const Vector3d &wfd)
{
  compute_internal_force_nodes(solid, in, ip, wf, wfd);
  //if (TL)
  //{
  //  grid->f.at(in) -= vol0PK1.at(ip)*wfd;

  //  if (domain->axisymmetric)
  //    grid->f.at(in)[0] -= vol0PK1.at(ip)(2, 2)*wf/x0.at(ip)[0];
  //}
  //else
  //{
  //  if (MLS)
  //    grid->f.at(in) -= vol.at(ip)*wf*sigma.at(ip)*Di*
  //                      (grid->x0.at(in) - (is_TL? x0: x).at(ip));
  //  else
  //    grid->f.at(in) -= vol.at(ip)*sigma.at(ip)*wfd;

  //  if (domain->axisymmetric)
  //    grid->f.at(in)[0] -= vol.at(ip)*sigma.at(ip)(2, 2)*wf/x.at(ip)[0];
  //}

  if (!solid.grid->rigid.at(in))
    solid.grid->mb.at(in) += wf*solid.mbp.at(ip);
}

void Method::compute_temperature_nodes(Solid &solid, int in, int ip, double wf)
{
  if (double node_mass = solid.grid->mass.at(in))
    solid.grid->T.at(in) += wf*solid.mass.at(ip)*solid.T.at(ip)/node_mass;
}

void Method::compute_temperature_driving_force_nodes(Solid &solid, int in, int ip, double wf, const Vector3d &wfd)
{
  if (solid.domain->axisymmetric)
  {
    error->one(FLERR, "Temperature and axisymmetric not yet supported.\n");
  }

  if (solid.grid->mass.at(in))
    solid.grid->Qext.at(in) += wf*solid.gamma.at(ip);

  solid.grid->Qint.at(in) += wfd.dot(solid.q.at(ip));
}
