/* ----------------------------------------------------------------------
 *
*                   ***       Karamelo       ***
*              Parallel Material Point Method Simulator
 *
*Copyright (2019) Alban de Vaucorbeil, alban.devaucorbeil@monash.edu
*Materials Science and Engineering, Monash University
*Clayton VIC 3800, Australia

*This software is distributed under the GNU General Public License.
 *
*----------------------------------------------------------------------- */

#include <solid.h>
#include <domain.h>
#include <error.h>
#include <input.h>
#include <material.h>
#include <memory.h>
#include <method.h>
#include <mpm.h>
#include <mpm_math.h>
#include <universe.h>
#include <update.h>
#include <var.h>
#include <matrix.h>
#include <math.h>
#include <mpi.h>
#include <string>
#include <vector>

using namespace std;

using namespace MPM_Math;


#define SQRT_3_OVER_2 1.224744871 // sqrt(3.0/2.0)
#define FOUR_THIRD 1.333333333333333333333333333333333333333



vector<string> split(string s, string delimiter)
{
  size_t pos_start = 0, pos_end, delim_len = delimiter.length();
  string token;
  vector<string> res;

  while ((pos_end = s.find(delimiter, pos_start)) != string::npos)
  {
    token = s.substr(pos_start, pos_end - pos_start);
    pos_start = pos_end + delim_len;
    res.push_back(token);
  }

  res.push_back(s.substr(pos_start));
  return res;
}

Solid::Solid(MPM *mpm, vector<string> args): Pointers(mpm)
{
  // Check that a method is available:
  if (update->method == nullptr)
  {
    error->all(FLERR, "Error: a method should be defined before creating a solid!\n");
  }

  if (args.size() < 2)
  {
    string error_str = "Error: solid command not enough arguments. ";
    for (auto &x : usage)
      error_str += x.second;
    error->all(FLERR, error_str);
  }

  id = args[0];

  if (universe->me == 0)
    cout << "Creating new solid with ID: " << id << endl;

  method_type = update->method_type;

  np = 0;

  if (update->method->is_CPDI)
  {
    nc = pow(2, domain->dimension);
  }
  else
    nc = 0;

  mat = nullptr;

  if (update->method->is_TL)
  {
    is_TL = true;
    grid = new Grid(mpm);
  }
  else
  {
    is_TL = false;
    grid = domain->grid;
  }

  if (update->sub_method_type == Update::SubMethodType::APIC ||
      update->sub_method_type == Update::SubMethodType::AFLIP ||
      update->sub_method_type == Update::SubMethodType::ASFLIP)
  {
    apic = true;
  }
  else
  {
    apic = false;
  }

  dtCFL = 1.0e22;
  max_p_wave_speed = 0;
  vtot = 0;
  mtot = 0;
  comm_n = 50; // Number of double to pack for particle exchange between CPUs.


  if (args[1].compare("restart") == 0)
  {
// If the keyword restart, we are expecting to have read_restart()
// launched right after.
    return;
  }

  if (usage.find(args[1]) == usage.end())
  {
    string error_str = "Error, keyword \033[1;31m" + args[1] + "\033[0m unknown!\n";
    for (auto &x : usage)
      error_str += x.second;
    error->all(FLERR, error_str);
  }

  if (args.size() < Nargs.find(args[1])->second)
  {
    error->all(FLERR, "Error: not enough arguments.\n"
               + usage.find(args[1])->second);
  }

  if (args[1].compare("region") == 0)
  {
    // Set material, cellsize, and initial temperature:
    options(&args, args.begin() + 4);

    // Create particles:
    populate(args);
  }
  else if (args[1].compare("mesh") == 0)
  {
    // Set material and cellsize and initial temperature:
    options(&args, args.begin() + 3);

    read_mesh(args[2]);
  }
  else if (args[1].compare("file") == 0)
  {
    // Set material and cellsize and initial temperature:
    options(&args, args.begin() + 3);

    read_file(args[2]);
  }

  if (update->method->temp)
    comm_n = 54; // Number of double to pack for particle exchange between CPUs.
  else
    comm_n = 49;
}

Solid::~Solid()
{
  if (is_TL) delete grid;
}

void Solid::init()
{
  if (universe->me == 0)
  {
    cout << "Bounds for " << id << ":\n";
    cout << "xlo xhi: " << solidlo[0] << " " << solidhi[0] << endl;
    cout << "ylo yhi: " << solidlo[1] << " " << solidhi[1] << endl;
    cout << "zlo zhi: " << solidlo[2] << " " << solidhi[2] << endl;
  }

  // Calculate total volume:
  double vtot_local = 0;
  double mtot_local = 0;
  for (int ip = 0; ip<np_local; ip++)
  {
    vtot_local += vol.at(ip);
    mtot_local += mass.at(ip);
  }

  MPI_Allreduce(&vtot_local, &vtot, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  MPI_Allreduce(&mtot_local, &mtot, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);

  if (universe->me == 0)
  {
    cout << "Solid " << id << " total volume = " << vtot << endl;
    cout << "Solid " << id << " total mass = " << mtot << endl;
  }

  if (grid->nnodes == 0) grid->init(solidlo, solidhi);

  if (np == 0)
  {
    error->one(FLERR, "Error: solid does not have any particles.\n");
  }
}

void Solid::options(vector<string> *args, vector<string>::iterator it)
{
  // cout << "In solid::options()" << endl;
  if (args->end() < it + 3)
  {
    error->all(FLERR, "Error: not enough arguments.\n");
  }
  if (args->end() > it)
  {
    int iMat = material->find_material(*it);

    if (iMat == -1)
    {
      cout << "Error: could not find material named " << *it << endl;
      error->all(FLERR, "\n");
    }

    mat = &material->materials[iMat]; // point mat to the right material

    it++;

    grid->setup(*it); // set the grid cellsize

    it++;
    T0 = input->parsev(*it); // set initial temperature

    it++;

    if (it != args->end())
    {
      string error_str = "Error: too many arguments\n";
      for (auto &x : usage)
        error_str += x.second;
    }
  }
}

void Solid::grow(int nparticles)
{
  ptag.resize(nparticles);
  x0.resize(nparticles);
  x.resize(nparticles);

  if (method_type.compare("tlcpdi") == 0
      || method_type.compare("ulcpdi") == 0)
  {

    if (update->method->style == 0)
    { // CPDI-R4
      rp0.resize(domain->dimension*nparticles);
      rp.resize(domain->dimension*nparticles);
    }
    if (update->method->style == 1)
    { // CPDI-Q4
      xpc0.resize(nc*nparticles);
      xpc.resize(nc*nparticles);
    }
  }

  if (method_type.compare("tlcpdi2") == 0
      || method_type.compare("ulcpdi2") == 0)
  {
    xpc0.resize(nparticles);
    xpc.resize(nparticles);
  }

  v.resize(nparticles);
  v_update.resize(nparticles);
  a.resize(nparticles);
  mbp.resize(nparticles);
  f.resize(nparticles);
  sigma.resize(nparticles);
  strain_el.resize(nparticles);
  vol0PK1.resize(nparticles);
  L.resize(nparticles);
  F.resize(nparticles);
  R.resize(nparticles);
  D.resize(nparticles);
  Finv.resize(nparticles);
  Fdot.resize(nparticles);
  vol0.resize(nparticles);
  vol.resize(nparticles);
  rho0.resize(nparticles);
  rho.resize(nparticles);
  mass.resize(nparticles);
  eff_plastic_strain.resize(nparticles);
  eff_plastic_strain_rate.resize(nparticles);
  damage.resize(nparticles);
  damage_init.resize(nparticles);
  ienergy.resize(nparticles);
  mask.resize(nparticles);
  J.resize(nparticles);

  if (mat->cp != 0)
  {
    T.resize(nparticles);
    gamma.resize(nparticles);
    q.resize(nparticles);
  }
}

void Solid::compute_mass_nodes(int in, int ip, double wf)
{
  if (!grid->rigid.at(in) || mat->rigid)
    grid->mass.at(in) += wf*mass.at(ip);
}

void Solid::compute_velocity_nodes(int in, int ip, double wf, bool APIC)
{
  if ((!grid->rigid.at(in) || mat->rigid) && grid->mass.at(in))
  {
    double wf_mass = wf*mass.at(ip);
    
    Vector3d vtemp = v.at(ip);

    if (APIC)
    {
      if (is_TL)
        vtemp += Fdot.at(ip)*(grid->x0.at(in) - x0.at(ip));
      else
        vtemp += L.at(ip)*(grid->x0.at(in) - x.at(ip));

      grid->v.at(in) += vtemp/grid->mass.at(in);
    }
    else
    {
      if (update->method->ge)
        vtemp += L.at(ip)*(grid->x0.at(in) - x.at(ip));

      grid->v.at(in) += wf_mass*vtemp/grid->mass.at(in);
    }

    if (grid->rigid.at(in))
      grid->mb.at(in) += wf_mass*v_update.at(ip)/grid->mass.at(in);
  }
}

void Solid::compute_force_nodes(int in, int ip, double wf, const Vector3d &wfd, bool TL, bool MLS)
{
  if (TL)
  {
    grid->f.at(in) -= vol0PK1.at(ip)*wfd;

    if (domain->axisymmetric)
      grid->f.at(in)[0] -= vol0PK1.at(ip)(2, 2)*wf/x0.at(ip)[0];
  }
  else
  {
    if (MLS)
      grid->f.at(in) -= vol.at(ip)*wf*sigma.at(ip)*Di*
                        (grid->x0.at(in) - (is_TL? x0: x).at(ip));
    else
      grid->f.at(in) -= vol.at(ip)*sigma.at(ip)*wfd;

    if (domain->axisymmetric)
      grid->f.at(in)[0] -= vol.at(ip)*sigma.at(ip)(2, 2)*wf/x.at(ip)[0];
  }

  if (!grid->rigid.at(in))
    grid->mb.at(in) += wf*mbp.at(ip);
}

void Solid::compute_temperature_nodes(int in, int ip, double wf)
{
  if (grid->mass.at(in))
  {
    grid->T.at(in) += wf*mass.at(ip)*T.at(ip)/grid->mass.at(in);
  }
}

void Solid::compute_temperature_driving_force_nodes(int in, int ip, double wf, const Vector3d &wfd)
{
  if (domain->axisymmetric)
  {
    error->one(FLERR, "Temperature and axisymmetric not yet supported.\n");
  }

  if (grid->mass.at(in))
  {
    grid->Qext.at(in) += wf*gamma.at(ip);
  }
  grid->Qint.at(in) += wfd.dot(q.at(ip));
}

void Solid::compute_velocity_acceleration(int in, int ip, double wf)
{
  v_update.at(ip) += wf*grid->v_update.at(in);

  if (!mat->rigid)
  {
    const Vector3d &delta_a = wf*(grid->v_update.at(in) - grid->v.at(in))/update->dt;
    a.at(ip) += delta_a;
    f.at(ip) += mass.at(ip)*delta_a;
  }
}

void Solid::compute_particle_temperature(int in, int ip, double wf)
{
  T.at(ip) += wf*grid->T_update.at(in);
}

void Solid::compute_heat_flux(int in, int ip, const Vector3d &wfd, bool doublemapping)
{
  const vector<double> &Tn = doublemapping? grid->T: grid->T_update;

  q.at(ip) -= wfd*Tn.at(in)*(is_TL? vol0: vol).at(ip)*mat->invcp*mat->kappa;
}

void Solid::compute_rate_deformation_gradient(int in, int ip, double wf, const Vector3d &wfd,
                                              bool doublemapping, bool TL, bool APIC)
{
  vector<Matrix3d> &gradients = TL? Fdot: L;
  const vector<Vector3d> &vn = doublemapping? grid->v: grid->v_update;

  if (APIC)
  {
    const Vector3d &dx = grid->x.at(in) - grid->x0.at(in);

    for (int j = 0; j < domain->dimension; j++)
      for (int k = 0; k < domain->dimension; k++)
        gradients.at(ip)(j, k) += vn.at(in)[j]*dx[k]*wf;
  }
  else
    for (int j = 0; j < domain->dimension; j++)
      for (int k = 0; k < domain->dimension; k++)
        gradients.at(ip)(j, k) += vn.at(in)[j]*wfd[k];

  if (domain->dimension == 2 && domain->axisymmetric)
    gradients.at(ip)(2, 2) += vn.at(in)[0]*wf/x0.at(ip)[0];
}

void Solid::compute_position_corners()
{
  for (int i = 0; i < neigh_n.size(); i++)
  {
    int in = neigh_n.at(i);
    int ip = neigh_p.at(i);
    
    if (update->method->style == 1)
      for (int ic = 0; ic < nc; ic++)
        xpc.at(nc*ip + ic) += update->dt*wf.at(nc*i + ic)*grid->v_update.at(in);
  }
}

void Solid::reset_velocity_acceleration()
{
  for (Vector3d &v_update: v_update)
    v_update = Vector3d();
  for (Vector3d &a: a)
    a = Vector3d();
  for (Vector3d &f: f)
    f = Vector3d();
}

void Solid::reset_heat_flux()
{
  for (Vector3d &q: q)
    q = Vector3d();
}

void Solid::reset_rate_deformation_gradient(bool TL)
{
  vector<Matrix3d> &gradients = TL? Fdot: L;

  for (Matrix3d &gradient: gradients)
    gradient = Matrix3d();
}

void Solid::update_position()
{
  for (int ip = 0; ip < np_local; ip++)
    x.at(ip) += update->dt*v_update.at(ip);

  if (!is_TL)
    for (int ip = 0; ip < np_local; ip++)
      if (!domain->inside(x.at(ip)))
      {
        cout << "Error: Particle " << ip << " left the domain ("
              << domain->boxlo[0] << "," << domain->boxhi[0] << ","
              << domain->boxlo[1] << "," << domain->boxhi[1] << ","
              << domain->boxlo[2] << "," << domain->boxhi[2] << "):\n"
              << x.at(ip) << endl;

        error->one(FLERR, "");
      }
}

void Solid::update_particle(double alpha, bool positions, bool velocities)
{
  for (int ip = 0; ip < np_local; ip++)
  {
    if (positions)
      x.at(ip) += update->dt*v.at(ip);
    if (velocities)
      v.at(ip) = (1 - alpha)*v_update.at(ip) + alpha*(v.at(ip) + update->dt*a.at(ip));
  }
}

void Solid::update_deformation_gradient()
{
  if (mat->rigid)
    return;

  bool status, nh, vol_cpdi;
  Matrix3d eye = Matrix3d::identity();

  if (mat->type == material->constitutive_model::NEO_HOOKEAN)
    nh = true;
  else
    nh = false;

  if ((method_type.compare("tlcpdi") == 0 ||
      method_type.compare("ulcpdi") == 0) &&
      (update->method->style == 1))
  {
    vol_cpdi = true;
  }
  else
    vol_cpdi = false;

  for (int ip = 0; ip < np_local; ip++)
  {

    if (is_TL)
      F.at(ip) += update->dt*Fdot.at(ip);
    else
      F.at(ip) = (eye + update->dt*L.at(ip))*F.at(ip);

    Finv.at(ip) = inverse(F.at(ip));

    if (vol_cpdi)
    {
      vol.at(ip) = 0.5*(xpc[nc*ip + 0][0]*xpc[nc*ip + 1][1] -
                       xpc[nc*ip + 1][0]*xpc[nc*ip + 0][1] +
                       xpc[nc*ip + 1][0]*xpc[nc*ip + 2][1] -
                       xpc[nc*ip + 2][0]*xpc[nc*ip + 1][1] +
                       xpc[nc*ip + 2][0]*xpc[nc*ip + 3][1] -
                       xpc[nc*ip + 3][0]*xpc[nc*ip + 2][1] +
                       xpc[nc*ip + 3][0]*xpc[nc*ip + 0][1] -
                       xpc[nc*ip + 0][0]*xpc[nc*ip + 3][1]);
      // rho.at(ip) = rho0.at(ip);
      J.at(ip) = vol.at(ip)/vol0.at(ip);
    }
    else
    {
      J.at(ip) = determinant(F.at(ip));
      vol.at(ip) = J.at(ip)*vol0.at(ip);
    }


    if (J.at(ip) <= 0.0 && damage.at(ip) < 1.0)
    {
      cout << "Error: J[" << ptag.at(ip) << "]<=0.0 == " << J.at(ip) << endl;
      cout << "F[" << ptag.at(ip) << "]:" << endl << F.at(ip) << endl;
      cout << "Fdot[" << ptag.at(ip) << "]:" << endl << Fdot.at(ip) << endl;
      cout << "damage[" << ptag.at(ip) << "]:" << endl << damage.at(ip) << endl;
      error->one(FLERR, "");
    }
    rho.at(ip) = rho0.at(ip)/J.at(ip);

    if (!nh)
    {
// Only done if not Neo-Hookean:

      if (is_TL)
      {
        status = PolDec(F.at(ip), R.at(ip)); // polar decomposition of the deformation
                                       // gradient, F = R*U

        // In TLMPM. L is computed from Fdot:
        L.at(ip) = Fdot.at(ip)*Finv.at(ip);
        D.at(ip) = 0.5*(R.at(ip).transpose()*(L.at(ip) + L.at(ip).transpose())*R.at(ip));

        if (!status)
        {
          cout << "Polar decomposition of deformation gradient failed for "
            "particle "
            << ip << ".\n";
          cout << "F:" << endl << F.at(ip) << endl;
          cout << "timestep" << endl << update->ntimestep << endl;
          error->one(FLERR, "");
        }

      }
      else
        D.at(ip) = 0.5*(L.at(ip) + L.at(ip).transpose());
    }

    // strain_increment.at(ip) = update->dt*D.at(ip);
  }
}

void Solid::update_stress()
{
  if (mat->rigid)
    return;

  max_p_wave_speed = 0;
  double flow_stress;
  Matrix3d eye, FinvT, PK1, strain_increment;
  bool lin, nh;

  if (mat->type == material->constitutive_model::LINEAR)
    lin = true;
  else
    lin = false;

  if (mat->type == material->constitutive_model::NEO_HOOKEAN)
    nh = true;
  else
    nh = false;

  eye = Matrix3d::identity();

  if (lin)
  {
    for (int ip = 0; ip < np_local; ip++)
    {
      strain_increment = update->dt*D.at(ip);
      strain_el.at(ip) += strain_increment;
      sigma.at(ip) += 2*mat->G*strain_increment +
        mat->lambda*strain_increment.trace()*eye;

      if (is_TL)
      {
        vol0PK1.at(ip) = vol0.at(ip)*J.at(ip) *
          (R.at(ip)*sigma.at(ip)*R.at(ip).transpose()) *
          Finv.at(ip).transpose();
      }
    }
  }
  else if (nh)
  {
    for (int ip = 0; ip < np_local; ip++)
    {
// Neo-Hookean material:
      FinvT = Finv.at(ip).transpose();
      PK1 = mat->G*(F.at(ip) - FinvT) + mat->lambda*log(J.at(ip))*FinvT;
      vol0PK1.at(ip) = vol0.at(ip)*PK1;
      sigma.at(ip) = 1.0/J.at(ip)*(F.at(ip)*PK1.transpose());

      strain_el.at(ip) =
        0.5*(F.at(ip).transpose()*F.at(ip) - eye); // update->dt*D.at(ip);
    }
  }
  else
  {

    vector<double> pH(np_local, 0);
    vector<double> plastic_strain_increment(np_local, 0);
    vector<Matrix3d> sigma_dev;
    sigma_dev.resize(np_local);
    double tav = 0;

    for (int ip = 0; ip < np_local; ip++)
    {

      if (mat->cp != 0)
      {
        mat->eos->compute_pressure(pH.at(ip), ienergy.at(ip), J.at(ip), rho.at(ip),
                                   damage.at(ip), D.at(ip), grid->cellsize, T.at(ip));
        pH.at(ip) += mat->temp->compute_thermal_pressure(T.at(ip));

        sigma_dev.at(ip) = mat->strength->update_deviatoric_stress(
          sigma.at(ip), D.at(ip), plastic_strain_increment.at(ip),
          eff_plastic_strain.at(ip), eff_plastic_strain_rate.at(ip), damage.at(ip),
          T.at(ip));
      }
      else
      {
        mat->eos->compute_pressure(pH.at(ip), ienergy.at(ip), J.at(ip), rho.at(ip),
                                   damage.at(ip), D.at(ip), grid->cellsize);
        sigma_dev.at(ip) = mat->strength->update_deviatoric_stress(
          sigma.at(ip), D.at(ip), plastic_strain_increment.at(ip),
          eff_plastic_strain.at(ip), eff_plastic_strain_rate.at(ip), damage.at(ip));
      }

      eff_plastic_strain.at(ip) += plastic_strain_increment.at(ip);

      // // compute a characteristic time over which to average the plastic
      // strain

      tav = 1000*grid->cellsize/mat->signal_velocity;

      eff_plastic_strain_rate.at(ip) -=
        eff_plastic_strain_rate.at(ip)*update->dt/tav;
      eff_plastic_strain_rate.at(ip) += plastic_strain_increment.at(ip)/tav;
      eff_plastic_strain_rate.at(ip) = MAX(0.0, eff_plastic_strain_rate.at(ip));

      if (mat->damage != nullptr)
      {
        if (update->method->temp)
        {
          mat->damage->compute_damage(damage_init.at(ip), damage.at(ip), pH.at(ip),
                                      sigma_dev.at(ip), eff_plastic_strain_rate.at(ip),
                                      plastic_strain_increment.at(ip), T.at(ip));
        }
        else
        {
          mat->damage->compute_damage(damage_init.at(ip), damage.at(ip), pH.at(ip),
                                      sigma_dev.at(ip), eff_plastic_strain_rate.at(ip),
                                      plastic_strain_increment.at(ip));
        }
      }

      if (mat->cp != 0)
      {
        flow_stress = SQRT_3_OVER_2*sigma_dev.at(ip).norm();
        mat->temp->compute_heat_source(T.at(ip), gamma.at(ip), flow_stress,
                                       eff_plastic_strain_rate.at(ip));
        if (is_TL)
          gamma.at(ip) *= vol0.at(ip)*mat->invcp;
        else
          gamma.at(ip) *= vol.at(ip)*mat->invcp;
      }

      if (damage.at(ip) == 0 || pH.at(ip) >= 0)
        sigma.at(ip) = -pH.at(ip)*eye + sigma_dev.at(ip);
      else
        sigma.at(ip) = -pH.at(ip)*(1.0 - damage.at(ip))* eye + sigma_dev.at(ip);

      if (damage.at(ip) > 1e-10)
      {
        strain_el.at(ip) =
          (update->dt*D.at(ip).trace() + strain_el.at(ip).trace())/3.0*eye +
          sigma_dev.at(ip)/(mat->G*(1 - damage.at(ip)));
      }
      else
      {
        strain_el.at(ip) =
          (update->dt*D.at(ip).trace() + strain_el.at(ip).trace())/3.0*eye +
          sigma_dev.at(ip)/mat->G;
      }

      if (is_TL)
      {
        vol0PK1.at(ip) = vol0.at(ip)*J.at(ip) *
          (R.at(ip)*sigma.at(ip)*R.at(ip).transpose()) *
          Finv.at(ip).transpose();
      }
    }
  }

  double min_h_ratio = 1.0;

  for (int ip = 0; ip < np_local; ip++)
  {
    if (damage.at(ip) >= 1.0)
      continue;

    max_p_wave_speed =
      MAX(max_p_wave_speed,
          sqrt((mat->K + FOUR_THIRD*mat->G)/rho.at(ip)) +
          MAX(MAX(fabs(v.at(ip)(0)), fabs(v.at(ip)(1))), fabs(v.at(ip)(2))));

    if (std::isnan(max_p_wave_speed))
    {
      cout << "Error: max_p_wave_speed is nan with ip=" << ip
        << ", ptag.at(ip)=" << ptag.at(ip) << ", rho0.at(ip)=" << rho0.at(ip)<< ", rho.at(ip)=" << rho.at(ip)
        << ", K=" << mat->K << ", G=" << mat->G << ", J.at(ip)=" << J.at(ip)
        << endl;
      error->one(FLERR, "");
    }
    else if (max_p_wave_speed < 0.0)
    {
      cout << "Error: max_p_wave_speed= " << max_p_wave_speed
        << " with ip=" << ip << ", rho.at(ip)=" << rho.at(ip) << ", K=" << mat->K
        << ", G=" << mat->G << endl;
      error->one(FLERR, "");
    }

    if (is_TL)
    {
      Matrix3d eigenvalues = F.at(ip);
      eigendecompose(eigenvalues);
      if (/*esF.info()!= Success*/false)
      { // EB: revisit
        min_h_ratio = MIN(min_h_ratio, 1.0);
      }
      else
      {
        min_h_ratio = MIN(min_h_ratio, fabs(eigenvalues(0, 0)));
        min_h_ratio = MIN(min_h_ratio, fabs(eigenvalues(1, 1)));
        min_h_ratio = MIN(min_h_ratio, fabs(eigenvalues(2, 2)));
      }

      if (min_h_ratio == 0)
      {
        cout << "min_h_ratio == 0 with ip=" << ip
          << "F=\n" <<  F.at(ip) << endl
          << "eigenvalues of F:" << eigenvalues(0, 0) << "\t" << eigenvalues(1, 1) << "\t" << eigenvalues(2, 2) << endl;
     //cout << "esF.info()=" << esF.info() << endl;
        error->one(FLERR, "");
      }

      // // dt should also be lower than the inverse of \dot{F}e_i.
      // EigenSolver<Matrix3d> esFdot(Fdot.at(ip), false);
      // if (esFdot.info()!= Success) {
      // 	double lambda = fabs(esFdot.eigenvalues()[0].real());
      // 	lambda = MAX(lambda, fabs(esFdot.eigenvalues()[1].real()));
      // 	lambda = MAX(lambda, fabs(esFdot.eigenvalues()[2].real()));
      // 	dtCFL = MIN(dtCFL, 0.5/lambda);
      // }
    }
  }

  dtCFL = MIN(dtCFL, grid->cellsize*min_h_ratio/max_p_wave_speed);

  if (std::isnan(dtCFL))
  {
    cout << "Error: dtCFL = " << dtCFL << "\n";
    cout << "max_p_wave_speed = " << max_p_wave_speed
      << ", grid->cellsize=" << grid->cellsize << endl;
    error->one(FLERR, "");
  }
}

void Solid::compute_inertia_tensor()
{
  Vector3d dx;

  vector<Vector3d> *pos;
  Matrix3d eye = Matrix3d::identity(), Dtemp;
  double cellsizeSqInv = 1.0/(grid->cellsize*grid->cellsize);

  if (update->shape_function == Update::ShapeFunctions::LINEAR)
  {
    if (is_TL)
      pos = &x0;
    else
    {
      error->all(FLERR, "Shape function not supported for APIC and ULMPM.\n");
    }
    if (np_per_cell == 1)
    {
      Di = 16.0/4.0*cellsizeSqInv*eye;
    }
    else if (np_per_cell == 2)
    {
      Di = 16.0/3.0*cellsizeSqInv*eye;
    }
    else
    {
      error->all(FLERR, "Number of particle per cell not supported with linear "
                 "shape functions and APIC.\n");
    }
  }
  else if (update->shape_function == Update::ShapeFunctions::CUBIC_SPLINE)
  {

    Di = 3.0*cellsizeSqInv*eye;
  }
  else if (update->shape_function ==
           Update::ShapeFunctions::QUADRATIC_SPLINE)
  {

    Di = 4.0*cellsizeSqInv*eye;
  }
  else
  {
    error->all(FLERR, "Shape function not supported for APIC.\n");
  }

  if (domain->dimension == 1)
  {
    Di(1, 1) = 1;
    Di(2, 2) = 1;
  }
  else if (domain->dimension == 2)
  {
    Di(2, 2) = 1;
  }

  // if (is_TL)
  //   pos = &x0;
  // else
  //   pos = &x;

  // if (domain->dimension == 2) {
  //   for (int ip = 0; ip < np_local; ip++) {
  //     Dtemp = Vector3d();
  //     for (int j = 0; j < numneigh_pn.at(ip); j++) {
  //       in = neigh_pn.at(ip)[j];
  //       dx = grid->x0[in] - (*pos).at(ip);
  //       Dtemp(0, 0) += wf_pn.at(ip)[j]*(dx[0]*dx[0]);
  //       Dtemp(0, 1) += wf_pn.at(ip)[j]*(dx[0]*dx[1]);
  //       Dtemp(1, 1) += wf_pn.at(ip)[j]*(dx[1]*dx[1]);
  //     }
  //     Dtemp(1, 0) = Dtemp(0, 1);
  //     Dtemp(2, 2) = 1;
  //     Di.at(ip) = Dtemp.inverse();
  //   }
  // } else if (domain->dimension == 3) {
  //   for (int ip = 0; ip < np_local; ip++) {
  //     Dtemp = Vector3d();
  //     for (int j = 0; j < numneigh_pn.at(ip); j++) {
  //       in = neigh_pn.at(ip)[j];
  //       dx = grid->x0[in] - (*pos).at(ip);
  //       Dtemp(0, 0) += wf_pn.at(ip)[j]*(dx[0]*dx[0]);
  //       Dtemp(0, 1) += wf_pn.at(ip)[j]*(dx[0]*dx[1]);
  //       Dtemp(0, 2) += wf_pn.at(ip)[j]*(dx[0]*dx[2]);
  //       Dtemp(1, 1) += wf_pn.at(ip)[j]*(dx[1]*dx[1]);
  //       Dtemp(1, 2) += wf_pn.at(ip)[j]*(dx[1]*dx[2]);
  //       Dtemp(2, 2) += wf_pn.at(ip)[j]*(dx[2]*dx[2]);
  //     }
  //     Dtemp(1, 0) = Dtemp(0, 1);
  //     Dtemp(2, 1) = Dtemp(1, 2);
  //     Dtemp(2, 0) = Dtemp(0, 2);
  //     Di.at(ip) = Dtemp.inverse();
  //   }
  // }

  // Matrix3d eye;
  // eye.setIdentity();

  // double cellsizeSqInv = 1.0/(grid->cellsize*grid->cellsize);

  // for (int ip = 0; ip < np_local; ip++)
  //   {
  //     if ( form_function.compare("linear") == 0)
  // 	{
  // 	  // If the form function is linear:
  // 	  if ((Di.at(ip)(0,0) != 16.0/4.0*cellsizeSqInv ) || (Di.at(ip)(1,1) != 16.0/4.0*cellsizeSqInv) || (Di.at(ip)(2,2) != 16.0/4.0*cellsizeSqInv ))
  // 	    cout << "2 - Di[" << ip << "]=\n" << Di.at(ip) << "\n and " << 4.0*cellsizeSqInv*eye << endl;
  // 	  // Di.at(ip) = 16.0/4.0*cellsizeSqInv*eye;
  // 	}
  //     else if (form_function.compare("quadratic-spline") == 0)
  // 	{
  // 	  // If the form function is a quadratic spline:
  // 	  if ((Di.at(ip)(0,0) != 4.0*cellsizeSqInv ) || (Di.at(ip)(1,1) != 4.0*cellsizeSqInv) || (Di.at(ip)(2,2) != 4.0*cellsizeSqInv ))
  // 	    cout << "2 - Di[" << ip << "]=\n" << Di.at(ip) << "\n and " << 4.0*cellsizeSqInv*eye << endl;
  // 	}
  //     else if (form_function.compare("cubic-spline") == 0)
  // 	{
  // 	  // If the form function is a quadratic spline:
  // 	  if ((Di.at(ip)(0,0) != 3.0*cellsizeSqInv ) || (Di.at(ip)(1,1) != 3.0*cellsizeSqInv) || (Di.at(ip)(2,2) != 3.0*cellsizeSqInv ))
  // 	    cout << "2 - Di[" << ip << "]=\n" << Di.at(ip) << "\n and " << 3.0*cellsizeSqInv*eye << endl;
  // 	  // If the form function is a cubic spline:
  // 	  // Di.at(ip) = 3.0*cellsizeSqInv*eye;
  // 	}
  //     else if (form_function.compare("Bernstein-quadratic") == 0)
  // 	  if ((Di.at(ip)(0,0) != 12.0*cellsizeSqInv ) || (Di.at(ip)(1,1) != 12.0*cellsizeSqInv) || (Di.at(ip)(2,2) != 12.0*cellsizeSqInv ))
  // 	    cout << "2 - Di[" << ip << "]=\n" << Di.at(ip) << "\n and " << 12.0*cellsizeSqInv*eye << endl;
  //     //Di.at(ip) = 12.0*cellsizeSqInv*eye;
  //     // cout << "2 - Di[" << ip << "]=\n" << Di.at(ip) << endl;
  //   }
}

void Solid::copy_particle(int i, int j)
{
  ptag[j] = ptag[i];
  x0[j] = x0[i];
  x[j] = x[i];
  v[j] = v[i];
  v_update[j] = v_update[i];
  a[j] = a[i];
  mbp[j] = mbp[i];
  f[j] = f[i];
  vol0[j] = vol0[i];
  vol[j] = vol[i];
  rho0[j] = rho0[i];
  rho[j] = rho[i];
  mass[j] = mass[i];
  eff_plastic_strain[j] = eff_plastic_strain[i];
  eff_plastic_strain_rate[j] = eff_plastic_strain_rate[i];
  damage[j] = damage[i];
  damage_init[j] = damage_init[i];
  if (update->method->temp)
  {
    T[j] = T[i];
    gamma[j] = gamma[i];
    q[j] = q[i];
  }
  ienergy[j] = ienergy[i];
  mask[j] = mask[i];
  sigma[j] = sigma[i];
  strain_el[j] = strain_el[i];
  vol0PK1[j] = vol0PK1[i];
  L[j] = L[i];
  F[j] = F[i];
  R[j] = R[i];
  D[j] = D[i];
  Finv[j] = Finv[i];
  Fdot[j] = Fdot[i];
  J[j] = J[i];
  // if (method_type.compare("tlcpdi") == 0 || method_type.compare("ulcpdi") == 0)
  //   {
  //     if (update->method->style == 0)
  // 	{ // CPDI-R4
  // 	  for (int id = 0; id < domain->dimension; id++)
  // 	    {
  // 	      rp0[domain->dimension*j + id] = rp0[domain->dimension*i + id];
  // 	      rp[domain->dimension*j + id]  = rp[domain->dimension*i + id];
  // 	    }
  // 	}

  //     if (update->method->style == 1)
  // 	{ // CPDI-Q4
  // 	  for (int ic = 0; ic < nc; ic++)
  // 	    {
  // 	      xpc0[nc*j + ic] = xpc0[nc*i + ic];
  // 	      xpc[nc*j + ic]  = xpc[nc*i + ic];
  // 	    }
  // 	}
  //   }
}

void Solid::pack_particle(int i, vector<double> &buf)
{
  buf.push_back(ptag[i]);

  buf.push_back(x[i](0));
  buf.push_back(x[i](1));
  buf.push_back(x[i](2));

  buf.push_back(x0[i](0));
  buf.push_back(x0[i](1));
  buf.push_back(x0[i](2));

  buf.push_back(v[i](0));
  buf.push_back(v[i](1));
  buf.push_back(v[i](2));

  buf.push_back(v_update[i](0));
  buf.push_back(v_update[i](1));
  buf.push_back(v_update[i](2));

  buf.push_back(a[i](0));
  buf.push_back(a[i](1));
  buf.push_back(a[i](2));

  buf.push_back(mbp[i](0));
  buf.push_back(mbp[i](1));
  buf.push_back(mbp[i](2));

  buf.push_back(f[i](0));
  buf.push_back(f[i](1));
  buf.push_back(f[i](2));

  buf.push_back(vol0[i]);
  buf.push_back(vol[i]);

  buf.push_back(rho0[i]);
  buf.push_back(rho[i]);

  buf.push_back(mass[i]);
  buf.push_back(eff_plastic_strain[i]);
  buf.push_back(eff_plastic_strain_rate[i]);
  buf.push_back(damage[i]);
  buf.push_back(damage_init[i]);
  if (update->method->temp)
  {
    buf.push_back(T[i]);
    buf.push_back(gamma[i]);
    buf.push_back(q[i][0]);
    buf.push_back(q[i][1]);
    buf.push_back(q[i][2]);
  }
  buf.push_back(ienergy[i]);
  buf.push_back(mask[i]);

  buf.push_back(sigma[i](0, 0));
  buf.push_back(sigma[i](1, 1));
  buf.push_back(sigma[i](2, 2));
  buf.push_back(sigma[i](0, 1));
  buf.push_back(sigma[i](0, 2));
  buf.push_back(sigma[i](1, 2));

  // buf.push_back(vol0PK1[i](0,0));
  // buf.push_back(vol0PK1[i](0,1));
  // buf.push_back(vol0PK1[i](0,2));
  // buf.push_back(vol0PK1[i](1,0));
  // buf.push_back(vol0PK1[i](1,1));
  // buf.push_back(vol0PK1[i](1,2));
  // buf.push_back(vol0PK1[i](2,0));
  // buf.push_back(vol0PK1[i](2,1));
  // buf.push_back(vol0PK1[i](2,2));

  // buf.push_back(L[i](0,0));
  // buf.push_back(L[i](0,1));
  // buf.push_back(L[i](0,2));
  // buf.push_back(L[i](1,0));
  // buf.push_back(L[i](1,1));
  // buf.push_back(L[i](1,2));
  // buf.push_back(L[i](2,0));
  // buf.push_back(L[i](2,1));
  // buf.push_back(L[i](2,2));

  buf.push_back(F[i](0, 0));
  buf.push_back(F[i](0, 1));
  buf.push_back(F[i](0, 2));
  buf.push_back(F[i](1, 0));
  buf.push_back(F[i](1, 1));
  buf.push_back(F[i](1, 2));
  buf.push_back(F[i](2, 0));
  buf.push_back(F[i](2, 1));
  buf.push_back(F[i](2, 2));

  buf.push_back(J[i]);
}

void Solid::unpack_particle(int &i, vector<int> list, vector<double> &buf)
{
  int m;
  for (auto j: list)
  {
    m = j;

    ptag[i] = (tagint)buf[m++];

    x[i](0) = buf[m++];
    x[i](1) = buf[m++];
    x[i](2) = buf[m++];

    x0[i](0) = buf[m++];
    x0[i](1) = buf[m++];
    x0[i](2) = buf[m++];

    v[i](0) = buf[m++];
    v[i](1) = buf[m++];
    v[i](2) = buf[m++];

    v_update[i](0) = buf[m++];
    v_update[i](1) = buf[m++];
    v_update[i](2) = buf[m++];

    a[i](0) = buf[m++];
    a[i](1) = buf[m++];
    a[i](2) = buf[m++];

    mbp[i](0) = buf[m++];
    mbp[i](1) = buf[m++];
    mbp[i](2) = buf[m++];

    f[i](0) = buf[m++];
    f[i](1) = buf[m++];
    f[i](2) = buf[m++];

    vol0[i] = buf[m++];
    vol[i] = buf[m++];

    rho0[i] = buf[m++];
    rho[i] = buf[m++];

    mass[i] = buf[m++];
    eff_plastic_strain[i] = buf[m++];
    eff_plastic_strain_rate[i] = buf[m++];
    damage[i] = buf[m++];
    damage_init[i] = buf[m++];
    if (update->method->temp)
    {
      T[i] = buf[m++];
      gamma[i] = buf[m++];

      q[i][0] = buf[m++];
      q[i][1] = buf[m++];
      q[i][2] = buf[m++];
    }

    ienergy[i] = buf[m++];
    mask[i] = buf[m++];

    sigma[i](0, 0) = buf[m++];
    sigma[i](1, 1) = buf[m++];
    sigma[i](2, 2) = buf[m++];
    sigma[i](0, 1) = buf[m++];
    sigma[i](0, 2) = buf[m++];
    sigma[i](1, 2) = buf[m++];
    sigma[i](1, 0) = sigma[i](0, 1);
    sigma[i](2, 0) = sigma[i](0, 2);
    sigma[i](2, 1) = sigma[i](1, 2);

    // vol0PK1[i](0,0) = buf[m++];
    // vol0PK1[i](0,1) = buf[m++];
    // vol0PK1[i](0,2) = buf[m++];
    // vol0PK1[i](1,0) = buf[m++];
    // vol0PK1[i](1,1) = buf[m++];
    // vol0PK1[i](1,2) = buf[m++];
    // vol0PK1[i](2,0) = buf[m++];
    // vol0PK1[i](2,1) = buf[m++];
    // vol0PK1[i](2,2) = buf[m++];

    // L[i](0,0) = buf[m++];
    // L[i](0,1) = buf[m++];
    // L[i](0,2) = buf[m++];
    // L[i](1,0) = buf[m++];
    // L[i](1,1) = buf[m++];
    // L[i](1,2) = buf[m++];
    // L[i](2,0) = buf[m++];
    // L[i](2,1) = buf[m++];
    // L[i](2,2) = buf[m++];

    F[i](0, 0) = buf[m++];
    F[i](0, 1) = buf[m++];
    F[i](0, 2) = buf[m++];
    F[i](1, 0) = buf[m++];
    F[i](1, 1) = buf[m++];
    F[i](1, 2) = buf[m++];
    F[i](2, 0) = buf[m++];
    F[i](2, 1) = buf[m++];
    F[i](2, 2) = buf[m++];

    J[i] = buf[m++];
    i++;
  }
}

void Solid::populate(vector<string> args)
{
  if (universe->me == 0)
  {
    cout << "Solid delimitated by region ID: " << args[2] << endl;
  }

  // Look for region ID:
  int iregion = domain->find_region(args[2]);
  if (iregion == -1)
  {
    error->all(FLERR, "Error: region ID " + args[2] + " not does not exist.\n");
  }

  if (domain->created==false)
  {
    error->all(FLERR, "The domain must be created before any solids can (create_domain(...)).");
  }

  double *sublo = domain->sublo;
  double *subhi = domain->subhi;

  vector<double> limits = domain->regions[iregion]->limits();

  solidlo[0] = limits[0];
  solidhi[0] = limits[1];
  solidlo[1] = limits[2];
  solidhi[1] = limits[3];
  solidlo[2] = limits[4];
  solidhi[2] = limits[5];

  solidsublo[0] = MAX(solidlo[0], sublo[0]);
  solidsublo[1] = MAX(solidlo[1], sublo[1]);
  solidsublo[2] = MAX(solidlo[2], sublo[2]);

  solidsubhi[0] = MIN(solidhi[0], subhi[0]);
  solidsubhi[1] = MIN(solidhi[1], subhi[1]);
  solidsubhi[2] = MIN(solidhi[2], subhi[2]);

#ifdef DEBUG
  cout << "proc " << universe->me
    << "\tsolidsublo=[" << solidsublo[0] << "," << solidsublo[1] << "," << solidsublo[2]
    << "]\t solidsubhi=["<< solidsubhi[0] << "," << solidsubhi[1] << "," << solidsubhi[2]
    << "]\n";

//   std::vector<double> x2plot, y2plot;
#endif

  // Calculate total number of particles np_local:
  int nsubx, nsuby, nsubz;
  double delta;
  //double hdelta;
  //double Lsubx, Lsuby, Lsubz;

  double *boundlo, *boundhi;

  delta = grid->cellsize;

  //////////////////////////////////////////////////////////

  // if (is_TL)
  //   {
  //   // The grid will be ajusted to the solid's domain (good for TLMPM).
  //   // and we need to create the corresponding grid:
  //     grid->init(solidlo, solidhi);

  //     boundlo = solidlo;

  //     Lsubx = solidsubhi[0] - solidsublo[0];
  //     if (domain->dimension >= 2) Lsuby = solidsubhi[1] - solidsublo[1];
  //     if (domain->dimension == 3) Lsubz = solidsubhi[2] - solidsublo[2];
  //   }
  // else
  //   {
  //     // The grid is most likely bigger than the solid's domain (good for ULMPM),
  //     // so all particles created won't lie in the region, they will need to be
  //     // checked:

  //     boundlo = domain->sublo;

  //     Lsubx = domain->subhi[0] - domain->sublo[0];
  //     if (domain->dimension >= 2)
  // 	Lsuby = domain->subhi[1] - domain->sublo[1];
  //     if (domain->dimension == 3)
  // 	Lsubz = domain->subhi[2] - domain->sublo[2];
  //   }

  // nsubx = (int) (Lsubx/delta);
  // while (nsubx*delta <= Lsubx - 0.1*delta) nsubx++;
  // nsubx++;

  // if (domain->dimension >= 2)
  //   {
  //     nsuby = (int) (Lsuby/delta);
  //   while (nsuby*delta <= Lsuby - 0.1*delta) nsuby++;
  //   nsuby++;
  //   }
  // else
  //   {
  //     nsuby = 1;
  //   }

  // if (domain->dimension == 3)
  //   {
  //     nsubz = (int) (Lsubz/delta);
  //     while (nsubz*delta <= Lsubz - 0.1*delta) nsubz++;
  //     nsubz++;
  //   }
  // else
  //   {
  //     nsubz = 1;
  //   }


  // cout << "1--proc " << universe->me << "\t nsub=["<< nsubx << "," << nsuby << "," << nsubz << "]\n";
  //////////////////////////////////////////////////////////

  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (is_TL)
  {
    grid->init(solidlo, solidhi);
    boundlo = solidlo;
    boundhi = solidhi;
  }
  else
  {
    boundlo = domain->boxlo;
    boundhi = domain->boxhi;
  }

  double Loffsetlo[3] = {MAX(0.0, sublo[0] - boundlo[0]),
    MAX(0.0, sublo[1] - boundlo[1]),
    MAX(0.0, sublo[2] - boundlo[2])};
  double Loffsethi[3] = {MAX(0.0, MIN(subhi[0], boundhi[0]) - boundlo[0]),
    MAX(0.0, MIN(subhi[1], boundhi[1]) - boundlo[1]),
    MAX(0.0, MIN(subhi[2], boundhi[2]) - boundlo[2])};

  int noffsetlo[3] = {(int)floor(Loffsetlo[0]/delta),
    (int)floor(Loffsetlo[1]/delta),
    (int)floor(Loffsetlo[2]/delta)};

  int noffsethi[3] = {(int)ceil(Loffsethi[0]/delta),
    (int)ceil(Loffsethi[1]/delta),
    (int)ceil(Loffsethi[2]/delta)};

  cout << "1--- proc " << universe->me << " noffsetlo=[" << noffsetlo[0]
    << "," << noffsetlo[1] << "," << noffsetlo[2] << "]\n";
  cout << "1--- proc " << universe->me << " noffsethi=[" << noffsethi[0]
    << "," << noffsethi[1] << "," << noffsethi[2] << "]\n";

// cout << "abs=" << abs(boundlo[0] + noffsethi[0]*delta - subhi[0])<< "]\n";
  if (universe->procneigh[0][1] >= 0 &&
      abs(boundlo[0] + noffsethi[0]*delta - subhi[0]) < 1.0e-12)
  {
    noffsethi[0]++;
  }
  if (domain->dimension >= 2 && universe->procneigh[1][1] >= 0 &&
      abs(boundlo[1] + noffsethi[1]*delta - subhi[1]) < 1.0e-12)
  {
    noffsethi[1]++;
  }
  if (domain->dimension == 3 && universe->procneigh[2][1] >= 0 &&
      abs(boundlo[2] + noffsethi[2]*delta - subhi[2]) < 1.0e-12)
  {
    noffsethi[2]++;
  }

  cout << "2--- proc " << universe->me << " noffsethi=[" << noffsethi[0]
    << "," << noffsethi[1] << "," << noffsethi[2] << "]\n";

  nsubx = MAX(0, noffsethi[0] - noffsetlo[0]);
  if (domain->dimension >= 2)
  {
    nsuby = MAX(0, noffsethi[1] - noffsetlo[1]);
  }
  else
  {
    nsuby = 1;
  }
  if (domain->dimension >= 3)
  {
    nsubz = MAX(0, noffsethi[2] - noffsetlo[2]);
  }
  else
  {
    nsubz = 1;
  }

  if (universe->procneigh[0][1] == -1)
  {
    while (boundlo[0] + delta*(noffsetlo[0] + nsubx - 0.5) <
           MIN(subhi[0], boundhi[0]))
      nsubx++;
  }
  if (universe->procneigh[1][1] == -1)
  {
    while (boundlo[1] + delta*(noffsetlo[1] + nsuby - 0.5) <
           MIN(subhi[1], boundhi[1]))
      nsuby++;
  }
  if (universe->procneigh[2][1] == -1)
  {
    while (boundlo[2] + delta*(noffsetlo[2] + nsubz - 0.5) <
           MIN(subhi[2], boundhi[2]))
      nsubz++;
  }
  cout << "2--proc " << universe->me << "\t nsub=[" << nsubx << "," << nsuby
    << "," << nsubz << "]\n";
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

// #ifdef DEBUG
//   cout << "proc " << universe->me << "\tLsub=[" << Lsubx << "," << Lsuby << "," << Lsubz << "]\t nsub=["<< nsubx << "," << nsuby << "," << nsubz << "]\n";
// #endif

  np_local = nsubx*nsuby*nsubz;

  // Create particles:

  if (universe->me == 0)
    cout << "delta = " << delta << endl;

  int l = 0;
  double vol_;

  if (domain->dimension == 1)
    vol_ = delta;
  else if (domain->dimension == 2)
    vol_ = delta*delta;
  else
    vol_ = delta*delta*delta;

  double mass_;
  if (mat->rigid)
    mass_ = 1;
  else
    mass_ = mat->rho0*vol_;

  np_per_cell = (int)input->parsev(args[3]);
  double xi = 0.5;
  double lp = delta;
  int nip = 1;
  vector<double> intpoints;

  if (np_per_cell == 1)
  {
    // One particle per cell at the center:

    xi = 0.5;
    lp *= 0.5;
    nip = 1;

    intpoints = {0, 0, 0};
  }
  else if (np_per_cell == 2)
  {
    // Quadratic elements:

    if (domain->dimension == 1) nip = 2;
    else if (domain->dimension == 2) nip = 4;
    else                             nip = 8;

    // if (is_TL && nc == 0)
    //   xi = 0.5/sqrt(3.0);
    // else
    xi = 0.25;

    lp *= 0.25;

    intpoints = {-xi, -xi, -xi, -xi, xi, -xi, xi, -xi, -xi, xi, xi, -xi,
      -xi, -xi, xi, -xi, xi, xi, xi, -xi, xi, xi, xi, xi};
  }
  else if (np_per_cell == 3)
  {
    // Berstein elements:

    if (nc == 0)
      xi = 0.7746/2;
    else
      xi = 1.0/3.0;

    lp *= 1.0/6.0;
    nip = 27;

    if (domain->dimension == 1) nip = 3;
    else if (domain->dimension == 2) nip = 9;
    else nip = 27;

    intpoints = {-xi, -xi, -xi,
      -xi, 0, -xi,
      -xi, xi, -xi,
      0, -xi, -xi,
      0, 0, -xi,
      0, xi, -xi,
      xi, -xi, -xi,
      xi, 0, -xi,
      xi, xi, -xi,
      -xi, -xi, 0,
      -xi, 0, 0,
      -xi, xi, 0,
      0, -xi, 0,
      0, 0, 0,
      0, xi, 0,
      xi, -xi, 0,
      xi, 0, 0,
      xi, xi, 0,
      -xi, -xi, xi,
      -xi, 0, xi,
      -xi, xi, xi,
      0, -xi, xi,
      0, 0, xi,
      0, xi, xi,
      xi, -xi, xi,
      xi, 0, xi,
      xi, xi, xi};

  }
  else
  {
    lp *= 1.0/(2*np_per_cell);

    if (domain->dimension == 1)
    {
      nip = np_per_cell;
    }
    else if (domain->dimension == 2)
    {
      nip = np_per_cell*np_per_cell;
    }
    else
    {
      nip = np_per_cell*np_per_cell*np_per_cell;
    }

    double d = 1.0/np_per_cell;

    for (int k = 0; k < np_per_cell; k++)
    {
      for (int i = 0; i < np_per_cell; i++)
      {
        for (int j = 0; j < np_per_cell; j++)
        {
          intpoints.push_back((i + 0.5)*d - 0.5);
          intpoints.push_back((j + 0.5)*d - 0.5);
          intpoints.push_back((k + 0.5)*d - 0.5);
        }
      }
    }
  }

  np_local *= nip;
// #ifdef DEBUG
//   cout << "proc " << universe->me << "\tnp_local=" << np_local << endl;
// #endif
  mass_ /= (double)nip;
  vol_ /= (double)nip;

  // Allocate the space in the vectors for np particles:
  grow(np_local);

  int dim = domain->dimension;
  bool r4 = false;
  bool q4 = false;

  if (method_type.compare("tlcpdi") == 0 || method_type.compare("ulcpdi") == 0)
  {
    if (update->method->style == 0)
    { // CPDI-R4
      r4 = true;
    }
    if (update->method->style == 1)
    { // CPDI-Q4
      q4 = true;
    }
  }

  for (int i = 0; i < nsubx; i++)
  {
    for (int j = 0; j < nsuby; j++)
    {
      for (int k = 0; k < nsubz; k++)
      {
        for (int ip = 0; ip < nip; ip++)
        {

          if (l >= np_local)
          {
            cout << "Error in Solid::populate(), exceeding the allocated number of particles.\n";
            cout << "l = " << l << endl;
            cout << "np_local = " << np_local << endl;
            error->all(FLERR, "");
          }

          x0[l][0] = x[l][0] =
            boundlo[0] + delta*(noffsetlo[0] + i + 0.5 + intpoints[3*ip+0]);
          x0[l][1] = x[l][1] =
            boundlo[1] + delta*(noffsetlo[1] + j + 0.5 + intpoints[3*ip+1]);
          if (dim == 3)
            x0[l][2] = x[l][2] =
            boundlo[2] + delta*(noffsetlo[2] + k + 0.5 + intpoints[3*ip+2]);
          else
            x0[l][2] = x[l][2] = 0;

          // Check if the particle is inside the region:
          if (domain->inside_subdomain(x0[l][0], x0[l][1], x0[l][2]) && domain->regions[iregion]->inside(x0[l][0], x0[l][1], x0[l][2]) == 1)
          {
            // cout << "Inside\n";

            if (r4)
            { // CPDI-R4
              rp0[dim*l][0] = rp[dim*l][0] = lp;
              rp0[dim*l][1] = rp[dim*l][1] = 0;
              rp0[dim*l][2] = rp[dim*l][2] = 0;

              if (dim >= 2)
              {
                rp0[dim*l + 1][0] = rp[dim*l + 1][0] = 0;
                rp0[dim*l + 1][1] = rp[dim*l + 1][1] = lp;
                rp0[dim*l + 1][2] = rp[dim*l + 1][2] = 0;

                if (dim == 3)
                {
                  rp0[dim*l + 2][0] = rp[dim*l + 1][0] = 0;
                  rp0[dim*l + 2][1] = rp[dim*l + 1][1] = 0;
                  rp0[dim*l + 2][0] = rp[dim*l + 1][0] = lp;
                }
              }
            }
            if (q4)
            { // CPDI-Q4
              xpc0[nc*l][0] = xpc[nc*l][0] = x0[l][0] - lp;

              xpc0[nc*l + 1][0] = xpc[nc*l + 1][0] = x0[l][0] + lp;

              if (dim >= 2)
              {
                xpc0[nc*l][1] = xpc[nc*l][1] = x0[l][1] - lp;
                xpc0[nc*l + 1][1] = xpc[nc*l + 1][1] = x0[l][1] - lp;

                xpc0[nc*l + 2][0] = xpc[nc*l + 2][0] = x0[l][0] + lp;
                xpc0[nc*l + 2][1] = xpc[nc*l + 2][1] = x0[l][1] + lp;

                xpc0[nc*l + 3][0] = xpc[nc*l + 3][0] = x0[l][0] - lp;
                xpc0[nc*l + 3][1] = xpc[nc*l + 3][1] = x0[l][1] + lp;
              }

              if (dim == 3)
              {
                xpc0[nc*l][2] = xpc[nc*l][2] = x0[l][2] - lp;
                xpc0[nc*l + 1][2] = xpc[nc*l + 1][2] = x0[l][2] - lp;
                xpc0[nc*l + 2][2] = xpc[nc*l + 2][2] = x0[l][2] - lp;
                xpc0[nc*l + 3][2] = xpc[nc*l + 3][2] = x0[l][2] - lp;

                xpc0[nc*l + 4][0] = xpc[nc*l + 4][0] = x0[l][0] - lp;
                xpc0[nc*l + 4][1] = xpc[nc*l + 4][1] = x0[l][1] - lp;
                xpc0[nc*l + 4][2] = xpc[nc*l + 4][2] = x0[l][2] + lp;

                xpc0[nc*l + 5][0] = xpc[nc*l + 5][0] = x0[l][0] + lp;
                xpc0[nc*l + 5][1] = xpc[nc*l + 5][1] = x0[l][1] - lp;
                xpc0[nc*l + 5][2] = xpc[nc*l + 5][2] = x0[l][2] + lp;

                xpc0[nc*l + 6][0] = xpc[nc*l + 6][0] = x0[l][0] + lp;
                xpc0[nc*l + 6][1] = xpc[nc*l + 6][1] = x0[l][1] + lp;
                xpc0[nc*l + 6][2] = xpc[nc*l + 6][2] = x0[l][2] + lp;

                xpc0[nc*l + 7][0] = xpc[nc*l + 7][0] = x0[l][0] - lp;
                xpc0[nc*l + 7][1] = xpc[nc*l + 7][1] = x0[l][1] + lp;
                xpc0[nc*l + 7][2] = xpc[nc*l + 7][2] = x0[l][2] + lp;
              }
            }
            l++;
          }
        }
      }
    }
  }

  if (np_local > l)
  {
    grow(l);
  }
  np_local = l; // Adjust np_local to account for the particles outside the domain

  tagint ptag0 = 0;

  for (int proc = 0; proc<universe->nprocs; proc++)
  {
    int np_local_bcast;
    if (proc == universe->me)
    {
// Send np_local
      np_local_bcast = np_local;
    }
    else
    {
   // Receive np_local
      np_local_bcast = 0;
    }
    MPI_Bcast(&np_local_bcast, 1, MPI_INT, proc, universe->uworld);
    if (universe->me > proc) ptag0 += np_local_bcast;
  }

// #ifdef DEBUG
//   cout << "proc " << universe->me << "\tptag0 = " << ptag0 << endl;
// #endif
  np_local = l; // Adjust np to account for the particles outside the domain
  cout << "np_local=" << np_local << endl;

  for (int i = 0; i < np_local; i++)
  {
    a[i] = Vector3d();
    v[i] = Vector3d();
    f[i] = Vector3d();
    mbp[i] = Vector3d();
    v_update[i] = Vector3d();
    rho0[i] = rho[i] = mat->rho0;

    if (domain->axisymmetric == true)
    {
      mass[i] = mass_*x0[i][0];
      vol0[i] = vol[i] = mass[i]/rho0[i];
    }
    else
    {
      mass[i] = mass_;
      vol0[i] = vol[i] = vol_;
    }

    eff_plastic_strain[i] = 0;
    eff_plastic_strain_rate[i] = 0;
    damage[i] = 0;
    damage_init[i] = 0;
    if (update->method->temp)
    {
      T[i] = T0;
      gamma[i] = 0;
      q[i] = Vector3d();
    }
    ienergy[i] = 0;
    strain_el[i] = Matrix3d();
    sigma[i] = Matrix3d();
    vol0PK1[i] = Matrix3d();
    L[i] = Matrix3d();
    F[i] = Matrix3d::identity();
    R[i] = Matrix3d::identity();
    D[i] = Matrix3d();
    Finv[i] = Matrix3d();
    Fdot[i] = Matrix3d();
    J[i] = 1;
    mask[i] = 1;

    ptag[i] = ptag0 + i + 1 + domain->np_total;
  }

  if (l != np_local)
  {
    cout << "Error l=" << l << " != np_local=" << np_local << endl;
    error->one(FLERR, "");
  }

  int np_local_reduced;
  MPI_Allreduce(&np_local, &np_local_reduced, 1, MPI_INT, MPI_SUM, universe->uworld);
  np += np_local_reduced;
  domain->np_total += np;
  domain->np_local += np_local;
}

void Solid::update_particle_domain()
{
  int dim = domain->dimension;

  if (update->method->style == 0)
  { // CPDI-R4
    for (int ip = 0; ip < np_local; ip++)
    {
      rp[dim*ip] = F.at(ip)*rp0[dim*ip];
      if (dim >= 2)
        rp[dim*ip + 1] = F.at(ip)*rp0[dim*ip + 1];
      if (dim == 3)
        rp[dim*ip + 2] = F.at(ip)*rp0[dim*ip + 2];
    }
  }
  // if (update->method->style == 1) { // CPDI-Q4
  //   for (int ip=0; ip<np; ip++) {
  //     for(int ic=0; ic<nc; ic++) {
  // 	xpc[nc*ip + ic] += update->dt*v_update.at(ip);
  //     }
  //   }
  // }
}

void Solid::read_file(string fileName)
{
}

void Solid::read_mesh(string fileName)
{

  string line;

  int id;
  int length;
  int elemType;
  vector<string> splitLine;
  array<double, 3> xn;

  int nodeCount; // Count the number of nodes
  vector<array<double, 3>> nodes;

// #ifdef DEBUG
//   std::vector<double> x2plot, y2plot;
//   std::vector<double> xcplot, ycplot;
// #endif

  ifstream file(fileName, std::ios::in);

  if (!file)
  {
    cout << "Error: unable to open mesh file " << fileName << endl;
    error->one(FLERR, "");
  }

  if (universe->me == 0)
    cout << "Reading Gmsh mesh file ...\n";

  while (getline(file, line))
  {

    if (line.compare("$MeshFormat") == 0)
    {
      // Read mesh format informations:
      double version;
      file >> version;

      if (version >= 3.0)
      {
        cout << "Gmsh mesh file version >=3.0 not supported.\n";
        error->one(FLERR, "");
      }

      getline(file, line);
      if (line.compare("$EndMeshFormat") == 0)
      {

        if (universe->me == 0)
          cout << "Reading format...done!\n";
        break;
      }
      else
        if (universe->me == 0)
          cout << "Unexpected line: " << line << ". $EndMeshFormat expected!!\n";
    }

    if (line.compare("$Nodes") == 0)
    {
      if (universe->me == 0)
        cout << "Reading nodes...\n";
            // Read mesh node informations:
      file >> nodeCount;

      nodes.resize(nodeCount);

      for (int in = 0; in < nodeCount; in++)
      {
        file >> id >> xn[0] >> xn[1] >> xn[2];
        if (abs(xn[0]) < 1.0e-12)
          xn[0] = 0;
        if (abs(xn[1]) < 1.0e-12)
          xn[1] = 0;
        if (abs(xn[2]) < 1.0e-12)
          xn[2] = 0;

        if (domain->dimension == 1 && (xn[1] != 0. || xn[2] != 0.))
        {
          cout << "Error: node " << id << " has non 1D component.\n";
          error->one(FLERR, "");
        }

        if (domain->dimension == 2 && xn[2] != 0.)
        {
          cout << "Error: node " << id << " has non zero z component.\n";
          error->one(FLERR, "");
        }

        nodes[in] = xn;

        // Adjust solid bounds:
        if (xn[0] < solidlo[0])
          solidlo[0] = xn[0];
        if (xn[1] < solidlo[1])
          solidlo[1] = xn[1];
        if (xn[2] < solidlo[2])
          solidlo[2] = xn[2];

        if (xn[0] > solidhi[0])
          solidhi[0] = xn[0];
        if (xn[1] > solidhi[1])
          solidhi[1] = xn[1];
        if (xn[2] > solidhi[2])
          solidhi[2] = xn[2];
      }

      getline(file, line);
      if (line.compare("$EndNodes") == 0)
      {
        if (universe->me == 0)
          cout << "Reading nodes...done!\n";
      }
      else
        if (universe->me == 0)
          cout << "Unexpected line: " << line << ". $EndNodes expected!!\n";
    }

    if (line.compare("$Elements") == 0)
    {
      if (universe->me == 0)
        cout << "Reading elements...\n";
      file >> np; // Number of elements
      getline(file, line);

      // Allocate the space in the vectors for np particles:
      grow(np);

      for (int ie = 0; ie < np; ie++)
      {
        getline(file, line);

        splitLine = split(line, "\t ");

        length = splitLine.size();

        elemType = stoi(splitLine[1]);
        // elemType == 1: 2-node   line element
        // elemType == 3: 4-node   quadrangle
        // elemType == 4: 4-node   tetrahedra

        if (elemType == 1)
        {
          int no1 = stoi(splitLine[5]) - 1;
          int no2 = stoi(splitLine[6]) - 1;

          xpc0[nc*ie][0] = xpc[nc*ie][0] = nodes[no1][0];
          xpc0[nc*ie][1] = xpc[nc*ie][1] = nodes[no1][1];
          xpc0[nc*ie][2] = xpc[nc*ie][2] = nodes[no1][2];

          xpc0[nc*ie + 1][0] = xpc[nc*ie + 1][0] = nodes[no2][0];
          xpc0[nc*ie + 1][1] = xpc[nc*ie + 1][1] = nodes[no2][1];
          xpc0[nc*ie + 1][2] = xpc[nc*ie + 1][2] = nodes[no2][2];

          x0[ie][0] = x[ie][0] = 0.5*(nodes[no1][0] + nodes[no2][0]);
          x0[ie][1] = x[ie][1] = 0.5*(nodes[no1][1] + nodes[no2][1]);
          x0[ie][2] = x[ie][2] = 0.5*(nodes[no1][2] + nodes[no2][2]);
        }
        else if (elemType == 3)
        {

// #ifdef DEBUG
//           xcplot.clear();
//           xcplot.resize(5, 0);
//           ycplot.clear();
//           ycplot.resize(5, 0);
// #endif

          int no1 = stoi(splitLine[5]) - 1;
          int no2 = stoi(splitLine[6]) - 1;
          int no3 = stoi(splitLine[7]) - 1;
          int no4 = stoi(splitLine[8]) - 1;

          if (method_type.compare("tlcpdi") == 0 ||
              method_type.compare("ulcpdi") == 0)
          {

            xpc0[nc*ie][0] = xpc[nc*ie][0] = nodes[no1][0];
            xpc0[nc*ie][1] = xpc[nc*ie][1] = nodes[no1][1];
            xpc0[nc*ie][2] = xpc[nc*ie][2] = nodes[no1][2];

            xpc0[nc*ie + 1][0] = xpc[nc*ie + 1][0] = nodes[no2][0];
            xpc0[nc*ie + 1][1] = xpc[nc*ie + 1][1] = nodes[no2][1];
            xpc0[nc*ie + 1][2] = xpc[nc*ie + 1][2] = nodes[no2][2];

            xpc0[nc*ie + 2][0] = xpc[nc*ie + 2][0] = nodes[no3][0];
            xpc0[nc*ie + 2][1] = xpc[nc*ie + 2][1] = nodes[no3][1];
            xpc0[nc*ie + 2][2] = xpc[nc*ie + 2][2] = nodes[no3][2];

            xpc0[nc*ie + 3][0] = xpc[nc*ie + 3][0] = nodes[no4][0];
            xpc0[nc*ie + 3][1] = xpc[nc*ie + 3][1] = nodes[no4][1];
            xpc0[nc*ie + 3][2] = xpc[nc*ie + 3][2] = nodes[no4][2];
          }

          x0[ie][0] = x[ie][0] = 0.25*(nodes[no1][0] + nodes[no2][0] +
                                         nodes[no3][0] + nodes[no4][0]);
          x0[ie][1] = x[ie][1] = 0.25*(nodes[no1][1] + nodes[no2][1] +
                                         nodes[no3][1] + nodes[no4][1]);
          x0[ie][2] = x[ie][2] = 0.25*(nodes[no1][2] + nodes[no2][2] +
                                         nodes[no3][2] + nodes[no4][2]);

          // vol0[ie] = vol[ie] = 0.5*(xpc[nc*ie+0][0]*xpc[nc*ie+1][1] -
          // xpc[nc*ie+1][0]*xpc[nc*ie+0][1]
          // 	    + xpc[nc*ie+1][0]*xpc[nc*ie+2][1] -
          // xpc[nc*ie+2][0]*xpc[nc*ie+1][1]
          // 	    + xpc[nc*ie+2][0]*xpc[nc*ie+3][1] -
          // xpc[nc*ie+3][0]*xpc[nc*ie+2][1]
          // 	    + xpc[nc*ie+3][0]*xpc[nc*ie+0][1] -
          // xpc[nc*ie+0][0]*xpc[nc*ie+3][1]);

          vol0[ie] = vol[ie] =
            0.5 *
            (nodes[no1][0]*nodes[no2][1] - nodes[no2][0]*nodes[no1][1] +
             nodes[no2][0]*nodes[no3][1] - nodes[no3][0]*nodes[no2][1] +
             nodes[no3][0]*nodes[no4][1] - nodes[no4][0]*nodes[no3][1] +
             nodes[no4][0]*nodes[no1][1] - nodes[no1][0]*nodes[no4][1]);
        }
        else if (elemType == 4)
        {

          int no1 = stoi(splitLine[5]) - 1;
          int no2 = stoi(splitLine[6]) - 1;
          int no3 = stoi(splitLine[7]) - 1;
          int no4 = stoi(splitLine[8]) - 1;

          double x1 = nodes[no1][0];
          double y1 = nodes[no1][1];
          double z1 = nodes[no1][2];
          double x2 = nodes[no2][0];
          double y2 = nodes[no2][1];
          double z2 = nodes[no2][2];
          double x3 = nodes[no3][0];
          double y3 = nodes[no3][1];
          double z3 = nodes[no3][2];
          double x4 = nodes[no4][0];
          double y4 = nodes[no4][1];
          double z4 = nodes[no4][2];

          double x21 = x2 - x1;
          double x32 = x3 - x2;
          double x43 = x4 - x3;
          double x42 = x4 - x2;
          double y23 = y2 - y3;
          double y34 = y3 - y4;
          double y12 = y1 - y2;
          double y42 = y4 - y2;
          double z34 = z3 - z4;
          double z23 = z2 - z3;
          double z12 = z1 - z2;
          double z42 = z4 - z2;

          x0[ie][0] = x[ie][0] = 0.25*(x1 + x2 + x3 + x4);
          x0[ie][1] = x[ie][1] = 0.25*(y1 + y2 + y3 + y4);
          x0[ie][2] = x[ie][2] = 0.25*(z1 + z2 + z3 + z4);

          vol0[ie] = vol[ie] = (1/6)*(x21*(y23*z34 - y34*z23) +
                                          x32*(y34*z12 - y12*z34) +
                                          x43*(y12*z23 - y23*z12));
        }
        else
        {
          cout << "Element type " << elemType << " not supported!!\n";
          error->one(FLERR, "");
        }
      }

      getline(file, line);
      if (line.compare("$EndElements") == 0)
      {
        if (universe->me == 0)
          cout << "Reading elements...done!\n";
        break;
      }
      else
        if (universe->me == 0)
          cout << "Unexpected line: " << line << ". $EndElements expected!!\n";
    }
  }

  if (universe->me == 0)
    cout << "np=" << np << endl;

  for (int i = 0; i < np; i++)
  {
    a[i] = Vector3d();
    v[i] = Vector3d();
    f[i] = Vector3d();
    mbp[i] = Vector3d();
    v_update[i] = Vector3d();
    rho0[i] = rho[i] = mat->rho0;
    mass[i] = mat->rho0*vol0[i];
    eff_plastic_strain[i] = 0;
    eff_plastic_strain_rate[i] = 0;
    damage[i] = 0;
    damage_init[i] = 0;
    if (update->method->temp)
    {
      T[i] = T0;
      gamma[i] = 0;
      q[i] = Vector3d();
    }
    ienergy[i] = 0;
    strain_el[i] = Matrix3d();
    sigma[i] = Matrix3d();
    vol0PK1[i] = Matrix3d();
    L[i] = Matrix3d();
    F[i] = Matrix3d::identity();
    R[i] = Matrix3d::identity();
    D[i] = Matrix3d();
    Finv[i] = Matrix3d();
    Fdot[i] = Matrix3d();
    J[i] = 1;

    if (x0[i][0] < solidlo[0])
      solidlo[0] = x0[i][0];
    if (x0[i][1] < solidlo[1])
      solidlo[1] = x0[i][1];
    if (x0[i][2] < solidlo[2])
      solidlo[2] = x0[i][2];

    if (x0[i][0] > solidhi[0])
      solidhi[0] = x0[i][0];
    if (x0[i][1] > solidhi[1])
      solidhi[1] = x0[i][1];
    if (x0[i][2] > solidhi[2])
      solidhi[2] = x0[i][2];
  }

  if (grid->nnodes == 0)
  {
    grid->init(solidlo, solidhi);
  }
}

void Solid::write_restart(ofstream *of)
{
// Write solid bounds:
  of->write(reinterpret_cast<const char *>(&solidlo[0]), 3*sizeof(double));
  of->write(reinterpret_cast<const char *>(&solidhi[0]), 3*sizeof(double));
  of->write(reinterpret_cast<const char *>(&solidsublo[0]), 3*sizeof(double));
  of->write(reinterpret_cast<const char *>(&solidsubhi[0]), 3*sizeof(double));

  // Write number of particles:
  of->write(reinterpret_cast<const char *>(&np), sizeof(bigint));
  of->write(reinterpret_cast<const char *>(&np_local), sizeof(int));
  of->write(reinterpret_cast<const char *>(&nc), sizeof(int));

  // Write material's info:
  int iMat = material->find_material(mat->id);
  of->write(reinterpret_cast<const char *>(&iMat), sizeof(int));

  // Write cellsize:
  of->write(reinterpret_cast<const char *>(&grid->cellsize), sizeof(double));


  // Write particle's attributes:
  // cout << x[0](0) << ", " << x[0](1) << ", " << x[0](2) << endl;
  for (int ip = 0; ip < np_local; ip++)
  {
    of->write(reinterpret_cast<const char *>(&ptag.at(ip)), sizeof(tagint));
    of->write(reinterpret_cast<const char *>(&x0.at(ip)), sizeof(Vector3d));
    of->write(reinterpret_cast<const char *>(&x.at(ip)), sizeof(Vector3d));
    of->write(reinterpret_cast<const char *>(&v.at(ip)), sizeof(Vector3d));
    of->write(reinterpret_cast<const char *>(&sigma.at(ip)), sizeof(Matrix3d));
    of->write(reinterpret_cast<const char *>(&strain_el.at(ip)), sizeof(Matrix3d));
    if (is_TL)
    {
      of->write(reinterpret_cast<const char *>(&vol0PK1.at(ip)), sizeof(Matrix3d));
    }
    of->write(reinterpret_cast<const char *>(&F.at(ip)), sizeof(Matrix3d));
    of->write(reinterpret_cast<const char *>(&J.at(ip)), sizeof(double));
    of->write(reinterpret_cast<const char *>(&vol0.at(ip)), sizeof(double));
    of->write(reinterpret_cast<const char *>(&rho0.at(ip)), sizeof(double));
    of->write(reinterpret_cast<const char *>(&eff_plastic_strain.at(ip)), sizeof(double));
    of->write(reinterpret_cast<const char *>(&eff_plastic_strain_rate.at(ip)), sizeof(double));
    of->write(reinterpret_cast<const char *>(&damage.at(ip)), sizeof(double));
    of->write(reinterpret_cast<const char *>(&damage_init.at(ip)), sizeof(double));
    of->write(reinterpret_cast<const char *>(&T.at(ip)), sizeof(double));
    of->write(reinterpret_cast<const char *>(&ienergy.at(ip)), sizeof(double));
    of->write(reinterpret_cast<const char *>(&mask.at(ip)), sizeof(int));
  }
}

void Solid::read_restart(ifstream *ifr)
{
// Read solid bounds:
  ifr->read(reinterpret_cast<char *>(&solidlo[0]), 3*sizeof(double));
  ifr->read(reinterpret_cast<char *>(&solidhi[0]), 3*sizeof(double));
  ifr->read(reinterpret_cast<char *>(&solidsublo[0]), 3*sizeof(double));
  ifr->read(reinterpret_cast<char *>(&solidsubhi[0]), 3*sizeof(double));
  // cout << "solidlo=[" << solidlo[0] << "," << solidlo[1] << "," << solidlo[2] << endl;
  // cout << "solidhi=[" << solidhi[0] << "," << solidhi[1] << "," << solidhi[2] << endl;
  // cout << "solidsublo=[" << solidsublo[0] << "," << solidsublo[1] << "," << solidsublo[2] << endl;
  // cout << "solidsubhi=[" << solidsubhi[0] << "," << solidsubhi[1] << "," << solidsubhi[2] << endl;

  //init();
  // Read number of particles:
  ifr->read(reinterpret_cast<char *>(&np), sizeof(bigint));
  ifr->read(reinterpret_cast<char *>(&np_local), sizeof(int));
  ifr->read(reinterpret_cast<char *>(&nc), sizeof(int));

  // Read material's info:
  int iMat = -1;
  ifr->read(reinterpret_cast<char *>(&iMat), sizeof(int));
  mat = &material->materials[iMat];

  // Read cellsize:
  ifr->read(reinterpret_cast<char *>(&grid->cellsize), sizeof(double));
  if (is_TL)
  {
    grid->init(solidlo, solidhi);
  }

  // Read particle's attributes:
  grow(np_local);

  for (int ip = 0; ip < np_local; ip++)
  {
    ifr->read(reinterpret_cast<char *>(&ptag.at(ip)), sizeof(tagint));
    ifr->read(reinterpret_cast<char *>(&x0.at(ip)), sizeof(Vector3d));
    ifr->read(reinterpret_cast<char *>(&x.at(ip)), sizeof(Vector3d));
    ifr->read(reinterpret_cast<char *>(&v.at(ip)), sizeof(Vector3d));
    v_update.at(ip) = Vector3d();
    a.at(ip) = Vector3d();
    mbp.at(ip) = Vector3d();
    f.at(ip) = Vector3d();
    ifr->read(reinterpret_cast<char *>(&sigma.at(ip)), sizeof(Matrix3d));
    ifr->read(reinterpret_cast<char *>(&strain_el.at(ip)), sizeof(Matrix3d));
    if (is_TL)
    {
      ifr->read(reinterpret_cast<char *>(&vol0PK1.at(ip)), sizeof(Matrix3d));
    }
    L.at(ip) = Matrix3d();
    ifr->read(reinterpret_cast<char *>(&F.at(ip)), sizeof(Matrix3d));
    R.at(ip) = Matrix3d();
    D.at(ip) = Matrix3d();
    Finv.at(ip) = Matrix3d();
    Fdot.at(ip) = Matrix3d();
    ifr->read(reinterpret_cast<char *>(&J.at(ip)), sizeof(double));
    ifr->read(reinterpret_cast<char *>(&vol0.at(ip)), sizeof(double));
    vol.at(ip) = J.at(ip)*vol0.at(ip);
    ifr->read(reinterpret_cast<char *>(&rho0.at(ip)), sizeof(double));
    rho.at(ip) = rho0.at(ip)/J.at(ip);
    mass.at(ip) = rho0.at(ip)*vol0.at(ip);
    ifr->read(reinterpret_cast<char *>(&eff_plastic_strain.at(ip)), sizeof(double));
    ifr->read(reinterpret_cast<char *>(&eff_plastic_strain_rate.at(ip)), sizeof(double));
    ifr->read(reinterpret_cast<char *>(&damage.at(ip)), sizeof(double));
    ifr->read(reinterpret_cast<char *>(&damage_init.at(ip)), sizeof(double));
    ifr->read(reinterpret_cast<char *>(&T.at(ip)), sizeof(double));
    ifr->read(reinterpret_cast<char *>(&ienergy.at(ip)), sizeof(double));
    ifr->read(reinterpret_cast<char *>(&mask.at(ip)), sizeof(int));
  }
  // cout << x[0](0) << ", " << x[0](1) << ", " << x[0](2) << endl;
}
