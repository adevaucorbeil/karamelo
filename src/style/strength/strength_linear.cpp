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

#include <strength_linear.h>
#include <domain.h>
#include <error.h>
#include <input.h>
#include <mpm_math.h>
#include <universe.h>
#include <update.h>
#include <var.h>
#include <matrix.h>
#include <iostream>

using namespace std;

using namespace MPM_Math;


StrengthLinear::StrengthLinear(MPM *mpm, vector<string> args) : Strength(mpm, args)
{
  if (universe->me == 0) {
    cout << "Initiate StrengthLinear" << endl;
  }

  if (args.size() < 3) {
    error->all(FLERR, "Error: too few arguments for the strength command.\n");
  }

  if (args[2].compare("restart") ==
      0) { // If the keyword restart, we are expecting to have read_restart()
           // launched right after.
    G_ = 0;
    return;
  }

  //options(&args, args.begin()+3);
  G_ = input->parsev(args[2]);
  if (universe->me == 0) {
    cout << "Linear strength model:\n";
    cout << "\tG: shear modulus " << G_ << endl;
  }
}

double StrengthLinear::G(){
  return G_;
}

void
StrengthLinear::update_deviatoric_stress(Solid &solid,
                                              Kokkos::View<double*> &plastic_strain_increment,
                                              Kokkos::View<Matrix3d*> &sigma_dev) const
{
  double G_ = this->G_;
  double dt = update->dt;

  Kokkos::View<Matrix3d*> ssigma = solid.sigma;
  Kokkos::View<Matrix3d*> sD = solid.D;
  Kokkos::View<double*> sdamage = solid.damage;

  Kokkos::parallel_for("EOSLinear::compute_pressure", solid.np_local,
  KOKKOS_LAMBDA (const int &ip)
  {
    sigma_dev[ip] = Deviator(ssigma[ip]) + dt*2*G_*(1 - sdamage[ip])*Deviator(sD[ip]);
  });
}

void StrengthLinear::write_restart(ofstream *of) {
  of->write(reinterpret_cast<const char *>(&G_), sizeof(double));
}

void StrengthLinear::read_restart(ifstream *ifr) {
  if (universe->me == 0) {
    cout << "Restart StrengthLinear" << endl;
  }
  ifr->read(reinterpret_cast<char *>(&G_), sizeof(double));
}

