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

#include <strength_swift.h>
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

#define SQRT_3_OVER_2 1.224744871 // sqrt(3.0/2.0)

StrengthSwift::StrengthSwift(MPM *mpm, vector<string> args)
    : Strength(mpm, args)
{
  if (universe->me == 0) {
    cout << "Initiate StrengthSwift" << endl;
  }

  if (args.size() < 3) {
    error->all(FLERR, "Error: too few arguments for the strength command.\n");
  }

  if (args[2].compare("restart") ==
      0) { // If the keyword restart, we are expecting to have read_restart()
           // launched right after.
    G_ = A = B = C = n = 0;
    return;
  }

  if (args.size() < Nargs) {
    error->all(FLERR,
               "Error: too few arguments for the strength command.\n" + usage);
  }

  // options(&args, args.begin()+3);
  G_      = input->parsev(args[2]);
  A       = input->parsev(args[3]);
  B       = input->parsev(args[4]);
  C       = input->parsev(args[5]);
  n       = input->parsev(args[6]);

  if (universe->me == 0) {
    cout << "Swift material strength model:\n";
    cout << "\tG: shear modulus " << G_ << endl;
    cout << "\tA: initial yield stress " << A << endl;
    cout << "\tB: proportionality factor for plastic strain dependency " << B
         << endl;
    cout << "\tC: plastic offset " << C << endl;
    cout << "\tn: exponent for plastic strain dependency " << n << endl;
    cout << "\tFlow stress: sigma_f = [" << A << " + " << B << " * (eps_p - "
         << C << ")^" << n << "]";
    cout << "(1 - D)\n";
  }
}

double StrengthSwift::G() { return G_; }

void
StrengthSwift::update_deviatoric_stress(Solid &solid,
                                        Kokkos::View<double*> &plastic_strain_increment,
                                        Kokkos::View<Matrix3d*> &sigma_dev) const
{
  double C = this->C;
  double A = this->A;
  double B = this->B;
  double n = this->n;
  double G_ = this->G_;
  double dt = update->dt;

  Kokkos::parallel_for("EOSLinear::compute_pressure", solid.np_local,
  KOKKOS_LAMBDA (const int &ip)
  {
    if (solid.damage[ip] >= 1)
    {
      sigma_dev[ip] = Matrix3d();
      return;
    }

    Matrix3d sigmaFinal_dev, sigmaTrial, sigmaTrial_dev;
    double J2, Gd, yieldStress;

    if (solid.eff_plastic_strain[ip] > 1.0e-10 && solid.eff_plastic_strain[ip] > C)
      yieldStress = A + B*Kokkos::Experimental::pow(solid.eff_plastic_strain[ip] - C, n);
    else
      yieldStress = A;

    /*
     * deviatoric rate of unrotated stress
     */
    Gd               = G_;

    if (solid.damage[ip] > 0)
    {
      Gd *= (1 - solid.damage[ip]);
      yieldStress *= (1 - solid.damage[ip]);
    }

    // sigmaInitial_dev = Deviator(sigma);

    /*
     * perform a trial elastic update to the deviatoric stress
     */
    sigmaTrial = solid.sigma[ip] + dt*2*Gd*solid.D[ip];
    sigmaTrial_dev = Deviator(sigmaTrial); // increment stress deviator using deviatoric rate

    /*
     * check yield condition
     */
    J2             = SQRT_3_OVER_2 * sigmaTrial_dev.norm();
    sigmaFinal_dev = sigmaTrial_dev;

    if (J2 < yieldStress)
    {
      /*
       * no yielding has occured.
       * final deviatoric stress is trial deviatoric stress
       */
      plastic_strain_increment[ip] = 0.0;
    }
    else
    {
      // printf("yiedl\n");
      /*
       * yielding has occured
       */

      plastic_strain_increment[ip] = (J2 - yieldStress)/(3*Gd);
      /*
       * new deviatoric stress:
       * obtain by scaling the trial stress deviator
       */
      sigmaFinal_dev *= yieldStress/J2;
    }

    sigma_dev[ip] = sigmaFinal_dev;
  });
}

void StrengthSwift::write_restart(ofstream *of) {
  of->write(reinterpret_cast<const char *>(&G_), sizeof(double));
  of->write(reinterpret_cast<const char *>(&A), sizeof(double));
  of->write(reinterpret_cast<const char *>(&B), sizeof(double));
  of->write(reinterpret_cast<const char *>(&C), sizeof(double));
  of->write(reinterpret_cast<const char *>(&n), sizeof(double));
}

void StrengthSwift::read_restart(ifstream *ifr) {
  // cout << "Restart StrengthSwift" << endl;
  ifr->read(reinterpret_cast<char *>(&G_), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&A), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&B), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&C), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&n), sizeof(double));
}
