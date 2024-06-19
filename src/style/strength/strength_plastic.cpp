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

#include <strength_plastic.h>
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


StrengthPlastic::StrengthPlastic(MPM *mpm, vector<string> args) : Strength(mpm, args)
{
  if (universe->me == 0) {
    cout << "Initiate StrengthPlastic" << endl;
  }

  if (args.size() < 3) {
    error->all(FLERR, "Error: too few arguments for the strength command.\n");
  }

  if (args[2].compare("restart") ==
      0) { // If the keyword restart, we are expecting to have read_restart()
           // launched right after.
    G_ = yieldStress = 0;
    return;
  }

  //options(&args, args.begin()+3);
  G_ = input->parsev(args[2]);
  yieldStress = input->parsev(args[3]);
  if (universe->me == 0) {
    cout << "Linear plastic strength model:\n";
    cout << "\tG: shear modulus " << G_ << endl;

    cout << "\t yield stress: " << yieldStress << endl;
  }
}

float StrengthPlastic::G(){
  return G_;
}

void
StrengthPlastic::update_deviatoric_stress(Solid &solid,
                                          Kokkos::View<float*> &plastic_strain_increment,
                                          Kokkos::View<Matrix3d*> &sigma_dev) const
{
  float G_ = this->G_;
  float yieldStress = this->yieldStress;
  
  float dt = update->dt;

  Kokkos::View<Matrix3d*> ssigma = solid.sigma;
  Kokkos::View<Matrix3d*> sD = solid.D;
  Kokkos::View<float*> sdamage = solid.damage;

  Kokkos::parallel_for("EOSLinear::compute_pressure", solid.np_local,
  KOKKOS_LAMBDA (const int &ip)
  {
    Matrix3d sigmaInitial_dev, sigmaFinal_dev, sigmaTrial_dev, dev_rate;
    float J2, Gd, yieldStressD;

    /*
     * deviatoric rate of unrotated stress
     */
    Gd = G_*(1 - sdamage[ip]);
    yieldStressD = yieldStress*(1 - sdamage[ip]);
    dev_rate = 2*Gd*Deviator(sD[ip]);
    sigmaInitial_dev = Deviator(ssigma[ip]);

    /*
     * perform a trial elastic update to the deviatoric stress
     */
    sigmaTrial_dev = sigmaInitial_dev + dt*dev_rate; // increment stress deviator using deviatoric rate

    /*
     * check yield condition
     */
    J2 = Kokkos::sqrt(1.5)*sigmaTrial_dev.norm();
    sigmaFinal_dev = sigmaTrial_dev;

    if (J2 < yieldStressD) {
      /*
       * no yielding has occured.
       * final deviatoric stress is trial deviatoric stress
       */
      plastic_strain_increment[ip] = 0;

    } else {
      //printf("yiedl\n");
      /*
       * yielding has occured
       */

      plastic_strain_increment[ip] = (J2 - yieldStressD)/(3*Gd);
      /*
       * new deviatoric stress:
       * obtain by scaling the trial stress deviator
       */
      sigmaFinal_dev *= yieldStressD/J2;

    }

    sigma_dev[ip] = sigmaFinal_dev;
  });
}

void StrengthPlastic::write_restart(ofstream *of) {
  of->write(reinterpret_cast<const char *>(&G_), sizeof(float));
  of->write(reinterpret_cast<const char *>(&yieldStress), sizeof(float));
}

void StrengthPlastic::read_restart(ifstream *ifr) {
  if (universe->me == 0) {
    cout << "Restart StrengthPlastic" << endl;
  }
  ifr->read(reinterpret_cast<char *>(&G_), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&yieldStress), sizeof(float));
}
