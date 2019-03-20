#include <iostream>
#include "strength_plastic.h"
#include "input.h"
#include "domain.h"
#include "update.h"
#include "mpm_math.h"
#include <Eigen/Eigen>
#include "var.h"

using namespace std;
using namespace Eigen;
using namespace MPM_Math;


StrengthPlastic::StrengthPlastic(MPM *mpm, vector<string> args) : Strength(mpm, args)
{
  cout << "Initiate StrengthPlastic" << endl;

  if (args.size()<3) {
    cout << "Error: region command not enough arguments" << endl;
    exit(1);
  }
  //options(&args, args.begin()+3);
  G_ = input->parsev(args[2]);
  cout << "Set G to " << G_ << endl;

  yieldStress = input->parsev(args[3]);
  cout << "Set yieldStress to " << yieldStress << endl;
}

double StrengthPlastic::G(){
  return G_;
}

Matrix3d StrengthPlastic::update_deviatoric_stress(const Matrix3d sigma, const Matrix3d D, double &plastic_strain_increment)
{
  Matrix3d sigmaInitial_dev, sigmaFinal_dev, sigmaTrial_dev, dev_rate;
  double J2;

  /*
   * deviatoric rate of unrotated stress
   */
  dev_rate = 2.0 * G_ * Deviator(D);
  sigmaInitial_dev = Deviator(sigma);

  /*
   * perform a trial elastic update to the deviatoric stress
   */
  sigmaTrial_dev = sigmaInitial_dev + update->dt * dev_rate; // increment stress deviator using deviatoric rate

  /*
   * check yield condition
   */
  J2 = sqrt(3. / 2.) * sigmaTrial_dev.norm();
  sigmaFinal_dev = sigmaTrial_dev;

  if (J2 < yieldStress) {
    /*
     * no yielding has occured.
     * final deviatoric stress is trial deviatoric stress
     */
    plastic_strain_increment = 0.0;

  } else {
    //printf("yiedl\n");
    /*
     * yielding has occured
     */

    plastic_strain_increment = (J2 - yieldStress) / (3.0 * G_);
    /*
     * new deviatoric stress:
     * obtain by scaling the trial stress deviator
     */
    sigmaFinal_dev *= (yieldStress / J2);

  }

  return sigmaFinal_dev;
}
