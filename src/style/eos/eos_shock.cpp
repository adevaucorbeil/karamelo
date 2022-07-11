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

#include <eos_shock.h>
#include <domain.h>
#include <error.h>
#include <input.h>
#include <math_special.h>
#include <mpm_math.h>
#include <universe.h>
#include <var.h>
#include <matrix.h>
#include <iostream>

using namespace std;

using namespace MPM_Math;
using namespace MathSpecial;


EOSShock::EOSShock(MPM *mpm, vector<string> args) : EOS(mpm, args)
{
  if (universe->me == 0) {
    cout << "Initiate EOSShock" << endl;
  }

  if (args.size() < 3) {
    error->all(FLERR, "Error: not enough arguments.\n");
  }

  if (args[2].compare("restart") ==
      0) { // If the keyword restart, we are expecting to have read_restart()
           // launched right after.
    rho0_ = K_ = e0 = c0 = S = Gamma = Tr = cv = alpha = Q1 = Q2 = 0;
    artificial_viscosity = false;
    return;
  }

  if (args.size() < Nargs) {
    error->all(FLERR, "Error: not enough arguments.\n" + usage);
  }
  //options(&args, args.begin()+3);
  rho0_ = input->parsev(args[2]);
  K_ = input->parsev(args[3]);
  c0 = input->parsev(args[4]);
  S = input->parsev(args[5]);
  Gamma = input->parsev(args[6]);
  cv = input->parsev(args[7]);
  Tr = input->parsev(args[8]);
  Q1 = input->parsev(args[9]);
  Q2 = input->parsev(args[10]);

  if (universe->me == 0) {
    cout << "Set rho0 to " << rho0_ << endl;
    cout << "Set K to " << K_ << endl;
    cout << "Set c0 to " << c0 << endl;
    cout << "Set S to " << S << endl;
    cout << "Set gamma to " << Gamma << endl;
    cout << "Set cv to " << cv << endl;
    cout << "Set Tr to " << Tr << endl;
    cout << "Set the artificial vicosity coefficient Q1 to " << Q1 << endl;
    cout << "Set the artificial vicosity coefficient Q2 to " << Q2 << endl;
  }

  if (Q1 == 0 && Q2 == 0) {
    artificial_viscosity = false;
  } else {
    artificial_viscosity = true;
  }

  alpha = cv*rho0_;
  e0 = 0;
}


EOSShock::~EOSShock()
{

}

float EOSShock::rho0(){
  return rho0_;
}

float EOSShock::K(){
  return K_;
}

//pH, solid.ienergy[ip], solid.J[ip], solid.rho[ip], solid.damage[ip], solid.D[ip], solid.grid->cellsize, T
//float &pFinal, float &e, const float J, const float rho, const float damage, const Matrix3d D, const float cellsize, const float T


void EOSShock::compute_pressure(Solid &solid, Kokkos::View<float*> &pH) const
{
  float rho0_ = this->rho0_;
  float c0 = this->c0;
  float S = this->S;
  float Tr = this->Tr;
  float alpha = this->alpha;
  float Gamma = this->Gamma;
  float e0 = this->e0;
  bool artificial_viscosity = this->artificial_viscosity;
  float Q1 = this->Q1;
  float Q2 = this->Q2;

  float cp = solid.mat->cp;
  float cellsize = solid.grid->cellsize;

  Kokkos::View<float*> srho = solid.rho;
  Kokkos::View<float*> sJ = solid.J;
  Kokkos::View<Matrix3d*> sD = solid.D;
  Kokkos::View<float*> sT = solid.T;
  Kokkos::View<float*> sdamage = solid.damage;
  Kokkos::View<float*> sienergy = solid.ienergy;

  Kokkos::parallel_for("EOSShock::compute_pressure", solid.np_local,
  KOKKOS_LAMBDA (const int &ip)
  {
    float mu = srho[ip]/rho0_ - 1;
    float pH0 = rho0_*c0*c0*mu*(1 + mu)/(1 - (S - 1)*mu)/(1 - (S - 1)*mu);

    if (cp && sT[ip] > Tr)
      sienergy[ip] = alpha*(sT[ip] - Tr);
    else
       sienergy[ip] = 0;
    pH[ip] = pH0 + Gamma*(sienergy[ip] - e0);

    if (sdamage[ip] > 0 && pH[ip] < 0)
    {
      if (sdamage[ip] >= 1.0)
      {
	    // pFinal = rho0_ * Gamma * (e - e0);
	    pH[ip] = 0;
      }
      else
      {
	    // float mu_damaged = (1.0 - damage) * mu;
	    // float pH_damaged = rho0_ * (1.0 - damage) * square(c0) * mu_damaged * (1.0 + mu_damaged) / square(1.0 - (S - 1.0) * mu_damaged);
	    // pFinal = (pH_damaged + rho0_ * (1 + mu_damaged) * Gamma * (e - e0));;
	    pH[ip] *= 1 - sdamage[ip];
      }
    }

    if (artificial_viscosity)
    {
      float tr_eps = sD[ip].trace();

      if (tr_eps < 0)
        pH[ip] += srho[ip]*cellsize*tr_eps*(Q1*cellsize*tr_eps - Q2*c0*Kokkos::Experimental::sqrt(sJ[ip]));
    }
  });
}

void EOSShock::write_restart(ofstream *of) {
  of->write(reinterpret_cast<const char *>(&rho0_), sizeof(float));
  of->write(reinterpret_cast<const char *>(&K_), sizeof(float));
  of->write(reinterpret_cast<const char *>(&c0), sizeof(float));
  of->write(reinterpret_cast<const char *>(&S), sizeof(float));
  of->write(reinterpret_cast<const char *>(&Gamma), sizeof(float));
  of->write(reinterpret_cast<const char *>(&Tr), sizeof(float));
  of->write(reinterpret_cast<const char *>(&cv), sizeof(float));
  of->write(reinterpret_cast<const char *>(&Q1), sizeof(float));
  of->write(reinterpret_cast<const char *>(&Q2), sizeof(float));
}

void EOSShock::read_restart(ifstream *ifr) {
  if (universe->me == 0) {
    cout << "Restart EOSShock" << endl;
  }
  ifr->read(reinterpret_cast<char *>(&rho0_), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&K_), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&c0), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&S), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&Gamma), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&Tr), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&cv), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&Q1), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&Q2), sizeof(float));
  if (Q1 == 0 && Q2 == 0) {
    artificial_viscosity = false;
  } else {
    artificial_viscosity = true;
  }

  alpha = cv*rho0_;
  e0 = 0;
}

