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

double EOSShock::rho0(){
  return rho0_;
}

double EOSShock::K(){
  return K_;
}

//pH, solid.ienergy[ip], solid.J[ip], solid.rho[ip], solid.damage[ip], solid.D[ip], solid.grid->cellsize, T
//double &pFinal, double &e, const double J, const double rho, const double damage, const Matrix3d D, const double cellsize, const double T


void EOSShock::compute_pressure(Solid &solid, Kokkos::View<double*> &pH) const
{
  double rho0_ = this->rho0_;
  double c0 = this->c0;
  double S = this->S;
  double Tr = this->Tr;
  double alpha = this->alpha;
  double Gamma = this->Gamma;
  double e0 = this->e0;
  bool artificial_viscosity = this->artificial_viscosity;
  double Q1 = this->Q1;
  double Q2 = this->Q2;

  double cp = solid.mat->cp;
  double cellsize = solid.grid->cellsize;

  Kokkos::parallel_for("EOSShock::compute_pressure", solid.np_local,
  KOKKOS_LAMBDA (const int &ip)
  {
    double mu = solid.rho[ip]/rho0_ - 1;
    double pH0 = rho0_*c0*c0*mu*(1 + mu)/(1 - (S - 1)*mu)/(1 - (S - 1)*mu);

    if (cp && solid.T[ip] > Tr)
      solid.ienergy[ip] = alpha*(solid.T[ip] - Tr);
    else
       solid.ienergy[ip] = 0;
    pH[ip] = pH0 + Gamma*(solid.ienergy[ip] - e0);

    if (solid.damage[ip] > 0 && pH[ip] < 0)
    {
      if (solid.damage[ip] >= 1.0)
      {
	    // pFinal = rho0_ * Gamma * (e - e0);
	    pH[ip] = 0;
      }
      else
      {
	    // double mu_damaged = (1.0 - damage) * mu;
	    // double pH_damaged = rho0_ * (1.0 - damage) * square(c0) * mu_damaged * (1.0 + mu_damaged) / square(1.0 - (S - 1.0) * mu_damaged);
	    // pFinal = (pH_damaged + rho0_ * (1 + mu_damaged) * Gamma * (e - e0));;
	    pH[ip] *= 1 - solid.damage[ip];
      }
    }

    if (artificial_viscosity)
    {
      double tr_eps = solid.D[ip].trace();

      if (tr_eps < 0)
        pH[ip] += solid.rho[ip]*cellsize*tr_eps*(Q1*cellsize*tr_eps - Q2*c0*Kokkos::Experimental::sqrt(solid.J[ip]));
    }
  });
}

void EOSShock::write_restart(ofstream *of) {
  of->write(reinterpret_cast<const char *>(&rho0_), sizeof(double));
  of->write(reinterpret_cast<const char *>(&K_), sizeof(double));
  of->write(reinterpret_cast<const char *>(&c0), sizeof(double));
  of->write(reinterpret_cast<const char *>(&S), sizeof(double));
  of->write(reinterpret_cast<const char *>(&Gamma), sizeof(double));
  of->write(reinterpret_cast<const char *>(&Tr), sizeof(double));
  of->write(reinterpret_cast<const char *>(&cv), sizeof(double));
  of->write(reinterpret_cast<const char *>(&Q1), sizeof(double));
  of->write(reinterpret_cast<const char *>(&Q2), sizeof(double));
}

void EOSShock::read_restart(ifstream *ifr) {
  if (universe->me == 0) {
    cout << "Restart EOSShock" << endl;
  }
  ifr->read(reinterpret_cast<char *>(&rho0_), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&K_), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&c0), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&S), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&Gamma), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&Tr), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&cv), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&Q1), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&Q2), sizeof(double));
  if (Q1 == 0 && Q2 == 0) {
    artificial_viscosity = false;
  } else {
    artificial_viscosity = true;
  }

  alpha = cv*rho0_;
  e0 = 0;
}

