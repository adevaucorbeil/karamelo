#include <eos_fluid.h>
#include <domain.h>
#include <error.h>
#include <input.h>
#include <math_special.h>
#include <mpm_math.h>
#include <universe.h>
#include <var.h>
#include <matrix.h>
#include <iostream>
#include <math.h>

using namespace std;

using namespace MPM_Math;
using namespace MathSpecial;


EOSFluid::EOSFluid(MPM *mpm, vector<string> args) : EOS(mpm, args)
{
  if (universe->me == 0) {
    cout << "Initiate EOSFluid" << endl;
  }

  if (args.size() < 3) {
    error->all(FLERR, "Error: not enough arguments.\n");
  }

  if (args[2].compare("restart") ==
      0) { // If the keyword restart, we are expecting to have read_restart()
           // launched right after.
    rho0_ = 0;
    K_ = 0;
    Gamma = 0;
    return;
  }

  if (args.size() < 5) {
    error->all(FLERR, "Error: not enough arguments.\n");
  }

  //options(&args, args.begin()+3);
  rho0_ = input->parsev(args[2]);
  K_ = input->parsev(args[3]);
  Gamma = input->parsev(args[4]);

  if (universe->me == 0) {
    cout << "Set rho0 to " << rho0_ << endl;
    cout << "Set K to " << K_ << endl;
    cout << "Set gamma to " << Gamma << endl;
  }
}


EOSFluid::~EOSFluid()
{

}

float EOSFluid::rho0(){
  return rho0_;
}

float EOSFluid::K(){
  return K_;
}

void EOSFluid::compute_pressure(Solid &solid, Kokkos::View<float*> &pH) const
{
  float rho0_ = this->rho0_;
  float K_ = this->K_;
  float Gamma = this->Gamma;

  Kokkos::View<float*> sienergy = solid.ienergy;
  Kokkos::View<float*> srho = solid.rho;

  Kokkos::parallel_for("EOSFluid::compute_pressure", solid.np_local,
  KOKKOS_LAMBDA (const int &ip)
  {
    float mu = srho[ip]/rho0_;
    pH[ip] = K_*(pow(mu, Gamma) - 1);

    sienergy[ip] = 0;
  });
}


void EOSFluid::write_restart(ofstream *of) {
  of->write(reinterpret_cast<const char *>(&rho0_), sizeof(float));
  of->write(reinterpret_cast<const char *>(&K_), sizeof(float));
  of->write(reinterpret_cast<const char *>(&Gamma), sizeof(float));
}

void EOSFluid::read_restart(ifstream *ifr) {
  if (universe->me == 0) {
    cout << "Restart EOSFluid" << endl;
  }
  ifr->read(reinterpret_cast<char *>(&rho0_), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&K_), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&Gamma), sizeof(float));
}

