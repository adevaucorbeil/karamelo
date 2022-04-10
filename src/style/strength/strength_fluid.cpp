#include <strength_fluid.h>
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


StrengthFluid::StrengthFluid(MPM *mpm, vector<string> args) : Strength(mpm, args)
{
  if (universe->me == 0) {
    cout << "Initiate StrengthFluid" << endl;
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
    cout << "Fluid strength model:\n";
    cout << "\tmu: shear viscosity " << G_ << endl;
  }
}

double StrengthFluid::G(){
  return G_;
}


void
StrengthFluid::update_deviatoric_stress(Solid &solid,
                                        Kokkos::View<double*, MemorySpace> &plastic_strain_increment,
                                        Kokkos::View<Matrix3d*, MemorySpace> &sigma_dev) const
{
  double G_ = this->G_;

  Kokkos::parallel_for("EOSLinear::compute_pressure", solid.np_local,
  KOKKOS_LAMBDA (const int &ip)
  {
    sigma_dev[ip] = 2*G_*Deviator(solid.D[ip]);
  });
}

void StrengthFluid::write_restart(ofstream *of) {
  of->write(reinterpret_cast<const char *>(&G_), sizeof(double));
}

void StrengthFluid::read_restart(ifstream *ifr) {
  if (universe->me == 0) {
    cout << "Restart StrengthFluid" << endl;
  }
  ifr->read(reinterpret_cast<char *>(&G_), sizeof(double));
}
