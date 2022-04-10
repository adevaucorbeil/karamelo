/* -*- c++ -*- ----------------------------------------------------------*/


#ifdef STRENGTH_CLASS

StrengthStyle(fluid,StrengthFluid)

#else

#ifndef MPM_STRENGTH_FLUID_H
#define MPM_STRENGTH_FLUID_H

#include <strength.h>
#include <matrix.h>

class StrengthFluid : public Strength {

public:
  StrengthFluid(class MPM *, vector<string>);
  ~StrengthFluid() {};

  double G();

  void write_restart(ofstream *);
  void read_restart(ifstream *);

  void update_deviatoric_stress(Solid &solid,
                                Kokkos::View<double*, MemorySpace> &plastic_strain_increment,
                                Kokkos::View<Matrix3d*, MemorySpace> &sigma_dev) const override;
  
protected:
  double G_;
};

#endif
#endif
