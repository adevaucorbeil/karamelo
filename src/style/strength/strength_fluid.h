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

  float G();

  void write_restart(ofstream *);
  void read_restart(ifstream *);

  void update_deviatoric_stress(Solid &solid,
                                Kokkos::View<float*> &plastic_strain_increment,
                                Kokkos::View<Matrix3d*> &sigma_dev) const override;
  
protected:
  float G_;
};

#endif
#endif
