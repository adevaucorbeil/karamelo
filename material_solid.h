/* -*- c++ -*- ----------------------------------------------------------*/


#ifdef MAT_CLASS

MatStyle(solid,MatSol)

#else

#ifndef MPM_MAT_SOLID_H
#define MPM_MAT_SOLID_H

#include "mat.h"
#include "eos.h"

class MatSol : public Mat {

 public:
  // general material properties
  double rho0;                     // Reference density
  double c0;                       // Reference speed of sound EOS
  double cp;                       // Heat capacity
  double K0;                       // Reference bulk modulus
  double G0;                       // Reference shear modulus
  double thermal_conductivity;     // thermal conductivity
  double thermal_diffusivity;      // thermal diffusivity

  class EOS *eos;                  // user-defined Equation-Of-State
  //class Strength *strength;        // user-defined strength
  //class Damage *damage;            // user-defined damage

  MatSol(class MPM *, vector<string>);
  ~MatSol();

};

#endif
#endif
