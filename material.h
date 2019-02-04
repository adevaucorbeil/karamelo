/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_MATERIAL_H
#define MPM_MATERIAL_H

#include "pointers.h"
#include <vector>

class Material : protected Pointers {
 public:
  string id;

  // general material properties
  double rho0;                     // Reference density
  double c0;                       // Reference speed of sound EOS
  double cp;                       // Heat capacity
  double K0;                       // Reference bulk modulus
  double G0;                       // Reference shear modulus
  double thermal_conductivity;     // thermal conductivity
  double thermal_diffusivity;      // thermal diffusivity
  
  Material(class MPM *, vector<string>);
  virtual ~Material();
  void options(vector<string> *, vector<string>::iterator);
};

#endif
