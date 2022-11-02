/* -*- c++ -*- ----------------------------------------------------------
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

#ifdef METHOD_CLASS

MethodStyle(ulmpm/contacts,ULMPM_CONTACTS)

#else

#ifndef LMP_ULMPM_CONTACTS_H
#define LMP_ULMPM_CONTACTS_H

#include <method.h>
#include <vector>
#include <matrix.h>


class ULMPM_CONTACTS : public Method {
 public:
  ULMPM_CONTACTS(class MPM *);
};

// float linear_basis_function(float, int);
// float derivative_linear_basis_function(float, int, float);
// float cubic_spline_basis_function(float, int);
// float derivative_cubic_spline_basis_function(float, int, float);
// float bernstein_quadratic_basis_function(float, int);
// float derivative_bernstein_quadratic_basis_function(float, int, float);


#endif
#endif
