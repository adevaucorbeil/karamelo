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


#ifdef EOS_CLASS

EOSStyle(shock,EOSShock)

#else

#ifndef MPM_EOS_SHOCK_H
#define MPM_EOS_SHOCK_H

#include "eos.h"
#include <Eigen/Eigen>

/*! \ingroup eos eoslinear eos_linear

\section Syntax Syntax
\code
eos(eoss, shock, rho0, K, c0, S, Gamma0, cv, Tr, Q1, Q2)
\endcode
Defines a Mie-Gruneisen equation of state (EOSShock) called 'eoss' with a reference bulk density rho, a bulk modulus K, a sound celerity c0,
a Hugoniot slope coefficient S, a Grüneisen Gamma at the reference state Gamma, a heat capacity at constant
volume cv, temperature of reference Tr, and artificial viscosity coefficients Q1 and Q2.

<ul>
<li>eos-ID: name of the eos to be created.</li>
<li>rho0: reference bulk density.</li>
<li>K: bulk modulus. </li>
<li>c0: the bulk speed of sound.</li>
<li>S: the Hugoniot slope coefficient</li>
<li>Gamma0: the Grüneisen Gamma</li>
<li>cv: heat capacity at constant volume</li>
<li>Tr: temperature of refrerence.</li>
<li>Q1: linear artificial velocity coefficient.</li>
<li>Q2: quadratic artificial velocity coefficient.</li>
</ul>

\section Examples Examples

\code
eos(eoss, shock, rho0, K, c0, S, Gamma, cv, Tr, Q1, Q2)
\endcode
Defines a Mie-Gruneisen equation of state (EOSShock) called 'eoss' with a reference bulk density rho, a bulk modulus K, a bulk speed of sound c0,
a Hugoniot slope coefficient S, a Grüneisen Gamma at the reference state Gamma, a heat capacity at constant
volume cv, temperature of reference Tr, and artificial viscosity coefficients Q1 and Q2.

\section Description Description

This command defines a Mie-Grüneisen equation of state. The hydrostatic pressure is calculated according to:\n
\f[ p = \dfrac{\rho_0c_0^2(\eta -1)\bigg[\eta - \dfrac{\Gamma_0}{2}(\eta - 1)\bigg]}{[\eta - S(\eta-1)]^2} + \Gamma_0e\f]
where \f$\eta = \dfrac{\rho}{\rho_0}\f$ and \f$e\f$ the internal energy.

\section Class Class description
*/

class EOSShock : public EOS {

public:
  EOSShock(class MPM *, vector<string>);
  ~EOSShock();

  double rho0();
  double K();
  double G();
  void compute_pressure(double &, double &, const double, const double, const double, const Eigen::Matrix3d, const double, const double T = 0);

protected:
  double rho0_, K_, e0, c0, S, Gamma, Tr, cv, alpha, Q1, Q2;
  string usage = "Usage: eos(eos-ID, shock, rho, K, c0, S, Gamma, cv, Tr, Q1, Q2)\n";
  int Nargs = 11;

private:
  bool artificial_viscosity;
};

#endif
#endif
