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

EOSStyle(linear,EOSLinear)

#else

#ifndef MPM_EOS_LINEAR_H
#define MPM_EOS_LINEAR_H

#include <eos.h>
#include <matrix.h>


/*! \ingroup eos eoslinear eos_linear

\section Syntax Syntax
\code
eos(eos-ID, linear, rho0, K)
\endcode

<ul>
<li>eos-ID: name of the eos to be created.</li>
<li>rho0: reference bulk density</li>
<li>K: bulk modulus. </li>
</ul>

\section Examples Examples
\code
eos(eosl, linear, rho, K)
\endcode
Defines a linear EOS called 'eosl' with a reference bulk density rho, a bulk modulus K

\section Description Description

This command defines a linear equation of state. The hydrostatic pressure is calculated according to:\n
\f[ p = K \left(\frac{\rho}{\rho_0} - 1\right)\f]

\section Class Class description
*/


class EOSLinear : public EOS {

public:
  EOSLinear(class MPM *, vector<string>);
  ~EOSLinear();

  double rho0();
  double K();
  double G();
  void compute_pressure(Solid &solid, Kokkos::View<double*> &pH) const override;
  void write_restart(ofstream *);
  void read_restart(ifstream *);

protected:
  double rho0_, K_;
};

#endif
#endif
