/* -*- c++ -*- ----------------------------------------------------------*/


#ifdef EOS_CLASS

EOSStyle(fluid,EOSFluid)

#else

#ifndef MPM_EOS_FLUID_H
#define MPM_EOS_FLUID_H

#include <eos.h>
#include <matrix.h>

/*! \ingroup eos eosfluid eos_fluid

\section Syntax Syntax
\code
eos(eos-ID, fluid, rho0, K, gamma)
\endcode

<ul>
<li>eos-ID: name of the eos to be created.</li>
<li>rho0: reference bulk density</li>
<li>K: bulk modulus. </li>
<li>gamma: Tait exponent. </li>
</ul>

\section Examples Examples
\code
eos(eosf, fluid, rho, K, gamma)
\endcode
Defines a non-linear EOS called 'eosf' typically used for the simulation of fluids in MPM 
with a reference bulk density rho, a bulk modulus K and the Tait coefficient gamma.

\section Description Description

This command defines a non-linear equation of state. The hydrostatic pressure is calculated according to:\n
\f[ p = K \left[\left(\frac{\rho}{\rho_0}\right)^{\gamma} - 1\right]\f]

Note that a typical value is \f$\gamma\f$ = 7 for hydraulic simulations of water.

\section Class Class description
*/

class EOSFluid : public EOS {

public:
  EOSFluid(class MPM *, vector<string>);
  ~EOSFluid();

  float rho0();
  float K();
  float G();
  void compute_pressure(Solid &solid, Kokkos::View<float*> &pH) const override;
  void write_restart(ofstream *);
  void read_restart(ifstream *);

protected:
  float rho0_, K_, Gamma;
};

#endif
#endif
