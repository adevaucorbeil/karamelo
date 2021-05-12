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

#ifndef MPM_MATERIAL_H
#define MPM_MATERIAL_H

#include "damage.h"
#include "eos.h"
#include "pointers.h"
#include "strength.h"
#include "temperature.h"
#include <vector>

/*! This class stores the information concerning a given material.
 * 
 * The simplest material is rigid for which Mat::rigid == true. 
 * The other pointers are set to NULL.
 * The corresponding user command is: material(material-ID, rigid).\n\n
 * The second type of material is linear which corresponds to linear elasticity.\n
 * The corresponding user command is: material(material-ID, linear, rho, E, nu, optional: damage-ID)
 * where rho is the reference density, E the Young's modulus, nu the Poisson's ratio, and damage-ID.
 * the ID of the optional damage law.\n\n
 * The third type of material is neo-hookean which corresponds to neo-hookean elasticity.
 * The corresponding user command is similar to that for linear elasticiy materials: \n
 * material(material-ID, neo-hookean, rho, E, nu, optional: damage-ID).\n\n
 * The fourth and last type of material is eos-strength which corresponds to an elasto-plastic material
 * the hydrostatic pressure is determined using an equation of state (EOS) while the deviatoric stress is
 * determined using the shear modulus and a flow stress model.\n
 * The corresponding user command is: material(material-ID, eos-strength, eos-ID, strength-ID, optional: 
 * damage-ID, optional: temperature-ID)
 * where eos-ID, strength-ID, damage-ID and temperature-ID are the ID of the EOS, flow stress model,
 * damage law and temperature law, respectively.\n
 * The damage and temperature laws are optional.
 */
class Mat {
public:
  string id;                                         ///< Identification name of the material
  int type;                                          ///< Either RIGID, LINEAR, NEO_HOOKEAN, or SHOCK with values from  Material::constitutive_model
  bool rigid = false;                                ///< True if the material is rigid, false otherwise
  class EOS *eos = NULL;                             ///< Pointer to the EOS
  class Strength *strength = NULL;                   ///< Pointer to the Strength (flow stress rule)
  class Damage *damage = NULL;                       ///< Pointer to the Damage law
  class Temperature *temp = NULL;                    ///< Pointer to the Temperature law
  double rho0;                                       ///< Density in the reference state \f$\rho_0\f$
  double E;                                          ///< Young's modulus
  double nu;                                         ///< Poisson's ratio \f$\nu\f$
  double G;                                          ///< Shear modulus
  double K;                                          ///< Bulk modulus
  double lambda;                                     ///< 1st Lame parameter \f$\lambda\f$
  double signal_velocity;                            ///< Signal velocity in the reference state \f$c_0 = \sqrt{\frac{\lambda+2G}{\rho_0}}\f$
  double cp;                                         ///< Heat capacity at constant pressure
  double invcp;                                      ///< Inverse of the heat capacity at constant pressure
  double kappa;                                      ///< Thermal conductivity
  Mat(string, int, class EOS *, class Strength *, class Damage *,
      class Temperature *);                          ///< Creates an elasto-plastic material
  Mat(string, int, double, double, double, double, double); ///< Creates a linear or Neo-Hookean material
  Mat(string, int);                                  ///< Creates a rigid material
};

/*! Stores all the user defined Equations of State, elasto-plastic, damage, and temperature
 *  laws as well as the different materials which are a combination of the formers.
 */
class Material : protected Pointers {
 public:
  vector<Mat> materials;                        ///< List of defined materials
  vector<class EOS *> EOSs;                     ///< List of defined Equations of State
  vector<class Strength *> strengths;           ///< List of defined Strengths
  vector<class Damage *> damages;               ///< List of defined Damage laws
  vector<class Temperature *> temperatures;     ///< List of defined Temperature laws
  
  Material(class MPM *);                        ///< Creator
  virtual ~Material();                          ///< Destructor

  void add_strength(vector<string>);            ///< Create a new Strength
  int find_strength(string);                    ///< Finds the ID of a Strength
  void add_EOS(vector<string>);                 ///< Create a new EOS
  int find_EOS(string);                         ///< Finds the ID of an EOS
  void add_material(vector<string>);            ///< Create a new Material
  int find_material(string);                    ///< Finds the ID of a Material
  void add_damage(vector<string>);              ///< Create a new Damage
  int find_damage(string);                      ///< Finds the ID of a Damage
  void add_temperature(vector<string>);         ///< Create a new Temperature
  int find_temperature(string);                 ///< Finds the ID of a Temperature

  void write_restart(ofstream*);                ///< Write restart
  void read_restart(ifstream*);                 ///< Read restart

  typedef Strength *(*StrengthCreator)(MPM *,vector<string>);
  typedef map<string,StrengthCreator> StrengthCreatorMap;
  StrengthCreatorMap *strength_map;

  typedef EOS *(*EOSCreator)(MPM *,vector<string>);
  typedef map<string,EOSCreator> EOSCreatorMap;
  EOSCreatorMap *EOS_map;

  typedef Damage *(*DamageCreator)(MPM *,vector<string>);
  typedef map<string,DamageCreator> DamageCreatorMap;
  DamageCreatorMap *damage_map;

  typedef Temperature *(*TemperatureCreator)(MPM *,vector<string>);
  typedef map<string,TemperatureCreator> TemperatureCreatorMap;
  TemperatureCreatorMap *temperature_map;

  enum constitutive_model {
    RIGID = 0,       ///< Rigid material
    LINEAR = 1,      ///< Linear elasticity
    NEO_HOOKEAN = 2, ///< Neo-Hookean model
    SHOCK = 3,       ///< EOS + Strength or fluids
  };

private:
  template <typename T> static Strength *strength_creator(MPM *,vector<string>);
  template <typename T> static EOS *EOS_creator(MPM *,vector<string>);
  template <typename T> static Damage *damage_creator(MPM *,vector<string>);
  template <typename T> static Temperature *temperature_creator(MPM *,vector<string>);

  const map<string, string> usage = {{"rigid",        "Usage: material(material-ID, \033[1;32mrigid\033[0m)\n"},
				     {"linear",       "Usage: material(material-ID, \033[1;32mlinear\033[0m, rho, E, nu, optional: cp, optional: kappa, optional: damage-ID)\n"},
				     {"neo-hookean",  "Usage: material(material-ID, \033[1;32mneo-hookean\033[0m, rho, E, nu, optional: cp, optional: kappa, optional: damage-ID)\n"},
				     {"eos-strength", "Usage: material(material-ID, \033[1;32meos-strength\033[0m, eos-ID, strength-ID, optional: damage-ID, optional: temperature-ID)\n"}};
  const map<string, int>    Nargs = {{"rigid",        2},
				     {"linear",       5},
				     {"neo-hookean",  5},
				     {"eos-strength", 4}};
};


#endif

/*! \defgroup material material

\section Syntax Syntax
\code
material(material-ID, rigid)
material(material-ID, linear, rho, E, nu, optional: damage-ID)
material(material-ID, neo-hookean, rho, E, nu, optional: damage-ID)
material(material-ID, eos-strength, eos-ID, strength-ID, optional: damage-ID, optional: temperature-ID)
\endcode

<ul>
<li>material-ID: name of the material to be created.</li>
<li>rho: density of the material.</li>
<li>E: Young's modulus.</li>
<li>nu: Poisson's ration.</li>
<li>damage-ID: name of the damage law to be used (see Damage).</li>
<li>temperature-ID: name of the temperature law to be used (see Temperature).</li>
<li>eos-ID: name of EOS (Equation of State) to be used.</li>
<li>strength-ID: name of hardening law to be used (see Strength).</li>
</ul>

\section Examples Examples
\code
E        = 115
nu       = 0.31
rho      = 8.94e-06
material(mat1, linear, rho, E, nu)
\endcode
Defines a linear elastic material called 'mat1'. Its density is 8.94e-06, Poisson's ration 0.31 and Young's modulus 115.

\section Description Description

Four different kinds of materials are supported:
<ul>
<li>Rigid. The material does not deform.</li>
<li>Linear, defined by their density, Young's modulus \f$E\f$ and Poisson's ratio \f$\nu\f$. In this case, the stress is directly related to the strain increment:\n
\f{equation}{
  \boldsymbol{\sigma} = 2G \Delta \boldsymbol{\varepsilon} + \lambda \text{tr}(\Delta \boldsymbol{\varepsilon}) \boldsymbol{I}
\f}
where the shear modulus \f$G\f$ and Lame parameter \f$\lambda\f$ are obtained as follows:
\f{equation}{
G = \frac{E}{2 (1 + \nu)}\\
\lambda = \frac{E  \nu}{(1 + \nu)  (1 - 2 \nu)}
\f}</li>
<li>Neo-Hookean, defined by their density, Young's modulus \f$E\f$ and Poisson's ratio \f$\nu\f$. In this model, the first Piola-Kirchhoff stress tensor is given by:
\f{equation}{
  \boldsymbol{P} = G (\boldsymbol{F} - \boldsymbol{F}^{-T}) + \lambda \left( \ln J \right)\boldsymbol{F}^{-T}
\f} where \f$J\f$ is the determinant of \f$\boldsymbol{F}\f$. </li>
<li>Elasto-plastic: defined by an equation of state (EOS) and a hardening law (Strength).</li>

For all non-rigid materials, damage (Damage) is supported, to add allow for the use of a damage law, its ID, damage-ID, needs to be given.

Finally, for elasto-plastic materials, the increase of temperature coming from the plastic work can be talken into account using a temperature law (Temperature). 

*/
