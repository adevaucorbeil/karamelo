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

#ifndef MPM_SOLID_H
#define MPM_SOLID_H

#include <pointers.h>
#include <material.h>
#include <grid.h>
#include <vector>
#include <deque>
#include <matrix.h>




/*! This class represents a given solid.
 * 
 * In this class is stored the variables related to the solid boundaries as well as all those related to the particles it contains.
 * Moreover, it holds the evaluation of the shape functions and its derivatives.
 * The functions member of this class are related to the Particle to Grid (P2G) and Grid to Particle (G2P) steps of the MPM algorithm.
 * Another member is a function that updates the particles' stresses.
 * This class holds a pointer of type Grid that points to either the local grid, as used in the Total Lagrangian MPM, or the global grid,
 * as used in the Updated Lagrangian MPM.
 * This class has Pointers as parent in order to be accessible from anywhere within the MPM class.
 */
class Solid : protected Pointers {
 public:
  string id;                                ///< Solid id

  double solidlo[3], solidhi[3];            ///< Solid global bounds
  double solidsublo[3], solidsubhi[3];      ///< Solid local bounds

  bigint np;                                ///< Total number of particles in the domain
  int np_local;                             ///< Number of local particles (in this CPU)
  int np_per_cell;                          ///< Number of particles per cell (at the beginning)
  int comm_n;                               ///< Number of double to pack for particle exchange between CPU
  double vtot;                              ///< Total volume
  double mtot;                              ///< Total mass

  vector<tagint> ptag;                      ///< Unique identifier for particles in the system

  vector<Vector3d> x;                ///< Particles' current position
  vector<Vector3d> x0;               ///< Particles' reference position

  
  vector<Vector3d> rp;               ///< Current domain vector (CPDI1)
  vector<Vector3d> rp0;              ///< Reference domain vector (CPDI1)
  vector<Vector3d> xpc;              ///< Current position of the corners of the particles' domain (CPDI2o)
  vector<Vector3d> xpc0;             ///< Reference position of the corners of the particles' domain (CPDI2)
  int nc;                                   ///< Number of corners per particles: \f$2^{dimension}\f$
  
  vector<Vector3d> v;                ///< Particles' current velocity
  vector<Vector3d> v_update;         ///< Particles' velocity at time t+dt

  vector<Vector3d> a;                ///< Particles' acceleration

  vector<Vector3d> mbp;              ///< Particles' external forces times mass
  vector<Vector3d> f;                ///< Particles' internal forces

  vector<Matrix3d> sigma;            ///< Stress matrix
  vector<Matrix3d> strain_el;        ///< Elastic strain matrix
  vector<Matrix3d> vol0PK1;          ///< Transpose of the 1st Piola-Kirchhoff matrix times vol0
  vector<Matrix3d> L;                ///< Velocity gradient matrix
  vector<Matrix3d> F;                ///< Deformation gradient matrix
  vector<Matrix3d> R;                ///< Rotation matrix
  vector<Matrix3d> D;                ///< Symmetric part of L
  vector<Matrix3d> Finv;             ///< Inverse of the deformation gradient matrix
  vector<Matrix3d> Fdot;             ///< Rate of deformation gradient matrix
  Matrix3d Di;                       ///< Inertia tensor
  // vector<Matrix3d> BDinv;            ///< APIC B*Dinv tensor

  vector<double> J;                         ///< Determinant of the deformation matrix
  vector<double> vol0;                      ///< Particles' reference volume
  vector<double> vol;                       ///< Particles' current volume
  vector<double> rho0;                      ///< Particles' reference density
  vector<double> rho;                       ///< Particles' current density
  vector<double> mass;                      ///< Particles' current mass
  vector<double> eff_plastic_strain;        ///< Particles' effective plastic strain
  vector<double> eff_plastic_strain_rate;   ///< Particles' effective plastic strain rate
  vector<double> damage;                    ///< Particles' damage variable
  vector<double> damage_init;               ///< Particles' damage initiation variable
  vector<double> ienergy;                   ///< Particles' internal energy
  vector<int> mask;                         ///< Particles' group mask

  vector<double> T;                         ///< Particles' current temperature
  vector<double> gamma;                     ///< Particles' heat source
  vector<Vector3d> q;                ///< Particles' heat flux

  double max_p_wave_speed;                  ///< Maximum of the particle wave speed
  double dtCFL;

  deque<int> neigh_p;             ///< List of particle indices for neighbour pairs
  deque<int> neigh_n;             ///< List of node indices for neighbour pairs
  deque<double> wf;               ///< List of the weight functions \f$\Phi\f$ for neighbour pairs
  deque<double> wf_corners;               ///< List of the weight functions \f$\Phi\f$ for neighbour pairs evaluated at the corners of the particle's domain (used in CPDI)
  deque<Vector3d> wfd;               ///< List ofthe derivative of the weight functions \f$\partial \Phi/ \partial x\f$ for neighbour pairs

  class Mat *mat;                          ///< Pointer to the material

  class Grid *grid;                         ///< Pointer to the background grid

  string method_type;                       ///< Either tlmpm, tlcpdi, tlcpdi2, ulmpm, ulcpdi, or ulcpdi2 (all kinds of MPM supported)

  Solid(class MPM *, vector<string>);       ///< The main task of the constructor is to read the input arguments and launch the creation of particles.
  virtual ~Solid();

  void init();                              ///< Launch the initialization of the grid.
  void options(vector<string> *, vector<string>::iterator); ///< Determines the material and temperature schemes used.
  void grow(int);                           ///< Allocate memory for the vectors used for particles or resize them.

  void compute_mass_nodes(int in, int ip, double wf);                    ///< Compute nodal mass step of the Particle to Grid step of the MPM algorithm.
  void compute_velocity_nodes(int in, int ip, double wf, bool APIC);                ///< Compute nodal velocity (via momentum) step of the Particle to Grid step of the MPM algorithm.
  void compute_force_nodes(int in, int ip, double wf, const Vector3d &wfd, bool TL, bool MLS);      ///< Compute both external and internal forces step of the Particle to Grid step of the updated Lagrangian MPM algorithm.
  void compute_particle(bool positions, bool velocities, bool accelerations); ///< Compute the particles' temporary velocities and position, part of the Grid to Particles step of the MPM algorithm.
  void update_particle(double alpha, bool positions, bool velocities);          ///< Update the particles' velocities based on either PIC and/or FLIP and update the positions using the updated velocities.
                                                    ///< The argument is the ratio \f$\alpha\f$ used between PIC and FLIP.
                                                    ///< \f$\alpha = 0\f$ for pure PIC, \f$\alpha = 1\f$ for pure FLIP.
  void compute_rate_deformation_gradient(bool doublemapping, bool TL, bool APIC);      ///< Compute the time derivative of the deformation matrix for TLMPM, when APIC is not used.
  void update_deformation_gradient();               ///< Update the deformation gradient, volume, density, and the necessary strain matrices
  void update_stress();                             ///< Calculate the stress, damage and temperature at each particle, and determine the maximum allowed time step.
  void compute_inertia_tensor();                    ///< Compute the inertia tensor necessary for the Affice PIC.
  void compute_deformation_gradient();              ///< Compute the deformation gradient directly from the grid nodes' positions
  void update_particle_domain();                    ///< Update the particle domain. Used with CPDI

  void copy_particle(int, int);                     ///< Copy particle i attribute and copy it to particle j.
                                                    ///< This function is used to re-order the memory arrangment of particles.
                                                    ///< Usually this is done when particle j is deleted.
  void pack_particle(int, vector<double> &);        ///< Pack particles attributes into a buffer (used for generating a restart).
  void unpack_particle(int &, vector<int>, vector<double> &); ///< Unpack particles attributes from a buffer (used when reading a restart).

  void write_restart(ofstream*);                    ///< Write solid information in the restart file
  void read_restart(ifstream*);                     ///< Read solid information from the restart file

  void compute_temperature_nodes(bool);             ///< Compute nodal temperature step of the particle
  void compute_external_temperature_driving_forces_nodes(bool); ///< Compute external temperature driving forces
  void compute_internal_temperature_driving_forces_nodes(); ///< Compute internal forces step of the Particle to Grid step of the total Lagrangian MPM algorithm.
  void update_particle_temperature();               ///< Update the particles' temperature
  void update_heat_flux(bool);                      ///< Update the particles' heat source and fluxes

private:
  void populate(vector<string>);
  void read_mesh(string);
  void read_file(string);

  const map<string, string> usage ={
    {"region", "Usage: solid(solid-ID, \033[1;32mregion\033[0m, region-ID, N_ppc1D, material-ID, cell-size, T0)\n"},
    {"mesh",    "Usage: solid(solid-ID, \033[1;32mmesh\033[0m, meshfile, material-ID, h, T0)\n"},
    {"file",    "Usage: solid(solid-ID, \033[1;32mfile\033[0m, filename, material-ID, h, T0)\n"}
  };
  const map<string, int>    Nargs = {
             {"region", 7},
	     {"mesh",   6},
             {"file",   6}
           };
  
  double T0;                     ///< Initial temperature
  bool is_TL, apic;              ///< Boolean variables that are true if using total Lagrangian MPM, and APIC, respectively
};

#endif

/*! \defgroup solid solid

\section Syntax Syntax
\code
solid(solid-ID, region, region-ID, N_ppc1D, material-ID, h, T0)
solid(solid-ID, mesh, meshfile, material-ID, h, T0)
solid(solid-ID, file, filename, material-ID, h, T0)
\endcode

<ul>
<li>solid-ID: name of the solid to be created.</li>
<li>region-ID: name of the region that occupies the solid.</li>
<li>N_ppc1D: equivalent number of particles per cell in 1D (along one axis).</li>
<li>material-ID: name of the material to be used by the solid.</li>
<li>h: cell-size of the background grid.</li>
<li>T0: temperature of the solid at the start of the simulation.</li>
<li>meshfile: name of the file containing the mesh.</li>
<li>filename: name of the file containing the list of particles as well as their positions.</li>
</ul>

\section Examples Examples
\code
E        = 115
nu       = 0.31
rho      = 8.94e-06
ppc      = 1
R        = 0.2
xc       = R/2
yc       = R/2
zc       = 0
region(rBall, sphere, xc, yc, zc, R)
material(mat1, linear, rho, E, nu)
solid(sBall, region, rBall, ppc, mat1, R / 10, 25)
\endcode
Defines a solid called 'sBall' delimited by region 'rBall', a sphere (Sphere) of radius 
\f$R=0.2\f$ centered around (xc, yc, zc) which constitutive law (Material) is linear elastic.

\section Description Description

A solid is a collection of particles. These particles are created according to either:
<ol>
<li>a region: the region delimits the extends of the solid, and the particles are created within it,</li>
<li>a mesh: the particles correspond to the integration points of the mesh,</li>
<li>a file: the particles taken from the list found in the file.</li>
</ol>

When the solid is delimited by a region, the number of particles occupying a given background grid cell is given by the N_ppc1D parameter.
*/
