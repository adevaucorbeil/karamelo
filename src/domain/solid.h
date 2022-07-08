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
class Solid : public Pointers {
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

  Kokkos::View<tagint*> ptag;               ///< Unique identifier for particles in the system

  Kokkos::View<Vector3d*> x;                ///< Particles' current position
  Kokkos::View<Vector3d*> x0;               ///< Particles' reference position

  
  Kokkos::View<Vector3d*> rp;               ///< Current domain vector (CPDI1)
  Kokkos::View<Vector3d*> rp0;              ///< Reference domain vector (CPDI1)
  Kokkos::View<Vector3d*> xpc;              ///< Current position of the corners of the particles' domain (CPDI2o)
  Kokkos::View<Vector3d*> xpc0;             ///< Reference position of the corners of the particles' domain (CPDI2)
  int nc;                                   ///< Number of corners per particles: \f$2^{dimension}\f$
  
  Kokkos::View<Vector3d*> v;                ///< Particles' current velocity
  Kokkos::View<Vector3d*> v_update;         ///< Particles' velocity at time t+dt

  Kokkos::View<Vector3d*> a;                ///< Particles' acceleration

  Kokkos::View<Vector3d*> mbp;              ///< Particles' external forces times mass
  Kokkos::View<Vector3d*> f;                ///< Particles' internal forces

  Kokkos::View<Matrix3d*> sigma;            ///< Stress matrix
  Kokkos::View<Matrix3d*> strain_el;        ///< Elastic strain matrix
  Kokkos::View<Matrix3d*> vol0PK1;          ///< Transpose of the 1st Piola-Kirchhoff matrix times vol0
  Kokkos::View<Matrix3d*> L;                ///< Velocity gradient matrix
  Kokkos::View<Matrix3d*> F;                ///< Deformation gradient matrix
  Kokkos::View<Matrix3d*> R;                ///< Rotation matrix
  Kokkos::View<Matrix3d*> D;                ///< Symmetric part of L
  Kokkos::View<Matrix3d*> Finv;             ///< Inverse of the deformation gradient matrix
  Kokkos::View<Matrix3d*> Fdot;             ///< Rate of deformation gradient matrix
  Matrix3d Di;                       ///< Inertia tensor
  // vector<Matrix3d> BDinv;            ///< APIC B*Dinv tensor

  Kokkos::View<double*> J;                         ///< Determinant of the deformation matrix
  Kokkos::View<double*> vol0;                      ///< Particles' reference volume
  Kokkos::View<double*> vol;                       ///< Particles' current volume
  Kokkos::View<double*> rho0;                      ///< Particles' reference density
  Kokkos::View<double*> rho;                       ///< Particles' current density
  Kokkos::View<double*> mass;                      ///< Particles' current mass
  Kokkos::View<double*> eff_plastic_strain;        ///< Particles' effective plastic strain
  Kokkos::View<double*> eff_plastic_strain_rate;   ///< Particles' effective plastic strain rate
  Kokkos::View<double*> damage;                    ///< Particles' damage variable
  Kokkos::View<double*> damage_init;               ///< Particles' damage initiation variable
  Kokkos::View<double*> ienergy;                   ///< Particles' internal energy
  Kokkos::View<int*> mask;                         ///< Particles' group mask

  Kokkos::View<double*> T;                         ///< Particles' current temperature
  Kokkos::View<double*> gamma;                     ///< Particles' heat source
  Kokkos::View<Vector3d*> q;                       ///< Particles' heat flux
  
  Kokkos::View<double*> dtCFL;

  Kokkos::View<int**> neigh_n;              ///< Particles' node neighbors
  Kokkos::View<double**> wf;                ///< Particles' node neighbors' weight functions \f$\Phi\f$
  Kokkos::View<double***> wf_corners;       ///< Particles' node neighbors' weight functions \f$\Phi\f$ evaluated at the corners of the particle's domain (used in CPDI)
  Kokkos::View<Vector3d**> wfd;             ///< Particles' node neighbors' weight function derivatives \f$\partial \Phi/ \partial x\f$

  Kokkos::View<int*> error_flag;            ///< Error codes

  Kokkos::MDRangePolicy<Kokkos::Rank<2>> neigh_policy;

  class Mat *mat;                          ///< Pointer to the material

  class Grid *grid;                         ///< Pointer to the background grid

  string method_type;                       ///< Either tlmpm, tlcpdi, tlcpdi2, ulmpm, ulcpdi, or ulcpdi2 (all kinds of MPM supported)

  Solid(class MPM *, vector<string>);       ///< The main task of the constructor is to read the input arguments and launch the creation of particles.
  virtual ~Solid();

  void init();                              ///< Launch the initialization of the grid.
  void options(vector<string> *, vector<string>::iterator); ///< Determines the material and temperature schemes used.
  void grow(int);                           ///< Allocate memory for the vectors used for particles or resize them.
  
  void compute_position_corners();
  
  void compute_inertia_tensor();                    ///< Compute the inertia tensor necessary for the Affice PIC.
  void update_particle_domain();                    ///< Update the particle domain. Used with CPDI

  void copy_particle(int, int);                     ///< Copy particle i attribute and copy it to particle j.
                                                    ///< This function is used to re-order the memory arrangment of particles.
                                                    ///< Usually this is done when particle j is deleted.
  void pack_particle(int, vector<double> &);        ///< Pack particles attributes into a buffer (used for generating a restart).
  void unpack_particle(int &, vector<int>, vector<double> &); ///< Unpack particles attributes from a buffer (used when reading a restart).

  void write_restart(ofstream*);                    ///< Write solid information in the restart file
  void read_restart(ifstream*);                     ///< Read solid information from the restart file

  void populate(vector<string>);

private:
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
  bool update_Di;

public:
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
