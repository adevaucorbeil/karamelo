/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_SOLID_H
#define MPM_SOLID_H

#include "pointers.h"
#include "eos.h"
#include "grid.h"
#include <vector>
#include <Eigen/Eigen>


using namespace Eigen;

class Solid : protected Pointers {
 public:
  string id;                 // solid id
  bigint np;                 // number of particles
  double solidlo[3], solidhi[3]; // solid bounds

  Eigen::Vector3d *x;        // particles' current position
  Eigen::Vector3d *x0;       // particles' reference position

  Eigen::Vector3d *v;        // particles' current position
  Eigen::Vector3d *v_update; // particles' velocity at time t+dt

  Eigen::Vector3d *a;        // particles' acceleration

  Eigen::Vector3d *b;        // particles' external forces
  Eigen::Vector3d *f;        // particles' internal forces

  Eigen::Matrix3d *sigma;            // stress matrix
  Eigen::Matrix3d *PK1;              // 1st Piola-Kirchhoff matrix
  Eigen::Matrix3d *L;                // velocity gradient matrix
  Eigen::Matrix3d *F;                // deformation gradient matrix
  Eigen::Matrix3d *R;                // Rotation matrix
  Eigen::Matrix3d *U;
  Eigen::Matrix3d *Finv;             // inverse of the deformation gradient matrix
  Eigen::Matrix3d *Fdot;             // rate of deformation gradient matrix
  Eigen::Matrix3d *strain_increment; // strain increment matrix

  double *J;                 // determinant of the deformation matrix
  double *vol0;              // particles' reference volume
  double *vol;               // particles' current volume
  double *rho0;              // particles' reference density
  double *rho;               // particles' current density
  double *mass;              // particles' current mass
  int *mask;                 // particles' group mask

  
  int *numneigh_pn;   // number of nodes neighbouring a given particle
  int *numneigh_np;   // number of nodes neighbouring a given node
  vector<int> *neigh_pn;     // List of the nodes neighbouring a given particle
  vector<int> *neigh_np;     // List of the particles neighbouring a given node

  vector< double > *wf_pn;
  vector< double > *wf_np;
  vector< array<double,3> > *wfd_pn;
  vector< array<double,3> > *wfd_np;

  class EOS *eos;                     // Equation-of-State

  class Grid *grid;                   // background grid

  Solid(class MPM *, vector<string>);
  virtual ~Solid();

  void init();
  void options(vector<string> *, vector<string>::iterator);
  void grow(int);

  void compute_mass_nodes();
  void compute_velocity_nodes();
  void compute_external_forces_nodes();
  void compute_internal_forces_nodes();
  void compute_particle_velocities();
  void compute_particle_acceleration();
  void update_particle_position();
  void update_particle_velocities(double);
  void compute_rate_deformation_gradient();
  void update_deformation_gradient();
  void update_stress();
};

#endif
