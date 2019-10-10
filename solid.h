/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_SOLID_H
#define MPM_SOLID_H

#include "pointers.h"
#include "material.h"
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

  
  Eigen::Vector3d *rp;       // current domain vector (CPDI1)
  Eigen::Vector3d *rp0;      // reference domain vector (CPDI1)
  Eigen::Vector3d *xpc;      // current position of the corners of the particles' domain (CPDI2)
  Eigen::Vector3d *xpc0;     // reference position of the corners of the particles' domain (CPDI2)
  int nc;                    // number of corners per particles: 2^dimension
  
  Eigen::Vector3d *v;        // particles' current velocity
  Eigen::Vector3d *v_update; // particles' velocity at time t+dt

  Eigen::Vector3d *a;        // particles' acceleration

  Eigen::Vector3d *mb;        // particles' external forces times mass
  Eigen::Vector3d *f;        // particles' internal forces

  Eigen::Matrix3d *sigma;            // stress matrix
  Eigen::Matrix3d *strain_el;        // elastic strain matrix
  // Eigen::Matrix3d *PK1;              // 1st Piola-Kirchhoff matrix
  Eigen::Matrix3d *vol0PK1;          // Transpose of the 1st Piola-Kirchhoff matrix times vol0
  Eigen::Matrix3d *L;                // velocity gradient matrix
  Eigen::Matrix3d *F;                // deformation gradient matrix
  Eigen::Matrix3d *R;                // Rotation matrix
  Eigen::Matrix3d *U;
  Eigen::Matrix3d *D;                // symmetric part of L
  Eigen::Matrix3d *Finv;             // inverse of the deformation gradient matrix
  Eigen::Matrix3d *Fdot;             // rate of deformation gradient matrix
  Eigen::Matrix3d *Di;               // inertia tensor

  double *J;                         // determinant of the deformation matrix
  double *vol0;                      // particles' reference volume
  double *vol;                       // particles' current volume
  double vtot;                       // total volume
  double *rho0;                      // particles' reference density
  double *rho;                       // particles' current density
  double *mass;                      // particles' current mass
  double *eff_plastic_strain;        // particles' effective plastic strain
  double *eff_plastic_strain_rate;   // particles' effective plastic strain rate
  double *damage;                    // particles' damage variable
  double *damage_init;               // particles' damage initiation variable
  int *mask;                         // particles' group mask

  double min_inv_p_wave_speed;   // minimum of the inverse of the particle wave speed
  double dtCFL;
  
  int *numneigh_pn;   // number of nodes neighbouring a given particle
  int *numneigh_np;   // number of nodes neighbouring a given node
  vector<int> *neigh_pn;     // List of the nodes neighbouring a given particle
  vector<int> *neigh_np;     // List of the particles neighbouring a given node

  vector< double > *wf_pn;
  vector< double > *wf_np;
  vector< double > *wf_pn_corners;

  vector< Eigen::Vector3d > *wfd_pn;
  vector< Eigen::Vector3d > *wfd_np;


  struct Mat *mat;                     // Material

  class Grid *grid;                   // background grid

  string method_type;

  Solid(class MPM *, vector<string>);
  virtual ~Solid();

  void init();
  void options(vector<string> *, vector<string>::iterator);
  void grow(int);

  void compute_mass_nodes(bool);
  void compute_velocity_nodes(bool);
  void compute_velocity_nodes_APIC(bool);
  void compute_external_forces_nodes(bool);
  void compute_internal_forces_nodes_TL();
  void compute_internal_forces_nodes_UL(bool);
  void compute_particle_velocities();
  void compute_particle_acceleration();
  void update_particle_position();
  void update_particle_velocities(double);
  void compute_rate_deformation_gradient_TL();
  void compute_rate_deformation_gradient_TL_APIC();
  void compute_rate_deformation_gradient_UL_USL();
  void compute_rate_deformation_gradient_UL_MUSL();
  void compute_rate_deformation_gradient_UL_APIC();
  void update_deformation_gradient();
  void update_stress();
  void compute_inertia_tensor(string);
  void compute_deformation_gradient();
  void copy_particle(int, int);
  void update_particle_domain();

private:
  void populate(vector<string>);
};

#endif
