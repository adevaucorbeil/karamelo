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
  Eigen::Matrix3d *Fdot;             // rate of deformation gradient matrix
  Eigen::Matrix3d *strain_increment; // strain increment matrix

  double *J;                 // determinant of the deformation matrix
  double *vol0;              // particles' reference volume
  double *vol;               // particles' current volume
  double *mass;              // particles' current mass
  int *mask;                 // particles' group mask

  class EOS *eos;                     // Equation-of-State

  class Grid *grid;                   // background grid

  Solid(class MPM *, vector<string>);
  virtual ~Solid();
  virtual void init();
  void options(vector<string> *, vector<string>::iterator);
  void grow(int);

};

#endif
