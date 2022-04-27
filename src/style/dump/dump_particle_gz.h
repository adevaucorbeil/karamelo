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

#ifdef DUMP_CLASS

DumpStyle(particle/gz,DumpParticleGz)

#else

#ifndef MPM_DUMP_PARTICLE_GZ_H
#define MPM_DUMP_PARTICLE_GZ_H

#include <dump.h>
#include <grid.h>

class DumpParticleGz : public Dump {
  Kokkos::View<tagint*>::HostMirror ptag;               ///< Unique identifier for particles in the system

  Kokkos::View<Vector3d*>::HostMirror x;                ///< Particles' current position
  Kokkos::View<Vector3d*>::HostMirror x0;               ///< Particles' reference position
  
  Kokkos::View<Vector3d*>::HostMirror v;                ///< Particles' current velocity

  Kokkos::View<Vector3d*>::HostMirror mbp;              ///< Particles' external forces times mass

  Kokkos::View<Matrix3d*>::HostMirror sigma;            ///< Stress matrix
  Kokkos::View<Matrix3d*>::HostMirror strain_el;        ///< Elastic strain matrix
  Kokkos::View<Matrix3d*>::HostMirror R;                ///< Rotation matrix

  Kokkos::View<double*>::HostMirror vol;                       ///< Particles' current volume
  Kokkos::View<double*>::HostMirror mass;                      ///< Particles' current mass
  Kokkos::View<double*>::HostMirror eff_plastic_strain;        ///< Particles' effective plastic strain
  Kokkos::View<double*>::HostMirror eff_plastic_strain_rate;   ///< Particles' effective plastic strain rate
  Kokkos::View<double*>::HostMirror damage;                    ///< Particles' damage variable
  Kokkos::View<double*>::HostMirror damage_init;               ///< Particles' damage initiation variable
  Kokkos::View<double*>::HostMirror ienergy;                   ///< Particles' internal energy

  Kokkos::View<double*>::HostMirror T;                         ///< Particles' current temperature
  Kokkos::View<double*>::HostMirror gamma;                     ///< Particles' heat source

 public:
  DumpParticleGz(MPM *, vector<string>);
  ~DumpParticleGz();

  void write();
  protected:
  vector<string> known_var = {"x", "y", "z",
			      "x0", "y0", "z0",
			      "vx", "vy", "vz",
			      "s11", "s22", "s33",
			      "s12", "s13", "s23",
			      "e11", "e22", "e33",
			      "e12", "e13", "e23",
			      "seq", "volume", "mass",
			      "damage", "damage_init",
			      "bx", "by", "bz",
			      "ep", "epdot", "T",
			      "ienergy", "gamma"};
};

#endif
#endif
