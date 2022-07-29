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

#include <deque>
#include <dump.h>
#include <thread>
#include <utility>

class DumpParticleGz : public Dump {
  deque<Kokkos::View<tagint*>::HostMirror> ptag;               ///< Unique identifier for particles in the system

  deque<Kokkos::View<Vector3d*>::HostMirror> x;                ///< Particles' current position
  deque<Kokkos::View<Vector3d*>::HostMirror> x0;               ///< Particles' reference position
  
  deque<Kokkos::View<Vector3d*>::HostMirror> v;                ///< Particles' current velocity

  deque<Kokkos::View<Vector3d*>::HostMirror> mbp;              ///< Particles' external forces times mass

  deque<Kokkos::View<Matrix3d*>::HostMirror> sigma;            ///< Stress matrix
  deque<Kokkos::View<Matrix3d*>::HostMirror> strain_el;        ///< Elastic strain matrix
  deque<Kokkos::View<Matrix3d*>::HostMirror> R;                ///< Rotation matrix

  deque<Kokkos::View<float*>::HostMirror> vol;                       ///< Particles' current volume
  deque<Kokkos::View<float*>::HostMirror> mass;                      ///< Particles' current mass
  deque<Kokkos::View<float*>::HostMirror> eff_plastic_strain;        ///< Particles' effective plastic strain
  deque<Kokkos::View<float*>::HostMirror> eff_plastic_strain_rate;   ///< Particles' effective plastic strain rate
  deque<Kokkos::View<float*>::HostMirror> damage;                    ///< Particles' damage variable
  deque<Kokkos::View<float*>::HostMirror> damage_init;               ///< Particles' damage initiation variable
  deque<Kokkos::View<float*>::HostMirror> ienergy;                   ///< Particles' internal energy

  deque<Kokkos::View<float*>::HostMirror> T;                         ///< Particles' current temperature
  deque<Kokkos::View<float*>::HostMirror> gamma;                     ///< Particles' heat source

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
private:
  deque<pair<thread, vector<float>>>  threads;        ///< Pair storing the threads and the buffer

  void write_to_file(bigint, string, bigint, bigint);
};

#endif
#endif
