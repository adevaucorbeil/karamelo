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
#include <grid.h>
#include <thread>
#include <unordered_map>

class DumpParticleGz : public Dump {
  deque<Kokkos::View<tagint*>::HostMirror> ptag;               ///< Unique identifier for particles in the system

  deque<Kokkos::View<Vector3d*>::HostMirror> x;                ///< Particles' current position
  deque<Kokkos::View<Vector3d*>::HostMirror> x0;               ///< Particles' reference position
  
  deque<Kokkos::View<Vector3d*>::HostMirror> v;                ///< Particles' current velocity

  deque<Kokkos::View<Vector3d*>::HostMirror> mbp;              ///< Particles' external forces times mass

  deque<Kokkos::View<Matrix3d*>::HostMirror> sigma;            ///< Stress matrix
  deque<Kokkos::View<Matrix3d*>::HostMirror> strain_el;        ///< Elastic strain matrix
  deque<Kokkos::View<Matrix3d*>::HostMirror> R;                ///< Rotation matrix

  deque<Kokkos::View<double*>::HostMirror> vol;                       ///< Particles' current volume
  deque<Kokkos::View<double*>::HostMirror> mass;                      ///< Particles' current mass
  deque<Kokkos::View<double*>::HostMirror> eff_plastic_strain;        ///< Particles' effective plastic strain
  deque<Kokkos::View<double*>::HostMirror> eff_plastic_strain_rate;   ///< Particles' effective plastic strain rate
  deque<Kokkos::View<double*>::HostMirror> damage;                    ///< Particles' damage variable
  deque<Kokkos::View<double*>::HostMirror> damage_init;               ///< Particles' damage initiation variable
  deque<Kokkos::View<double*>::HostMirror> ienergy;                   ///< Particles' internal energy

  deque<Kokkos::View<double*>::HostMirror> T;                         ///< Particles' current temperature
  deque<Kokkos::View<double*>::HostMirror> gamma;                     ///< Particles' heat source

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
  unordered_map<bigint, thread>        threads;        ///< Map storing the threads for each write_to_file
  unordered_map<bigint, bool>          thread_status;  ///< true when thread is running
  unordered_map<bigint, deque<double>> bufs;           ///< Buffers containing the values to be written in dump

  void write_to_file(bigint, string, bigint);
};

#endif
#endif
