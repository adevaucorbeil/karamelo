#pragma once

#include <Mat2.h>
#include <svd.h>
#include <graphics.h>

#include <chrono>
#include <iostream>
#include <limits>

using namespace std;
using namespace Kokkos;
using namespace std::chrono;

constexpr double sqrt1(double x, double curr, double prev)
{
    return curr == prev? curr: sqrt1(x, 0.5*(curr + x/curr), curr);
}

void time_step(time_point<steady_clock> &time_now, const string &message) {
  fence();
  time_point<steady_clock> new_time_now = steady_clock::now();
  if (!message.empty())
  cout << message << ": " << duration_cast<milliseconds>(new_time_now - time_now).count() << "ms" << endl;
  time_now = new_time_now;
}

template<int n_particles>
void mls_mpm() {
  time_point<steady_clock> time_now = steady_clock::now(); 

  constexpr int quality = 1;  // Use a larger value for higher-res simulations
  // constexpr int n_particles = 9000*quality*quality;
  constexpr int n_grid = sqrt1(n_particles, n_particles, 0)*3/2;

  constexpr double dx = 1.0/n_grid;
  constexpr double inv_dx = n_grid;
  constexpr double dt = 1e-4/quality;
  constexpr double p_vol = (dx*0.5)*(dx*0.5);
  constexpr double p_rho = 1;
  constexpr double p_mass = p_vol*p_rho;
  constexpr double E = 5e3; // Young's modulus
  constexpr double nu = 0.2; // Poisson's ratio
  constexpr double mu_0 = E/(2*(1 + nu)); // Lame mu
  constexpr double lambda_0 = E*nu/((1 + nu)*(1 - 2*nu)); // Lame lambda

  View<Vec2[n_particles]> x("x"); // position
  View<Vec2[n_particles]> v("v"); // velocity
  View<Mat2[n_particles]> C("C"); // affine velocity field
  View<Mat2[n_particles]> F("F"); // deformation gradient
  View<int[n_particles]> material("material"); // material id
  View<double[n_particles]> Jp("Jp"); // plastic deformation

  constexpr int np = n_particles;
  typename View<Vec2[np]>::HostMirror x_host = create_mirror_view(x); // position mirror
  typename View<int[np]>::HostMirror material_host = create_mirror_view(material); // material mirror

  View<Vec2[n_grid][n_grid]> grid_v("grid_v"); // gird node momentum/velocity
  View<Vec2[n_grid][n_grid], MemoryTraits<Atomic>> grid_v_atomic = grid_v; // gird node momentum/velocity
  View<double[n_grid][n_grid]> grid_m("grid_m"); // grid node mass
  View<double[n_grid][n_grid], MemoryTraits<Atomic>> grid_m_atomic = grid_m; // grid node mass
  Vec2 gravity{};
  double attractor_strength = 0;
  Vec2 attractor_pos{};

  auto substep = [&]() {
    fence();
    Profiling::pushRegion("Grid reset");
    parallel_for(MDRangePolicy<Rank<2>>({ 0, 0 }, { n_grid, n_grid }), KOKKOS_LAMBDA(int i, int j) {
      grid_v(i, j) = Vec2();
      grid_m(i, j) = 0;
    });
    fence();
    Profiling::popRegion();

    //fence();
     //time_step(time_now, /*"Reset grid"*/"");

    // Particle state update and scatter to grid (P2G)
    fence();
    Profiling::pushRegion("P2G");
    parallel_for(n_particles, KOKKOS_LAMBDA(int p) {
      Vec2 base = x(p)*inv_dx - 0.5;
      base.x = (int)base.x;
      base.y = (int)base.y;
      Vec2 fx = x(p)*inv_dx - base;

      Vec2 w[3] = { 0.5*(1.5 - fx).square(), 0.75 - (fx - 1).square(), 0.5*(fx -  0.5).square() };

      F(p) = (Mat2() + dt*C(p))*F(p); // deformation gradient update
      double h = max(0.1, min(5.0, exp(10*(1 - Jp(p))))); // Hardening coefficient: snow gets harder when compressed
      if (material(p) == 1) // jelly, make it softer
        h = 0.3;
      double mu = mu_0*h, la = lambda_0*h;
//      if (material(p) == 0) // liquid
        //mu = 0.0;
      Mat2 U, sig, V;
      svd(F(p), U, sig, V);
      double J = 1;
      //for (int d = 0; d < 2; d++) {
      //  double new_sig = sig(d, d);
      //  if (material(p) == 2) { // Snow
      //    new_sig = min(max(sig(d, d), 1 - 2.5e-2), 1 + 4.5e-3); // Plasticity
      //    Jp(p) *= sig(d, d)/new_sig;
      //    sig(d, d) = new_sig;
      //    J *= new_sig;
      //  }
      //}
//      if (material(p) == 0) // Reset deformation gradient to avoid numerical instability
//             F(p) = Mat2()*sqrt(J);
      //else if (material(p) == 2)
      //  F(p) = U*sig*V.transpose(); // Reconstruct elastic deformation gradient after plasticity
      Mat2 stress = 2*mu*(F(p) - U*V.transpose())*F(p).transpose() + Mat2()*la*J*(J - 1);
      stress = (-dt*p_vol*4*inv_dx*inv_dx)*stress;
      Mat2 affine = stress + p_mass * C(p);
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) { // Loop over 3x3 grid node neighborhood
          Vec2 offset(i, j);
          Vec2 dpos = (offset - fx)*dx;
          double weight = w[i].x*w[j].y;
          Vec2 index = base + offset;
          grid_v_atomic((int)index.x, (int)index.y) += weight*(p_mass*v(p) + affine*dpos);
          grid_m_atomic((int)index.x, (int)index.y) += weight*p_mass;
        }
    });
    fence();
    Profiling::popRegion();

    //fence();
     //time_step(time_now, "P2G");
    fence();
    Profiling::pushRegion("Grid update");
    parallel_for(MDRangePolicy<Rank<2>>({ 0, 0 }, { n_grid, n_grid }), KOKKOS_LAMBDA(int i, int j) {
      if (grid_m(i, j) > 0) { // No need for epsilon here
        grid_v(i, j) = (1/grid_m(i, j))*grid_v(i, j); // Momentum to velocity
        grid_v(i, j) += dt*gravity*30; // gravity
        Vec2 dist = attractor_pos - dx*Vec2(i, j);
        grid_v(i, j) += dist/(0.01 + dist.norm()) * attractor_strength*dt*100;
        if (i < 3 && grid_v(i, j).x < 0)
          grid_v(i, j).x = 0; // Boundary conditions
        if (i > n_grid - 3 && grid_v(i, j).x > 0)
          grid_v(i, j).x = 0;
        if (j < 3 && grid_v(i, j).y < 0)
          grid_v(i, j).y = 0;
        if (j > n_grid - 3 && grid_v(i, j).y > 0)
          grid_v(i, j).y = 0;
      }
    });
    fence();
    Profiling::popRegion();

  //   time_step(time_now, "Grid update");
  // time_now = steady_clock::now();
    fence();
    Profiling::pushRegion("G2P");
    parallel_for(n_particles, KOKKOS_LAMBDA(int p) { // grid to particle(G2P)
      Vec2 base = x(p)*inv_dx - 0.5;
      base.x = (int)base.x;
      base.y = (int)base.y;
      Vec2 fx = x(p)*inv_dx - base;

      Vec2 w[3] = { 0.5*(1.5 - fx).square(), 0.75 - (fx - 1).square(), 0.5*(fx -  0.5).square() };
      v(p).x = 0; v(p).y = 0;
      C(p).m00 = 0; C(p).m01 = 0; C(p).m10 = 0; C(p).m11 = 0;
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) { // loop over 3x3 grid node neighborhood
          Vec2 dpos = Vec2(i, j) - fx;
          Vec2 offset = base + Vec2(i, j);
          Vec2 g_v = grid_v((int)offset.x, (int)offset.y);
          double weight = w[i].x*w[j].y;
          v(p) += weight*g_v;
          C(p) += 4*inv_dx*weight*Mat2(g_v.x*dpos.x, g_v.x*dpos.y, g_v.y*dpos.x, g_v.y*dpos.y);
        }

      x(p) += dt*v(p); // advection
    });
    fence();
    Profiling::popRegion();

    //time_step(time_now, "G2P");
  };

  auto reset = [&]() {
    int group_size = n_particles/3;
    parallel_for(n_particles, KOKKOS_LAMBDA(int i) {
      // terrible way to generate random numbers
      x(i) = Vec2((((69471233*i + 787723)*(73491213*i + 7919123))%1046527/1046527.0/2 + 1)*0.2 + 0.3 + 0.1*(int)(i/group_size),
                  (((67331233*i + 729723)*(73331123*i + 7351123))%1676903/1676903.0/2 + 1)*0.2 + 0.05 + 0.29*(int)(i/group_size));

      material[i] = i/group_size; // 0: fluid 1: jelly 2: snow
  //    material[i] = 1;
      v(i) = Vec2();
      F(i) = Mat2();
      Jp(i) = 1;
      C(i) = Mat2(0, 0, 0, 0);
    });
  };

  // cout << "[Hint] Use WSAD/arrow keys to control gravity. Use left/right mouse bottons to attract/repel. Press R to reset.";

  reset();
  gravity = Vec2(0, -1);

  vector<float> vertices(5*n_particles);

#if 0
  draw(n_particles, &vertices.front(),
       [&](bool r,
           bool left, bool right, bool up, bool down,
           bool mouse_left, bool mouse_right,
           double xpos, double ypos) {
//      time_step(time_now, "Draw");

      if (r)
        reset();

      if (left || right || up || down)
        gravity = Vec2();
      if (left)
        gravity.x -= 1;
      if (right)
        gravity.x += 1;
      if (up)
        gravity.y += 1;
      if (down)
        gravity.y -= 1;

      attractor_pos = Vec2(xpos, ypos);
      attractor_strength = 0;
      if (mouse_left)
        attractor_strength++;
      if (mouse_right)
        attractor_strength--;

//      time_step(time_now, "Process input");

      for (int s = 0; s < max(1, (int)(2e-3/dt)); s++) {
        substep();
      }

      deep_copy(Cuda(), x_host, x);
      deep_copy(Cuda(), material_host, material);

//      time_step(time_now, "Deep copy");

      for (int i = 0; i < n_particles; i++) {
        vertices.at(5*i    ) = 2*x_host(i).x - 1;
        vertices.at(5*i + 1) = 2*x_host(i).y - 1;
        vertices.at(5*i + 2) = material_host(i) == 0;
        vertices.at(5*i + 3) = material_host(i) == 1;
        vertices.at(5*i + 4) = material_host(i) == 2;
      }

//      time_step(time_now, "Vertex update");
    });
#else
  for (int i = 0; i < 5000; i++) {
    substep();
    //cout << i << endl;
    //time_step(time_now, "step");
  }
  cout << n_particles << ", " << n_grid;
  time_step(time_now, ", ");
#endif
}