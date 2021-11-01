#pragma once

#include <svd.h>
#include <graphics.h>

#include <chrono>
#include <iostream>
#include <limits>

using namespace std;
using namespace Kokkos;
using namespace std::chrono;

void mls_mpm(int quality) {
  int n_particles = 3*4*4*quality*quality;
  int n_grid = 9*quality;

  double dx = 1.0/n_grid;
  double inv_dx = n_grid;
  double dt = 1e-4/quality;
  double p_vol = (dx*0.5)*(dx*0.5);
  constexpr double p_rho = 1;
  double p_mass = p_vol*p_rho;
  constexpr double E = 5e3; // Young's modulus
  constexpr double nu = 0.2; // Poisson's ratio
  constexpr double mu_0 = E/(2*(1 + nu)); // Lame mu
  constexpr double lambda_0 = E*nu/((1 + nu)*(1 - 2*nu)); // Lame lambda

  View<Vec2*> x("x", n_particles); // position
  View<Vec2*> v("v", n_particles); // velocity
  View<Mat2*> C("C", n_particles); // affine velocity field
  View<Mat2*> F("F", n_particles); // deformation gradient
  View<int*> material("material", n_particles); // material id
  View<double*> Jp("Jp", n_particles); // plastic deformation

  View<Vec2**> grid_v("grid_v", n_grid, n_grid); // gird node momentum/velocity
  View<double**> grid_m("grid_m", n_grid, n_grid); // grid node mass

  Vec2 gravity(0, -1);

  auto reset = [&]() {
    int group_size = n_particles/3;
    parallel_for(n_particles, KOKKOS_LAMBDA(int i) {
      int i_in_group = i%group_size;
      int dim_size = 4*quality;
      double row = i_in_group%dim_size;
      double col = i_in_group/dim_size;
      x(i) = Vec2(row/dim_size*0.2 + 0.3 + 0.1*(int)(i/group_size),
                  col/dim_size*0.2 + 0.05 + 0.29*(int)(i/group_size));
      material[i] = 1;//i/group_size; // 0: fluid 1: jelly 2: snow
      v(i) = Vec2();
      F(i) = Mat2();
      Jp(i) = 1;
      C(i) = Mat2(0, 0, 0, 0);
    });
  };

  auto substep = [&]() {
    // reset grid
    parallel_for(MDRangePolicy<Rank<2>>({ 0, 0 }, { n_grid, n_grid }), KOKKOS_LAMBDA(int i, int j) {
      grid_v(i, j) = Vec2();
      grid_m(i, j) = 0;
    });

    // P2G
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
      if (material(p) == 0) // liquid
        mu = 0.0;
      Mat2 U, sig, V;
      svd(F(p), U, sig, V);
      double J = 1;
      for (int d = 0; d < 2; d++) {
        double new_sig = sig(d, d);
        if (material(p) == 2) { // Snow
          new_sig = min(max(sig(d, d), 1 - 2.5e-2), 1 + 4.5e-3); // Plasticity
        }
        Jp(p) *= sig(d, d)/new_sig;
        sig(d, d) = new_sig;
        J *= new_sig;
      }
      if (material(p) == 0) // Reset deformation gradient to avoid numerical instability
        F(p) = Mat2()*sqrt(J);
      else if (material(p) == 2)
        F(p) = U*sig*V.transpose(); // Reconstruct elastic deformation gradient after plasticity
      Mat2 stress = 2*mu*(F(p) - U*V.transpose())*F(p).transpose() + Mat2()*la*J*(J - 1);
      stress = (-dt*p_vol*4*inv_dx*inv_dx)*stress;
      Mat2 affine = stress + p_mass * C(p);
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) { // Loop over 3x3 grid node neighborhood
          Vec2 offset(i, j);
          Vec2 dpos = (offset - fx)*dx;
          double weight = w[i].x*w[j].y;
          Vec2 index = base + offset;
          const Vec2 &dv = weight*(p_mass*v(p) + affine*dpos);
          atomic_add(&grid_v((int)index.x, (int)index.y).x, dv.x);
          atomic_add(&grid_v((int)index.x, (int)index.y).y, dv.y);
          atomic_add(&grid_m((int)index.x, (int)index.y), weight*p_mass);
        }
    });

    // grid update
    parallel_for(MDRangePolicy<Rank<2>>({ 0, 0 }, { n_grid, n_grid }), KOKKOS_LAMBDA(int i, int j) {
      if (grid_m(i, j) > 0) { // No need for epsilon here
        grid_v(i, j) = (1/grid_m(i, j))*grid_v(i, j); // Momentum to velocity
        grid_v(i, j) += dt*gravity*30; // gravity
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

    // G2P
    parallel_for(n_particles, KOKKOS_LAMBDA(int p) {
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
  };

  // cout << n_particles << ", " << n_grid << ": ";
  // cout.flush();

#if 1

  time_point<steady_clock> start_time = steady_clock::now();

  reset();

  for (int i = 0; i < 5000; i++) {
    substep();
  }

  fence();

  cout << duration_cast<milliseconds>(steady_clock::now() - start_time).count() << endl;
#else
  View<Vec2*>::HostMirror x_host = create_mirror_view(x); // position mirror
  View<int*>::HostMirror material_host = create_mirror_view(material); // material mirror

  vector<float> vertices(5*n_particles);

  reset();

  draw(n_particles, &vertices.front(),
       [&](bool r, bool left, bool right, bool up, bool down,
           bool mouse_left, bool mouse_right,
           double xpos, double ypos) {

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

      for (int s = 0; s < max(1, (int)(2e-3/dt)); s++) {
        substep();
      }

      deep_copy(x_host, x);
      deep_copy(material_host, material);

      for (int i = 0; i < n_particles; i++) {
        vertices.at(5*i    ) = 2*x_host(i).x - 1;
        vertices.at(5*i + 1) = 2*x_host(i).y - 1;
        vertices.at(5*i + 2) = material_host(i) == 0;
        vertices.at(5*i + 3) = material_host(i) == 1;
        vertices.at(5*i + 4) = material_host(i) == 2;
      }
    });
#endif
}
