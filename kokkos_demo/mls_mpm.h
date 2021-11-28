#pragma once

#include <svd.h>
#include <graphics.h>
#include <LazyReference.h>

#include <chrono>
#include <iostream>
#include <limits>

using namespace std;
using namespace Kokkos;
using namespace std::chrono;

using FloatingType = float;
using Vector2 = Vector<FloatingType, 2>;
using Matrix2 = Matrix<FloatingType, 2, 2>;
using LazyVector2 = Vector<LazyReference<FloatingType>, 2>;

void mls_mpm(int quality) {
  const int n_particles = 3*4*4*quality*quality;
  const int n_grid = 9*quality;

  const FloatingType dx = 1.0/n_grid;
  const FloatingType inv_dx = n_grid;
  const FloatingType dt = 1e-4/quality;
  const FloatingType p_vol = (dx*0.5)*(dx*0.5);
  constexpr FloatingType p_rho = 1;
  const FloatingType p_mass = p_vol*p_rho;
  constexpr FloatingType E = 5e3; // Young's modulus
  constexpr FloatingType nu = 0.2; // Poisson's ratio
  constexpr FloatingType mu_0 = E/(2*(1 + nu)); // Lame mu
  constexpr FloatingType lambda_0 = E*nu/((1 + nu)*(1 - 2*nu)); // Lame lambda

  View<FloatingType*> xx("x", n_particles);
  View<FloatingType*> xy("x", n_particles);
  View<FloatingType*> vx("v", n_particles);
  View<FloatingType*> vy("v", n_particles);
  View<FloatingType*> C00("C", n_particles);
  View<FloatingType*> C01("C", n_particles);
  View<FloatingType*> C10("C", n_particles);
  View<FloatingType*> C11("C", n_particles);
  View<FloatingType*> F00("F", n_particles);
  View<FloatingType*> F01("F", n_particles);
  View<FloatingType*> F10("F", n_particles);
  View<FloatingType*> F11("F", n_particles);
  View<int*> material("material", n_particles);
  View<FloatingType*> Jp("Jp", n_particles);

  View<FloatingType**> grid_vx("grid_v", n_grid, n_grid);
  View<FloatingType**> grid_vy("grid_v", n_grid, n_grid);
  View<FloatingType**> grid_m("grid_m", n_grid, n_grid);

  Vector2 gravity(0, -1);

  auto reset = [&]() {
    int group_size = n_particles/3;
    parallel_for(n_particles, KOKKOS_LAMBDA(int i) {
      int i_in_group = i%group_size;
      int dim_size = 4*quality;
      FloatingType row = i_in_group%dim_size;
      FloatingType col = i_in_group/dim_size;

      xx(i) = row/dim_size*0.2 + 0.3  + 0.1 *(i/group_size);
      xy(i) = col/dim_size*0.2 + 0.05 + 0.29*(i/group_size);
      vx(i) = 0;
      vy(i) = 0;
      C00(i) = 0;
      C01(i) = 0;
      C10(i) = 0;
      C11(i) = 0;
      F00(i) = 1;
      F01(i) = 0;
      F10(i) = 0;
      F11(i) = 1;
      material(i) = 1;
      Jp(i) = 1;
    });
  };

  auto substep = [&]() {
    // reset node
    parallel_for(MDRangePolicy<Rank<2>>({ 0, 0 }, { n_grid, n_grid }), KOKKOS_LAMBDA(int i, int j) {
      grid_vx(i, j) = 0;
      grid_vy(i, j) = 0;
      grid_m(i, j) = 0;
    });

    // P2G
    parallel_for(n_particles, KOKKOS_LAMBDA(int i) {
      Vector2 x(xx(i), xy(i));
      Vector2 v(vx(i), vy(i));
      Matrix2 C(C00(i), C01(i), C10(i), C11(i));
      Matrix2 F(F00(i), F01(i), F10(i), F11(i));
      int mat = material(i);
      FloatingType jp = Jp(i);

      Vector<int, 2> base = x*inv_dx - Vector2(0.5, 0.5);
      Vector2 fx = x*inv_dx - base;

      Vector2 w[]{ 0.5*(Vector2(1.5, 1.5) - fx).square(),
        Vector2(0.75, 0.75) - (fx - Vector2(1, 1)).square(), 
        0.5*(fx - Vector2(0.5, 0.5)).square() };

      F = (Matrix2(1, 0, 0, 1) + dt*C)*F; // deformation gradient update
      FloatingType h = max((FloatingType)0.1, min((FloatingType)5.0, exp(10*(1 - jp)))); // Hardening coefficient: snow gets harder when compressed
      if (mat == 1) // jelly, make it softer
        h = 0.3;
      FloatingType mu = mu_0*h, la = lambda_0*h;
      if (mat == 0) // liquid
        mu = 0.0;
      Matrix2 U, sig, V;
      svd(F, U, sig, V);
      FloatingType J = 1;
      for (int d = 0; d < 2; d++) {
        FloatingType new_sig = sig(d, d);
        if (mat == 2) { // Snow
          new_sig = min(max(sig(d, d), (FloatingType)(1 - 2.5e-2)), (FloatingType)(1 + 4.5e-3)); // Plasticity
        }
        jp *= sig(d, d)/new_sig;
        sig(d, d) = new_sig;
        J *= new_sig;
      }
      if (mat == 0) // Reset deformation gradient to avoid numerical instability
        F = Matrix2(1, 0, 0, 1)*sqrt(J);
      else if (mat == 2)
        F = U*sig*V.transpose(); // Reconstruct elastic deformation gradient after plasticity
      Matrix2 stress = 2*mu*(F - U*V.transpose())*F.transpose() + Matrix2(1, 0, 0, 1)*la*J*(J - 1);
      stress = (-dt*p_vol*4*inv_dx*inv_dx)*stress;
      Matrix2 affine = stress + p_mass*C;
      Vector2 dpos;
      Vector<int, 2> offset, index;
      FloatingType weight;
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) { // Loop over 3x3 node node neighborhood
          offset.x() = i;
          offset.y() = j;
                
          weight = w[i].x()*w[j].y();
          index = base + offset;

          atomic_add(&grid_m((int)index.x(), (int)index.y()), weight*p_mass);

          dpos = (offset - fx)*dx;
          const Vector2 &dv = weight*(p_mass*v + affine*dpos);
          atomic_add(&grid_vx(index.x(), index.y()), dv.x());
          atomic_add(&grid_vy(index.x(), index.y()), dv.y());
        }

      F00(i) = F(0, 0);
      F01(i) = F(0, 1);
      F10(i) = F(1, 0);
      F11(i) = F(1, 1);
      Jp(i) = jp;
    });

    // node update
    parallel_for(MDRangePolicy<Rank<2>>({ 0, 0 }, { n_grid, n_grid }), KOKKOS_LAMBDA(int i, int j) {
      Vector2 v(grid_vx(i, j), grid_vy(i, j));
      FloatingType m = grid_m(i, j);

      if (m > 0) { // No need for epsilon here
        v /= m; // Momentum to velocity
        v += dt*gravity*30; // gravity
        if (i < 3 && v.x() < 0)
          v.x() = 0; // Boundary conditions
        if (i > n_grid - 3 && v.x() > 0)
          v.x() = 0;
        if (j < 3 && v.y() < 0)
          v.y() = 0;
        if (j > n_grid - 3 && v.y() > 0)
          v.y() = 0;
      }
      grid_vx(i, j) = v.x();
      grid_vy(i, j) = v.y();
    });

    // G2P
    parallel_for(n_particles, KOKKOS_LAMBDA(int i) {
      Vector2 x(xx(i), xy(i));
      Vector2 v(vx(i), vy(i));
      Matrix2 C(C00(i), C01(i), C10(i), C11(i));

      Vector<int, 2> base = x*inv_dx - Vector2(0.5, 0.5);
      Vector2 fx = x*inv_dx - base;

      Vector2 w[]{ 0.5*(Vector2(1.5, 1.5) - fx).square(),
        Vector2(0.75, 0.75) - (fx - Vector2(1, 1)).square(), 
        0.5*(fx - Vector2(0.5, 0.5)).square() };

      v.x() = 0; v.y() = 0;
      C(0, 0) = 0; C(0, 1) = 0; C(1, 0) = 0; C(1, 1) = 0;
      
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) { // loop over 3x3 node node neighborhood
          Vector2 dpos = Vector2(i, j) - fx;
          Vector2 offset = base + Vector2(i, j);
          Vector2 g_v;
          g_v.x() = grid_vx((int)offset.x(), (int)offset.y());
          g_v.y() = grid_vy((int)offset.x(), (int)offset.y());
          FloatingType weight = w[i].x()*w[j].y();
          v += weight*g_v;
          C += 4*inv_dx*weight*Matrix2(g_v.x()*dpos.x(), g_v.x()*dpos.y(), g_v.y()*dpos.x(), g_v.y()*dpos.y());
        }

      x += dt*v; // advection
      xx(i) = x.x();
      xy(i) = x.y();
      vx(i) = v.x();
      vy(i) = v.y();
      C00(i) = C(0, 0);
      C01(i) = C(0, 1);
      C10(i) = C(1, 0);
      C11(i) = C(1, 1);
    });
  };

#if 1
  time_point<steady_clock> start_time = steady_clock::now();

  reset();

  for (int i = 0; i < 10000; i++) {
    substep();
  }

  fence();

  cout //<< n_particles << " " << n_grid << " "
       << duration_cast<milliseconds>(steady_clock::now() - start_time).count() << endl;
#else
  View<FloatingType*>::HostMirror xx_host = create_mirror_view(xx);
  View<FloatingType*>::HostMirror xy_host = create_mirror_view(xy);
  View<int*>::HostMirror material_host = create_mirror_view(material);

  vector<float> vertices(5*n_particles);

  reset();

  draw(n_particles, &vertices.front(),
       [&](bool r, bool left, bool right, bool up, bool down,
           bool mouse_left, bool mouse_right,
           FloatingType xpos, FloatingType ypos) {

      if (r)
        reset();

      if (left || right || up || down)
        gravity = Vector2();
      if (left)
        gravity.x() -= 1;
      if (right)
        gravity.x() += 1;
      if (up)
        gravity.y() += 1;
      if (down)
        gravity.y() -= 1;

      for (int s = 0; s < max(1, (int)(2e-3/dt)); s++) {
        substep();
      }

      deep_copy(xx_host, xx);
      deep_copy(xy_host, xy);
      deep_copy(material_host, material);

      for (int i = 0; i < n_particles; i++) {
        vertices.at(5*i    ) = 2*xx_host(i) - 1;
        vertices.at(5*i + 1) = 2*xy_host(i) - 1;
        vertices.at(5*i + 2) = material_host(i) == 0;
        vertices.at(5*i + 3) = material_host(i) == 1;
        vertices.at(5*i + 4) = material_host(i) == 2;
      }
    });
#endif
}
