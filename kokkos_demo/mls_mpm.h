#pragma once

#include <mat2.h>
#include <graphics.h>

using namespace std;
using namespace Kokkos;

KOKKOS_INLINE_FUNCTION
void svd(const mat2 &A,
         mat2 &U,
         mat2 &sig,
         mat2 &V);

void mls_mpm() {
  constexpr int quality = 1;  // Use a larger value for higher-res simulations
  constexpr int n_particles = 1800*quality*quality;
  constexpr int n_grid = 15*quality;

  constexpr double dx = 1.0/n_grid;
  constexpr double inv_dx = n_grid;
  constexpr double dt = 1e-3/quality;
  constexpr double p_vol = (dx*0.5)*(dx*0.5);
  constexpr double p_rho = 1;
  constexpr double p_mass = p_vol*p_rho;
  constexpr double E = 5e3; // Young's modulus
  constexpr double nu = 0.2; // Poisson's ratio
  constexpr double mu_0 = E/(2*(1 + nu)); // Lame mu
  constexpr double lambda_0 = E*nu/((1 + nu)*(1 - 2*nu)); // Lame lambda

  View<vec2[n_particles]> x("x"); // position
  View<vec2[n_particles]> v("v"); // velocity
  View<mat2[n_particles]> C("C"); // affine velocity field
  View<mat2[n_particles]> F("F"); // deformation gradient
  View<int[n_particles]> material("material"); // material id
  View<double[n_particles]> Jp("Jp"); // plastic deformation

  View<vec2[n_particles]>::HostMirror x_host = create_mirror_view(x); // position mirror
  View<int[n_particles]>::HostMirror material_host = create_mirror_view(material); // material mirror

  View<vec2[n_particles][n_particles]> grid_v("grid_v"); // gird node momentum/velocity
  View<vec2[n_particles][n_particles], MemoryTraits<Atomic>> grid_v_atomic = grid_v; // gird node momentum/velocity
  View<double[n_particles][n_particles]> grid_m("grid_m"); // grid node mass
  View<double[n_particles][n_particles], MemoryTraits<Atomic>> grid_m_atomic = grid_m; // grid node mass
  vec2 gravity{};
  double attractor_strength;
  vec2 attractor_pos{};

  auto substep = [&]() {
    parallel_for(MDRangePolicy<Rank<2>>({ 0, 0 }, { n_particles, n_particles }), KOKKOS_LAMBDA(int i, int j) {
      grid_v(i, j) = vec2();
      grid_m(i, j) = 0;
    });

    // Particle state update and scatter to grid (P2G)
    parallel_for(n_particles, KOKKOS_LAMBDA(int p) {
      vec2 base = x(p)*inv_dx - 0.5;
      base.x = (int)base.x;
      base.y = (int)base.y;
      vec2 fx = x(p)*inv_dx - base;

      vec2 w[3] = { 0.5*(1.5 - fx).square(), 0.75 - (fx - 1).square(), 0.5*(fx -  0.5).square() };

      F(p) = (mat2() + dt*C(p))*F(p); // deformation gradient update
      double h = max(0.1, min(5.0, exp(10*(1 - Jp(p))))); // Hardening coefficient: snow gets harder when compressed
      if (material(p) == 1) // jelly, make it softer
        h = 0.3;
      double mu = mu_0*h, la = lambda_0*h;
      if (material(p) == 0) // liquid
        mu = 0.0;
      mat2 U, sig, V;
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
      if (material(p) == 0) // Reset deformation gradient to avoid numerical instability
        F(p) = mat2()*sqrt(J);
      //else if (material(p) == 2)
      //  F(p) = U*sig*V.transpose(); // Reconstruct elastic deformation gradient after plasticity
      mat2 stress = 2*mu*(F(p) - U*V.transpose())*F(p).transpose() + mat2()*la*J*(J - 1);
      stress = (-dt*p_vol*4*inv_dx*inv_dx)*stress;
      mat2 affine = stress + p_mass * C(p);
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) { // Loop over 3x3 grid node neighborhood
          vec2 offset(i, j);
          vec2 dpos = (offset - fx)*dx;
          double weight = w[i].x*w[j].y;
          vec2 index = base + offset;
          grid_v_atomic((int)index.x, (int)index.y) += weight*(p_mass*v(p) + affine*dpos);
          grid_m_atomic((int)index.x, (int)index.y) += weight*p_mass;
        }
    });
    parallel_for(MDRangePolicy<Rank<2>>({ 0, 0 }, { n_particles, n_particles }), KOKKOS_LAMBDA(int i, int j) {
      if (grid_m(i, j) > 0) { // No need for epsilon here
        grid_v(i, j) = (1/grid_m(i, j))*grid_v(i, j); // Momentum to velocity
        grid_v(i, j) += dt*gravity*30; // gravity
        vec2 dist = attractor_pos - dx*vec2(i, j);
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
    parallel_for(n_particles, KOKKOS_LAMBDA(int p) { // grid to particle(G2P)
      vec2 base = x(p)*inv_dx - 0.5;
      base.x = (int)base.x;
      base.y = (int)base.y;
      vec2 fx = x(p)*inv_dx - base;

      vec2 w[3] = { 0.5*(1.5 - fx).square(), 0.75 - (fx - 1).square(), 0.5*(fx -  0.5).square() };
      vec2 new_v;
      mat2 new_C(0, 0, 0, 0);
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) { // loop over 3x3 grid node neighborhood
          vec2 dpos = vec2(i, j) - fx;
          vec2 offset = base + vec2(i, j);
          vec2 g_v = grid_v((int)offset.x, (int)offset.y);
          double weight = w[i].x * w[j].y;
          new_v += weight*g_v;
          new_C += 4*inv_dx*weight*mat2(g_v.x*dpos.x, g_v.x*dpos.y, g_v.y*dpos.x, g_v.y*dpos.y);
        }

      v(p) = new_v;
      C(p) = new_C;
      x(p) += dt*v(p); // advection
    });
  };

  auto reset = [&]() {
    int group_size = n_particles/3;
    parallel_for(n_particles, KOKKOS_LAMBDA(int i) {
      // terrible way to generate random numbers
      x(i) = vec2((((69471233*i + 787723)*(73491213*i + 7919123))%1046527/1046527.0/2 + 1)*0.2 + 0.3 + 0.1*(int)(i/group_size),
                  (((67331233*i + 729723)*(73331123*i + 7351123))%1676903/1676903.0/2 + 1)*0.2 + 0.05 + 0.32*(int)(i/group_size));

      material[i] = i/group_size; // 0: fluid 1: jelly 2: snow
      v(i) = vec2();
      F(i) = mat2();
      Jp(i) = 1;
      C(i) = mat2(0, 0, 0, 0);
    });
  };

  cout << "[Hint] Use WSAD/arrow keys to control gravity. Use left/right mouse bottons to attract/repel. Press R to reset.";

  reset();
  gravity = vec2(0, -1);

  vector<float> vertices(5*n_particles);

  draw(n_particles, &vertices.front(),
       [&](bool r,
           bool left, bool right, bool up, bool down,
           bool mouse_left, bool mouse_right,
           double xpos, double ypos) {
      if (r)
        reset();

      if (left || right || up || down)
        gravity = vec2();
      if (left)
        gravity.x -= 1;
      if (right)
        gravity.x += 1;
      if (up)
        gravity.y += 1;
      if (down)
        gravity.y -= 1;

      attractor_pos = vec2(xpos, ypos);
      attractor_strength = 0;
      if (mouse_left)
        attractor_strength++;
      if (mouse_right)
        attractor_strength--;

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
}

KOKKOS_INLINE_FUNCTION
void svd(const mat2 &A,
         mat2 &U,
         mat2 &sig,
         mat2 &V) {
  double x = A(0, 0) + A(1, 1), y = A(1, 0) - A(0, 1);
  double scale = 1/sqrt(x*x + y*y);
  double c = x*scale;
  double s = y*scale;
  mat2 R(c, -s, s, c);
  mat2 S = R.transpose()*A;
  double s1 = 0, s2 = 0;
  if (abs(S(0, 1)) < 1e-5) {
    c = 1;
    s = 0;
  }
  else {
   double tao = 0.5*(S(0, 0) - S(1, 1));
   double w = sqrt(tao*tao + S(0, 1)* S(0, 1));
   double t = 0;
   if (tao > 0)
     t = S(0, 1)/(tao + w);
   else
     t = S(0, 1)/(tao - w);
   c = 1/sqrt(t*t + 1);
   s = -t*c;
   s1 = c*c*S(0, 0) - 2*c*s*S(0, 1) + s*s*S(1, 1);
   s2 = s*s*S(0, 0) + 2*c*s*S(0, 1) + c*c*S(1, 1);
  }
       
  if (s1 < s2) {
    double tmp = s1;
    s1 = s2;
    s2 = tmp;
    V = mat2(-s, c, -c, -s);
  }
  else
    V = mat2(c, s, -s, c);
  U = R*V;
  sig = mat2(s1, 0, 0, s2);
}