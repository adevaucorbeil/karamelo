#pragma once

#include <svd.h>
#include <nvToolsExt.h>
// #include <graphics.h>

#include <vector>
#include <chrono>
#include <iostream>
#include <limits>

using namespace std;
using namespace std::chrono;

struct Particle {
  Vec2 x;
  Vec2 v;
  Mat2 C;
  Mat2 F;
  int material;
  double Jp;
};

struct Node {
  Vec2 v;
  double m;
};

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

  vector<Particle> particles(n_particles);

  vector<vector<Node>> nodes(n_grid, vector<Node>(n_grid));

  Vec2 gravity(0, -1);

  auto reset = [&]() {

  };

  auto substep = [&]() {

    // nvtxRangePop();
  };

#if 1

  time_point<steady_clock> start_time = steady_clock::now();

  int group_size = n_particles/3;
  #pragma acc parallel loop
  for (int i = 0; i < n_particles; i++) {
    Particle particle = particles.at(i);

    int i_in_group = i%group_size;
    int dim_size = 4*quality;
    double row = i_in_group%dim_size;
    double col = i_in_group/dim_size;
    particle.x = Vec2(row/dim_size*0.2 + 0.3 + 0.1*(int)(i/group_size),
                      col/dim_size*0.2 + 0.05 + 0.29*(int)(i/group_size));
    particle.material = 1;//i/group_size; // 0: fluid 1: jelly 2: snow
    particle.v = Vec2();
    particle.F = Mat2();
    particle.Jp = 1;
    particle.C = Mat2(0, 0, 0, 0);

    particles.at(i) = particle;
  }

  for (int i = 0; i < 5000; i++) {
    nvtxRangePush("Reset_node");
    #pragma acc kernels copyout(nodes)
    {
      for (int i = 0; i < n_grid; i++) {
        for (int j = 0; j < n_grid; j++) {
            Node node = nodes.at(i).at(j);

            node.v = Vec2();
            node.m = 0;

            nodes.at(i).at(j) = node;
        }
      }
    }
    nvtxRangePop();

    // P2G
    nvtxRangePush("P2G");
    //#pragma acc parallel loop
    for (int p = 0; p < n_particles; p++) {
      Particle particle = particles.at(p);

      Vec2 base = particle.x*inv_dx - 0.5;
      base.x = (int)base.x;
      base.y = (int)base.y;
      Vec2 fx = particle.x*inv_dx - base;

      Vec2 w[3] = { 0.5*(1.5 - fx).square(), 0.75 - (fx - 1).square(), 0.5*(fx -  0.5).square() };

      particle.F = (Mat2() + dt*particle.C)*particle.F; // deformation gradient update
      double h = max(0.1, min(5.0, exp(10*(1 - particle.Jp)))); // Hardening coefficient: snow gets harder when compressed
      if (particle.material == 1) // jelly, make it softer
        h = 0.3;
      double mu = mu_0*h, la = lambda_0*h;
      if (particle.material == 0) // liquid
        mu = 0.0;
      Mat2 U, sig, V;
      svd(particle.F, U, sig, V);
      double J = 1;
      for (int d = 0; d < 2; d++) {
        double new_sig = sig(d, d);
        if (particle.material == 2) { // Snow
          new_sig = min(max(sig(d, d), 1 - 2.5e-2), 1 + 4.5e-3); // Plasticity
        }
        particle.Jp *= sig(d, d)/new_sig;
        sig(d, d) = new_sig;
        J *= new_sig;
      }
      if (particle.material == 0) // Reset deformation gradient to avoid numerical instability
        particle.F = Mat2()*sqrt(J);
      else if (particle.material == 2)
        particle.F = U*sig*V.transpose(); // Reconstruct elastic deformation gradient after plasticity
      Mat2 stress = 2*mu*(particle.F - U*V.transpose())*particle.F.transpose() + Mat2()*la*J*(J - 1);
      stress = (-dt*p_vol*4*inv_dx*inv_dx)*stress;
      Mat2 affine = stress + p_mass * particle.C;
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) { // Loop over 3x3 node node neighborhood
          Vec2 offset(i, j);
          Vec2 dpos = (offset - fx)*dx;
          double weight = w[i].x*w[j].y;
          Vec2 index = base + offset;
          const Vec2 &dv = weight*(p_mass*particle.v + affine*dpos);
          #pragma acc atomic update
          nodes.at((int)index.x).at((int)index.y).v.x += dv.x;
          #pragma acc atomic update
          nodes.at((int)index.x).at((int)index.y).v.y += dv.y;
          #pragma acc atomic update
          nodes.at((int)index.x).at((int)index.y).m += weight*p_mass;
        }

      particles.at(p) = particle;
    }
    nvtxRangePop();

    // node update
    nvtxRangePush("Grid_update");
    //#pragma acc parallel loop
    for (int i = 0; i < n_grid; i++) {
      for (int j = 0; j < n_grid; j++) {
        Node node = nodes.at(i).at(j);

        if (node.m > 0) { // No need for epsilon here
          node.v = (1/node.m)*node.v; // Momentum to velocity
          node.v += dt*gravity*30; // gravity
          if (i < 3 && node.v.x < 0)
            node.v.x = 0; // Boundary conditions
          if (i > n_grid - 3 && node.v.x > 0)
            node.v.x = 0;
          if (j < 3 && node.v.y < 0)
            node.v.y = 0;
          if (j > n_grid - 3 && node.v.y > 0)
            node.v.y = 0;
        }

        nodes.at(i).at(j) = node;
      }
    }
    nvtxRangePop();

    // G2P
    nvtxRangePush("G2P");
    //#pragma acc parallel loop
    for (int p = 0; p < n_particles; p++) {
      Particle particle = particles.at(p);

      Vec2 base = particle.x*inv_dx - 0.5;
      base.x = (int)base.x;
      base.y = (int)base.y;
      Vec2 fx = particle.x*inv_dx - base;

      Vec2 w[3] = { 0.5*(1.5 - fx).square(), 0.75 - (fx - 1).square(), 0.5*(fx -  0.5).square() };
      particle.v.x = 0; particle.v.y = 0;
      particle.C.m00 = 0; particle.C.m01 = 0; particle.C.m10 = 0; particle.C.m11 = 0;
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) { // loop over 3x3 node node neighborhood
          Vec2 dpos = Vec2(i, j) - fx;
          Vec2 offset = base + Vec2(i, j);
          Vec2 g_v = nodes.at((int)offset.x).at((int)offset.y).v;
          double weight = w[i].x*w[j].y;
          particle.v += weight*g_v;
          particle.C += 4*inv_dx*weight*Mat2(g_v.x*dpos.x, g_v.x*dpos.y, g_v.y*dpos.x, g_v.y*dpos.y);
        }

      particle.x += dt*particle.v; // advection

      particles.at(p) = particle;
    }
    nvtxRangePop();
  }

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
