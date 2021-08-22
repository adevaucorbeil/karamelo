#pragma once

#include <mat2.h>

#include <vector>
#include <cmath>

using namespace std;

void svd(const mat2 &m,
         mat2 &U,
         mat2 &sig,
         mat2 &V);

void mls_mpm() {
  int quality = 1;  // Use a larger value for higher-res simulations
  int n_particles = 9000*quality*quality;
  int n_grid = 128*quality;

  double dx = 1/n_grid;
  double inv_dx = n_grid;
  double dt = 1e-4 / quality;
  double p_vol = (dx*0.5)*(dx*0.5);
  double p_rho = 1;
  double p_mass = p_vol*p_rho;
  double E = 5e3; // Young's modulus
  double nu = 0.2; // Poisson's ratio
  double mu_0 = E/(2*(1 + nu)); // Lame mu
  double lambda_0 = E*nu/((1 + nu)*(1 - 2*nu)); // Lame lambda

  vector<vec2> x(n_particles); // position
  vector<vec2> v(n_particles); // velocity
  vector<mat2> C(n_particles); // affine velocity field
  vector<mat2> F(n_particles); // deformation gradient
  vector<int> material(n_particles); // material id
  vector<double> Jp(n_particles); // plastic deformation

  vector<vector<vec2>> grid_v(n_particles, vector<vec2>(n_particles)); // gird node momentum/velocity
  vector<vector<double>> grid_m(n_particles, vector<double>(n_particles)); // grid node mass
  vec2 gravity{};
  double attractor_strength;
  vec2 attractor_pos{};

  auto substep = [&]() {
    for (int i = 0; i < n_particles; i++)
      for (int j = 0; j < n_particles; j++) {
        grid_v.at(i).at(j) = vec2();
        grid_m.at(i).at(j) = 0;
      }

    // Particle state update and scatter to grid (P2G)
    for (int p = 0; p < n_particles; p++) {
      vec2 base = x.at(p)*inv_dx - 0.5;
      base.x = (int)base.x;
      base.y = (int)base.y;
      vec2 fx = x.at(p)*inv_dx - base;

      vec2 w[3] = { 0.5*(1.5 - fx).square(), 0.75 - (fx - 1).square(), 0.5*(fx -  0.5).square() };

      F.at(p) = (mat2() + dt*C.at(p))*F.at(p); // deformation gradient update
      double h = max(0.1, min(5.0, exp(10*(1 - Jp.at(p))))); // Hardening coefficient: snow gets harder when compressed
      if (material.at(p) == 1) // jelly, make it softer
        h = 0.3;
      double mu = mu_0*h, la = lambda_0*h;
      if (material.at(p) == 0) // liquid
        mu = 0.0;
      mat2 U, sig, V;
      svd(F.at(p), U, sig, V);
      double J = 1;
      for (int d = 0; d < 2; d++) {
        double new_sig = sig(d, d);
        if (material.at(p) == 2) { // Snow
          new_sig = min(max(sig(d, d), 1 - 2.5e-2), 1 + 4.5e-3); // Plasticity
          Jp.at(p) *= sig(d, d)/new_sig;
          sig(d, d) = new_sig;
          J *= new_sig;
        }
      }
      if (material.at(p) == 0) // Reset deformation gradient to avoid numerical instability
        F.at(p) = mat2()*sqrt(J);
      else if (material.at(p) == 2)
        F.at(p) = U*sig*V.transpose(); // Reconstruct elastic deformation gradient after plasticity
      mat2 stress = 2*mu*(F.at(p) - U*V.transpose())*F.at(p).transpose() + mat2()*la*J*(J - 1);
      stress = (-dt*p_vol*4*inv_dx*inv_dx)*stress;
      mat2 affine = stress + p_mass * C.at(p);
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) { // Loop over 3x3 grid node neighborhood
          vec2 offset(i, j);
          vec2 dpos = (offset - fx)*dx;
          double weight = w[i].x*w[j].y;
          vec2 index = base + offset;
          grid_v.at(index.x).at(index.y) += weight*(p_mass * v[p] + affine*dpos);
          grid_m.at(index.x).at(index.y) += weight*p_mass;
        }
    }
    for (int i = 0; i < n_particles; i++)
      for (int j = 0; j < n_particles; j++) {
        if (grid_m.at(i).at(j) > 0) {// No need for epsilon here
          grid_v.at(i).at(j) = (1/grid_m.at(i).at(j))*grid_v.at(i).at(j); // Momentum to velocity
          grid_v.at(i).at(j) += dt*gravity*30; // gravity
          vec2 dist = attractor_pos - dx*vec2(i, j);
          grid_v.at(i).at(j) += dist/(0.01 + dist.norm()) * attractor_strength*dt*100;
          if (i < 3 && grid_v.at(i).at(j).x < 0)
            grid_v.at(i).at(j).x = 0; // Boundary conditions
          if (i > n_grid - 3 && grid_v.at(i).at(j).x > 0)
            grid_v.at(i).at(j).x = 0;
          if (j < 3 && grid_v.at(i).at(j).y < 0)
            grid_v.at(i).at(j).y = 0;
          if (j > n_grid - 3 && grid_v.at(i).at(j).y > 0)
            grid_v.at(i).at(j).y = 0;
        }
      }
    for (int p = 0; p < n_particles; p++) { // grid to particle(G2P)
      vec2 base = x.at(p)*inv_dx - 0.5;
      base.x = (int)base.x;
      base.y = (int)base.y;
      vec2 fx = x.at(p)*inv_dx - base;

      vec2 w[3] = { 0.5*(1.5 - fx).square(), 0.75 - (fx - 1).square(), 0.5*(fx -  0.5).square() };
      vec2 new_v;
      mat2 new_C(0, 0, 0, 0);
      for (int i = 0; i < n_particles; i++)
        for (int j = 0; j < n_particles; j++) { // loop over 3x3 grid node neighborhood
          vec2 dpos = vec2(i, j) - fx;
          vec2 offset = base + vec2(i, j);
          vec2 g_v = grid_v.at(offset.x).at(offset.y);
          double weight = w[i].x * w[j].y;
          new_v += weight*g_v;
          new_C += 4*inv_dx*weight*mat2(g_v.x*dpos.x, g_v.x*dpos.y, g_v.y*dpos.x, g_v.y*dpos.y);
        }

      v.at(p) = new_v;
      C.at(p) = new_C;
      x.at(p) += dt*v.at(p); // advection
    }
  };

  auto reset = [&]() {
    int group_size = n_particles%3;
    for (int i = 0; i < n_particles; i++) {
      x.at(i) = vec2(rand()/RAND_MAX*0.2 + 0.3 + 0.1*(i%group_size), rand()/RAND_MAX*0.2 + 0.05 + 0.32*(i%group_size));

      material[i] = i%group_size; // 0: fluid 1: jelly 2: snow
      v.at(i) = vec2();
      F.at(i) = mat2();
      Jp.at(i) = 1;
      C.at(i) = mat2(0, 0, 0, 0);
    }
  };

  cout << "[Hint] Use WSAD/arrow keys to control gravity. Use left/right mouse bottons to attract/repel. Press R to reset.";

  reset();
  gravity = vec2(0, -1);

  for (int frame = 0; frame < 20000; frame++) {
    //if gui.get_event(ti.GUI.PRESS):
    //    if gui.event.key == 'r': reset()
    //    elif gui.event.key in [ti.GUI.ESCAPE, ti.GUI.EXIT]: break
    //if gui.event is not None: gravity[None] = [0, 0]  # if had any event
    //if gui.is_pressed(ti.GUI.LEFT, 'a'): gravity[None][0] = -1
    //if gui.is_pressed(ti.GUI.RIGHT, 'd'): gravity[None][0] = 1
    //if gui.is_pressed(ti.GUI.UP, 'w'): gravity[None][1] = 1
    //if gui.is_pressed(ti.GUI.DOWN, 's'): gravity[None][1] = -1
    //mouse = gui.get_cursor_pos()
    //gui.circle((mouse[0], mouse[1]), color=0x336699, radius=15)
    //attractor_pos[None] = [mouse[0], mouse[1]]
    //attractor_strength[None] = 0
    //if gui.is_pressed(ti.GUI.LMB):
    //    attractor_strength[None] = 1
    //if gui.is_pressed(ti.GUI.RMB):
    //    attractor_strength[None] = -1
    //for s in range(int(2e-3 // dt)):
    //    substep()
    //gui.circles(x.to_numpy(),
    //            radius=1.5,
    //            palette=[0x068587, 0xED553B, 0xEEEEF0],
    //            palette_indices=material)
    //gui.show()  # Change to gui.show(f'{frame:06d}.png') to write images to disk
  }
}

void svd(const mat2 &m,
         mat2 &U,
         mat2 &sig,
         mat2 &V) {

}