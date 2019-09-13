/* -*- c++ -*- ----------------------------------------------------------*/

#ifdef FIX_CLASS

FixStyle(indent,FixIndent)

#else

#ifndef MPM_FIX_INDENT_H
#define MPM_FIX_INDENT_H

#include "fix.h"
#include "var.h"
#include <vector>

class FixIndent : public Fix {
 public:
  FixIndent(class MPM *, vector<string>);
  ~FixIndent();
  void setmask();
  void init();
  void setup();
  
  void initial_integrate();
  void post_particles_to_grid() {};
  void post_update_grid_state() {};
  void post_grid_to_point() {};
  void post_advance_particles() {};
  void post_velocities_to_grid() {};
  void final_integrate() {};

private:
  string type; // sphere
  int Kpos, xpos, ypos, zpos, Rpos; // Positions of K, the position of the sphere, and its radius in the argument list (args)
};

#endif
#endif

