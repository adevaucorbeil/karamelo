/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_MATERIAL_H
#define MPM_MATERIAL_H

#include "pointers.h"
#include <vector>
#include "strength.h"
#include "eos.h"
#include "damage.h"

class Mat{
public:
  string id;
  int type;
  class EOS* eos;
  class Strength* strength;
  class Damage* damage;
  double rho0, E, nu, G, K, lambda, signal_velocity;
  Mat(string, int, class EOS*, class Strength*, class Damage* = NULL);
  Mat(string, int, double, double, double);
};

class Material : protected Pointers {
 public:
  string id;
  vector<Mat> materials;                 // list of materials
  vector<class EOS *> EOSs;              // list of defined Equations of State
  vector<class Strength *> strengths;    // list of defined Strengths
  vector<class Damage *> damages;        // list of defined Damage laws
  
  Material(class MPM *);
  virtual ~Material();
  void options(vector<string> *, vector<string>::iterator);

  void add_strength(vector<string>);
  int find_strength(string);
  void add_EOS(vector<string>);
  void set_EOS(vector<string>);
  int find_EOS(string);
  void add_material(vector<string>);
  int find_material(string);
  void add_damage(vector<string>);
  int find_damage(string);
  
  typedef Strength *(*StrengthCreator)(MPM *,vector<string>);
  typedef map<string,StrengthCreator> StrengthCreatorMap;
  StrengthCreatorMap *strength_map;

  typedef EOS *(*EOSCreator)(MPM *,vector<string>);
  typedef map<string,EOSCreator> EOSCreatorMap;
  EOSCreatorMap *EOS_map;

  typedef Damage *(*DamageCreator)(MPM *,vector<string>);
  typedef map<string,DamageCreator> DamageCreatorMap;
  DamageCreatorMap *damage_map;

  enum constitutive_model {
    LINEAR = 0,      // Linear elasticity
    NEO_HOOKEAN = 1, // Neo-Hookean model
    SHOCK = 2        // EOS + Strength
  };

private:
  template <typename T> static Strength *strength_creator(MPM *,vector<string>);
  template <typename T> static EOS *EOS_creator(MPM *,vector<string>);
  template <typename T> static Damage *damage_creator(MPM *,vector<string>);
};


#endif
