/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_MATERIAL_H
#define MPM_MATERIAL_H

#include "pointers.h"
#include <vector>
#include "strength.h"
#include "eos.h"
#include "damage.h"
#include "temperature.h"

class Mat{
public:
  string id;
  int type;
  bool rigid = false;
  class EOS* eos;
  class Strength* strength;
  class Damage* damage;
  vector<class Temperature*> temp;
  double rho0, E, nu, G, K, lambda, signal_velocity;
  Mat(string, int, class EOS*, class Strength*, class Damage* = NULL);
  Mat(string, int, double, double, double);
  Mat(string, bool);
};

class Material : protected Pointers {
 public:
  string id;
  vector<Mat> materials;                        // list of materials
  vector<class EOS *> EOSs;                     // list of defined Equations of State
  vector<class Strength *> strengths;           // list of defined Strengths
  vector<class Damage *> damages;               // list of defined Damage laws
  vector<class Temperature *> temperatures;     // list of defined Temperature laws
  
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
  void add_temperature(vector<string>);
  int find_temperature(string);
  
  typedef Strength *(*StrengthCreator)(MPM *,vector<string>);
  typedef map<string,StrengthCreator> StrengthCreatorMap;
  StrengthCreatorMap *strength_map;

  typedef EOS *(*EOSCreator)(MPM *,vector<string>);
  typedef map<string,EOSCreator> EOSCreatorMap;
  EOSCreatorMap *EOS_map;

  typedef Damage *(*DamageCreator)(MPM *,vector<string>);
  typedef map<string,DamageCreator> DamageCreatorMap;
  DamageCreatorMap *damage_map;

  typedef Temperature *(*TemperatureCreator)(MPM *,vector<string>);
  typedef map<string,TemperatureCreator> TemperatureCreatorMap;
  TemperatureCreatorMap *temperature_map;

  enum constitutive_model {
    LINEAR = 0,      // Linear elasticity
    NEO_HOOKEAN = 1, // Neo-Hookean model
    SHOCK = 2       // EOS + Strength or fluids
  };

private:
  template <typename T> static Strength *strength_creator(MPM *,vector<string>);
  template <typename T> static EOS *EOS_creator(MPM *,vector<string>);
  template <typename T> static Damage *damage_creator(MPM *,vector<string>);
  template <typename T> static Temperature *temperature_creator(MPM *,vector<string>);
  
  const map<string, string> usage = {{"rigid",        "Usage: material(material-ID, rigid)\n"},
				     {"linear",       "Usage: material(material-ID, linear, rho, E, nu, optional: damage-ID)\n"},
				     {"neo-hookean",  "Usage: material(material-ID, neo-hookean, rho, E, nu, optional: damage-ID)\n"},
				     {"eos-strength", "Usage: material(material-ID, eos-strength, eos-ID, strength-ID, optional: damage-ID)\n"}};
  const map<string, int>    Nargs = {{"rigid",        2},
				     {"linear",       5},
				     {"neo-hookean",  5},
				     {"eos-strength", 4}};
};


#endif
