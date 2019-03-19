/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_MATERIAL_H
#define MPM_MATERIAL_H

#include "pointers.h"
#include <vector>
#include "strength.h"
#include "eos.h"

struct Mat{
  string id;
  class EOS* eos;
  class Strength* strength;
};

class Material : protected Pointers {
 public:
  string id;
  vector<Mat> materials;                 // list of materials
  vector<class EOS *> EOSs;              // list of defined Equations of State
  vector<class Strength *> strengths;    // list of defined Strengths
  
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
  
  typedef Strength *(*StrengthCreator)(MPM *,vector<string>);
  typedef map<string,StrengthCreator> StrengthCreatorMap;
  StrengthCreatorMap *strength_map;

  typedef EOS *(*EOSCreator)(MPM *,vector<string>);
  typedef map<string,EOSCreator> EOSCreatorMap;
  EOSCreatorMap *EOS_map;

private:
  template <typename T> static Strength *strength_creator(MPM *,vector<string>);
  template <typename T> static EOS *EOS_creator(MPM *,vector<string>);
};


#endif
