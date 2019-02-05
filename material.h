/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_MATERIAL_H
#define MPM_MATERIAL_H

#include "pointers.h"
#include <vector>
#include "mat.h"
#include "eos.h"

class Material : protected Pointers {
 public:
  string id;
  vector<class Mat *> materials;    // list of defined Materials
  vector<class EOS *> EOSs;    // list of defined Materials
  
  Material(class MPM *);
  virtual ~Material();
  void options(vector<string> *, vector<string>::iterator);

  void add_mat(vector<string>);
  int find_mat(string);
  void add_EOS(vector<string>);
  void set_EOS(vector<string>);
  int find_EOS(string);
  
  typedef Mat *(*MatCreator)(MPM *,vector<string>);
  typedef map<string,MatCreator> MatCreatorMap;
  MatCreatorMap *mat_map;

  typedef EOS *(*EOSCreator)(MPM *,vector<string>);
  typedef map<string,EOSCreator> EOSCreatorMap;
  EOSCreatorMap *EOS_map;

private:
  template <typename T> static Mat *mat_creator(MPM *,vector<string>);
  template <typename T> static EOS *EOS_creator(MPM *,vector<string>);
};

#endif
