#include "mpm.h"
#include "material.h"
#include "style_mat.h"
#include "style_eos.h"
#include <vector>

using namespace std;


Material::Material(MPM *mpm) : Pointers(mpm)
{
  mat_map = new MatCreatorMap;
  EOS_map = new EOSCreatorMap;

#define MAT_CLASS
#define MatStyle(key,Class) \
  (*ma_map)[#key] = &mat_creator<Class>;
#include "style_mat.h"
#undef MatStyle
#undef MAT_CLASS

#define EOS_CLASS
#define EOSStyle(key,Class) \
  (*EOS_map)[#key] = &EOS_creator<Class>;
#include "style_eos.h"
#undef EOSStyle
#undef EOS_CLASS
}

Material::~Material()
{
  delete mat_map;
  delete EOS_map;
}


void Material::options(vector<string> *args, vector<string>::iterator it)
{
  cout << "In material::options()" << endl;
  if (args->end() < it) {
    cout << "Error: not enough arguments" << endl;
    exit(1);
  }
  if (args->end() > it) {
    cout << "Ignoring optional arguments: ";
    for (it; it != args->end(); ++it){
      cout << *it << "\t";
    }
    cout << endl;
  }
}

/* ----------------------------------------------------------------------
   create the EOS
------------------------------------------------------------------------- */

void Material::add_EOS(vector<string> args){
  cout << "In add_EOS" << endl;

  if (find_EOS(args[0]) >= 0) {
    cout << "Error: reuse of EOS ID" << endl;
    exit(1);
  }

    // create the EOS

  if (EOS_map->find(args[1]) != EOS_map->end()) {
    cout << "Create EOS\n";
    EOSCreator EOS_creator = (*EOS_map)[args[1]];
    EOSs.push_back(EOS_creator(mpm, args));
    EOSs.back()->init();
  }
  else {
    cout << "Unknown EOS style " << args[1] << endl;
    exit(1);
  }
  
}

void Material::set_EOS(vector<string> args){
  cout << "In Material::set_EOS" << endl;
}

int Material::find_EOS(string name)
{
  cout << "In find_EOS\n";
  for (int iEOS = 0; iEOS < EOSs.size(); iEOS++) {
    cout << "EOSs["<< iEOS <<"]->id=" << EOSs[iEOS]->id << endl;
    if (name.compare(EOSs[iEOS]->id) == 0) return iEOS;
  }
  return -1;
}

/* ----------------------------------------------------------------------
   create a new material
------------------------------------------------------------------------- */

void Material::add_mat(vector<string> args){
  cout << "In add_mat" << endl;

  if (find_mat(args[0]) >= 0) {
    cout << "Error: reuse of material ID" << endl;
    exit(1);
  }

    // create the Material

  string *estyle = &args[1];

  if (mat_map->find(*estyle) != mat_map->end()) {
    MatCreator mat_creator = (*mat_map)[*estyle];
    materials.push_back(mat_creator(mpm, args));
    //materials.back()->init();
  }
  else {
    cout << "Unknown mat style " << *estyle << endl;
    exit(1);
  }
  
}

int Material::find_mat(string name)
{
  for (int imat = 0; imat < materials.size(); imat++)
    if (name.compare(materials[imat]->id) == 0) return imat;
  return -1;
}

/* ----------------------------------------------------------------------
   one instance per mat style in style_mat.h
------------------------------------------------------------------------- */

template <typename T>
Mat *Material::mat_creator(MPM *mpm, vector<string> args)
{
  return new T(mpm, args);
}

/* ----------------------------------------------------------------------
   one instance per EOS style in style_eos.h
------------------------------------------------------------------------- */

template <typename T>
EOS *Material::EOS_creator(MPM *mpm, vector<string> args)
{
  return new T(mpm, args);
}
