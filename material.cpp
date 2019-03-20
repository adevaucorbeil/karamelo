#include "mpm.h"
#include "material.h"
#include "style_strength.h"
#include "style_eos.h"
#include <vector>

using namespace std;


Material::Material(MPM *mpm) : Pointers(mpm)
{
  strength_map = new StrengthCreatorMap();
  EOS_map = new EOSCreatorMap();

#define STRENGTH_CLASS
#define StrengthStyle(key,Class) \
  (*strength_map)[#key] = &strength_creator<Class>;
#include "style_strength.h"
#undef StrengthStyle
#undef STRENGTH_CLASS

#define EOS_CLASS
#define EOSStyle(key,Class) \
  (*EOS_map)[#key] = &EOS_creator<Class>;
#include "style_eos.h"
#undef EOSStyle
#undef EOS_CLASS
}

Material::~Material()
{
  for (int i = 0; i < strengths.size(); i++) delete strengths[i];
  for (int i = 0; i < EOSs.size(); i++) delete EOSs[i];

  delete strength_map;
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
   create a new strength
------------------------------------------------------------------------- */

void Material::add_strength(vector<string> args){
  cout << "In add_strength" << endl;

  if (find_strength(args[0]) >= 0) {
    cout << "Error: reuse of strength ID" << endl;
    exit(1);
  }

    // create the Strength

  string *estyle = &args[1];

  if (strength_map->find(*estyle) != strength_map->end()) {
    StrengthCreator strength_creator = (*strength_map)[*estyle];
    strengths.push_back(strength_creator(mpm, args));
    //materials.back()->init();
  }
  else {
    cout << "Unknown strength style " << *estyle << endl;
    exit(1);
  }
}

int Material::find_strength(string name)
{
  for (int istrength = 0; istrength < strengths.size(); istrength++)
    if (name.compare(strengths[istrength]->id) == 0) return istrength;
  return -1;
}

/* ----------------------------------------------------------------------
   create a new material
------------------------------------------------------------------------- */

void Material::add_material(vector<string> args){
  cout << "In add_material" << endl;

  if (args.size()<3) {
    cout << "Error: material command not enough arguments" << endl;
    exit(1);
  }

  if (find_material(args[0]) >= 0) {
    cout << "Error: reuse of material ID" << endl;
    exit(1);
  }

  // create the Material
  int iEOS = material->find_EOS(args[1]);

  if (iEOS == -1) {
    cout << "Error: could not find EOS named: " << args[1] << endl;
    exit(1);
  }

  int iStrength = material->find_strength(args[2]);
  if (iStrength == -1) {
    cout << "Error: could not find strength named: " << args[2] << endl;
    exit(1);
  }
  
  cout << "Creating new mat with ID: " << args[0] << endl;
  Mat new_material = {args[0], EOSs[iEOS], strengths[iStrength]};
  materials.push_back(new_material);
}

int Material::find_material(string name)
{
  for (int imaterial = 0; imaterial < materials.size(); imaterial++)
    if (name.compare(materials[imaterial].id) == 0) return imaterial;
  return -1;
}

/* ----------------------------------------------------------------------
   one instance per strength style in style_strength.h
------------------------------------------------------------------------- */

template <typename T>
Strength *Material::strength_creator(MPM *mpm, vector<string> args)
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
