#include "mpm.h"
#include "material.h"
#include "input.h"
#include "var.h"
#include "style_strength.h"
#include "style_eos.h"
#include "style_damage.h"
#include "style_temperature.h"
#include <vector>

using namespace std;


Material::Material(MPM *mpm) : Pointers(mpm)
{
  strength_map = new StrengthCreatorMap();
  EOS_map = new EOSCreatorMap();
  damage_map = new DamageCreatorMap();
  temperature_map = new TemperatureCreatorMap();

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

#define DAMAGE_CLASS
#define DamageStyle(key,Class) \
  (*damage_map)[#key] = &damage_creator<Class>;
#include "style_damage.h"
#undef DamageStyle
#undef DAMAGE_CLASS

#define TEMPERATURE_CLASS
#define TemperatureStyle(key,Class) \
  (*temperature_map)[#key] = &temperature_creator<Class>;
#include "style_temperature.h"
#undef TemperatureStyle
#undef TEMPERATURE_CLASS
}

Material::~Material()
{
  for (int i = 0; i < strengths.size(); i++) delete strengths[i];
  for (int i = 0; i < EOSs.size(); i++) delete EOSs[i];
  for (int i = 0; i < damages.size(); i++) delete damages[i];
  for (int i = 0; i < temperatures.size(); i++) delete temperatures[i];

  delete strength_map;
  delete EOS_map;
  delete damage_map;
  delete temperature_map;
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
   create a new damage
------------------------------------------------------------------------- */

void Material::add_damage(vector<string> args){
  cout << "In add_damage" << endl;

  if (find_damage(args[0]) >= 0) {
    cout << "Error: reuse of damage ID" << endl;
    exit(1);
  }

    // create the Damage

  string *estyle = &args[1];

  if (damage_map->find(*estyle) != damage_map->end()) {
    DamageCreator damage_creator = (*damage_map)[*estyle];
    damages.push_back(damage_creator(mpm, args));
    //materials.back()->init();
  }
  else {
    cout << "Unknown damage style " << *estyle << endl;
    exit(1);
  }
}

int Material::find_damage(string name)
{
  for (int idamage = 0; idamage < damages.size(); idamage++)
    if (name.compare(damages[idamage]->id) == 0) return idamage;
  return -1;
}

/* ----------------------------------------------------------------------
   create a new temperature
------------------------------------------------------------------------- */

void Material::add_temperature(vector<string> args){
  cout << "In add_temperature" << endl;

  if (find_temperature(args[0]) >= 0) {
    cout << "Error: reuse of temperature ID" << endl;
    exit(1);
  }

    // create the Temperature

  string *estyle = &args[1];

  if (temperature_map->find(*estyle) != temperature_map->end()) {
    TemperatureCreator temperature_creator = (*temperature_map)[*estyle];
    temperatures.push_back(temperature_creator(mpm, args));
    //materials.back()->init();
  }
  else {
    cout << "Unknown temperature style " << *estyle << endl;
    exit(1);
  }
}

int Material::find_temperature(string name)
{
  for (int itemperature = 0; itemperature < temperatures.size(); itemperature++)
    if (name.compare(temperatures[itemperature]->id) == 0) return itemperature;
  return -1;
}

/* ----------------------------------------------------------------------
   create a new material
------------------------------------------------------------------------- */

void Material::add_material(vector<string> args){
  // cout << "In add_material" << endl;

  if (args.size()<2) {
    cout << "Error: material command not enough arguments" << endl;
    for (auto& x: usage) cout << x.second;
    exit(1);
  }

  if (find_material(args[0]) >= 0) {
    cout << "Error: reuse of material ID" << endl;
    exit(1);
  }

  if (usage.find(args[1]) == usage.end()) {
    cout << "Error, keyword \033[1;31m" << args[1] << "\033[0m unknown!\n";
    for (auto& x: usage) cout << x.second;
    exit(1);
  }

  if (args.size() < Nargs.find(args[1])->second) {
    cout << "Error: not enough arguments.\n";
    cout << usage.find(args[1])->second;
    exit(1);
  }

  // if (args.size() > Nargs.find(args[1])->second) {
  //   cout << "Error: too many arguments.\n";
  //   cout << usage.find(args[1])->second;
  //   exit(1);
  // }

  if (args[1].compare("linear")==0 || args[1].compare("neo-hookean")==0) {
    int type;
    if (args[1].compare("linear")==0) type = LINEAR;
    else type = NEO_HOOKEAN;

    Mat new_material(args[0], type, input->parsev(args[2]), input->parsev(args[3]), input->parsev(args[4]));
    materials.push_back(new_material);
    
  } else if (args[1].compare("rigid")==0) {
    Mat new_material(args[0], true);
    materials.push_back(new_material);
  } else {
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

    if (args.size() > 3) {
      int iDamage = material->find_damage(args[3]);
      if (iDamage == -1) {
	cout << "Error: could not find damage named: " << args[3] << endl;
	exit(1);
      }
      Mat new_material(args[0], SHOCK, EOSs[iEOS], strengths[iStrength], damages[iDamage]);
      materials.push_back(new_material);
    } else {
      Mat new_material(args[0], SHOCK, EOSs[iEOS], strengths[iStrength]);
      materials.push_back(new_material);
    }
  }

  int iTemp;
  for(int i=Nargs.find(args[1])->second; i<args.size(); i++) {
    iTemp = material->find_temperature(args[i]);
    if (iTemp == -1) {
      cout << "Error: could not find temperature named: " << args[i] << endl;
      exit(1);
    }
    materials.back().temp.push_back(temperatures[iTemp]);
  }

  cout << "Creating new mat with ID: " << args[0] << endl;
  
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

/* ----------------------------------------------------------------------
   one instance per damage style in style_damage.h
------------------------------------------------------------------------- */

template <typename T>
Damage *Material::damage_creator(MPM *mpm, vector<string> args)
{
  return new T(mpm, args);
}

/* ----------------------------------------------------------------------
   one instance per temperature style in style_temperature.h
------------------------------------------------------------------------- */

template <typename T>
Temperature *Material::temperature_creator(MPM *mpm, vector<string> args)
{
  return new T(mpm, args);
}


Mat::Mat(string id_, int type_, class EOS* eos_, class Strength* strength_, class Damage* damage_){
  id = id_;
  type = type_;
  eos = eos_;
  strength = strength_;
  damage = damage_;
  rho0 = eos->rho0();
  K = eos->K();
  G = strength->G();
  E = 9*K*G/(3*K+G);
  nu = (3*K-2*G)/(2*(3*K+G));
  lambda = K - 2*G/3;
  signal_velocity = sqrt((lambda+2*G)/rho0);

  cout << "Properties for material " << id << endl;
  cout << "\tReference density: " << rho0 << endl;
  cout << "\tYoung\'s modulus: " << E << endl;
  cout << "\tPoisson\'s ratio: " << nu << endl;
  cout << "\tShear modulus: " << G << endl;
  cout << "\tBulk modulus: " << K << endl;
  cout << "\tLame first parameter (Lambda): " << lambda << endl;
  cout << "\tSignal velocity: " << signal_velocity << endl;
}

Mat::Mat(string id_, int type_, double rho0_, double E_, double nu_){
  id = id_;
  type = type_;
  eos = NULL;
  strength = NULL;
  damage = NULL;
  rho0 = rho0_;
  E = E_;
  nu = nu_;
  G = E/(2*(1+nu));
  lambda = E*nu/((1+nu)*(1-2*nu));
  K = E/(3*(1-2*nu));
  signal_velocity = sqrt(K/rho0);
}

Mat::Mat(string id_, bool rigid_){
  id = id_;
  rigid = rigid_;
  if (!rigid_) {
    cout << "Error in Mat::Mat(string id_, bool rigid_), rigid_==false.\n";
    exit(1);
  }
}
