/* -*- c++ -*- ----------------------------------------------------------
 *
 *                    ***       Karamelo       ***
 *               Parallel Material Point Method Simulator
 * 
 * Copyright (2019) Alban de Vaucorbeil, alban.devaucorbeil@monash.edu
 * Materials Science and Engineering, Monash University
 * Clayton VIC 3800, Australia

 * This software is distributed under the GNU General Public License.
 *
 * ----------------------------------------------------------------------- */

#include "material.h"
#include "error.h"
#include "input.h"
#include "mpm.h"
#include "style_damage.h"
#include "style_eos.h"
#include "style_strength.h"
#include "style_temperature.h"
#include "var.h"
#include <vector>

using namespace std;

/*! Creates the maps of the different known EOS, Strength, Damage, and Temperature.
 * The list of these different laws are fetched in headers style_eos.h, style_strength.h, 
 * style_damage.h, style_temperature.h
 */
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

/*! Destroys the vectors and maps of the known EOS, Strength, Damage, and Temperature.
 */
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


/* ----------------------------------------------------------------------
   create the EOS
------------------------------------------------------------------------- */

/*! This function is the C++ equivalent to the eos() user function.\n
 * Syntax: eos(name, type, type specific arguments)\n
 * This function checks if the EOS name was already used.\n
 * If not, it creates an entry in the vector Material::EOSs and calls EOS::EOS()
 */
void Material::add_EOS(vector<string> args){
  cout << "In add_EOS" << endl;

  if (find_EOS(args[0]) >= 0) {
    error->all(FLERR, "Error: reuse of EOS ID.\n");
  }

    // create the EOS

  if (EOS_map->find(args[1]) != EOS_map->end()) {
    cout << "Create EOS\n";
    EOSCreator EOS_creator = (*EOS_map)[args[1]];
    EOSs.push_back(EOS_creator(mpm, args));
    EOSs.back()->init();
  }
  else {
    error->all(FLERR, "Unknown EOS style " + args[1] + "\n");
  } 
}

/*! This function checks if 'name' is already used for an EOS.\n
 * If an EOS named 'name' exists, it returns its ID. It returns -1 otherwise.
 */
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

/*! This function is the C++ equivalent to the strength() user function.\n
 * Syntax: strength(name, type, type specific arguments)\n
 * This function checks if the Strength name was already used.\n
 * If not, it creates an entry in the vector Material::strengths and calls Strength::Strength()
 */
void Material::add_strength(vector<string> args){
  cout << "In add_strength" << endl;

  if (find_strength(args[0]) >= 0) {
    error->all(FLERR, "Error: reuse of strength ID.\n");
  }

    // create the Strength

  string *estyle = &args[1];

  if (strength_map->find(*estyle) != strength_map->end()) {
    StrengthCreator strength_creator = (*strength_map)[*estyle];
    strengths.push_back(strength_creator(mpm, args));
    //materials.back()->init();
  }
  else {
    error->all(FLERR, "Unknown strength style " + *estyle + ".\n");
  }
}

/*! This function checks if 'name' is already used for a Strength.\n
 * If a Strength named 'name' exists, it returns its ID. It returns -1 otherwise.
 */
int Material::find_strength(string name)
{
  for (int istrength = 0; istrength < strengths.size(); istrength++)
    if (name.compare(strengths[istrength]->id) == 0) return istrength;
  return -1;
}


/* ----------------------------------------------------------------------
   create a new damage
------------------------------------------------------------------------- */

/*! This function is the C++ equivalent to the damage() user function.\n
 * Syntax: damage(name, type, type specific arguments)\n
 * This function checks if the damage name was already used.\n
 * If not, it creates an entry in the vector Material::damages and calls Damage::Damage()
 */
void Material::add_damage(vector<string> args){
  cout << "In add_damage" << endl;

  if (find_damage(args[0]) >= 0) {
    error->all(FLERR, "Error: reuse of damage ID.\n");
  }

    // create the Damage

  string *estyle = &args[1];

  if (damage_map->find(*estyle) != damage_map->end()) {
    DamageCreator damage_creator = (*damage_map)[*estyle];
    damages.push_back(damage_creator(mpm, args));
    //materials.back()->init();
  }
  else {
    error->all(FLERR,"Unknown damage style " + *estyle + "\n");
  }
}

/*! This function checks if 'name' is already used for a Damage.\n
 * If a Damage named 'name' exists, it returns its ID. It returns -1 otherwise.
 */
int Material::find_damage(string name)
{
  for (int idamage = 0; idamage < damages.size(); idamage++)
    if (name.compare(damages[idamage]->id) == 0) return idamage;
  return -1;
}

/* ----------------------------------------------------------------------
   create a new temperature
------------------------------------------------------------------------- */

/*! This function is the C++ equivalent to the temperature() user function.\n
 * Syntax: temperature(name, type, type specific arguments)\n
 * This function checks if the Temperature name was already used.\n
 * If not, it creates an entry in the vector Material::temperatures and calls Temperature::Temperature()
 */
void Material::add_temperature(vector<string> args){
  cout << "In add_temperature" << endl;

  if (find_temperature(args[0]) >= 0) {
    error->all(FLERR, "Error: reuse of temperature ID.\n");
  }

    // create the Temperature

  string *estyle = &args[1];

  if (temperature_map->find(*estyle) != temperature_map->end()) {
    TemperatureCreator temperature_creator = (*temperature_map)[*estyle];
    temperatures.push_back(temperature_creator(mpm, args));
    //materials.back()->init();
  }
  else {
    error->all(FLERR, "Error: unknown temperature style " + *estyle + "\n.");
  }
}

/*! This function checks if 'name' is already used for an EOS.\n
 * If an EOS named 'name' exists, it returns its ID. It returns -1 otherwise.
 */
int Material::find_temperature(string name)
{
  for (int itemperature = 0; itemperature < temperatures.size(); itemperature++)
    if (name.compare(temperatures[itemperature]->id) == 0) return itemperature;
  return -1;
}

/* ----------------------------------------------------------------------
   create a new material
------------------------------------------------------------------------- */

/*! This function is the C++ equivalent to the eos() user function.\n
 * Syntax: eos(name, type, type specific arguments)\n
 * This function checks if the EOS name was already used.\n
 * If not, it creates an entry in the vector Material::EOSs and calls EOS::EOS()
 */
void Material::add_material(vector<string> args) {
  // cout << "In add_material" << endl;

  if (args.size() < 2) {
    string error_str = "Error: material command not enough arguments\n";
    for (auto &x : usage)
      error_str += x.second;
    error->all(FLERR, error_str);
  }

  if (find_material(args[0]) >= 0) {
    error->all(FLERR, "Error: reuse of material ID.\n");
  }

  if (usage.find(args[1]) == usage.end()) {
    string error_str =
        "Error, keyword \033[1;31m" + args[1] + "\033[0m unknown!\n";
    for (auto &x : usage)
      error_str += x.second;
    error->all(FLERR, error_str);
  }

  if (args.size() < Nargs.find(args[1])->second) {
    error->all(FLERR,
               "Error: not enough arguments.\n" + usage.find(args[1])->second);
  }

  if (args[1].compare("linear") == 0 || args[1].compare("neo-hookean") == 0) {
    int type;
    if (args[1].compare("linear") == 0)
      type = LINEAR;
    else
      type = NEO_HOOKEAN;

    materials.push_back(Mat{args[0], type, input->parsev(args[2]),
                            input->parsev(args[3]), input->parsev(args[4])});

  } else if (args[1].compare("rigid") == 0) {
    materials.push_back(Mat{args[0], RIGID});
  } else {
    // create the Material
    int iEOS = material->find_EOS(args[2]);

    if (iEOS == -1) {
      error->all(FLERR, "Error: could not find EOS named: " + args[2] + ".\n");
    }

    int iStrength = material->find_strength(args[3]);
    if (iStrength == -1) {
      error->all(FLERR,
                 "Error: could not find strength named: " + args[3] + ".\n");
    }

    Damage *damage_ = NULL;
    Temperature *temp_ = NULL;

    if (args.size() > Nargs.find(args[1])->second) {
      int iDamage, iTemp;

      iDamage = material->find_damage(args[4]);

      if (iDamage == -1) {
        // args[4] does not correspond to any damage law, maybe it is a
        // temperature law.

        iTemp = material->find_temperature(args[4]);

        if (iTemp == -1) {
          // It is not a temperature law either => error
          error->all(FLERR,
                     "Error: could not find damage named: " + args[4] + ".\n");
          error->all(FLERR, "Error: could not find temperature named: " +
                                args[4] + ".\n");
        } else {
          // It is a temperature law!
          temp_ = temperatures[iTemp];
	  if (args.size() > Nargs.find(args[1])->second + 1) {
	    string error_str = "Error: the last argument of the material command should be the temperature\n";
	    error_str += usage.find(args[1])->second;
	    error->all(FLERR, error_str);
	  }
        }

      } else {
        // args[4] corresponds to a damage law
        damage_ = damages[iDamage];

        if (args.size() == Nargs.find(args[1])->second + 2) {

          iTemp = material->find_temperature(args[5]);

          if (iTemp == -1) {
            error->all(FLERR, "Error: could not find temperature named: " +
                                  args[5] + ".\n");
          } else {
            temp_ = temperatures[iTemp];
          }
        }
      }
    }

    materials.push_back(
        Mat{args[0], SHOCK, EOSs[iEOS], strengths[iStrength], damage_, temp_});
  }
  cout << "Creating new mat with ID: " << args[0] << endl;
}

/*! This function checks if 'name' is already used for a Material.\n
 * If a Material named 'name' exists, it returns its ID. It returns -1 otherwise.
 */
int Material::find_material(string name)
{
  for (int imaterial = 0; imaterial < materials.size(); imaterial++)
    if (name.compare(materials[imaterial].id) == 0) return imaterial;
  return -1;
}

void Material::write_restart(ofstream *of) {
  // Save EOSs:
  size_t N = EOSs.size();
  of->write(reinterpret_cast<const char *>(&N), sizeof(int));

  for (int i = 0; i < N; i++) {
    size_t Neos = EOSs[i]->id.size();
    of->write(reinterpret_cast<const char *>(&Neos), sizeof(size_t));
    of->write(reinterpret_cast<const char *>(EOSs[i]->id.c_str()), Neos);
    cout << "id = " << EOSs[i]->id << endl;

    Neos = EOSs[i]->style.size();
    of->write(reinterpret_cast<const char *>(&Neos), sizeof(size_t));
    of->write(reinterpret_cast<const char *>(EOSs[i]->style.c_str()), Neos);
    EOSs[i]->write_restart(of);
    cout << "style = " << EOSs[i]->style << endl;
  }

  // Save strengths:
  N = strengths.size();
  of->write(reinterpret_cast<const char *>(&N), sizeof(int));

  for (int i = 0; i < N; i++) {
    size_t Nstrengths = strengths[i]->id.size();
    of->write(reinterpret_cast<const char *>(&Nstrengths), sizeof(size_t));
    of->write(reinterpret_cast<const char *>(strengths[i]->id.c_str()), Nstrengths);
    cout << "id = " << strengths[i]->id << endl;

    Nstrengths = strengths[i]->style.size();
    of->write(reinterpret_cast<const char *>(&Nstrengths), sizeof(size_t));
    of->write(reinterpret_cast<const char *>(strengths[i]->style.c_str()), Nstrengths);
    strengths[i]->write_restart(of);
    cout << "style = " << strengths[i]->style << endl;
  }

  // Save damages:
  N = damages.size();
  of->write(reinterpret_cast<const char *>(&N), sizeof(int));

  for (int i = 0; i < N; i++) {
    size_t Ndamages = damages[i]->id.size();
    of->write(reinterpret_cast<const char *>(&Ndamages), sizeof(size_t));
    of->write(reinterpret_cast<const char *>(damages[i]->id.c_str()), Ndamages);
    cout << "id = " << damages[i]->id << endl;

    Ndamages = damages[i]->style.size();
    of->write(reinterpret_cast<const char *>(&Ndamages), sizeof(size_t));
    of->write(reinterpret_cast<const char *>(damages[i]->style.c_str()), Ndamages);
    damages[i]->write_restart(of);
    cout << "style = " << damages[i]->style << endl;
  }

  // Save temperatures:
  N = temperatures.size();
  of->write(reinterpret_cast<const char *>(&N), sizeof(int));

  for (int i = 0; i < N; i++) {
    size_t Ntemperatures = temperatures[i]->id.size();
    of->write(reinterpret_cast<const char *>(&Ntemperatures), sizeof(size_t));
    of->write(reinterpret_cast<const char *>(temperatures[i]->id.c_str()), Ntemperatures);
    cout << "id = " << temperatures[i]->id << endl;

    Ntemperatures = temperatures[i]->style.size();
    of->write(reinterpret_cast<const char *>(&Ntemperatures), sizeof(size_t));
    of->write(reinterpret_cast<const char *>(temperatures[i]->style.c_str()), Ntemperatures);
    temperatures[i]->write_restart(of);
    cout << "style = " << temperatures[i]->style << endl;
  }

  // Save materials:
  N = materials.size();
  of->write(reinterpret_cast<const char *>(&N), sizeof(int));

  for (int i = 0; i < N; i++) {
    size_t Nmats = materials[i].id.size();
    of->write(reinterpret_cast<const char *>(&Nmats), sizeof(size_t));
    of->write(reinterpret_cast<const char *>(materials[i].id.c_str()), Nmats);
    cout << "id = " << materials[i].id << endl;

    of->write(reinterpret_cast<const char *>(&materials[i].type), sizeof(int));
    cout << "type = " << materials[i].type << endl;

    if (materials[i].type == SHOCK) {
      int iEOS = find_EOS(materials[i].eos->id);
      of->write(reinterpret_cast<const char *>(&iEOS), sizeof(int));

      int iStrength = find_strength(materials[i].strength->id);
      of->write(reinterpret_cast<const char *>(&iStrength), sizeof(int));

      int iDamage = -1;
      if (materials[i].damage != NULL) {
	iDamage = find_damage(materials[i].damage->id);
      }
      of->write(reinterpret_cast<const char *>(&iDamage), sizeof(int));

      int iTemperature = -1;
      if (materials[i].temp != NULL) {
	iTemperature = find_temperature(materials[i].temp->id);
      }
      of->write(reinterpret_cast<const char *>(&iTemperature), sizeof(int));
    } else if (materials[i].type == LINEAR || materials[i].type == NEO_HOOKEAN) {
      of->write(reinterpret_cast<const char *>(&materials[i].rho0), sizeof(double));
      of->write(reinterpret_cast<const char *>(&materials[i].E), sizeof(double));
      of->write(reinterpret_cast<const char *>(&materials[i].nu), sizeof(double));
    }
  }

  // -2 flag signals end of Material::write_restart()
  int flag = -2;
  of->write(reinterpret_cast<const char *>(&flag), sizeof(int));
}


void Material::read_restart(ifstream *ifr) {
  cout << "In Material::read_restart" << endl;

  // Pull EOSs:
  size_t N = 0;
  ifr->read(reinterpret_cast<char *>(&N), sizeof(int));
  EOSs.resize(N);

  for (int i = 0; i < N; i++) {
    size_t Neos = 0;
    string id = "";

    ifr->read(reinterpret_cast<char *>(&Neos), sizeof(size_t));
    id.resize(Neos);

    ifr->read(reinterpret_cast<char *>(&id[0]), Neos);
    cout << "id = " << id << endl;

    string style = "";
    ifr->read(reinterpret_cast<char *>(&Neos), sizeof(size_t));
    style.resize(Neos);

    ifr->read(reinterpret_cast<char *>(&style[0]), Neos);
    cout << "style = " << style << endl;
    EOSCreator EOS_creator = (*EOS_map)[style];
    EOSs[i] = EOS_creator(mpm, vector<string>{id, style, "restart"});
    EOSs[i]->read_restart(ifr);
    EOSs[i]->init();
  }

  // Pull Strengths:
  N = 0;
  ifr->read(reinterpret_cast<char *>(&N), sizeof(int));
  strengths.resize(N);

  for (int i = 0; i < N; i++) {
    size_t Nstrengths = 0;
    string id = "";

    ifr->read(reinterpret_cast<char *>(&Nstrengths), sizeof(size_t));
    id.resize(Nstrengths);

    ifr->read(reinterpret_cast<char *>(&id[0]), Nstrengths);
    cout << "id = " << id << endl;

    string style = "";
    ifr->read(reinterpret_cast<char *>(&Nstrengths), sizeof(size_t));
    style.resize(Nstrengths);

    ifr->read(reinterpret_cast<char *>(&style[0]), Nstrengths);
    cout << "style = " << style << endl;
    StrengthCreator strength_creator = (*strength_map)[style];
    strengths[i] = strength_creator(mpm, vector<string>{id, style, "restart"});
    strengths[i]->read_restart(ifr);
    strengths[i]->init();
  }

  // Pull Damages:
  N = 0;
  ifr->read(reinterpret_cast<char *>(&N), sizeof(int));
  damages.resize(N);

  for (int i = 0; i < N; i++) {
    size_t Ndamages = 0;
    string id = "";

    ifr->read(reinterpret_cast<char *>(&Ndamages), sizeof(size_t));
    id.resize(Ndamages);

    ifr->read(reinterpret_cast<char *>(&id[0]), Ndamages);
    cout << "id = " << id << endl;

    string style = "";
    ifr->read(reinterpret_cast<char *>(&Ndamages), sizeof(size_t));
    style.resize(Ndamages);

    ifr->read(reinterpret_cast<char *>(&style[0]), Ndamages);
    cout << "style = " << style << endl;
    DamageCreator damage_creator = (*damage_map)[style];
    damages[i] = damage_creator(mpm, vector<string>{id, style, "restart"});
    damages[i]->read_restart(ifr);
    damages[i]->init();
  }

  // Pull Temperatures:
  N = 0;
  ifr->read(reinterpret_cast<char *>(&N), sizeof(int));
  temperatures.resize(N);

  for (int i = 0; i < N; i++) {
    size_t Ntemperatures = 0;
    string id = "";

    ifr->read(reinterpret_cast<char *>(&Ntemperatures), sizeof(size_t));
    id.resize(Ntemperatures);

    ifr->read(reinterpret_cast<char *>(&id[0]), Ntemperatures);
    cout << "id = " << id << endl;

    string style = "";
    ifr->read(reinterpret_cast<char *>(&Ntemperatures), sizeof(size_t));
    style.resize(Ntemperatures);

    ifr->read(reinterpret_cast<char *>(&style[0]), Ntemperatures);
    cout << "style = " << style << endl;
    TemperatureCreator temperature_creator = (*temperature_map)[style];
    temperatures[i] = temperature_creator(mpm, vector<string>{id, style, "restart"});
    temperatures[i]->read_restart(ifr);
    temperatures[i]->init();
  }

  // Pull Materials:
  N = 0;
  ifr->read(reinterpret_cast<char *>(&N), sizeof(int));

  for (int i = 0; i < N; i++) {
    size_t Nmaterials = 0;
    string id = "";

    ifr->read(reinterpret_cast<char *>(&Nmaterials), sizeof(size_t));
    id.resize(Nmaterials);

    ifr->read(reinterpret_cast<char *>(&id[0]), Nmaterials);
    cout << "id = " << id << endl;

    int type = 0;
    ifr->read(reinterpret_cast<char *>(&type), sizeof(int));

    if (type == SHOCK) {
      int iEOS = -1;
      int iStrength = -1;
      int iDamage = -1;
      int iTemp = -1;

      ifr->read(reinterpret_cast<char *>(&iEOS), sizeof(int));
      ifr->read(reinterpret_cast<char *>(&iStrength), sizeof(int));
      ifr->read(reinterpret_cast<char *>(&iDamage), sizeof(int));
      ifr->read(reinterpret_cast<char *>(&iTemp), sizeof(int));

      Damage *damage_ = NULL;
      Temperature *temp_ = NULL;

      if (iDamage != -1) {
        damage_ = damages[iDamage];
      }
      if (iTemp != -1) {
        temp_ = temperatures[iTemp];
      }

      materials.push_back(
          Mat{id, SHOCK, EOSs[iEOS], strengths[iStrength], damage_, temp_});
    } else if (type == LINEAR || type == NEO_HOOKEAN) {
      double rho0, E, nu = 0;
      ifr->read(reinterpret_cast<char *>(&rho0), sizeof(double));
      ifr->read(reinterpret_cast<char *>(&E), sizeof(double));
      ifr->read(reinterpret_cast<char *>(&nu), sizeof(double));

      materials.push_back(Mat{id, type, rho0, E, nu});
    } else if (type == RIGID) {
      materials.push_back(Mat{id, RIGID});
    } else {
      error->one(FLERR, "Error: unkown material type" + to_string(type) + ".\n");
    }
  }

  // -2 flag signals end of Material::read_restart()
  int flag = 0;
  ifr->read(reinterpret_cast<char *>(&flag), sizeof(int));

  if (flag != -2) {
      error->one(FLERR, "Error: unexpected end to Material::read_restart(): flag = " + to_string(flag) + ". Number of read entities unexpected.\n");    
  }
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

/*! The arguments are: material ID, material type (see Material::constitutive_model)
 * a pointer to an EOS, a pointer to a Strength, a pointer to a Damage, and pointer to a Temperature.
 * The last two can be NULL.
 */
Mat::Mat(string id_, int type_, class EOS* eos_, class Strength* strength_, class Damage* damage_, class Temperature* temp_){
  id = id_;
  type = type_;
  eos = eos_;
  strength = strength_;
  damage = damage_;
  temp = temp_;
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

/*! The arguments are: material ID, material type (see Material::constitutive_model)
 * the density in the reference state, the Young's modulus and Poisson's ratio.
 */
Mat::Mat(string id_, int type_, double rho0_, double E_, double nu_) {
  id = id_;
  type = type_;
  eos = NULL;
  strength = NULL;
  damage = NULL;
  temp = NULL;
  rho0 = rho0_;
  E = E_;
  nu = nu_;
  G = E / (2 * (1 + nu));
  lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
  K = E / (3 * (1 - 2 * nu));
  signal_velocity = sqrt(K / rho0);
}

/*! The arguments are: material ID, material type (see Material::constitutive_model)
 */
Mat::Mat(string id_, int type_) {
  type = type_;
  id = id_;
  rigid = true;
}
