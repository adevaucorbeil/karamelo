#include "mpm.h"
#include "input.h"
#include "output.h"
#include "style_command.h"
#include "domain.h"
#include "material.h"
#include "group.h"
#include "update.h"
#include "variable.h"
#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
#include <stack>
#include "var.h"

#define DELTALINE 256
#define DELTA 4

using namespace std;


Input::Input(MPM *mpm, int argc, char **argv) : Pointers(mpm)
{
  maxline = maxcopy = 0;
  maxarg = 0;
  arg = NULL;

    // fill map with commands listed in style_command.h

  command_map = new CommandCreatorMap();

#define COMMAND_CLASS
#define CommandStyle(key,Class)				\
  (*command_map)[#key] = &command_creator<Class>;
#include "style_command.h"
#undef CommandStyle
#undef COMMAND_CLASS
}


Input::~Input()
{
  delete command_map;
}

/* ----------------------------------------------------------------------
   process all input from infile
   infile = stdin or file if command-line arg "-in" was used
------------------------------------------------------------------------- */

void Input::file()
{
  // Test to create a Var:

  Var v(mpm);

  bool ignore = false;

  istream is(infile);
  while (is) {
    char c = char(is.get());
    if (c != '\n') {

      if (c == '#') ignore = true; // ignore everything after #
      if (c == '\377') {
	cout << line << endl;
	cout << "gives: " << parsev(line).result() << endl;
      } else {
	if (!ignore) line.append(&c,1);
      }

    } else {

      ignore = false;
      cout << line << endl;
      cout << "gives: " << parsev(line).result() << endl;
      line.clear();

    }
  }
}

// Function to find precedence of  
// operators. 
double Input::precedence(char op){ 
    if(op == '+'||op == '-') return 1;
    if(op == '*'||op == '/') return 2;
    if(op == '^') return 3;
    if(op == 'e'|| op == 'E') return 4;
    return 0;
} 
  
// Function to perform arithmetic operations. 
double Input::applyOp(double a, double b, char op){ 
    switch(op){ 
    case '+': return a + b;
    case '-': return a - b;
    case '*': return a * b;
    case '/': return a / b;
    case '^': return pow(a,b);
    case 'e': return a*pow(10,b);
    case 'E': return a*pow(10,b);
    case '(':
      printf("Error: unmatched parenthesis (\n");
      exit(1);
    default:
      printf("Error: unknown operator %c\n", op);
      exit(1);
    } 
}


bool Input::is_operator(char op){
  if (op=='+') return true;
  if (op=='-') return true;
  if (op=='*') return true;
  if (op=='/') return true;
  if (op=='^') return true;
  return false;
}

// check if op is either of +-/*()
bool Input::is_math_char(char op){
  if (op=='+') return true;
  if (op=='-') return true;
  if (op=='*') return true;
  if (op=='/') return true;
  if (op=='^') return true;
  if (op=='(') return true;
  if (op==')') return true;
  if (op=='=') return true;
  return false;
}

// evaluate function func with argument arg:
Var Input::evaluate_function(string func, string arg){
}

// remove white spaces from string
string Input::remove_whitespace(string str){
  string str_;

  for(int i=0; i<str.length(); i++){
    if( str[i] != ' ') str_.append(&str[i],1); // Add the non-whitespace character to str_
  }
  return str_;
}

double Input::parse(string str){
  cout << "Error: deprecated function" << endl;
  exit(1);
}

Var Input::parsev(string str)
{
}


int Input::dimension(vector<string> args){

  if (args.size()==0) {
    cout << "Error: dimension did not receive enough arguments: 1 required" << endl;
    exit(1);
  }

  if (args.size() > 1) {
    cout << "Error: dimension received too many arguments: 1 required" << endl;
    exit(1);
  }

  int dim = (int) parse(args[0]);


  if (dim != 2 && dim != 3) {
    cout << "Error: dimension argument: " << dim << endl;
    exit(1);
  }
  else domain->dimension = dim;

  cout << "Set dimension to " << dim << endl;
  return 0;
}


int Input::region(vector<string> args){
  domain->add_region(args);
  return 0;
}

int Input::solid(vector<string> args){
  domain->add_solid(args);
  return 0;
}

int Input::EOS(vector<string> args){
  material->add_EOS(args);
  return 0;
}

int Input::dump(vector<string> args){
  output->add_dump(args);
  return 0;
}

int Input::group_command(vector<string> args){
  group->assign(args);
  return 0;
}

int Input::log(vector<string> args){
  output->set_log(args);
  return 0;
}

int Input::method_modify(vector<string> args){
  update->modify_method(args);
  return 0;
}

/* ----------------------------------------------------------------------
   one instance per command in style_command.h
------------------------------------------------------------------------- */

template <typename T>
void Input::command_creator(MPM *mpm, vector<string> args)
{
  T cmd(mpm);
  cmd.command(args);
}
