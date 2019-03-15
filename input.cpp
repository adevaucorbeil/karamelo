#include "mpm.h"
#include "input.h"
#include "output.h"
#include "style_command.h"
#include "domain.h"
#include "material.h"
#include "group.h"
#include "update.h"
#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
#include <stack>
#include "var.h"
#include <map>
#include "modify.h"

#define DELTALINE 256
#define DELTA 4

using namespace std;


Input::Input(MPM *mpm, int argc, char **argv) : Pointers(mpm)
{
  maxline = maxcopy = 0;
  maxarg = 0;
  arg = NULL;
  vars = new map<string, Var>;

  (*vars)["time"] = Var("time", 0); // test


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
  delete vars;
}

/* ----------------------------------------------------------------------
   process all input from infile
   infile = stdin or file if command-line arg "-in" was used
------------------------------------------------------------------------- */

void Input::file()
{
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
  // cout << "t.constant()" << (*vars)["t"].is_constant() << endl;
  // (*vars)["time"] = Var("time", 2, false);
  // cout << "t=" << (*vars)["t"].result() << endl;
  // cout << "t=" << (*vars)["t"].eq() << endl;
  // cout << "t=" << (*vars)["t"].result(mpm) << endl;
  
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
Var Input::applyOp(Var a, Var b, char op){
  cout << "in applyOp with a=" << a.eq() << "=" << a.result() << " b=" << b.eq() << "=" << b.result() << " op=" << op << endl; 
  switch(op){
  case '+': return a + b;
  case '-': return a - b;
  case '*': return a * b;
  case '/': return a / b;
  case '^': return a ^ b;
  case 'e': return a*powv(10,b);
  case 'E': return a*powv(10,b);
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
  cout << "Evaluate function " << func << " with argument: " << arg << endl;

  // Separate arguments:
  vector<string> args;

  int j = 0;
  int start = 0;
  for (int i=0; i<arg.length(); i++) {
    // Locate comas.
    if (arg[i] == ',' || i==arg.length()-1)  {
      if (i==start && i!=arg.length()-1) {
	cout << "Error: missing argument" << endl;
	exit(1);
      }

      args.resize(args.size()+1);
      args.back().append(&arg[start], i - start + (i==arg.length()-1) );
      //cout << "Received argument " << j+1 << " :" << args.back() << endl;
      start = i+1;
      j++;
    }
  }

  if (func.compare("dimension") == 0) return Var(dimension(args));
  if (func.compare("region") == 0) return Var(region(args));
  if (func.compare("solid") == 0) return Var(solid(args));
  if (func.compare("eos") == 0) return Var(EOS(args));
  if (func.compare("dump") == 0) return Var(dump(args));
  if (func.compare("group") == 0) return Var(group_command(args));
  if (func.compare("log") == 0) return Var(log(args));
  if (func.compare("method_modify") == 0) return Var(method_modify(args));
  if (func.compare("fix") == 0) return Var(fix(args));

  // invoke commands added via style_command.h

  if (command_map->find(func) != command_map->end()) {
    CommandCreator command_creator = (*command_map)[func];
    command_creator(mpm,args);
    return Var(0);
  }

  else if (func.compare("exp") == 0) return expv(parsev(arg));
  cout << "Error: Unknown function " << func << endl;
  exit(1);
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
  cout << "Error: Input::parse deprecated function" << endl;
  exit(1);
}

Var Input::parsev(string str)
{
  // stack to store integer values.
  stack <Var> values;

  // stack to store operators.
  stack <char> ops;

  // stack to store functions.
  stack <string> funcs;

  string returnvar;

  str = remove_whitespace(str);
  //cout << "New string: " << str << endl;

  bool negative = false; // flag that indicates that the next value will have to be multiplied by -1

  for (int i=0; i<str.length(); i++) {
    if (isdigit(str[i])) {

      string number;

      // The number could be longer than one digit:
      int j;
      for (j=1; j<str.length()-i;j++) {
	if (!isdigit(str[i+j]) && str[i+j]!='.') break;
      }

      number.append(&str[i],j);
      i += j-1;

      if (negative) {
	Var result((-1)*stof(number));
	values.push(result);
	negative = false;
      } else {
	values.push(Var (stof(number)));
      }

      //cout << "Pushed number: " << values.top() << " to values stack" << endl;
    }

    // If the first character in the righ-hand-side is (, push it to ops
    else if (str[i] == '(') {
      //cout << "Pushed \'(\' to ops stack" << endl;
      ops.push(str[i]);

      if (i+1 < str.length() && str[i+1] == '-') {
	negative = true;
	i++;
      }
    }

    // Closing brace encountered, solve
    // entire brace.
    else if (str[i] == ')') {
      //cout << "Found \')\'" << endl;

      if (ops.empty() && values.size() < 2) {
	printf("Error, unmatched parenthesis )\n");
	exit(1);
      }

      while(ops.top() != '(')
	{
	  Var val2 = values.top();
	  values.pop();

	  Var val1 = values.top();
	  values.pop();

	  char op = ops.top();
	  ops.pop();
	  //cout << val1 << " " << val2 << " " << op << endl;

	  values.push(applyOp(val1, val2, op));

	  //cout << "Pushed number: " << values.top() << " to values stack" << endl;
	  if (ops.empty()) {
	    printf("Error, unmatched parenthesis )\n");
	    exit(1);
	  }
	}

      // pop opening brace.
      ops.pop();
    }

    else if (is_operator(str[i])){
      char new_op = str[i];
      //printf("found operator %c\n", new_op);

      if (values.empty()) {
	if (new_op == '-') {
	  negative = true;
	  continue;
	}
      }

      if (i+1 >= str.length()) {
	printf("Error: end-of-line character detected after operator %c\n", new_op);
	exit(1);
      }

      if (i+1 < str.length() && str[i+1] == '\0') {
	printf("Error: end-of-line character detected after operator %c\n", new_op);
	exit(1);
      }

      else if (i+1 < str.length() && str[i+1] == '*') {
	new_op = '^';
	i++;
      }

      else if (i+1 < str.length() && is_operator(str[i+1]) && str[i+1] != '-') {
	printf("Error: unknown operator sequence %c%c\n", new_op, str[i+1]);
	exit(1);
      }


      // While top of 'ops' has same or greater
      // precedence to current token, which
      // is an operator. Apply operator on top
      // of 'ops' to top two elements in values stack.
      while(!ops.empty() &&
	    precedence(ops.top()) >= precedence(new_op)){
	Var val2 = values.top();
	values.pop();

	Var val1 = values.top();
	values.pop();

	char op = ops.top();
	ops.pop();
	//cout << val1 << " " << val2 << " " << op << endl;

	values.push(applyOp(val1, val2, op));
      }

      // Push current token to 'ops'.
      ops.push(new_op);

      if (i+1 < str.length() && str[i+1] == '-') {
	negative = true;
	i++;
      }
    }


    else {
      string word;

      // Check how long is the word
      int j;
      for (j=1; j<str.length()-i;j++) {
	if (is_math_char(str[i+j])) break;
      }

      word.append(&str[i],j);
      i += j-1;

      //cout << "Found keyword: " << word << endl;

      if (word == "E" || word == "e") { // E or e have to be followed by + or - to indicate that it is 10^+xx
      	if (!values.empty() && isdigit(str[i-1]) && i+1 < str.length() && (str[i+1] == '+' || str[i+1] == '-')) {
	  // Push current token to 'ops'.
	  ops.push('e');

	  if (str[i+1] == '-') {
	    negative = true;
	    i++;
	  }

	  if (str[i+1] == '+') i++;
	  continue;
	}
      }

      // Check if word is a variable:
      map<string, Var>::iterator it;
      it = vars->find(word);

      if (it != vars->end()){
	// word is a variable
	//cout << "word is a variable\n";
	if (i+1 < str.length() && str[i+1] == '=') {

	  if (!values.empty() || !ops.empty() ) {
	    printf("Error: I do not understand when '=' is located in the middle of an expression\n");
	    exit(1);
	  }

	  else {
	    returnvar = word;
	    //cout << "The computed value will be stored in " <<  returnvar << endl;
	    i++;
	  }
	}

	else if (i+1 >= str.length() && values.empty() && ops.empty() ) {
	  if (negative) {
	    return -(*vars)[word];
	  }
	  else {
	    return (*vars)[word];
	  }
	}

	else {
	  if (negative) {
	    values.push(-(*vars)[word]);
	    negative = false;
	  } else {
	    values.push((*vars)[word]);
	  }
	  //cout << "push " << word << "=" << values.top() << " to values\n";
	}
      }


      else if (i+1 < str.length() && str[i+1]=='=' && values.empty() && ops.empty()){
	// Check if there is an '=':
	//cout << "Check if there is an =\n";
	returnvar = word;
	//cout << "The computed value will be stored in " <<  returnvar << endl;
	i++;
      }

      else if (i+1 < str.length() && str[i+1]=='(') {
	// It's a function:
	i+=2;

	// Extract the argument:
	int k = 0;
	int nLparenthesis = 0;
	int nRparenthesis = 0;
	while(str[i+k]!=')' || nLparenthesis != nRparenthesis){
	  if (str[i+k]=='(') nLparenthesis++;
	  if (str[i+k]==')') nRparenthesis++;
	  k++;
	  if (i+k > str.length()) {
	    cout << "Error: Unbalanced parenthesis '('" << endl;
	    exit(1);
	  }
	}
	string arg;
	arg.append(&str[i],k);
	//cout << "Found function " << word << " with argument: " << arg << endl;
	i += k;
	values.push(evaluate_function(word, arg));
      }

      else {
	cout << "Error: " << word << " is unknown\n";
	exit(1);
      }
    }
  }

  // Entire expression has been parsed at this
  // point, apply remaining ops to remaining
  // values.
  while(!ops.empty()){
    Var val2 = values.top();
    values.pop();

    Var val1 = values.top();
    values.pop();

    char op = ops.top();
    ops.pop();

    values.push(applyOp(val1, val2, op));
  }

  // Top of 'values' contains result, return it.
  if (values.empty()) {
    return Var(-1);
  }
  else {
    //cout << "value = " << values.top() << endl;
    if (!returnvar.empty()) {
      (*vars)[returnvar] = values.top();
    }
    return values.top();
  }
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

  int dim = (int) parsev(args[0]);


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

int Input::fix(vector<string> args){
  modify->add_fix(args);
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
