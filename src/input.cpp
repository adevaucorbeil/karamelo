/* ----------------------------------------------------------------------
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

#include "input.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "log.h"
#include "material.h"
#include "modify.h"
#include "mpi_wrappers.h"
#include "mpm.h"
#include "output.h"
#include "scheme.h"
#include "style_command.h"
#include "universe.h"
#include "update.h"
#include "var.h"
#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <mpi.h>
#include <stack>
#include <string.h>

#define DELTALINE 256
#define DELTA 4

using namespace std;


Input::Input(MPM *mpm, int argc, char **argv) : Pointers(mpm)
{
  MPI_Comm_rank(universe->uworld,&me);

  line_number = 0;
  maxline = maxcopy = 0;
  maxarg = 0;
  arg = NULL;
  vars = new map<string, Var>;

  (*vars)["time"]     = Var("time", 0);
  (*vars)["timestep"] = Var("timestep", 0);
  (*vars)["dt"]       = Var("dt", 0);
  (*vars)["x"]        = Var("x", 0);
  (*vars)["y"]        = Var("y", 0);
  (*vars)["z"]        = Var("z", 0);
  (*vars)["x0"]       = Var("x0", 0);
  (*vars)["y0"]       = Var("y0", 0);
  (*vars)["z0"]       = Var("z0", 0);
  (*vars)["PI"]       = Var(M_PI);


  // fill map with commands listed in style_command.h

  command_map = new CommandCreatorMap();

#define COMMAND_CLASS
#define CommandStyle(key,Class)				\
  (*command_map)[#key] = &command_creator<Class>;
#include "style_command.h"
#undef CommandStyle
#undef COMMAND_CLASS

  // Protected variables:
  string s = "x";
  protected_vars.push_back(s);
  s = "y";
  protected_vars.push_back(s);
  s = "z";
  protected_vars.push_back(s);
  s = "time";
  protected_vars.push_back(s);
  s = "dt";
  protected_vars.push_back(s);
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
  cout << "In Input::file()\n";
  bool ignore = false;
  int end = 0;

  istream is(infile);

  while(1) {

    line_number++;

    if (me == 0) {
      while(1) {
	char c = char(is.get());
	if (c != '\n') {
	  if (c == '#') ignore = true; // ignore everything after #
	  if (c == '\377') {
	    end = 1;
	    break;
	  } else {
	    if (!ignore) line.append(&c,1);
	  }

	} else {
	  ignore = false;
	  break;
	}
      }
    }

    MPI_string_bcast(line, MPI_CHAR, 0, universe->uworld);
    if (line.compare("quit")==0) break;
    parsev(line).result();

    line.clear();

    MPI_Bcast(&end,1,MPI_INT,0, universe->uworld);
    if (end) break;
  }

}

/*! The higher the returned number, the higher the precedence of the operator.\n
 * Precedence == 5 for the E operator for the nth power of 10.\n
 * Precedence == 4 for the square operators '**' or '^'.\n
 * Precedence == 3 for the multipication and division operators.\n
 * Precedence == 2 for the addition and subtraction operators.\n
 * Precedence == 1 for the other known operators.\n
 * The function returns 0 if the operator is not known.*/
double Input::precedence(const string op){
  if (op[0] == '>') return 1;
  // if (op == ">=") return 1;
  if (op[0] == '<') return 1;
  //if (op == "<=") return 1;
  if (op[0] == '=') return 1;
  //if (op == "==") return 1;
  if (op[0] == '!') return 1;  
  //if (op == "!=") return 1;

  if (op[0] == '+') return 2;
  if (op[0] == '-') return 2;

  if (op[0] == '/') return 3;

  if (op[0] == '^') return 4;

  if (op[0] == 'e') return 5;
  if (op[0] == 'E') return 5;

  if (op == "*") return 3;
  if (op == "**") return 4;
  return 0;
}

/*! This function takes performs the following operation: a 'op' b.\n
 * For instance, if op == '+', it performs a + b.\n
 * Known operators are: '+', '-', '*', '\', '**' and '^' for power, 
 * 'e' and 'E' for the power of 10, and all the ordering operator like '<', '<=', '!=', '==', ...\n
 * It returns the value of the operation, or generates an error if op is a parenthesis.  
 */
Var Input::applyOp(Var a, const string op, Var b){
  int lop = op.length();
  auto it = op.begin();
  if (lop == 1) {
    if (*it=='+') return a + b;
    else if (*it=='-') return a - b;
    else if (*it=='*') return a * b;
    else if (*it=='/') return a / b;
    else if (*it=='^') return a ^ b;
    else if (*it=='e') return a*powv(10,b);
    else if (*it=='>') return a > b;
    else if (*it=='<') return a < b;
    else if (*it=='(') {
      error->all(FLERR, "Error: unmatched parenthesis (\n");
    } else {
      error->all(FLERR, "Error: unknown operator " + op + "\n");
    }
  } else {
    auto it1 = it + 1;
    if (*it=='>' && *it1=='=') return a >= b;
    else if (*it=='<' && *it1=='=') return a <= b;
    else if (*it=='=' && *it1=='=') return a == b;
    else if (*it=='!' && *it1=='=') return a != b;
    else if (*it=='*' && *it1=='*') return a ^ b;
    else {
      error->all(FLERR, "Error: unknown operator " + op + "\n");
    }
  }
}

/*! This function checks if op is a known on-character operator: '+', '-', '*', '/', '^', '<', '>', '!'.\n
 * It returns true if so, false if not.
 */
bool Input::is_operator(char op){
  if (op=='+') return true;
  if (op=='-') return true;
  if (op=='*') return true;
  if (op=='/') return true;
  if (op=='^') return true;
  if (op=='>') return true;
  if (op=='<') return true;
  if (op=='!') return true;
  return false;
}

/*! Checks if the argument is either of: '+', '-', '/', '*', '(', ')'.\n
 * It returns true if so, and false otherwise.
 */
bool Input::is_math_char(char op){
  if (is_operator(op)) return true;
  if (op=='(') return true;
  if (op==')') return true;
  if (op=='=') return true;
  return false;
}

/*! Evaluates the user function func with argument arg.\n
 * An example of function is dimension() that sets the domain dimension, or run() which runs the code.
 */
Var Input::evaluate_function(string func, string arg){
  // cout << "Evaluate function " << func << " with argument: " << arg << endl;

  // Separate arguments:
  vector<string> args;

  int j = 0;
  int start = 0;
  for (int i=0; i<arg.length(); i++) {
    // Locate comas.
    if (arg[i] == ',' || i==arg.length()-1)  {
      if (i==start && i!=arg.length()-1) {
	error->all(FLERR, "Error: missing argument.\n");
      }

      args.resize(args.size()+1);
      args.back().append(&arg[start], i - start + (i==arg.length()-1) );
      //cout << "Received argument " << j+1 << " :" << args.back() << endl;
      start = i+1;
      j++;
    }
  }


  if (func.compare("exp") == 0)
    return expv(parsev(arg));
  if (func.compare("sqrt") == 0)
    return sqrtv(parsev(arg));
  if (func.compare("cos") == 0)
    return cosv(parsev(arg));
  if (func.compare("sin") == 0)
    return sinv(parsev(arg));
  if (func.compare("tan") == 0)
    return tanv(parsev(arg));
  if (func.compare("atan2") == 0) {
    if ((args.size() < 2) || (args.size() > 2)) {
      error->all(FLERR, "Error: atan2 takes exactly two positional arguments.\n");
    }
    return atan2v(parsev(args[0]), parsev(args[1]));
  }
  if (func.compare("log") == 0)
    return logv(parsev(arg));
  if (func.compare("dimension") == 0)
    return Var(dimension(args));
  if (func.compare("axisymmetric") == 0)
    return Var(axisymmetric(args));
  if (func.compare("region") == 0)
    return Var(region(args));
  if (func.compare("solid") == 0)
    return Var(solid(args));
  if (func.compare("eos") == 0)
    return Var(add_EOS(args));
  if (func.compare("strength") == 0)
    return Var(add_strength(args));
  if (func.compare("material") == 0)
    return Var(add_material(args));
  if (func.compare("damage") == 0)
    return Var(add_damage(args));
  if (func.compare("temperature") == 0)
    return Var(add_temperature(args));
  if (func.compare("dump") == 0)
    return Var(dump(args));
  if (func.compare("group") == 0)
    return Var(group_command(args));
  if (func.compare("set_output") == 0)
    return Var(set_output(args));
  if (func.compare("log_modify") == 0)
    return Var(log_modify(args));
  if (func.compare("method") == 0)
    return Var(method(args));
  if (func.compare("scheme") == 0)
    return Var(scheme(args));
  if (func.compare("fix") == 0)
    return Var(fix(args));
  if (func.compare("delete_fix") == 0)
    return Var(delete_fix(args));
  if (func.compare("compute") == 0)
    return Var(compute(args));
  if (func.compare("delete_compute") == 0)
    return Var(delete_compute(args));
  if (func.compare("dt_factor") == 0)
    return Var(set_dt_factor(args));
  if (func.compare("set_dt") == 0)
    return Var(set_dt(args));
  if (func.compare("value") == 0)
    return value(args);
  if (func.compare("plot") == 0)
    return Var(plot(args));
  if (func.compare("create_domain") == 0)
    return Var(create_domain(args));
  if (func.compare("save_plot") == 0)
    return Var(save_plot(args));

  // invoke commands added via style_command.h

  if (command_map->find(func) != command_map->end()) {
    CommandCreator command_creator = (*command_map)[func];
    return command_creator(mpm,args);
  }
  else if (func.compare("evaluate") == 0)
    return Var(parsev(arg).result(mpm));
  else if (func.compare("print") == 0)
    return Var(print(args));
  error->all(FLERR, "Error: Unknown function " + func + "\n");
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
  error->all(FLERR, "Error: Input::parse deprecated function.\n");
}

/*! This function translates the input file syntax, either mathematical expressions or functions, 
 * to C++ and evaluate or execute them.
 */
Var Input::parsev(string str)
{
  // stack to store integer values.
  stack <Var> values;

  // stack to store operators.
  stack <string> ops;

  // stack to store functions.
  stack <string> funcs;

  string returnvar;

  str = remove_whitespace(str);
  // cout << "New string: " << str << endl;

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
      string bracket;
      bracket.push_back(str[i]);
      ops.push(bracket);

      if (i+1 < str.length() && str[i+1] == '-' && !(i+2 < str.length() && str[i+2] == '(')) {
	negative = true;
	i++;
      }
    }

    // Closing brace encountered, solve
    // entire brace.
    else if (str[i] == ')') {
      //cout << "Found \')\'" << endl;

      if (ops.empty() && values.size() < 2) {
	error->all(FLERR, "Error, unmatched parenthesis )\n");
      }

      while(ops.top() != "(")
	{
	  if (values.empty()) {
	    error->all(FLERR, "Error: Ops is not empty with top element being " + ops.top() + ", while values is.\n");
	  } else if (values.size() == 1) {
	    if (ops.top() == "-") {
	      Var val1 = values.top();
	      values.pop();
	      ops.pop();
	      values.push(-val1);
	    } else if (ops.top() == "!") {
	      Var val1 = values.top();
	      values.pop();
	      ops.pop();
	      values.push(!val1);
	    } else if (ops.top() == "+") {
	      ops.pop();
	    } else {
	      Var val1 = values.top();
	      error->all(FLERR, "Error: do not know how to apply " + ops.top() + ", to " + val1.eq() + ".\n");
	    }
	  } else {
	    Var val2 = values.top();
	    values.pop();

	    Var val1 = values.top();
	    values.pop();

	    string op = ops.top();
	    ops.pop();
	    //cout << val1 << " " << val2 << " " << op << endl;

	    values.push(applyOp(val1, op, val2));

	    //cout << "Pushed number: " << values.top() << " to values stack" << endl;
	    if (ops.empty()) {
	      error->all(FLERR, "Error, unmatched parenthesis )\n");
	    }
	  }
	}
      // pop opening brace.
      ops.pop();
    }

    else if (is_operator(str[i]) || (str[i]=='=' && i+1 < str.length() && str[i+1] == '=')){
      string new_op;
      new_op.push_back(str[i]);
      //printf("found operator %c\n", new_op);

      if (values.empty() && !(i+1 < str.length() && str[i+1] == '(')) {
	if (new_op == "-") {
	  negative = true;
	  continue;
	}
      }

      if (i+1 >= str.length()) {
	error->all(FLERR, "Error: end-of-line character detected after operator " + new_op + ".\n");
      }

      if (i+1 < str.length() && str[i+1] == '\0') {
	error->all(FLERR, "Error: end-of-line character detected after operator " + new_op + ".\n");
      }

      else if (i+1 < str.length() && str[i+1] == '*') {
	new_op.push_back('*');
	i++;
      }

      else if (i+1 < str.length() && str[i+1] == '=') {
	new_op.push_back('=');
	i++;
      }

      else if (i+1 < str.length() && is_operator(str[i+1]) && str[i+1] != '-') {
	error->all(FLERR, "Error: unknown operator sequence " + new_op + str[i+1] + ".\n");
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

	string op = ops.top();
	ops.pop();
	//cout << val1 << " " << val2 << " " << op << endl;

	values.push(applyOp(val1, op, val2));
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
	  ops.push("e");

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
	if (i+1 < str.length() && str[i+1] == '=' && str[i+2] != '=') {

	  if (!values.empty() || !ops.empty() ) {
	    error->all(FLERR, "Error: I do not understand when '=' is located in the middle of an expression\n");
	  }

	  else {
	    returnvar = word;
	    //cout << "The computed value will be stored in " <<  returnvar << endl;
	    i++;
	  }
	}

	else if (i+1 >= str.length() && values.empty() && ops.empty() ) {
	  if (negative) {
	    if (!returnvar.empty()) {
	      (*vars)[returnvar] = -(*vars)[word];
	      cout << returnvar << " = " << (*vars)[returnvar].result() << endl;
	    }
	    return -(*vars)[word];
	  }
	  else {
	    if (!returnvar.empty()) {
	      (*vars)[returnvar] = (*vars)[word];
	      cout << returnvar << " = " << (*vars)[returnvar].result() << endl;
	    }
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
	if (protected_variable(returnvar)) {
	  error->all(FLERR, "Error: " + returnvar + " is a protected variable: it cannot be modified!\n");
	}
	//cout << "The computed value will be stored in " <<  returnvar << endl;
	i++;
      }

      else if (i + 1 < str.length() && str[i + 1] == '(')
      {
        // It's a function:
        i += 2;

        // Extract the argument:
        int k             = 0;
        int nLparenthesis = 0;
        int nRparenthesis = 0;
        while (str[i + k] != ')' || nLparenthesis != nRparenthesis)
        {
          if (str[i + k] == '(')
            nLparenthesis++;
          if (str[i + k] == ')')
            nRparenthesis++;
          k++;
          if (i + k > str.length())
          {
            cout << "Error: Unbalanced parenthesis '('" << endl;
            exit(1);
          }
        }
        string arg;
        arg.append(&str[i], k);
        // cout << "Found function " << word << " with argument: " << arg <<
        // endl;
        i += k;
        if (negative) values.push(-evaluate_function(word, arg));
	else values.push(evaluate_function(word, arg));
      }

      else {
	error->all(FLERR, "Error: " + word + " is unknown.\n");
      }
    }
  }

  // Entire expression has been parsed at this
  // point, apply remaining ops to remaining
  // values.
  while(!ops.empty()){
    if (values.empty()) {
      error->all(FLERR, "Error: Ops is not empty with top element being " + ops.top() + ", while values is.\n");
    } else if (values.size() == 1) {
      if (ops.top() == "-") {
	Var val1 = values.top();
	values.pop();
	ops.pop();
	values.push(-val1);
      } else if (ops.top() == "!") {
	Var val1 = values.top();
	values.pop();
	ops.pop();
	values.push(!val1);
      } else if (ops.top() == "+") {
	ops.pop();
      } else {
	Var val1 = values.top();
	error->all(FLERR, "Error: do not know how to apply " + ops.top() + ", to " + val1.eq() + ".\n");
      }
    } else {
      Var val2 = values.top();
      values.pop();

      Var val1 = values.top();
      values.pop();

      string op = ops.top();
      ops.pop();

      values.push(applyOp(val1, op, val2));
    }
  }

  // Top of 'values' contains result, return it.
  if (values.empty()) {
    if (!returnvar.empty()) {
      (*vars)[returnvar] = -1;
    }
    return Var(-1);
  }
  else {
    if (!returnvar.empty()) {
      (*vars)[returnvar] = values.top();
      cout << returnvar << " = " << (*vars)[returnvar].result(mpm) << endl;
    }
    return values.top();
  }
}

int Input::dimension(vector<string> args) {
  // Check that a method is available:
  if (update->method == NULL) {
    error->all(
        FLERR,
        "Error: a method should be defined before calling dimension()!\n");
  }

  domain->set_dimension(args);
  return 0;
}

int Input::axisymmetric(vector<string> args) {
  domain->set_axisymmetric(args);
  return 0;
}

int Input::region(vector<string> args) {
  domain->add_region(args);
  return 0;
}

int Input::solid(vector<string> args) {
  domain->add_solid(args);
  return 0;
}

int Input::add_EOS(vector<string> args) {
  material->add_EOS(args);
  return 0;
}

int Input::add_strength(vector<string> args) {
  material->add_strength(args);
  return 0;
}

int Input::add_damage(vector<string> args) {
  material->add_damage(args);
  return 0;
}

int Input::add_temperature(vector<string> args) {
  material->add_temperature(args);
  return 0;
}

int Input::add_material(vector<string> args) {
  material->add_material(args);
  return 0;
}

int Input::dump(vector<string> args) {
  output->add_dump(args);
  return 0;
}

int Input::group_command(vector<string> args) {
  group->assign(args);
  return 0;
}

int Input::set_output(vector<string> args) {
  output->set_log(args);
  return 0;
}

int Input::log_modify(vector<string> args) {
  output->log->modify(args);
  return 0;
}

int Input::method(vector<string> args) {
  update->create_method(args);
  return 0;
}

int Input::scheme(vector<string> args) {
  update->create_scheme(args);
  return 0;
}

int Input::fix(vector<string> args) {
  modify->add_fix(args);
  return 0;
}

int Input::delete_fix(vector<string> args) {
  modify->delete_fix(args[0]);
  return 0;
}

int Input::compute(vector<string> args) {
  modify->add_compute(args);
  return 0;
}

int Input::delete_compute(vector<string> args) {
  modify->delete_compute(args[0]);
  return 0;
}

int Input::set_dt_factor(vector<string> args) {
  update->set_dt_factor(args);
  return 0;
}

/* This function overruns the CFL timestep.\n
 * By using the set_dt() function, users set their own timestep 
 * that will be fixed during the whole simulation.
 * Warning! This can create problems when the CFL condition is not respected.
 */
int Input::set_dt(vector<string> args) {
  update->set_dt(args);
  return 0;
}

/* The returned value is a constant user-variables that will no longer change.
 */
Var Input::value(vector<string> args) {
  if (args.size() < 1) {
    error->all(FLERR, "Error: too few arguments for command value().\n");
  } else if (args.size() > 1) {
    error->all(FLERR, "Error: too many arguments for command value().\n");
  }
  Var v = parsev(args[0]);
  v.make_constant(mpm);
  return v;
}

int Input::plot(vector<string> args) {
  output->add_plot(args);
  return 0;
}

int Input::save_plot(vector<string> args) {
  if (args.size() < 1) {
    error->all(FLERR, "Error: too few arguments for save_plot(). It should be save_plot(filename)\n");
  } else if (args.size() > 1) {
    error->all(FLERR, "Error: too many arguments for save_plot(). It should be save_plot(filename)\n");
  }
  output->save_plot = true;
  output->ofile_plot = args[0];
  return 0;
}

/* This command is for debugging purposes.
 */
int Input::print(vector<string> args) {
  if (args.size() < 1) {
    error->all(FLERR, "Error: too few arguments for command value().\n");
  } else if (args.size() > 1) {
    error->all(FLERR, "Error: too many arguments for command value().\n");
  }

  Var v = parsev(args[0]);
  cout << args[0] << " = {equation=\"" << v.eq()
       << "\", value=" << v.result(mpm) << ", constant=";

  if (v.is_constant())
    cout << "true";
  else
    cout << "false";

  cout << "}\n";

  return 0;
}
/* ----------------------------------------------------------------------
   one instance per command in style_command.h
------------------------------------------------------------------------- */

template <typename T>
Var Input::command_creator(MPM *mpm, vector<string> args) {
  T cmd(mpm);
  return cmd.command(args);
}

bool Input::protected_variable(string variable) {
  for (int i = 0; i < protected_vars.size(); i++) {
    if (variable.compare(protected_vars[i]) == 0)
      return 1;
  }
  return 0;
}

int Input::create_domain(vector<string> args) {
  domain->create_domain(args);
  return 0;
}
