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

#include <input.h>
#include <domain.h>
#include <error.h>
#include <group.h>
#include <log.h>
#include <material.h>
#include <modify.h>
#include <mpi_wrappers.h>
#include <mpm.h>
#include <output.h>
#include <scheme.h>
#include <style_command.h>
#include <universe.h>
#include <update.h>
#include <var.h>
#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <mpi.h>
#include <stack>
#include <string>
#include <regex>

#include <expression_function_parenthesis.h>
#include <expression_function_evaluate.h>
#include <expression_functions_basic.h>
#include <expression_functions_error_gamma.h>
#include <expression_functions_exponential.h>
#include <expression_functions_hyperbolic.h>
#include <expression_functions_nearest_integer.h>
#include <expression_functions_power.h>
#include <expression_functions_trigonometric.h>

#include <expression_operand_constant.h>
#include <expression_operand_expression.h>
#include <expression_operand_literal.h>
#include <expression_operand_index.h>
#include <expression_operand_vector.h>
#include <expression_operand_time.h>

#include <expression_operators_arithmetic.h>
#include <style_expression.h>

#define DELTALINE 256
#define DELTA 4

using namespace std;

const bool Input::DEBUG_EXPRESSIONS = false;

Input::Input(MPM *mpm, int argc, char **argv) : Pointers(mpm)
{
  MPI_Comm_rank(universe->uworld,&me);

  line_number = 0;
  maxline = maxcopy = 0;
  maxarg = 0;
  arg = nullptr;
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
#include <style_command.h>
#undef CommandStyle
#undef COMMAND_CLASS

  // Protected variables:
  string s = "x";
  protected_vars.push_back(s);
  s = "y";
  protected_vars.push_back(s);
  s = "z";
  protected_vars.push_back(s);
  s = "x0";
  protected_vars.push_back(s);
  s = "y0";
  protected_vars.push_back(s);
  s = "z0";
  protected_vars.push_back(s);
  s = "time";
  protected_vars.push_back(s);
  s = "dt";
  protected_vars.push_back(s);
  
  operation_factory.register_class<ExpressionFunctionParenthesis>("(");
  operation_factory.register_class<ExpressionFunctionEvaluate>("evaluate(");

  operation_factory.register_class<ExpressionFunctionAbs>("abs(");
  operation_factory.register_class<ExpressionFunctionRemainder>("remainder(");
  operation_factory.register_class<ExpressionFunctionMax>("max(");
  operation_factory.register_class<ExpressionFunctionMin>("min(");
  
  operation_factory.register_class<ExpressionFunctionErf>("erf(");
  operation_factory.register_class<ExpressionFunctionErfc>("erfc(");
  operation_factory.register_class<ExpressionFunctionTgamma>("tgamma(");
  operation_factory.register_class<ExpressionFunctionLgamma>("lgamma(");
  
  operation_factory.register_class<ExpressionFunctionExp>("exp(");
  operation_factory.register_class<ExpressionFunctionExp2>("exp2(");
  operation_factory.register_class<ExpressionFunctionExpm1>("expm1(");
  operation_factory.register_class<ExpressionFunctionLog>("log(");
  operation_factory.register_class<ExpressionFunctionLog10>("log10(");
  operation_factory.register_class<ExpressionFunctionLog2>("log2(");
  operation_factory.register_class<ExpressionFunctionLog1p>("log1p(");

  operation_factory.register_class<ExpressionFunctionSinh>("sinh(");
  operation_factory.register_class<ExpressionFunctionCosh>("cosh(");
  operation_factory.register_class<ExpressionFunctionTanh>("tanh(");
  operation_factory.register_class<ExpressionFunctionAsinh>("asinh(");
  operation_factory.register_class<ExpressionFunctionAcosh>("acosh(");
  operation_factory.register_class<ExpressionFunctionAtanh>("atanh(");

  operation_factory.register_class<ExpressionFunctionCeil>("ceil(");
  operation_factory.register_class<ExpressionFunctionFloor>("floor(");
  operation_factory.register_class<ExpressionFunctionTrunc>("trunc(");

  operation_factory.register_class<ExpressionFunctionPow>("pow(");
  operation_factory.register_class<ExpressionFunctionSqrt>("sqrt(");
  operation_factory.register_class<ExpressionFunctionCbrt>("cbrt(");
  operation_factory.register_class<ExpressionFunctionHypot>("hypot(");

  operation_factory.register_class<ExpressionFunctionSin>("sin(");
  operation_factory.register_class<ExpressionFunctionCos>("cos(");
  operation_factory.register_class<ExpressionFunctionTan>("tan(");
  operation_factory.register_class<ExpressionFunctionAsin>("asin(");
  operation_factory.register_class<ExpressionFunctionAcos>("acos(");
  operation_factory.register_class<ExpressionFunctionAtan>("atan(");
  operation_factory.register_class<ExpressionFunctionAtan2>("atan2(");

  infix_factory.register_class<ExpressionSum       >("+");
  infix_factory.register_class<ExpressionDifference>("-");
  infix_factory.register_class<ExpressionProduct   >("*");
  infix_factory.register_class<ExpressionQuotient  >("/");
  operation_factory.register_class<ExpressionNegation>("-");

  operation_factory.register_class<ExpressionOperandConstant<Kokkos::Experimental::    pi_v<float>>>("PI"    );
  operation_factory.register_class<ExpressionOperandConstant<Kokkos::Experimental::     e_v<float>>>("E"     );
  operation_factory.register_class<ExpressionOperandConstant<Kokkos::Experimental::egamma_v<float>>>("EGAMMA");
  operation_factory.register_class<ExpressionOperandConstant<Kokkos::Experimental::   phi_v<float>>>("PHI"   );

  operation_factory.register_class<ExpressionOperandIndex>("i");

  operation_factory.register_class<ExpressionOperandVector<&Solid::x, &Grid::x, 0>>("x");
  operation_factory.register_class<ExpressionOperandVector<&Solid::x, &Grid::x, 1>>("y");
  operation_factory.register_class<ExpressionOperandVector<&Solid::x, &Grid::x, 2>>("z");

  operation_factory.register_class<ExpressionOperandVector<&Solid::x0, &Grid::x0, 0>>("x0");
  operation_factory.register_class<ExpressionOperandVector<&Solid::x0, &Grid::x0, 1>>("y0");
  operation_factory.register_class<ExpressionOperandVector<&Solid::x0, &Grid::x0, 2>>("z0");

  operation_factory.register_class<ExpressionOperandTime>("time");

#define EXPRESSION_CLASS
#define ExpressionStyle(key, Class) operation_factory.register_class<Class>(#key"(");
#include <style_expression.h>
#undef ExpressionStyle
#undef EXPRESSION_CLASS

  if (DEBUG_EXPRESSIONS)
    for (const string &str: {
    "1/2/3",
    "1 - 2 - 3",
    "sin(cos(2)) + 1",
    "cos(2)",
    "cos(PI)",
    "-1",
    "-sqrt(14.386*-sin(4))",
    "-PI*-1",
    "1 - -PI"
    })
    {
      cout << "----- Parsing \"" << str << "\" -----" << endl;
      parsev(str);
      cout << "-------- Evaluating --------" << endl;
      Expression &expression = expressions.at(str);
      expression.evaluate();
      cout << str << " = " << expression.registers(0, 0) << endl;
      cout << "----------------------------" << endl << endl;
    }
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
  // cout << "In Input::file()\n";
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
float Input::precedence(const string op){
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
  return Var();
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


  if (func == "exp")
    return expv(parsev(arg));
  if (func == "sqrt")
    return sqrtv(parsev(arg));
  if (func == "cos")
    return cosv(parsev(arg));
  if (func == "sin")
    return sinv(parsev(arg));
  if (func == "tan")
    return tanv(parsev(arg));
  if (func == "atan2") {
    if ((args.size() < 2) || (args.size() > 2)) {
      error->all(FLERR, "Error: atan2 takes exactly two positional arguments.\n");
    }
    return atan2v(parsev(args[0]), parsev(args[1]));
  }
  if (func == "if") {
    if ((args.size() < 3) || (args.size() > 3)) {
      error->all(FLERR, "Error: if takes exactly three positional arguments.\n");
    }
    return ifv(parsev(args[0]), parsev(args[1]), parsev(args[2]));
  }
  if (func == "log")
    return logv(parsev(arg));
  if (func == "dimension")
    return Var(dimension(args));
  if (func == "axisymmetric")
    return Var(axisymmetric(args));
  if (func == "region")
    return Var(region(args));
  if (func == "solid")
    return Var(solid(args));
  if (func == "eos")
    return Var(add_EOS(args));
  if (func == "strength")
    return Var(add_strength(args));
  if (func == "material")
    return Var(add_material(args));
  if (func == "damage")
    return Var(add_damage(args));
  if (func == "temperature")
    return Var(add_temperature(args));
  if (func == "dump")
    return Var(dump(args));
  if (func == "group")
    return Var(group_command(args));
  if (func == "set_output")
    return Var(set_output(args));
  if (func == "log_modify")
    return Var(log_modify(args));
  if (func == "method")
    return Var(method(args));
  if (func == "scheme")
    return Var(scheme(args));
  if (func == "fix")
    return Var(fix(args));
  if (func == "delete_fix")
    return Var(delete_fix(args));
  if (func == "compute")
    return Var(compute(args));
  if (func == "delete_compute")
    return Var(delete_compute(args));
  if (func == "dt_factor")
    return Var(set_dt_factor(args));
  if (func == "set_dt")
    return Var(set_dt(args));
  if (func == "value")
    return value(args);
  if (func == "plot")
    return Var(plot(args));
  if (func == "create_domain")
    return Var(create_domain(args));
  if (func == "save_plot")
    return Var(save_plot(args));

  // invoke commands added via style_command.h

  if (command_map->find(func) != command_map->end()) {
    CommandCreator command_creator = (*command_map)[func];
    return command_creator(mpm,args);
  }
  else if (func == "evaluate")
    return Var(parsev(arg).result(mpm));
  else if (func == "print")
    return Var(print(args));
  else if (func == "restart")
    return Var(restart(args));
  error->all(FLERR, "Error: Unknown function " + func + "\n");
  return Var();
}

// remove white spaces from string
string Input::remove_whitespace(string str){
  string str_;

  bool quote = false;

  for(int i=0; i<str.length(); i++){
    if (str[i] == '"')
      quote = !quote;
    else if (str[i] != ' ' || quote)
      str_.append(&str[i],1); // Add the non-whitespace character to str_
  }
  return str_;
}

float Input::parse(string str){
  error->all(FLERR, "Error: Input::parse deprecated function.\n");
  return nan("");
}

Expression &
Input::parsev(const string &name, float value)
{
    const pair<map<string, Expression>::iterator, bool> it = expressions.emplace(piecewise_construct, forward_as_tuple(name), tuple<>());
    
    Expression &expression = it.first->second;

    if (it.second)
      expression.operations.emplace_back(new ExpressionOperandLiteral(value));
    else
    {
      ExpressionOperandLiteral *literal;

      if (expression.operations.size() == 1)
        literal = dynamic_cast<ExpressionOperandLiteral *>(expression.operations.front().get());
      else
        literal = nullptr;
    
      if (literal)
        literal->value = value;
      else
        error->all(FLERR, name + " was not a literal expression.\n");
    }

    expression.registers = Kokkos::View<float**>("expression", 1, 1);

    return expression;
}

/*! This function translates the input file syntax, either mathematical expressions or functions, 
 * to C++ and evaluate or execute them.
 */
Var Input::parsev(string str)
{
  static int depth = -1;
  depth++;
  /////////////////////////////////////////// EXPRESSIONS ///////////////////////////////////////////////////////////
  if (!str.empty() && !depth)
  {
    smatch match;
    string name, expression_string;
    if (regex_search(str, match, regex("\\s*=\\s*")))
    {
      name = match.prefix();
      expression_string = match.suffix();
    }
    else
    {
      name = str;
      expression_string = str;
    }

    const pair<map<string, Expression>::iterator, bool> it = expressions.emplace(piecewise_construct, forward_as_tuple(name), tuple<>());

    if (it.second)
    {
      Expression &expression = it.first->second;

      deque<unique_ptr<Expression::Operation>> operator_stack;
  
      string current_token;

      enum class CharacterType: char
      {
        NUMBER,
        LETTER,
        OPEN_PARENTHESIS,
        CLOSE_PARENTHESIS,
        COMMA,
        SYMBOL,
        NONE
      } current_type = CharacterType::NONE;

      bool after_operand = false;

      for (string::const_iterator it = expression_string.cbegin(), end = expression_string.cend();; it++)
      {
        char current_character = it == end? ' ': *it;

        if (DEBUG_EXPRESSIONS)
          cout << "\"" << current_character << "\", ";

        CharacterType next_type = current_character == ' '?   CharacterType::NONE:
                                  current_character >= '0' &&
                                  current_character <= '9' ||
                                  current_character == '.'?   CharacterType::NUMBER:
                                  current_character >= 'a' &&
                                  current_character <= 'z' ||
                                  current_character >= 'A' &&
                                  current_character <= 'Z' ||
                                  current_character == '_'?   CharacterType::LETTER:
                                  current_character == '('?   CharacterType::OPEN_PARENTHESIS:
                                  current_character == ')'?   CharacterType::CLOSE_PARENTHESIS:
                                  current_character == ','?   CharacterType::COMMA:
                                                              CharacterType::SYMBOL;

        if (next_type == CharacterType::NUMBER && current_type == CharacterType::LETTER)
          next_type = CharacterType::LETTER;

        if (current_type == CharacterType::NUMBER && (current_character == 'e' ||
            current_token.back() == 'e' && current_character == '-'))
          next_type = CharacterType::NUMBER;

        if ((current_type != next_type || next_type == CharacterType::SYMBOL || next_type == CharacterType::OPEN_PARENTHESIS ||
             next_type == CharacterType::CLOSE_PARENTHESIS) &&
            (current_type != CharacterType::LETTER || next_type != CharacterType::OPEN_PARENTHESIS))
        {
          if (current_type == CharacterType::CLOSE_PARENTHESIS ||
              current_type == CharacterType::COMMA)
          {
            while (!operator_stack.back()->isFunction())
            {
              expression.operations.push_back(move(operator_stack.back()));
              operator_stack.pop_back();
            }

            if (after_operand = current_type == CharacterType::CLOSE_PARENTHESIS)
            {
              if (dynamic_cast<const ExpressionFunctionEvaluate *>(operator_stack.back().get()))
              {
                if (const ExpressionOperandExpression *operand = dynamic_cast<const ExpressionOperandExpression *>
                                                                 (expression.operations.back().get()))
                {
                  double evaluation = expressions.at(operand->name).getConstant();
                  expression.operations.pop_back();
                  expression.operations.emplace_back(new ExpressionOperandLiteral(evaluation));

                  cout << "EVALUATE: " << evaluation << endl;
                }
                else
                  error->all(FLERR, "Argument to evaluate was not an Expression");
              }
              else
              {
                expression.operations.push_back(move(operator_stack.back()));
                operator_stack.pop_back();
              }

              after_operand = true;
            }

            if (DEBUG_EXPRESSIONS)
              cout << endl << "CLOSE_PARENTHESIS OR COMMA" << endl;
          }
          else if (current_type != CharacterType::NONE)
          {
            Expression::Operation *new_operation = current_type == CharacterType::NUMBER? new ExpressionOperandLiteral(stod(current_token)): nullptr;
            if (current_type == CharacterType::NUMBER && DEBUG_EXPRESSIONS)
              cout << endl << "NUMBER " << stod(current_token) << endl;
            
            if (!new_operation)
            {
              const map<string,Expression>::iterator &it = current_token == name? expressions.end(): expressions.find(current_token);

              if (it == expressions.end())
              {
                if (DEBUG_EXPRESSIONS)
                  cout << endl << "(\"" << current_token << "\", " << (after_operand? "infix": "unary") << "), " << endl;
                new_operation = (after_operand? infix_factory: operation_factory).new_instance(current_token);
              }
              else
              {
                if (DEBUG_EXPRESSIONS)
                  cout << endl << "EXPRESSION FROM EXPRESSION \"" << current_token << "\"" << endl;
                expression.expression_dependencies.insert(&it->second);
                new_operation = new ExpressionOperandExpression(current_token);
              }
            }

            if (!new_operation)
            {
              cout << current_token << " does not appear to be an operation. Aborting." << endl;
              expressions.erase(name);
              goto end_of_expressions;
            }

            if (after_operand = new_operation->isOperand())
              expression.operations.emplace_back(new_operation);
            else
            {
              while (!operator_stack.empty() && operator_stack.back()->precedence() <= new_operation->precedence() &&
                     !operator_stack.back()->isFunction())
              {
                expression.operations.push_back(move(operator_stack.back()));
                operator_stack.pop_back();
              }
              operator_stack.emplace_back(new_operation);
            }
          }

          current_token.clear();
        }

        if (it == end)
          break;

        current_token += current_character;
        current_type = next_type;
      }

      while (!operator_stack.empty())
      {
        expression.operations.push_back(move(operator_stack.back()));
        operator_stack.pop_back();
      }

      int current_index = 0, max_index = 0;

      for (const unique_ptr<Expression::Operation> &operation: expression.operations)
      {
        max_index = max(max_index, current_index -= operation->arity() - 1);
      }

      if (!max_index)
        cout << "MAX INDEX CANNOT BE ZERO" << endl;
        
      expression.registers = Kokkos::View<float**>("expression", max_index, 1);
    }
  }

  end_of_expressions:

  /////////////////////////////////////////////// VARS //////////////////////////////////////////////////////////////

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
	      if (universe->me == 0) {
		cout << returnvar << " = " << (*vars)[returnvar].result() << endl;
	      } else {
		(*vars)[returnvar].result();
	      }
	    }
      depth--;
	    return -(*vars)[word];
	  }
	  else {
	    if (!returnvar.empty()) {
	      (*vars)[returnvar] = (*vars)[word];
	      if (universe->me == 0) {
		cout << returnvar << " = " << (*vars)[returnvar].result() << endl;
	      } else {
		(*vars)[returnvar].result();
	      }
	    }
      depth--;
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
	    error->all(FLERR, "Error: Unbalanced parenthesis '('.\n");
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
    depth--;
    return Var(-1);
  }
  else {
    if (!returnvar.empty()) {
      (*vars)[returnvar] = values.top();
      if (universe->me == 0) {
	cout << returnvar << " = " << (*vars)[returnvar].result(mpm) << endl;
      } else {
	(*vars)[returnvar].result(mpm);
      }
    }
    depth--;
    return values.top();
  }
}

int Input::dimension(vector<string> args) {
  // Check that a method is available:
  if (update->method == nullptr) {
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
  if (universe->me == 0) {
    cout << args[0] << " = {equation=\"" << v.eq()
	 << "\", value=" << v.result(mpm) << ", constant=";

    if (v.is_constant())
      cout << "true";
    else
      cout << "false";

    cout << "}\n";
  } else {
    v.result(mpm);
  }

  return 0;
}


int Input::restart(vector<string> args) {
  output->create_restart(args);
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
    if (variable == protected_vars[i])
      return 1;
  }
  return 0;
}

int Input::create_domain(vector<string> args) {
  domain->create_domain(args);
  return 0;
}
