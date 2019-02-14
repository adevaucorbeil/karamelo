#include <iostream>
#include <string>
#include <string.h>
#include <map>
#include <math.h>
#include <stack>
#include "variable.h"
#include <stdio.h>
#include <stdlib.h>
#include <vector>

using namespace std;


Variable::Variable(map<string, Variable> * ptr, double v, bool c)
{
  known_var = ptr;
  equation = to_string(v);
  value = v;
  constant = c;
}

Variable::Variable(map<string, Variable> * ptr, string eq, double v, bool c)
{
  known_var = ptr;
  equation = eq;
  value = v;
  constant = c;
}

bool Variable::is_constant() const {
  return constant;
}

void Variable::evaluate()
{
  if (constant) return;

  cout << "Evaluate equation: " <<  equation << endl;
  value = parse(known_var, equation).value;
  if (constant) equation = to_string(value);
}

double Variable::result()
{
  if (constant) return value;
  else {
    evaluate();
    return value;
  }
}

double Variable::result() const
{
  return value;
}

string Variable::str() const
{
  if (equation != "") return equation;
  else return to_string(value);
}

string Variable::eq() const
{
  return equation;
}



Variable Variable::operator+(const Variable& right)
{
  if (this->known_var != right.known_var) {
    cout << "Error: this->known_var != right.known_var, I don't know how to deal with this" << endl;
    exit(1);
  }
  if (this->constant && right.constant) {
    Variable result(this->known_var, this->value + right.value);
    return result;
  } else {
    Variable result(this->known_var, "(" + this->str() + "+" + right.str() + ")", this->value + right.value, this->constant && right.constant);
    return result;
  }
}

Variable Variable::operator-(const Variable& right)
{
  if (this->known_var != right.known_var) {
    cout << "Error: this->known_var != right.known_var, I don't know how to deal with this" << endl;
    exit(1);
  }

  if (this->constant && right.constant) {
    Variable result(this->known_var, this->value - right.value);
    return result;
  } else {
    Variable result(this->known_var, "(" + this->str() + "-" + right.str() + ")", this->value - right.value, this->constant && right.constant);
    return result;
  }
}

Variable Variable::operator-()
{
  if (this->constant) {
    Variable result(this->known_var, -this->value);
    return result;
  } else {
    Variable result(this->known_var, "(-" + this->str() + ")", -this->value, this->constant);
    return result;
  }
}

Variable Variable::operator*(const Variable& right)
{
  if (this->known_var != right.known_var) {
    cout << "Error: this->known_var != right.known_var, I don't know how to deal with this" << endl;
    exit(1);
  }

  if (this->constant && right.constant) {
    Variable result(this->known_var, this->value * right.value);
    return result;
  } else {
    Variable result(this->known_var, "(" + this->str() + "*" + right.str() + ")", this->value * right.value, this->constant && right.constant);
    return result;
  }
}


Variable Variable::operator/(const Variable& right)
{
  if (this->known_var != right.known_var) {
    cout << "Error: this->known_var != right.known_var, I don't know how to deal with this" << endl;
    exit(1);
  }

  if (this->constant && right.constant) {
    Variable result(this->known_var, this->value / right.value);
    return result;
  } else {
    Variable result(this->known_var, "(" + this->str() + "/" + right.str() + ")", this->value / right.value, this->constant && right.constant);
    return result;
  }
}


Variable Variable::operator^(const Variable& right)
{
  if (this->known_var != right.known_var) {
    cout << "Error: this->known_var != right.known_var, I don't know how to deal with this" << endl;
    exit(1);
  }

  if (this->constant && right.constant) {
    Variable result(this->known_var, pow(this->value, right.value));
    return result;
  } else {
    Variable result(this->known_var, "(" + this->str() + "^" + right.str() + ")", pow(this->value, right.value), this->constant && right.constant);
    return result;
  }
}

Variable operator*(int left, Variable right){
  if (right.is_constant()) {
    Variable result(right.known_var, left*right.result());
    return result;
  } else {
    Variable result(right.known_var, "(" + to_string(left) + "*" + right.str() + ")", left*right.result(), right.is_constant());
    return result;
  }
}

Variable operator*(Variable left, int right){
  if (left.is_constant()) {
    Variable result(left.known_var, left.result()*right);
    return result;
  } else {
    Variable result(left.known_var, "(" + left.str() + "*" + to_string(right) + ")", left.result()*right, left.is_constant());
    return result;
  }
}

Variable powv(int base, Variable p){
  if (p.is_constant()) {
    Variable result(p.known_var, pow(base, p.result()));
    return result;
  } else {
    Variable result(p.known_var, "pow(" + to_string(base) + "," + p.str() + ")", pow(base, p.result()), p.is_constant());
    return result;
  }
}

Variable powv(Variable base, Variable p){
  if (base.is_constant() && p.is_constant()) {
    Variable result(base.known_var, pow(base.result(), p.result()));
    return result;
  } else {
    Variable result(base.known_var, "pow(" + base.str() + "," + p.str() + ")", pow(base.result(), p.result()), base.is_constant() && p.is_constant());
    return result;
  }
}

Variable expv(Variable x){
  if (x.is_constant()) {
    Variable result(x.known_var, exp(x.result()));
    return result;
  } else {
    Variable result(x.known_var, "exp(" + x.str() + ")", exp(x.result()), x.is_constant());
    return result;
  }
}

// Function to find precedence of  
// operators. 
double precedence(char op){ 
    if(op == '+'||op == '-') return 1;
    if(op == '*'||op == '/') return 2;
    if(op == '^') return 3;
    if(op == 'e'|| op == 'E') return 4;
    return 0;
} 
  
// Function to perform arithmetic operations. 
double applyOp(double a, double b, char op){ 
  switch(op){ 
  case '+': return a + b;
  case '-': return a - b;
  case '*': return a * b;
  case '/': return a / b;
  case '^': return pow(a,b);
  case 'e': return a*pow((int) 10,b);
  case 'E': return a*pow((int) 10,b);
  case '(':
    printf("Error: unmatched parenthesis (\n");
    exit(1);
  default:
    printf("Error: unknown operator %c\n", op);
    exit(1);
  } 
}

// Function to perform arithmetic operations. 
Variable applyOp(Variable a, Variable b, char op){
  switch(op){ 
  case '+': return a + b;
  case '-': return a - b;
  case '*': return a * b;
  case '/': return a / b;
  case '^': return powv(a,b);
  case 'e': return powv(10,b);
  case 'E': return powv(10,b);
  case '(':
    printf("Error: unmatched parenthesis (\n");
    exit(1);
  default:
    printf("Error: unknown operator %c\n", op);
    exit(1);
  } 
}

bool is_operator(char op){
  if (op=='+') return true;
  if (op=='-') return true;
  if (op=='*') return true;
  if (op=='/') return true;
  if (op=='^') return true;
  return false;
}

// check if op is either of +-/*()
bool is_math_char(char op){
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
Variable evaluate_function(map<string, Variable> *known_var, string func, string arg){
  //cout << "Evaluate function " << func << " with argument: " << arg << endl;

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

  // if (func.compare("dimension") == 0) return dimension(known_var, args);
  // if (func.compare("region") == 0) return (double) region(args);
  // if (func.compare("solid") == 0) return (double) solid(args);
  // if (func.compare("eos") == 0) return (double) EOS(args);
  // if (func.compare("dump") == 0) return (double) dump(args);
  // if (func.compare("group") == 0) return (double) group_command(args);
  // if (func.compare("log") == 0) return (double) log(args);
  // if (func.compare("method_modify") == 0) return (double) method_modify(args);
  // if (func.compare("fix") == 0) return (double) fix(args);

  // invoke commands added via style_command.h

  // if (command_map->find(func) != command_map->end()) {
  //   CommandCreator command_creator = (*command_map)[func];
  //   command_creator(mpm,args);
  //   return 0;
  // }



  //else
    if (func.compare("exp") == 0) return expv(parse(known_var, arg));
  cout << "Error: Unknown function " << func << endl;
  exit(1);
}

// remove white spaces from string
string remove_whitespace(string str){
  string str_;

  for(int i=0; i<str.length(); i++){
    if( str[i] != ' ') str_.append(&str[i],1); // Add the non-whitespace character to str_
  }
  return str_;
}


Variable parse(map<string, Variable> *known_var, string str)
{
  // stack to store integer values. 
  stack <Variable> values; 
      
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
	Variable result(known_var, (-1)*stof(number));
      	values.push(result);
      	negative = false;
      } else {
	Variable result(known_var, stof(number));
      	values.push(result);
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
	  Variable val2 = values.top();
	  values.pop();

	  Variable val1 = values.top();
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
	Variable val2 = values.top();
	values.pop();

	Variable val1 = values.top();
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

      	  if (str[i+1] == '+') {
      	    i++;
      	  }
      	  continue;
      	}
      }

      // Check if word is a variable:
      map<string, Variable>::iterator it;
      it = known_var->find(word);

      if (it != known_var->end()){
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
	    Variable result = -(*known_var)[word];
	    return result;
	  }
	  else {
	    return (*known_var)[word];
	  }
	}
	
	else {
	  if (negative) {
	    values.push((-1)*(*known_var)[word]);
	    negative = false;
	  } else {
	    values.push((*known_var)[word]);
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
	values.push(evaluate_function(known_var, word, arg));
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
    Variable val2 = values.top();
    values.pop();

    Variable val1 = values.top();
    values.pop();

    char op = ops.top();
    ops.pop();
                  
    values.push(applyOp(val1, val2, op)); 
  }

  // Top of 'values' contains result, return it.
  if (values.empty()) {
    return Variable(known_var, -1);
  }
  else {
    //cout << "value = " << values.top() << endl;
    if (!returnvar.empty()) {
      (*known_var)[returnvar] = values.top();
    }
    return values.top();
  }
}

