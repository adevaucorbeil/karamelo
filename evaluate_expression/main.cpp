#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <bits/stdc++.h> 

using namespace std; 

map<string, double> variables;

// Function to find precedence of  
// operators. 
double precedence(char op){ 
    if(op == '+'||op == '-') return 1;
    if(op == '*'||op == '/') return 2;
    if(op == '^') return 3;
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
    case '(':
      printf("Error: unmatched parenthesis (\n");
      exit(1);
    default:
      printf("Error: unknown operator %c\n", op);
      exit(1);
    } 
}

void extract_lhs(char *str, char *&lhs, char *&rhs)
{
  size_t size;
  size_t space_pos;

  char operators[] = "+-*/()=^";
  
  size =  strcspn(str, operators);
  //if (size == strlen(str)) return double()
  lhs = NULL;
  if (size > 0) lhs = (char*) malloc(size+1);
  rhs = strpbrk(str, operators);
  if (size > 0) {
    strncpy(lhs,str,size);
    space_pos = strcspn(lhs," ");
    if (space_pos < size) lhs[(int) space_pos] = '\0';
    else lhs[(int) size] = '\0'; // Manually add null character at the end of the string
  }
  //printf("extract_lhs: str=%s, lhs=%s, rhs=%s\n", str, lhs, rhs);
}

bool is_operator(char op){
  if (op=='+') return true;
  if (op=='-') return true;
  if (op=='*') return true;
  if (op=='/') return true;
  if (op=='^') return true;
  return false;
}

double evaluate_function(char *func, double value){
  if (strcmp(func, "exp") == 0) return (double) exp(value);
  printf("Error: Unknown function %s\n", func);
  exit(1);
}


char* remove_whitespace(char *&str){
  size_t space_pos = strcspn(str," ");
  size_t size_str = strlen(str);

  char * str_;
  str_ = (char*) malloc(strlen(str)+1);
  int k = 0;

  for(int i=0; i<size_str; i++) {
    if (str[i] != ' ') {
      str_[k] = str[i];
      k++;
    }
  }
  str_[k++] = '\0';

  return str_;
}


double evaluate(char *str)
{
  // stack to store integer values. 
  stack <double> values; 
      
  // stack to store operators. 
  stack <char> ops;

  // stack to store functions.
  stack <char*> funcs;

  string returnvar;
  
  char *e, *lhs, *rhs;
  e = NULL;
  lhs = NULL;
  rhs = remove_whitespace(str);

  bool negative = false; // flag that indicates that the next value will have to be multiplied by -1

  while (rhs!=NULL){
    if(isdigit(rhs[0])){
      e = rhs;
      rhs = NULL;
      extract_lhs(e, lhs, rhs);

      if (negative) {
	values.push((-1)*atof(lhs));
	negative = false;
      }

      else values.push(atof(lhs));
      //printf("push %f to values\n", values.top());
    }

    // If the first character in the righ-hand-side is (, push it to ops
    else if(rhs[0] == '('){
      //printf("found (\n");
      ops.push(rhs[0]);
      rhs++;
      if (rhs[0] == '-') {
	negative = true;
	rhs++;
      }
      //printf("( lhs=%s, rhs=%s\n", lhs, rhs);
    }

    // Closing brace encountered, solve  
    // entire brace. 
    else if(rhs[0] == ')') {
      //printf("found )\n");
      while(ops.top() != '(') 
	{
	  double val2 = values.top();
	  values.pop();

	  double val1 = values.top();
	  values.pop();

	  char op = ops.top();
	  ops.pop();
	  //printf("%f %f %c\n", val1, val2, op);

	  if (ops.empty()) {
	    printf("Error, unmatched parenthesis )\n");
	    exit(1);
	  }
	  values.push(applyOp(val1, val2, op));
	}

      // pop opening brace. 
      ops.pop();

      // evaluate function:
      if (!funcs.empty()) {
	double value_func = evaluate_function(funcs.top(), values.top());
	//printf("%s of %f = %f\n", funcs.top(), values.top(), value_func);
	funcs.pop();
	values.pop();
	values.push(value_func);
      }

      rhs++;
      if (rhs[0] == '\0') rhs = NULL;
      //printf(") lhs=%s, rhs=%s\n", lhs, rhs);
    }

    else if (is_operator(rhs[0])){
      char new_op = rhs[0];
      //printf("found operator %c\n", new_op);

      if (values.empty()) {
	if (new_op == '-') {
	  negative = true;
	  rhs++;
	  continue;
	}
      }


      if (rhs[1] == '\0') {
	printf("Error: end-of-line character detected after operator %c\n", new_op);
	exit(1);
      }

      else if (rhs[1] == '*') {
	new_op = '^';
	rhs++;
      }

      else if (is_operator(rhs[1]) && rhs[1] != '-') {
	printf("Error: unknown operator sequence %c%c\n", new_op, rhs[1]);
	exit(1);
      }


      // While top of 'ops' has same or greater
      // precedence to current token, which
      // is an operator. Apply operator on top
      // of 'ops' to top two elements in values stack.
      while(!ops.empty() &&
	    precedence(ops.top()) >= precedence(new_op)){
	double val2 = values.top();
	values.pop();

	double val1 = values.top();
	values.pop();

	char op = ops.top();
	ops.pop();
	//printf("%f %f %c\n", val1, val2, op);

	values.push(applyOp(val1, val2, op));
      } 
              
      // Push current token to 'ops'. 
      ops.push(new_op);
      rhs++;
      if (rhs[0] == '-') {
	negative = true;
	rhs++;
      }
    }

    else {
      e = rhs;
      rhs = NULL;
      extract_lhs(e, lhs, rhs);

      //printf("else lhs=%s, rhs=%s\n", lhs, rhs);

      // Check if lhs is a variable:
      map<string, double>::iterator it;
      it = variables.find((string) lhs);

      if (it != variables.end()){
	// lhs is a variable
	if (rhs != NULL && rhs[0] == '=') {

	  if (!values.empty() || !ops.empty() ) {
	    printf("Error: I do not understand when '=' is located in the middle of an expression\n");
	    exit(1);
	  }

	  else {
	    returnvar = string(lhs);
	    rhs++;
	  }
	}

	else if (rhs == NULL && values.empty() && ops.empty() ) {
	  if (negative) {
	    return (-1)*variables[(string) lhs];
	  }
	  else {
	    return variables[(string) lhs];
	  }
	}
	
	else {
	  if (negative) {
	    values.push((-1)*variables[(string) lhs]);
	    negative = false;
	  } else {
	    values.push(variables[(string) lhs]);
	  }
	  //printf("push %s=%f to values\n", lhs, values.top());
	}
      }

      else if (rhs[0]=='=' && values.empty() && ops.empty()){
	// Check if there is an '=':
	  returnvar = string(lhs);
	  rhs++;
	  //cout << "The computed value will be stored in " <<  returnvar << endl;
	}

      else if (rhs[0]=='(') {
	// It's a function:
	//printf("Found function %s\n", lhs);
	funcs.push(lhs);
      }

      else {
	printf("%s is unknown\n", lhs);
	exit(1);
      }
    }
  }

  // Entire expression has been parsed at this
  // point, apply remaining ops to remaining
  // values.
  while(!ops.empty()){
    double val2 = values.top();
    values.pop();

    double val1 = values.top();
    values.pop();

    char op = ops.top();
    ops.pop();
                  
    values.push(applyOp(val1, val2, op)); 
  } 

  free(lhs);
  free(rhs);
  // Top of 'values' contains result, return it.
  printf("value=%f\n", values.top());
  if (!returnvar.empty()) {
    variables[returnvar] = values.top();
  }
  return values.top(); 
}


int main(){
  variables["ab"] = 10;
  variables["p"] = 2.5;
  
  double result;
  char *expr = "x = 2 + (2*ab^2) + 10";

  // cout << "ab=" << variables["ab"] << endl;
  // printf("%s = %f\n", expr, evaluate(expr));
  // //printf("%s = %f\n", "x+2", evaluate("x+2"));
  // printf("%f\n", evaluate("1*-2"));
  // printf("%f\n", evaluate("1--2"));
  printf("%f\n", evaluate("1 2"));
  printf("%f\n", evaluate("p = -p+(-2+2*(3-4))**2"));
  cout << "p=" << variables["p"] << endl;
  printf("%f\n", evaluate("2*p"));
  cout << "p=" << variables["p"] << endl;
  printf("%f\n", evaluate("p = 2"));
  cout << "p=" << variables["p"] << endl;
}


