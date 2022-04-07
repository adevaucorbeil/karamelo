#include <expression_function_parenthesis.h>
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

#include <expression_operators_arithmetic.h>
#include <style_expression.h>

#include <cstdlib>

StyleFactory<Expression::Operation> Expression::operation_factory;
std::map<std::string, std::unique_ptr<Expression>> Expression::named_expressions;

void Expression::initialize()
{
  operation_factory.register_class<ExpressionFunctionParenthesis>("(");

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

  operation_factory.register_class<ExpressionSum       >("+");
  operation_factory.register_class<ExpressionDifference>("-");
  operation_factory.register_class<ExpressionProduct   >("*");
  operation_factory.register_class<ExpressionQuotient  >("/");

  operation_factory.register_class<ExpressionOperandConstant<Kokkos::Experimental::    pi_v<double>>>("PI"    );
  operation_factory.register_class<ExpressionOperandConstant<Kokkos::Experimental::     e_v<double>>>("E"     );
  operation_factory.register_class<ExpressionOperandConstant<Kokkos::Experimental::egamma_v<double>>>("EGAMMA");
  operation_factory.register_class<ExpressionOperandConstant<Kokkos::Experimental::   phi_v<double>>>("PHI"   );

  operation_factory.register_class<ExpressionOperandIndex>("i");

  operation_factory.register_class<ExpressionOperandVector<&Solid::x, &Grid::x, 0>>("x");
  operation_factory.register_class<ExpressionOperandVector<&Solid::x, &Grid::x, 1>>("y");
  operation_factory.register_class<ExpressionOperandVector<&Solid::x, &Grid::x, 2>>("z");

#define EXPRESSION_CLASS
#define ExpressionStyle(key, Class) operation_factory.register_class<Class>(#key"(");
#include <style_expression.h>
#undef ExpressionStyle
#undef EXPRESSION_CLASS
}

void Expression::finalize()
{
  named_expressions.clear();
}

Expression &
Expression::make_named_expression(const std::string &name,
                                  const std::string &expression)
{
    return *named_expressions.try_emplace(name, new Expression(expression)).first->second;
}

Expression::Expression(const std::string &expression)
{
  std::deque<std::unique_ptr<Operation>> operator_stack;
  
  std::string current_token;

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

  for (std::string::const_iterator it = expression.cbegin(), end = expression.cend();; it++)
  {
    char current_character = it == end? ' ': *it;

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

    if ((current_type != next_type || next_type == CharacterType::OPEN_PARENTHESIS ) &&
        (current_type != CharacterType::LETTER || next_type != CharacterType::OPEN_PARENTHESIS))
    {
      if (current_type == CharacterType::CLOSE_PARENTHESIS ||
          current_type == CharacterType::COMMA)
      {
        while (!operator_stack.back()->isFunction())
        {
          operations.push_back(std::move(operator_stack.back()));
          operator_stack.pop_back();
        }

        if (current_type == CharacterType::CLOSE_PARENTHESIS)
        {
          operations.push_back(std::move(operator_stack.back()));
          operator_stack.pop_back();
        }
      }
      else if (current_type != CharacterType::NONE)
      {
        Operation *new_operation = current_type == CharacterType::NUMBER? new ExpressionOperandLiteral(std::stod(current_token)): nullptr;
        if (!new_operation)
        {
          const std::map<std::string, std::unique_ptr<Expression>>::const_iterator &it = named_expressions.find(current_token);

          new_operation = it == named_expressions.cend()? operation_factory.new_instance(current_token): new ExpressionOperandExpression(*it->second);
        }


        if (new_operation->isOperand())
          operations.emplace_back(new_operation);
        else
        {
          while (!operator_stack.empty() && operator_stack.back()->precedence() < new_operation->precedence() &&
                 !operator_stack.back()->isFunction())
          {
            operations.push_back(std::move(operator_stack.back()));
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
    operations.push_back(std::move(operator_stack.back()));
    operator_stack.pop_back();
  }

  registers = Kokkos::View<double**>("expression", 10, 20);

  for (const std::unique_ptr<Operation> &operation: operations)
    operation->setExpression(this);
}

Expression::~Expression()
{}

void Expression::evaluate() const
{
  for (const std::unique_ptr<Operation> &operation: operations)
  {
    operation->apply();
  }
}

#include <iostream>

void Expression::print() const
{
  Kokkos::View<double**>::HostMirror mirror = Kokkos::create_mirror(registers);
  Kokkos::deep_copy(mirror, registers);
  for (int i = 0; i < registers.extent(1); i++)
  {
    std::cout << i << ": " << mirror(0, i) << std::endl;
  }
}
