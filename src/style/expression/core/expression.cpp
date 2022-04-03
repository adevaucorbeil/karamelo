#include <expression_operand_constant.h>
#include <expression_operand_index.h>
#include <expression_operators_arithmetic.h>
#include <expression_function_parenthesis.h>
#include <expression_function_trigonometry.h>
#include <style_expression.h>

#include <cstdlib>

StyleFactory<Expression::Operation> Expression::operation_factory;

void Expression::initialize()
{
  operation_factory.register_class<ExpressionSum       >("+");
  operation_factory.register_class<ExpressionDifference>("-");
  operation_factory.register_class<ExpressionProduct   >("*");
  operation_factory.register_class<ExpressionQuotient  >("/");

  operation_factory.register_class<ExpressionOperandIndex>("i");

  operation_factory.register_class<ExpressionFunctionParenthesis>("(");
  operation_factory.register_class<ExpressionFunctionSin>("sin(");
  operation_factory.register_class<ExpressionFunctionCos>("cos(");
  operation_factory.register_class<ExpressionFunctionTan>("tan(");
  operation_factory.register_class<ExpressionFunctionAsin>("asin(");
  operation_factory.register_class<ExpressionFunctionAcos>("acos(");
  operation_factory.register_class<ExpressionFunctionAtan>("atan(");
  operation_factory.register_class<ExpressionFunctionAtan2>("atan2(");

#define EXPRESSION_CLASS
#define ExpressionStyle(key, Class) operation_factory.register_class<Class>(#key"(");
#include <style_expression.h>
#undef ExpressionStyle
#undef EXPRESSION_CLASS
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
          if (!operator_stack.back()->isParenthesis())
          {
            operations.push_back(std::move(operator_stack.back()));
          }
          operator_stack.pop_back();
        }
      }
      else if (current_type != CharacterType::NONE)
      {
        Operation *new_operation = current_type == CharacterType::NUMBER? new ExpressionOperandConstant(std::stod(current_token)):
                                                                          operation_factory.new_instance(current_token);

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
    if (!operator_stack.back()->isParenthesis())
    {
      operations.push_back(std::move(operator_stack.back()));
    }
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
