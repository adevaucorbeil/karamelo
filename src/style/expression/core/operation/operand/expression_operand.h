#pragma once

#include <expression_operation.h>

template<typename DERIVED>
using ExpressionOperand = ExpressionOperation<DERIVED, true, false, 0, 0>;
