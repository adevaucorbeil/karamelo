#pragma once

#include <expression_operation.h>

template<typename DERIVED, int PRECEDENCE>
using ExpressionOperatorBinary = ExpressionOperation<DERIVED, false, false, PRECEDENCE, 2>;
