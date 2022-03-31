#pragma once

#include <expression_operation.h>

template<typename DERIVED, int ARITY>
using ExpressionFunction = ExpressionOperation<DERIVED, false, true, 0, ARITY>;
