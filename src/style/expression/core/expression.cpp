#include <expression_operation.h>

#include <solid.h>

Expression::~Expression()
{}

void Expression::evaluate()
{
  for (Expression *expression_dependency: expression_dependencies)
    expression_dependency->evaluate();

  index = 0;
  cout << "registers.extent(0)=" << registers.extent(0) << endl;
  registers = Kokkos::View<double**>("expression", 1, 1);

  for (const std::unique_ptr<Operation> &operation: operations)
    operation->apply(this);
}

void Expression::evaluate(Solid &solid)
{
  for (Expression *expression_dependency: expression_dependencies)
    expression_dependency->evaluate(solid);

  index = 0;
  registers = Kokkos::View<double**>("expression", registers.extent(0), solid.np_local);

  for (const std::unique_ptr<Operation> &operation: operations)
  {
    operation->set_solid(solid);
    operation->apply(this);
  }
}

void Expression::evaluate(Grid &grid)
{
  for (Expression *expression_dependency: expression_dependencies)
    expression_dependency->evaluate(grid);

  index = 0;
  registers = Kokkos::View<double**>("expression", registers.extent(0), grid.nnodes);

  for (const std::unique_ptr<Operation> &operation: operations)
  {
    operation->set_grid(grid);
    operation->apply(this);
  }
}
