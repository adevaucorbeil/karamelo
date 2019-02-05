#include "mpm.h"
#include "grid.h"
#include "material.h"
#include <vector>

using namespace std;


Grid::Grid(MPM *mpm) :
  Pointers(mpm)
{
  cout << "Creating new grid" << endl;
}

Grid::~Grid()
{
}

void Grid::grow(int nn){
  nnodes = nn;
  x.reserve(nnodes);
}
