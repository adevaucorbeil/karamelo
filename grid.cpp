#include "mpm.h"
#include "grid.h"
#include "material.h"
#include "input.h"
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

void Grid::init(string cs){
  cellsize = input->parse(cs);
  cout << "Set grid cellsize to " << cellsize << endl;
}

void Grid::grow(int nn){
  nnodes = nn;
  x.reserve(nnodes);
}
