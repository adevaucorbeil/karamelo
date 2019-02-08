#include "mpm.h"
#include "grid.h"
#include "material.h"
#include "input.h"
#include <vector>
#include "memory.h"

using namespace std;


Grid::Grid(MPM *mpm) :
  Pointers(mpm)
{
  cout << "Creating new grid" << endl;

  x= NULL;
  v = v_update = NULL;
  b = f = NULL;
}

Grid::~Grid()
{
  memory->destroy(x);
  memory->destroy(v);
  memory->destroy(v_update);
  memory->destroy(b);
  memory->destroy(f);
}

void Grid::init(string cs){
  cellsize = input->parse(cs);
  cout << "Set grid cellsize to " << cellsize << endl;
}

void Grid::grow(int nn){
  nnodes = nn;

  x = memory->grow(x, nn, 3, "grid:x");
  v = memory->grow(x, nn, 3, "grid:v");
  v_update = memory->grow(x, nn, 3, "grid:v_update");
  b = memory->grow(x, nn, 3, "grid:b");
  f = memory->grow(x, nn, 3, "grid:f");
}
