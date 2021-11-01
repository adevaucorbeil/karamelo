#include <Kokkos_Core.hpp>

#include <mls_mpm.h>

int main(int argc,
         char *argv[])
{
  Kokkos::initialize(argc, argv);

  for (int i = 0; i < 40; i++)
    mls_mpm(i + 1);

  Kokkos::finalize();
}
