#include <Kokkos_Core.hpp>

#include <mls_mpm.h>

int main(int argc,
         char *argv[])
{
  Kokkos::initialize(argc, argv);

  mls_mpm<1000>();

  Kokkos::finalize();
}
