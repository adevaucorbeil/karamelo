#pragma once

#include <Kokkos_Core.hpp>

#include <iostream>
#include <chrono>

using namespace std;
using namespace chrono;

void kokkos_test() {
  steady_clock::time_point t0 = high_resolution_clock::now();

  int result;
  Kokkos::parallel_reduce("test", 10000000, KOKKOS_LAMBDA(int64_t i, int& result) {
    result += 305175781 % (i + 1);
    result %= 12207031;
  }, result);

  cout << result << " calculated in " << duration<double, milli>(high_resolution_clock::now() - t0).count() << "ms" << endl;
}