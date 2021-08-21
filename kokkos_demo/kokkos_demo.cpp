#include <Kokkos_Core.hpp>

#include <iostream>
#include <chrono>

using namespace std;
using namespace chrono;
using namespace Kokkos;

int main(int argc, char *argv[])
{
    initialize(argc, argv);

    steady_clock::time_point t0 = high_resolution_clock::now();

    int result;
    parallel_reduce("test", 10000000, [=](int64_t i, int &result) {
        result += i;
        result %= 100;
    }, result);

    cout << result << " calculated in " << duration<double, std::milli>(high_resolution_clock::now() - t0).count() << "ms" << endl;

    finalize();
}
