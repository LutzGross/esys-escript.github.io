#include <cusp/dia_matrix.h>
#include <cusp/monitor.h>
#include <cusp/krylov/lsqr.h>
#include <cusp/gallery/poisson.h>

// where to perform the computation
//typedef cusp::device_memory MemorySpace;
typedef cusp::host_memory MemorySpace;

// which floating point type to use
typedef double ValueType;

int main(void)
{
    // create an empty sparse matrix structure (DIA format)
    cusp::dia_matrix<int, ValueType, MemorySpace> A;

    // initialize matrix
    cusp::gallery::poisson9pt(A, 100, 100);

    // allocate storage for solution (x) and right hand side (b)
    cusp::array1d<ValueType, MemorySpace> x(A.num_rows, 0);
    cusp::array1d<ValueType, MemorySpace> b(A.num_rows, 1);

    // set stopping criteria:
    //  iteration_limit    = 100
    //  relative_tolerance = 1e-3
    cusp::verbose_monitor<ValueType> monitor(b, 100, 1.e-3);

    // solve the linear system A x = b
    cusp::krylov::lsqr(A, x, b, cusp::krylov::lsqr_parameters<ValueType>(), monitor);
    return 0;
}
