#include "InverseSolvers.h"
using namespace NTPoly;

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "InverseSolvers_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
////////////////////////////////////////////////////////////////////////////////
void InverseSolvers::Invert(
    const DistributedSparseMatrix &Overlap, DistributedSparseMatrix &InverseMat,
    const IterativeSolverParameters &solver_parameters) {
  Invert_wrp(GetIH(Overlap), GetIH(InverseMat), GetIH(solver_parameters));
}
} // namespace NTPoly