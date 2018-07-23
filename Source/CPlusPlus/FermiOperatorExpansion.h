#ifndef FERMIOPERATOREXPANSION_h
#define FERMIOPERATOREXPANSION_h

#include "SolverBase.h"

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class SolverParameters;
class Matrix_ps;
//! A Class For Solving Chemistry Systems Based On Sparse Matrices.
class FermiOperatorExpansion : public SolverBase {
public:
  //! Compute the density matrix using the fermi operator expansion.
  //! This subroutine is called if you already know the chemical potential.
  //!\param Hamiltonian the matrix to compute the corresponding density from.
  //!\param InverseSquareRoot of the overlap matrix.
  //!\param nel the number of electrons.
  //!\param degree the degree of the polynomial to use.
  //!\param chemical_potential the chemical potential.
  //!\param Density the density matrix computed by this routine.
  //!\param solver_parameters parameters for the solver
  static void Compute(const Matrix_ps &Hamiltonian,
                      const Matrix_ps &InverseSquareRoot, int nel,
                      Matrix_ps &Density, int degree, double chemical_potential,
                      const SolverParameters &solver_parameters);
  //! Compute the density matrix using the fermi operator expansion.
  //! This subroutine is called if you don't know the chemical potential.
  //!\param Hamiltonian the matrix to compute the corresponding density from.
  //!\param InverseSquareRoot of the overlap matrix.
  //!\param nel the number of electrons.
  //!\param degree the degree of the polynomial to use.
  //!\param Density the density matrix computed by this routine.
  //!\param solver_parameters parameters for the solver
  static void Compute(const Matrix_ps &Hamiltonian,
                      const Matrix_ps &InverseSquareRoot, int nel,
                      Matrix_ps &Density, int degree,
                      const SolverParameters &solver_parameters);
};
} // namespace NTPoly
#endif