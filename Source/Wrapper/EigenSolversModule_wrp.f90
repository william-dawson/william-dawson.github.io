!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the eigensolveres module for calling from other languages.
MODULE EigenSolversModule_wrp
  USE EigenSolversModule, ONLY : EigenDecomposition
  USE DistributedSparseMatrixModule_wrp, ONLY : DistributedSparseMatrix_wrp
  USE IterativeSolversModule_wrp, ONLY : IterativeSolverParameters_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: EigenDecomposition_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the eigenvalues and eigenvectors of a matrix.
  SUBROUTINE EigenDecomposition_wrp(ih_this, ih_eigenvectors, ih_eigenvalues, &
       & ih_solver_parameters) bind(c,name="EigenDecomposition_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_eigenvectors(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_eigenvalues(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_this
    TYPE(DistributedSparseMatrix_wrp) :: h_eigenvectors
    TYPE(DistributedSparseMatrix_wrp) :: h_eigenvalues
    TYPE(IterativeSolverParameters_wrp) :: h_solver_parameters

    h_this = TRANSFER(ih_this,h_this)
    h_eigenvectors = TRANSFER(ih_eigenvectors,h_eigenvectors)
    h_eigenvalues = TRANSFER(ih_eigenvalues,h_eigenvalues)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL EigenDecomposition(h_this%data, h_eigenvectors%data, &
         & h_eigenvalues%data, h_solver_parameters%data)
  END SUBROUTINE EigenDecomposition_wrp
END MODULE EigenSolversModule_wrp