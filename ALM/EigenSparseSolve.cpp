#include "EigenSparseSolve.h"
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>
#include "MUMPSSupport"

namespace NuTo
{

template <typename TSolver>
Eigen::VectorXd SolveWithSolver(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& b)
{
    TSolver solver;
    solver.analyzePattern(A);
    solver.factorize(A);
    return solver.solve(b);
}

Eigen::VectorXd EigenSparseSolve(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& b, std::string solver)
{
    // direct solvers
    if (solver == "EigenSparseLU")
        return SolveWithSolver<Eigen::SparseLU<Eigen::SparseMatrix<double>>>(A, b);
    if (solver == "EigenSparseQR")
        return SolveWithSolver<Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>>(A, b);
    if (solver == "EigenSimplicialLLT")
        return SolveWithSolver<Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>>(A, b);
    if (solver == "EigenSimplicialLDLT")
        return SolveWithSolver<Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>>(A, b);

    // iterative solvers
    if (solver == "EigenConjugateGradient")
        return SolveWithSolver<Eigen::ConjugateGradient<Eigen::SparseMatrix<double>>>(A, b);
    if (solver == "EigenLeastSquaresConjugateGradient")
        return SolveWithSolver<Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>>>(A, b);
    if (solver == "EigenBiCGSTAB")
        return SolveWithSolver<Eigen::BiCGSTAB<Eigen::SparseMatrix<double>>>(A, b);

    if (solver == "MumpsLU")
        return SolveWithSolver<Eigen::MUMPSLU<Eigen::SparseMatrix<double>>>(A, b);
    if (solver == "MumpsLDLT")
        return SolveWithSolver<Eigen::MUMPSLDLT<Eigen::SparseMatrix<double>, Eigen::Upper>>(A, b);
}


EigenSparseSolver::EigenSparseSolver(std::string solver)
    : mSolver(solver)
{
}

Eigen::VectorXd EigenSparseSolver::Solve(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& b) const
{
    return EigenSparseSolve(A, b, mSolver);
}

} // namespace NuTo

