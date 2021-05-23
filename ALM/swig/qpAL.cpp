#include "mex.h"
#define EIGEN_USE_MKL_ALL
#include "MUMPSSupport"
#include <Eigen/Sparse>
//#include <Eigen/SparseCholesky>


Eigen::SparseVector<double> qpAL(const Eigen::SparseMatrix<double>& H, const Eigen::SparseVector<double>& f)
{
   Eigen::MUMPSLDLT<Eigen::SparseMatrix<double>, Eigen::Upper|Eigen::Lower> solver;;
   //Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
   //Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(H);
    Eigen::SparseVector<double>  x = solver.solve(std::move(-f));
    return x;
}

