#include "mex.h"
#include "MexEig"
#include "MUMPSSupport.h"
#include <Eigen/Sparse>
#include <type_traits>
#include <limits>

typedef Eigen::SparseMatrix<double,Eigen::ColMajor,std::make_signed<mwIndex>::type> MatlabSparse;

Eigen::Map<MatlabSparse > matlab_to_eigen_sparse(const mxArray * mat)
{
    mxAssert(mxGetClassID(mat) == mxDOUBLE_CLASS,
             "Type of the input matrix isn't double");
    mwSize     m = mxGetM(mat);
    mwSize     n = mxGetN(mat);
    mwSize    nz = mxGetNzmax(mat);
    /*Theoretically fails in very very large matrices*/
    mxAssert(nz <= std::numeric_limits< std::make_signed<mwIndex>::type>::max(),
             "Unsupported Data size.");
    double  * pr = mxGetPr(mat);
    MatlabSparse::StorageIndex* ir = reinterpret_cast<MatlabSparse::StorageIndex*>(mxGetIr (mat));
    MatlabSparse::StorageIndex* jc = reinterpret_cast<MatlabSparse::StorageIndex*>(mxGetJc (mat));
    Eigen::Map<MatlabSparse> result(m, n, nz, jc, ir, pr);
    return result;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])

{

//    if(!(mxIsSparse(prhs[0]))){
//        mexErrMsgIdAndTxt("MATLAB: ", "First input matrix is not sparse.\n");
//    }
//
//    if((mxIsSparse(prhs[1]))){
//        mexErrMsgIdAndTxt("MATLAB: ", "Second input matrix is not dense.\n");
//    }

// prhs[0]: A
    Eigen::Map<MatlabSparse> A = matlab_to_eigen_sparse(prhs[0]);
    mwSize N = mxGetM(prhs[0]); // size of BAUs for constructing the GGM

// prhs[1]: b
    Eigen::Map<MatlabSparse> b = matlab_to_eigen_sparse(prhs[1]);
    mwSize n = mxGetM(prhs[1]);


    //Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    //Eigen::MUMPSLDLT<MatlabSparse ,Eigen::Lower> solver;
    //Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    //Eigen::MUMPSLU<Eigen::SparseMatrix<double>> solver;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    Eigen::VectorXd xsol = solver.solve(b);

    plhs[0] = EigenToMxArray(xsol);

}



