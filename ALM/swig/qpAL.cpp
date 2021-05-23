#include "mex.h"
#include <Eigen/Sparse>
#include <fmt/core.h>
#include <iostream>
#include <type_traits>
#include <limits>

#define EIGEN_USE_MKL_ALL
#include "MUMPSSupport"
//#include <Eigen/Sparse>
#include<Eigen/SparseCholesky>

//https://stackoverflow.com/questions/49952275/passing-sparse-arrays-from-matlab-to-eigen-c-and-back-to-matlab
Eigen::Map<MatlabSparse > matlab_to_eigen_sparse(const mxArray * mat)
{
    mxAssert(mxGetClassID(mat) == mxDOUBLE_CLASS,
             "Type of the input matrix isn't double");
    mwSize     m = mxGetM (mat);
    mwSize     n = mxGetN (mat);
    mwSize    nz = mxGetNzmax (mat);
    /*Theoretically fails in very very large matrices*/
    mxAssert(nz <= std::numeric_limits< std::make_signed<mwIndex>::type>::max(),
             "Unsupported Data size."
    );
    double  * pr = mxGetPr (mat);
    MatlabSparse::StorageIndex* ir = reinterpret_cast<MatlabSparse::StorageIndex*>(mxGetIr (mat));
    MatlabSparse::StorageIndex* jc = reinterpret_cast<MatlabSparse::StorageIndex*>(mxGetJc (mat));
    Eigen::Map<MatlabSparse> result (m, n, nz, jc, ir, pr);
    return result;
}

mxArray* eigen_to_matlab_sparse(const Eigen::Ref<const Eigen::SparseVector<double>,Eigen::StandardCompressedFormat>& mat)
{
    mxArray * result = mxCreateSparse (mat.rows(), mat.cols(), mat.nonZeros(), mxREAL);
    //const MatlabSparse::StorageIndex* ir = mat.innerIndexPtr();
    //const MatlabSparse::StorageIndex* jc = mat.outerIndexPtr();
    const Eigen::SparseVector<double>::StorageIndex* ir = mat.innerIndexPtr();
    const Eigen::SparseVector<double>::StorageIndex* jc = mat.outerIndexPtr();
    const double* pr = mat.valuePtr();

    mwIndex * ir2 = mxGetIr (result);
    mwIndex * jc2 = mxGetJc (result);
    double  * pr2 = mxGetPr (result);

    for (mwIndex i = 0; i < mat.nonZeros(); i++) {
        pr2[i] = pr[i];
        ir2[i] = ir[i];
    }
    for (mwIndex i = 0; i < mat.cols() + 1; i++) {
        jc2[i] = jc[i];
    }
    return result;
}

// [x, ..] = almbound(H,f,Aeb,beq,lb,ub,x0,lambda,options)
// debug: save matrices and vectors in a *.mat file!
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    try {
        if (nrhs < 1)
            throw std::invalid_argument("required one input arg");
        if (nlhs != 1)
            throw std::invalid_argument("required one output arg");

        Eigen::MatrixXd H;
        Eigen::VectorXd f;

        MxArrayToEigen(H, prhs[0]);
        MxArrayToEigen(f, prhs[1]);
        Eigen::SparseMatrix<double> Hs   = H.sparseView();
        Eigen::SparseVector<double> fs   = f.sparseView();

        Eigen::MUMPSLDLT<Eigen::SparseMatrix<double>, Eigen::Upper|Eigen::Lower> solver;;
        //Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        //Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(Hs);
        auto infoH = solver.info();
        //if(!solver.info()) {
            // decomposition failed
        //    return;
        //}




        //Eigen::SparseVector<double> delta_x = solver.solve(std::move(-state.gradient));
        Eigen::SparseVector<double>  delta_x = solver.solve(std::move(-fs));
        auto infoX = solver.info();
//        if(!solver.info()) {
//            // solving failed
//            return;
//        }
        plhs[0] = EigenToMxArray(delta_x.toDense()); // change to eigen_to_matlab_sparse!
        //plhs[0] = eigen_to_matlab_sparse(x0);
        //plhs[0] = eigen_to_matlab_sparse(delta_x);
    }
    catch (std::exception& ex){
        mexErrMsgIdAndTxt("tmp::error", ex.what());
    }
}

