#include "mex.h"
#include "MexEig"
#include "Function.h"
#include "NewtonALM.h"
#include <Eigen/Sparse>
#include <iostream>
#include <type_traits>
#include <limits>

using FunctionXd   = cppoptlib::function::Function<double>;
using MatlabSparse =  Eigen::SparseMatrix<double,Eigen::ColMajor,std::make_signed<mwIndex>::type> ;

using FunctionXd = cppoptlib::function::Function<double>;

class Function : public FunctionXd {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using FunctionXd::hessian_t;
    using FunctionXd::vector_t;

    scalar_t operator()(const vector_t &x, const matrix_t &H, const vector_t &f, const matrix_t &Aeq, const vector_t &beq, const vector_t &lambda, const scalar_t c) const override {
        Eigen::SparseMatrix<double> xt = x.transpose();
        Eigen::SparseMatrix<double> ft = f.transpose();
        Eigen::SparseMatrix<double> lambdaT = lambda.transpose();
        Eigen::SparseMatrix<double> consT = (Aeq * x - beq).transpose();
        scalar_t c1 = (0.5 * xt * H * x).sum();
        scalar_t c2 = (ft * x).sum();
        scalar_t c3 = (lambdaT * ( Aeq* x - beq)).sum();
        scalar_t c4 = (c/2 * consT * (Aeq * x - beq)).sum();
        return ( c1 + c2 + c3 + c4 );
    }

    void Gradient(const vector_t &x,  const matrix_t &H, const vector_t &f, const matrix_t &Aeq, const vector_t &beq, const vector_t &lambda, const scalar_t c, vector_t *grad) const override {
        *grad = H * x + f + Aeq.transpose() * lambda + c * Aeq.transpose()* (Aeq*x - beq) ;
    }

    void Hessian(const vector_t &x,  const matrix_t &H, const vector_t &f, const matrix_t &Aeq, const vector_t &beq, const vector_t &lambda, const scalar_t c, hessian_t *hessian) const override {
        *hessian = H + c*Aeq.transpose()*Aeq;
    }
};

// copy methods from here!
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

mxArray* eigen_to_matlab_sparse(const Eigen::Ref<const MatlabSparse,Eigen::StandardCompressedFormat>& mat)
{
    mxArray * result = mxCreateSparse (mat.rows(), mat.cols(), mat.nonZeros(), mxREAL);
    const MatlabSparse::StorageIndex* ir = mat.innerIndexPtr();
    const MatlabSparse::StorageIndex* jc = mat.outerIndexPtr();
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

        Eigen::SparseMatrix<double> H   = matlab_to_eigen_sparse(prhs[0]);
        Eigen::SparseVector<double> f   = matlab_to_eigen_sparse(prhs[1]);
        Eigen::SparseMatrix<double> Aeq = matlab_to_eigen_sparse(prhs[2]);
        Eigen::SparseVector<double> beq = matlab_to_eigen_sparse(prhs[3]);
        Eigen::SparseVector<double> lb  = matlab_to_eigen_sparse(prhs[4]);
        Eigen::SparseVector<double> ub  = matlab_to_eigen_sparse(prhs[5]);
        Eigen::SparseVector<double> x0  = matlab_to_eigen_sparse(prhs[6]);
        Eigen::SparseVector<double> lambda = matlab_to_eigen_sparse(prhs[7]);

//        Aeq.setZero();
//        beq.setZero();
//        lb.setZero();
//        ub.setZero();
//        lambda.setZero();

        //plhs[0] = eigen_to_matlab_sparse(lambda);

//        Function fx;
//        Function::scalar_t c(10);
//        auto state = fx.Eval(x0, H, f, Aeq, beq, lb, ub, lambda, c);
//        std::cout << "this" << std::endl;
//
////        plhs[0] = EigenToMxArray(state.Aeq);
//
////        std::cout << fx(x0,lambda,c) << std::endl;
////        std::cout << state.gradient << std::endl;
////        std::cout << state.hessian << std::endl;
//
//        cppoptlib::solver::NewtonBound<Function> solver;
//
//        //double cons{0};
//        Eigen::Vector<double, Eigen::Dynamic> cons;
//        Eigen::Vector<double, Eigen::Dynamic> x;
//        double eta0{0.1258925};
//        double c0{10};
//        double epsilon0{1/c0};
//        double tau{10};
//        double alpha{0.1};
//        double beta{0.9};
//        double epsilonk = 1/c;
//        double etak = eta0 / pow(c,alpha);
//        double eta{1e-6};
//
//
//        // ganti stopping criteria buat ALM (liat buku N&W dan paper Andy)
//        while(state.gradient.template lpNorm<Eigen::Infinity>() > eta) {
//            solver.setStoppingCriteria(epsilonk);
//            auto[solution, solver_state] = solver.Minimize(fx, x0, H, f, Aeq, beq,  lb, ub, lambda, c); // think how to supply stopping criteria here!
//            state = fx.Eval(solution.x, solution.H, solution.f, solution.Aeq, solution.beq, solution.lb, solution.ub, solution.lambda, solution.c);
//
//            // compute constraint value
//            //cons = solution.x[0]*solution.x[0] + solution.x[1]*solution.x[1] - 1;
//            cons = solution.Aeq * solution.x - solution.beq;
//            if (cons.norm() <= etak){
//                lambda   = lambda + c*cons;
//                epsilonk = epsilonk/c;
//                etak     = etak / pow(c,beta);
//            } else {
//                c        = tau*c;
//                epsilonk = epsilon0/c;
//                etak     = eta0/pow(c,alpha);
//
//            }
//            x = solution.x;
//        }
//        plhs[0] = EigenToMxArray(x);
//
    }
    catch (std::exception& ex){
        mexErrMsgIdAndTxt("tmp::error", ex.what());
    }
}
