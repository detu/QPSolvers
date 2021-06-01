#include "mex.h"
#include "Function.h"
#include "NewtonALM.h"
#include <fmt/core.h>
#include <iostream>
#include <type_traits>
#include <limits>

using FunctionXd   = cppoptlib::function::Function<double>;
using MatlabSparse =  Eigen::SparseMatrix<double,Eigen::ColMajor,std::make_signed<mwIndex>::type> ;

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

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    try {
        if (nrhs < 1)
            throw std::invalid_argument("required one input arg");
        if (nlhs != 1)
            throw std::invalid_argument("required one output arg");

//        Eigen::SparseMatrix<double> H = matlab_to_eigen_sparse(prhs[0]);
//        Eigen::SparseVector<double> f = matlab_to_eigen_sparse(prhs[1]);
//        Eigen::SparseMatrix<double> Aeq = matlab_to_eigen_sparse(prhs[2]);
//        Eigen::SparseVector<double> beq = matlab_to_eigen_sparse(prhs[3]);
//        Eigen::SparseVector<double> lb = matlab_to_eigen_sparse(prhs[4]);
//        Eigen::SparseVector<double> ub = matlab_to_eigen_sparse(prhs[5]);
//        Eigen::SparseVector<double> x0 = matlab_to_eigen_sparse(prhs[6]);
//        Eigen::SparseVector<double> lambda = matlab_to_eigen_sparse(prhs[7]);

        Eigen::SparseMatrix<double> H = matlab_to_eigen_sparse(std::move(prhs[0]));
        Eigen::SparseVector<double> f = matlab_to_eigen_sparse(std::move(prhs[1]));
        Eigen::SparseMatrix<double> Aeq = matlab_to_eigen_sparse(std::move(prhs[2]));
        Eigen::SparseVector<double> beq = matlab_to_eigen_sparse(std::move(prhs[3]));
        Eigen::SparseVector<double> lb = matlab_to_eigen_sparse(std::move(prhs[4]));
        Eigen::SparseVector<double> ub = matlab_to_eigen_sparse(std::move(prhs[5]));
        Eigen::SparseVector<double> x0 = matlab_to_eigen_sparse(std::move(prhs[6]));
        Eigen::SparseVector<double> lambda = matlab_to_eigen_sparse(std::move(prhs[7]));

        int numX = Aeq.cols();
        int numC = Aeq.rows();

        Function fx;
        Function::scalar_t c(10);

        cppoptlib::solver::NewtonBound<Function> solver;

        Eigen::SparseVector<double> cons(numC);
        Eigen::SparseVector<double> x(numX);
        double eta0{0.1258925};
        double c0{10};
        double epsilon0{1 / c0};
        double tau{10};
        double alpha{0.1};
        double beta{0.9};
        double epsilonk = 1 / c;
        double etak = eta0 / pow(c, alpha);
        double eta{1e-2};
        int verbose = 1;

        if (verbose) {
            // print header of outer iteration
            std::cout << "---------------------------------------------------------------------------------------"
                      << std::endl;
            std::cout << "k    " << "\t";
            std::cout << "f(x_k)" << "\t" << "\t";
            std::cout << "gradf(x_k)" << "\t";
            std::cout << "cons(x_k)" << "\t";
            std::cout << "stepsize" << "\t";
            std::cout << "epsilonk" << std::endl;
            std::cout << "---------------------------------------------------------------------------------------"
                      << std::endl;
        }


        // need to check if x0 is feasible (look at CT Kelley code)

        // ganti stopping criteria buat ALM (liat buku N&W dan paper Andy)
        double normX{100};
        double grad{100};
        double jac;
        int k = 1;
        auto state = fx.Eval(x0, H, f, Aeq, beq, lb, ub, lambda, c);


        while (grad > eta && normX > eta) {
            //solver.setStoppingCriteria(epsilonk);
            //auto[solution, solver_state] = solver.Minimize(fx, x0s, Hs, fs, Aeqs, beqs,  lbs, ubs, lambdas, c);
            auto[solution, solver_state] = solver.Minimize(fx, state);

            // compute constraint value
            cons = solution.Aeq * solution.x - solution.beq;
            jac = cons.norm();
            if (jac <= etak) {
                lambda = lambda + c * cons;
                if (epsilonk > 0.01) {
                    epsilonk = epsilonk / c;
                }
                //epsilonk = epsilonk/c;
                etak = etak / pow(c, beta);
            } else {
                c = tau * c;
                if (epsilonk > 0.01) {
                    epsilonk = epsilon0 / c;
                }
                //epsilonk = epsilon0/c;
                etak = eta0 / pow(c, alpha);

            }
            normX = (solution.x - x0).norm();
            grad = solver_state.gradient_norm;
            x0 = solution.x;
            if (verbose) {
                std::cout << k << "\t";
                std::cout << fmt::format("{:.4e}", solution.value) << "\t";
                std::cout << fmt::format("{:.4e}", grad) << "\t";
                std::cout << fmt::format("{:.4e}", jac) << "\t";
                std::cout << fmt::format("{:.4e}", normX) << "\t";
                std::cout << fmt::format("{:.4e}", epsilonk) << std::endl;
            }
            auto state = solution;
            k++;

        }
        plhs[0] = eigen_to_matlab_sparse(std::move(x0));
    }
    catch (std::exception& ex){
        mexErrMsgIdAndTxt("tmp::error", ex.what());
    }
}

