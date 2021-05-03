#include "mex.h"
#include "MexEig"
#include "Function.h"
#include "NewtonALM.h"
#include <Eigen/Dense>
#include <iostream>

using FunctionXd = cppoptlib::function::Function<double>;

class Function : public FunctionXd {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using FunctionXd::hessian_t;
    using FunctionXd::vector_t;

    scalar_t operator()(const vector_t &x, const matrix_t &H, const vector_t &f, const matrix_t &Aeq, const vector_t &beq, const vector_t lambda, const scalar_t c) const override {
        scalar_t c1 = 0.5 * x.transpose() * H * x;
        scalar_t c2 = f.transpose() * x;
        scalar_t c3 = lambda.transpose() * ( Aeq* x - beq);
        scalar_t c4 = c/2 * (Aeq * x - beq).squaredNorm();
        return ( c1 + c2 + c3 + c4 );
    }

    void Gradient(const vector_t &x,  const matrix_t &H, const vector_t &f, const matrix_t &Aeq, const vector_t &beq, const vector_t lambda, const scalar_t c, vector_t *grad) const override {
        (*grad) = H * x + f.transpose() + Aeq * lambda.transpose() + c/2 * Aeq * Aeq.transpose();
    }

    void Hessian(const vector_t &x,  const matrix_t &H, const vector_t &f, const matrix_t &Aeq, const vector_t &beq, const vector_t lambda, const scalar_t c, hessian_t *hessian) const override {
        (*hessian) = H;
    }
};

// almbound(H,f,Aeb,beq,lb,ub,x0,options)
// debug: save matrices and vectors in a *.mat file!
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    try {
        if (nrhs != 1)
            throw std::invalid_argument("required one input arg");
        if (nlhs != 1)
            throw std::invalid_argument("required one output arg");

        Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> xi(0,0);
        MxArrayToEigen(xi, prhs[0]);
        plhs[0] = EigenToMxArray(xi);
    }
    catch (std::exception& ex){
        mexErrMsgIdAndTxt("tmp::error", ex.what());
    }
}

