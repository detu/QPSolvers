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

//    scalar_t operator()(const vector_t &x, const scalar_t lambda, const scalar_t c) const override {
//        return  -x[0] - x[1] + lambda * (x[0]*x[0] + x[1]*x[1] - 1) + c/2 * ( (x[0]*x[0] + x[1]*x[1] - 1) * (x[0]*x[0] + x[1]*x[1] - 1) );
//    }
    scalar_t operator()(const vector_t &x, const matrix_t &H, const vector_t &f, const matrix_t &Aeq, const vector_t &beq, const vector_t lambda, const scalar_t c) const override {
        //return 0.5*x.transpose()*H*x + f.transpose()*x + (lambda.transpose() * (Aeq * x - beq)) + (c/2 * ( (Aeq * x - beq).squaredNorm() ));
        return static_cast<scalar_t>(0.5 * x.transpose() * H * x + f.transpose() * x +
                                     lambda.transpose() * (Aeq * x - beq));
    }

    void Gradient(const vector_t &x, const scalar_t lambda, const scalar_t c, vector_t *grad) const override {
        (*grad)[0] = -1 + 2*lambda*x[0] + c * ( (x[0]*x[0] + x[1]*x[1] - 1) * 2 * x[0] );
        (*grad)[1] = -1 + 2*lambda*x[1] + c * ( (x[0]*x[0] + x[1]*x[1] - 1) * 2 * x[1] );
    }

    void Hessian(const vector_t &x, const scalar_t lambda, const scalar_t c, hessian_t *hessian) const override {
        (*hessian)(0, 0) = 2*lambda + 2*c*(3*x[0]*x[0] + x[1]*x[1] - 1);
        (*hessian)(0, 1) = 4*c*x[1]*x[0];
        (*hessian)(1, 0) = 4*c*x[0]*x[1];
        (*hessian)(1, 1) = 2*lambda + 2*c*(x[0]*x[0] + 3*x[1]*x[1] - 1);
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

