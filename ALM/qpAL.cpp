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

    scalar_t operator()(const vector_t &x, const matrix_t &H, const vector_t &f, const matrix_t &Aeq, const vector_t &beq, const vector_t &lambda, const scalar_t c) const override {
        scalar_t c1 = 0.5 * x.transpose() * H * x;
        scalar_t c2 = f.transpose() * x;
        scalar_t c3 = lambda.transpose() * ( Aeq* x - beq);
        scalar_t c4 = c/2 * (Aeq * x - beq).squaredNorm() * (Aeq * x - beq).squaredNorm();
        return ( c1 + c2 + c3 + c4 );
    }

    void Gradient(const vector_t &x,  const matrix_t &H, const vector_t &f, const matrix_t &Aeq, const vector_t &beq, const vector_t &lambda, const scalar_t c, vector_t *grad) const override {
        *grad = H * x + f + Aeq * lambda + c * (Aeq*x - beq);
    }

    void Hessian(const vector_t &x,  const matrix_t &H, const vector_t &f, const matrix_t &Aeq, const vector_t &beq, const vector_t &lambda, const scalar_t c, hessian_t *hessian) const override {
        *hessian = H;
    }
};

// copy methods from here!
//https://stackoverflow.com/questions/49952275/passing-sparse-arrays-from-matlab-to-eigen-c-and-back-to-matlab


// [x, ..] = almbound(H,f,Aeb,beq,lb,ub,x0,lambda,options)
// debug: save matrices and vectors in a *.mat file!
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    try {
        if (nrhs < 1)
            throw std::invalid_argument("required one input arg");
        if (nlhs != 1)
            throw std::invalid_argument("required one output arg");

        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> H;
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> f;
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Aeq;
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> beq;
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> lb;
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> ub;
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x0;
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> lambda;

        MxArrayToEigen(H, prhs[0]);
        MxArrayToEigen(f, prhs[1]);
        MxArrayToEigen(Aeq, prhs[2]);
        MxArrayToEigen(beq, prhs[3]);
        MxArrayToEigen(lb, prhs[4]);
        MxArrayToEigen(ub, prhs[5]);
        MxArrayToEigen(x0, prhs[6]);
        MxArrayToEigen(lambda, prhs[7]);

//        Aeq.setZero();
//        beq.setZero();
//        lb.setZero();
//        ub.setZero();
//        lambda.setZero();

//        plhs[0] = EigenToMxArray(lambda);

        Function fx;
        Function::scalar_t c(10);
        auto state = fx.Eval(x0, H, f, Aeq, beq, lb, ub, lambda, c);
        std::cout << "this" << std::endl;

//        plhs[0] = EigenToMxArray(state.Aeq);

//        std::cout << fx(x0,lambda,c) << std::endl;
//        std::cout << state.gradient << std::endl;
//        std::cout << state.hessian << std::endl;

        cppoptlib::solver::NewtonBound<Function> solver;

        //double cons{0};
        Eigen::Vector<double, Eigen::Dynamic> cons;
        Eigen::Vector<double, Eigen::Dynamic> x;
        double eta0{0.1258925};
        double c0{10};
        double epsilon0{1/c0};
        double tau{10};
        double alpha{0.1};
        double beta{0.9};
        double epsilonk = 1/c;
        double etak = eta0 / pow(c,alpha);
        double eta{1e-6};


        // ganti stopping criteria buat ALM (liat buku N&W dan paper Andy)
        while(state.gradient.template lpNorm<Eigen::Infinity>() > eta) {
            solver.setStoppingCriteria(epsilonk);
            auto[solution, solver_state] = solver.Minimize(fx, x0, H, f, Aeq, beq,  lb, ub, lambda, c); // think how to supply stopping criteria here!
            state = fx.Eval(solution.x, solution.H, solution.f, solution.Aeq, solution.beq, solution.lb, solution.ub, solution.lambda, solution.c);

            // compute constraint value
            //cons = solution.x[0]*solution.x[0] + solution.x[1]*solution.x[1] - 1;
            cons = solution.Aeq * solution.x - solution.beq;
            if (cons.norm() <= etak){
                lambda   = lambda + c*cons;
                epsilonk = epsilonk/c;
                etak     = etak / pow(c,beta);
            } else {
                c        = tau*c;
                epsilonk = epsilon0/c;
                etak     = eta0/pow(c,alpha);

            }
            x = solution.x;
        }
        plhs[0] = EigenToMxArray(x);

    }
    catch (std::exception& ex){
        mexErrMsgIdAndTxt("tmp::error", ex.what());
    }
}

