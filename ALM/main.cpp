#include "Function.h"
#include "NewtonALM.h"
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
        scalar_t c4 = c/2 * (Aeq * x - beq).transpose() * (Aeq * x - beq);
        return ( c1 + c2 + c3 + c4 );
    }

    void Gradient(const vector_t &x,  const matrix_t &H, const vector_t &f, const matrix_t &Aeq, const vector_t &beq, const vector_t &lambda, const scalar_t c, vector_t *grad) const override {
        *grad = H * x + f + Aeq * lambda + c * Aeq * (Aeq*x - beq) ;
    }

    void Hessian(const vector_t &x,  const matrix_t &H, const vector_t &f, const matrix_t &Aeq, const vector_t &beq, const vector_t &lambda, const scalar_t c, hessian_t *hessian) const override {
        *hessian = H + c*Aeq*Aeq;
    }
};

int main(){
    Eigen::Vector3d x0;
    Eigen::Matrix3d H;
    Eigen::Vector3d f;
    Eigen::Matrix3d Aeq;
    Eigen::Vector3d beq;
    Eigen::Vector3d lb,ub;
    Eigen::Vector3d lambda;

    x0 << 0,
          0,
          0;
    H << 5, -2, -1,
         -2, 4, 3,
         -1, 3, 5;
    f << 2,
         -35,
         -47;

//    Aeq.setZero();
//    beq.setZero();
//    lb.setZero();
//    ub.setZero();
    lambda.setZero();
    Aeq << 1, 0, 0,
           0, 1, 0,
           0, 0, 0;
    beq << 2,
           3,
           0;
    lb << 0,
          0,
          0;
    ub << 5,
          5,
          5;


    Function fx;
    Function::scalar_t c(10);
    auto state = fx.Eval(x0, H, f, Aeq, beq, lb, ub, lambda, c);
    std::cout << "this" << std::endl;

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
    double eta{1e-4};
    int verbose = 1;

    if (verbose){
        // print header of outer iteration
        std::cout << "--------------------------------------------------------" << std::endl;
        std::cout << " k "  << "\t";
        std::cout << "  f(x_k)   " << "\t";
        std::cout << "  ||gradf(x_k)||  "  << "\t";
        std::cout << "  ||constraint(x_k)||  "  << "\t";
        std::cout << "  stepsize    " << std::endl;
        std::cout << "--------------------------------------------------------" << std::endl;
    }


    // need to check if x0 is feasible (look at CT Kelley code)

    // ganti stopping criteria buat ALM (liat buku N&W dan paper Andy)
    double normX = 100;
    double grad  = 100;
    double jac   = 100;
    int k = 0;
    while(grad > eta && normX > eta) {
        solver.setStoppingCriteria(epsilonk);
        auto[solution, solver_state] = solver.Minimize(fx, x0, H, f, Aeq, beq,  lb, ub, lambda, c); // think how to supply stopping criteria here!
        state = fx.Eval(solution.x, solution.H, solution.f, solution.Aeq, solution.beq, solution.lb, solution.ub, solution.lambda, solution.c);

        // compute constraint value
        cons = solution.Aeq * solution.x - solution.beq;
        jac  = cons.norm();
        if (jac <= etak){
            lambda   = lambda + c*cons;
            epsilonk = epsilonk/c;
            etak     = etak / pow(c,beta);
        } else {
            c        = tau*c;
            epsilonk = epsilon0/c;
            etak     = eta0/pow(c,alpha);

        }
        //normX = (solution.x - x0).norm();
        normX = (solution.x - x0).lpNorm<Eigen::Infinity>();
        grad  = state.gradient.template lpNorm<Eigen::Infinity>();
        x0    = solution.x;
        if (verbose){
            std::cout <<  k      << "\t";
            std::cout <<  state.value << "\t";
            std::cout <<  grad << "\t";
            std::cout <<  jac << "\t";
            std::cout <<  normX << std::endl;
        }
        k++;
    }

    std::cout << "argmin " << x0.transpose() << std::endl;
    return 0;
}

