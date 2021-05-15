#include "Function.h"
#include "NewtonALM.h"
#include "MATio"

#include <iostream>
#include <iomanip>
#include <Eigen/Sparse>
#include <type_traits>
#include <limits>
#include <fmt/core.h>

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

int main(){

    matio::MatioFile file("qpData.mat");
    //matio::MatioFile file("qpSparse.mat");
    Eigen::MatrixXd H;
    Eigen::VectorXd f;
    Eigen::MatrixXd Aeq;
    Eigen::VectorXd beq;
    Eigen::VectorXd lb;
    Eigen::VectorXd ub;
    Eigen::VectorXd x0;
    Eigen::VectorXd lambda;


    if (file.read_mat("Hm", H)) {
        std::cout << "error: " << file.lasterr() << std::endl;
    }
    if (file.read_mat("f", f)) {
        std::cout << "error: " << file.lasterr() << std::endl;
    }
    if (file.read_mat("Aeq", Aeq)) {
        std::cout << "error: " << file.lasterr() << std::endl;
    }
    if (file.read_mat("beq", beq)) {
        std::cout << "error: " << file.lasterr() << std::endl;
    }
    if (file.read_mat("lb", lb)) {
        std::cout << "error: " << file.lasterr() << std::endl;
    }
    if (file.read_mat("ub", ub)) {
        std::cout << "error: " << file.lasterr() << std::endl;
    }


//    Eigen::Vector3d x0;
//    Eigen::Matrix3d H;
//    Eigen::Vector3d f;
//    Eigen::Matrix3d Aeq;
//    Eigen::Vector3d beq;
//    Eigen::Vector3d lb,ub;
//    Eigen::Vector3d lambda;
//
//    x0 << 0,
//          0,
//          0;
//    H << 5, -2, -1,
//         -2, 4, 3,
//         -1, 3, 5;
//    f << 2,
//         -35,
//         -47;
//
//    lambda.setZero();
//    Aeq << 0, 0, 0,
//           0, 0, 0,
//           0, 0, 0;
//    beq << 2,
//           3,
//           0;
//    lb << 0,
//          0,
//          0;
//    ub << 5,
//          5,
//          5;


    int numX = Aeq.cols();
    int numC = Aeq.rows();
    x0.setZero(numX);
    lambda.setZero(numC);

    // convert to sparse
    Eigen::SparseMatrix<double> Hs   = H.sparseView();
    Eigen::SparseVector<double> fs   = f.sparseView();
    Eigen::SparseMatrix<double> Aeqs = Aeq.sparseView();
    Eigen::SparseVector<double> beqs = beq.sparseView();
    Eigen::SparseVector<double> lbs  = lb.sparseView();
    Eigen::SparseVector<double> ubs  = ub.sparseView();
    Eigen::SparseVector<double> x0s  = x0.sparseView();
    Eigen::SparseVector<double> lambdas  = lambda.sparseView();
    Function fx;
    Function::scalar_t c(10);
    //Function::scalar_t c(1);
    //auto state = fx.Eval(x0, H, f, Aeq, beq, lb, ub, lambda, c);
    auto state = fx.Eval(x0s, Hs, fs, Aeqs, beqs, lbs, ubs, lambdas, c);

    cppoptlib::solver::NewtonBound<Function> solver;

    //double cons{0};
    //Eigen::Vector<double, Eigen::Dynamic> cons;
    Eigen::SparseVector<double> cons(numC);
    //Eigen::Vector<double, Eigen::Dynamic> x;
    Eigen::SparseVector<double> x(numX);
    double eta0{0.1258925};
    double c0{10};
    double epsilon0{1/c0};
    double tau{10};
    double alpha{0.1};
    double beta{0.9};
    double epsilonk = 1/c;
    double etak = eta0 / pow(c,alpha);
    double eta{1e-1};
    int verbose = 1;

    if (verbose){
        // print header of outer iteration
        std::cout << "-------------------------------------------------------------------------------------------" << std::endl;
        std::cout << "k    "     << "\t";
        std::cout << "f(x_k)"    << "\t" << "\t";
        std::cout << "gradf(x_k)"<< "\t";
        std::cout << "cons(x_k)" << "\t";
        std::cout << "stepsize"  << "\t";
	std::cout << "epsilonk"  << std::endl;
        std::cout << "-------------------------------------------------------------------------------------------" << std::endl;
    }


    // need to check if x0 is feasible (look at CT Kelley code)

    // ganti stopping criteria buat ALM (liat buku N&W dan paper Andy)
    double normX{100};
    double grad{100};
    double jac;
    int k = 0;
    while(grad > eta && normX > eta) {
        solver.setStoppingCriteria(epsilonk);
        //auto[solution, solver_state] = solver.Minimize(fx, x0, H, f, Aeq, beq,  lb, ub, lambda, c); // think how to supply stopping criteria here!
        auto[solution, solver_state] = solver.Minimize(fx, x0s, Hs, fs, Aeqs, beqs,  lbs, ubs, lambdas, c);
        state = fx.Eval(solution.x, solution.H, solution.f, solution.Aeq, solution.beq, solution.lb, solution.ub, solution.lambda, solution.c);

        // compute constraint value
        cons = solution.Aeq * solution.x - solution.beq;
        jac  = cons.norm();
        if (jac <= etak){
            lambda   = lambda + c*cons;
	    if (epsilonk > 0.01){
		epsilonk = epsilonk/c; 
	    }
            //epsilonk = epsilonk/c;
            etak     = etak / pow(c,beta);
        } else {
            c        = tau*c;
	    if (epsilonk > 0.01){
		epsilonk = epsilon0/c;
	    }
            //epsilonk = epsilon0/c;
            etak     = eta0/pow(c,alpha);

        }
        normX = (solution.x - x0).norm();
        //normX = (solution.x - x0).lpNorm<Eigen::Infinity>();
        //grad  = state.gradient.template lpNorm<Eigen::Infinity>();
        grad  = state.gradient.norm();
        x0    = solution.x;
        if (verbose){
            std::cout <<  k      << "\t" ;
            std::cout << fmt::format("{:.4e}", state.value) << "\t";
            std::cout << fmt::format("{:.4e}", grad)        << "\t" ;
            std::cout << fmt::format("{:.4e}", jac)         << "\t";
            std::cout << fmt::format("{:.4e}", normX)       << "\t";
	    std::cout << fmt::format("{:.4e}", epsilonk)    << std::endl;
        }
        k++;
    }

    //std::cout << "argmin " << x0.transpose() << std::endl;
    return 0;
}

