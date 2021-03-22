#include <Eigen/Dense>

struct GradProj{

    GradProj(Eigen::MatrixXd Q, Eigen::VectorXd b, Eigen::VectorXd x0, Eigen::VectorXd lb, Eigen::VectorXd ub,int maxIter=100):Q_(Q), b_(b), x0_(x0), lb_(lb), ub_(ub), maxIter_(maxIter){

    }

    Eigen::VectorXd solve() const noexcept{
        return x_;
    }

private:
    Eigen::MatrixXd Q_;
    Eigen::VectorXd b_;
    Eigen::VectorXd x_;
    Eigen::VectorXd lb_;
    Eigen::VectorXd ub_;
    Eigen::VectorXd x0_;
    int maxIter_;
};
