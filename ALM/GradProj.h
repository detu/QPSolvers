#include <Eigen/Dense>

struct GradProj{

    GradProj(Eigen::MatrixXd Q, Eigen::VectorXd b, Eigen::VectorXd x0, Eigen::VectorXd lb, Eigen::VectorXd ub, double tol = 1e-6, int maxIter=100):Q_(Q), b_(b), x0_(x0), lb_(lb), ub_(ub), tol_(tol), maxIter_(maxIter){
        auto ndim = ub_.size();
        Eigen::VectorXd xc = x0_;
        Eigen::VectorXd kku = Eigen::MatrixXd::Zero(ndim,1);
        Eigen::VectorXd kkl = Eigen::MatrixXd::Zero(ndim,1);
        kku = ub_;
        kkl = lb_;
        // There should be a checking that for every element ub_ > lb_ here !
        //double fc =
        Eigen::VectorXd gc  = Q_*x_ + b_;
        Eigen::VectorXd xt  = kkProj(xc - gc, kku, kkl);
        Eigen::VectorXd pgc = xc - xt;
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
    double tol_;
    int maxIter_;

    Eigen::VectorXd kkProj(Eigen::VectorXd &x, Eigen::VectorXd &kku, Eigen::VectorXd && kkl){
        auto ndim = x.size();
        Eigen::VectorXd px = Eigen::MatrixXd::Zero(ndim, 1);
        px = x.array().min(kku.array());
        px = px.array().max(kkl.array());
        return px;
    }
};
