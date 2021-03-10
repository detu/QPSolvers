// Class implementation of ALM
// [x, fval, exitflag, output, lambda] = ALMqp(H,f,A,b,Aeq,beq,lb,ub,x0,options)
// like matlab's quadprog

#include <Eigen/Dense>
#include <Eigen/Sparse>

struct ALMqp{

    ALMqp(Eigen::SparseMatrix<double> H,Eigen::SparseVector<double> f, Eigen::SparseMatrix<double> Aeq, Eigen::SparseVector<double> beq, Eigen::SparseVector<double> lb, Eigen::SparseVector<double> ub): H_{H}, f_(f), Aeq_(Aeq), beq_(beq), lb_(lb), ub_(ub)
    {
        
    }

    Eigen::VectorXd sol() const noexcept{
        return x_;
    }

private:
    Eigen::SparseMatrix<double> H_;
    Eigen::SparseVector<double> f_;
    Eigen::SparseMatrix<double> Aeq_;
    Eigen::SparseVector<double> beq_;
    Eigen::SparseVector<double> lb_;
    Eigen::SparseVector<double> ub_;
    Eigen::VectorXd x_;
};