//#include "DirectQP.h"
#include "NewtonLS.h"
#include <Eigen/Dense>
#include <iostream>

int main(){
    Eigen::Matrix4d Q;
    Eigen::Vector4d b;
    Q << 1,1,1,1,
         1,2,2,2,
         1,2,3,3,
         1,2,3,4;
    b << -4,
         -7,
         -9,
         -10;

    //auto qp = DirectQP(Q,b);
    Eigen::Vector4d x0{0,0,0,0};
    auto qp = NewtonLS(Q,b,x0);
    auto x  = qp.solve();
    auto fobj = (x.transpose()*Q*x + x.transpose()*b).value();
    std::cout << "The solution is:\n" << x << std::endl;
    std::cout << "The optimal objective function value:\n" << fobj << std::endl;

    return 0;
}

