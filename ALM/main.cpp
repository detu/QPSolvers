#include "DirectQP.h"
#include <Eigen/Dense>
#include <iostream>

int main(){
    Eigen::Matrix4f Q;
    Eigen::Vector4f b;
    Q << 1,1,1,1,
         1,2,2,2,
         1,2,3,3,
         1,2,3,4;
    b << -4,
         -7,
         -9,
         -10;

    auto qp = DirectQP(Q,b);
    auto x  = qp.solve();
    std::cout << "The solution is:\n" << x << std::endl;

    return 0;
}

