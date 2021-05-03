#include "mex.h"
#include "MexEig"
#include "Function.h"
#include "NewtonALM.h"
#include <Eigen/Dense>
#include <iostream>

// almbound(H,f,Aeb,beq,lb,ub,x0,options)
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

