#include "mex.h"
#include "sofa.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // 检查输入参数
    if (nrhs != 5) {
        mexErrMsgIdAndTxt("SOFA:inputError", "需要5个输入参数：uta,utb,tta,ttb,rnpb");
    }
    
    // 获取标量输入
    double uta = mxGetScalar(prhs[0]);
    double utb = mxGetScalar(prhs[1]);
    double tta = mxGetScalar(prhs[2]);
    double ttb = mxGetScalar(prhs[3]);
    
    // 获取3x3矩阵输入
    if (mxGetM(prhs[4]) != 3 || mxGetN(prhs[4]) != 3) {
        mexErrMsgIdAndTxt("SOFA:inputError", "rnpb必须是3x3矩阵");
    }
    double *rnpb_data = mxGetPr(prhs[4]);
    double rnpb[3][3];
    
    // 将MATLAB矩阵转换为C数组
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            rnpb[i][j] = rnpb_data[j*3 + i]; // MATLAB是列优先存储
        }
    }
    
    // 调用SOFA函数
    double result = iauGst06(uta, utb, tta, ttb, rnpb);
    
    // 返回结果
    plhs[0] = mxCreateDoubleScalar(result);
}