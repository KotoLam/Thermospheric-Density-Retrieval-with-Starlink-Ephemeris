#include "mex.h"
#include "sofa.h"
// By Deepseek
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // 检查输入参数
    if (nrhs != 3) {
        mexErrMsgIdAndTxt("SOFA:inputError", "需要3个输入参数：xp, yp, sp");
    }
    
    // 获取输入参数
    double xp = mxGetScalar(prhs[0]);
    double yp = mxGetScalar(prhs[1]);
    double sp = mxGetScalar(prhs[2]);
    
    // 准备输出矩阵
    double rpom[3][3];
    
    // 调用SOFA函数
    iauPom00(xp, yp, sp, rpom);
    
    // 创建MATLAB输出矩阵
    plhs[0] = mxCreateDoubleMatrix(3, 3, mxREAL);
    double *out = mxGetPr(plhs[0]);
    
    // 复制数据到MATLAB矩阵（注意MATLAB是列优先）
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            out[j*3 + i] = rpom[i][j];
        }
    }
}