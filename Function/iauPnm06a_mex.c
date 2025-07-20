#include "mex.h"
#include "sofa.h"
// By Deepseek
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // 检查输入参数
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("SOFA:inputError", "需要2个输入参数：date1, date2");
    }
    
    // 获取输入参数
    double date1 = mxGetScalar(prhs[0]);
    double date2 = mxGetScalar(prhs[1]);
    
    // 准备输出矩阵
    double rbpn[3][3];
    
    // 调用SOFA函数
    iauPnm06a(date1, date2, rbpn);
    
    // 创建MATLAB输出矩阵
    plhs[0] = mxCreateDoubleMatrix(3, 3, mxREAL);
    double *out = mxGetPr(plhs[0]);
    
    // 复制数据到MATLAB矩阵（注意MATLAB是列优先）
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            out[j*3 + i] = rbpn[i][j];
        }
    }
}