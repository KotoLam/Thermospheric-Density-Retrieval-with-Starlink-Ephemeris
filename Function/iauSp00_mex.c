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
    
    // 调用SOFA函数
    double result = iauSp00(date1, date2);
    
    // 返回结果
    plhs[0] = mxCreateDoubleScalar(result);
}