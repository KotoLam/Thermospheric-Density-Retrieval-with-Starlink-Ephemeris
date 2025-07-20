#include "mex.h"
#include "sofa.h"
// By Deepseek
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /* 检查输入参数 */
    if (nrhs != 4) {
        mexErrMsgIdAndTxt("SOFA:inputError", "需要4个参数：uta, utb, tta, ttb");
    }
    
    /* 获取输入参数 */
    double uta = mxGetScalar(prhs[0]);
    double utb = mxGetScalar(prhs[1]);
    double tta = mxGetScalar(prhs[2]);
    double ttb = mxGetScalar(prhs[3]);
    
    /* 调用SOFA函数 */
    double gmst = iauGmst06(uta, utb, tta, ttb);
    
    /* 返回结果 */
    plhs[0] = mxCreateDoubleScalar(gmst);
}