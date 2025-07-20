#include "mex.h"
#include <math.h>
#include <stdlib.h>
// By Deepseek
void Legendre(int n, int m, double fi, double **pnm, double **dpnm);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // 检查输入参数
    if (nrhs != 3) {
        mexErrMsgIdAndTxt("Legendre:inputError", "需要3个输入参数：n, m, fi");
    }
    
    // 获取输入参数
    int n = (int)mxGetScalar(prhs[0]);
    int m = (int)mxGetScalar(prhs[1]);
    double fi = mxGetScalar(prhs[2]);
    
    // 检查m <= n
    if (m > n) {
        mexErrMsgIdAndTxt("Legendre:inputError", "m必须小于等于n");
    }
    
    // 为pnm和dpnm分配内存
    double **pnm = (double **)malloc((n+1)*sizeof(double *));
    double **dpnm = (double **)malloc((n+1)*sizeof(double *));
    for (int i = 0; i <= n; i++) {
        pnm[i] = (double *)malloc((m+1)*sizeof(double));
        dpnm[i] = (double *)malloc((m+1)*sizeof(double));
    }
    
    // 调用Legendre函数
    Legendre(n, m, fi, pnm, dpnm);
    
    // 创建MATLAB输出矩阵
    plhs[0] = mxCreateDoubleMatrix(n+1, m+1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(n+1, m+1, mxREAL);
    double *pnm_out = mxGetPr(plhs[0]);
    double *dpnm_out = mxGetPr(plhs[1]);
    
    // 复制数据到MATLAB矩阵（注意MATLAB是列优先）
    for (int j = 0; j <= m; j++) {
        for (int i = 0; i <= n; i++) {
            pnm_out[j*(n+1)+i] = pnm[i][j];
            dpnm_out[j*(n+1)+i] = dpnm[i][j];
        }
    }
    
    // 释放内存
    for (int i = 0; i <= n; i++) {
        free(pnm[i]);
        free(dpnm[i]);
    }
    free(pnm);
    free(dpnm);
}