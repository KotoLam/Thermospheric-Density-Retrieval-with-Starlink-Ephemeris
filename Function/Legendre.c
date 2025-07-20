#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// translate from the following matlab function, by Deepseek
// %--------------------------------------------------------------------------
// % Inputs:
// % n         maximum degree
// % m         maximum order
// % fi        angle [rad]
// %
// % Outputs:
// % pnm       normalized Legendre polynomial values
// % dpnm      normalized Legendre polynomial first derivative values
// %
// % Last modified:   2018/01/27   Meysam Mahooti
// %--------------------------------------------------------------------------
void Legendre(int n, int m, double fi, double **pnm, double **dpnm) {
    // Allocate memory for pnm and dpnm (n+1 x m+1 matrices)
    for (int i = 0; i <= n; i++) {
        pnm[i] = (double *)calloc(m + 1, sizeof(double));
        dpnm[i] = (double *)calloc(m + 1, sizeof(double));
    }

    pnm[0][0] = 1.0;
    dpnm[0][0] = 0.0;
    
    if (n >= 1 && m >= 1) {
        pnm[1][1] = sqrt(3.0) * cos(fi);
        dpnm[1][1] = -sqrt(3.0) * sin(fi);
    }

    // Diagonal coefficients
    for (int i = 2; i <= n; i++) {
        double factor = sqrt((2.0 * i + 1.0) / (2.0 * i));
        pnm[i][i] = factor * cos(fi) * pnm[i-1][i-1];
        dpnm[i][i] = factor * (cos(fi) * dpnm[i-1][i-1] - sin(fi) * pnm[i-1][i-1]);
    }

    // Horizontal first step coefficients
    for (int i = 1; i <= n; i++) {
        pnm[i][i-1] = sqrt(2.0 * i + 1.0) * sin(fi) * pnm[i-1][i-1];
        dpnm[i][i-1] = sqrt(2.0 * i + 1.0) * (cos(fi) * pnm[i-1][i-1] + sin(fi) * dpnm[i-1][i-1]);
    }

    // Horizontal second step coefficients
    int j = 0;
    int k = 2;
    while (1) {
        for (int i = k; i <= n; i++) {
            double factor = sqrt((2.0 * i + 1.0) / ((i - j) * (i + j)));
            double term1 = sqrt(2.0 * i - 1.0) * sin(fi) * pnm[i-1][j];
            double term2 = sqrt(((i + j - 1.0) * (i - j - 1.0)) / (2.0 * i - 3.0)) * pnm[i-2][j];
            pnm[i][j] = factor * (term1 - term2);

            double dterm1 = sqrt(2.0 * i - 1.0) * sin(fi) * dpnm[i-1][j];
            double dterm2 = sqrt(2.0 * i - 1.0) * cos(fi) * pnm[i-1][j];
            double dterm3 = sqrt(((i + j - 1.0) * (i - j - 1.0)) / (2.0 * i - 3.0)) * dpnm[i-2][j];
            dpnm[i][j] = factor * (dterm1 + dterm2 - dterm3);
        }
        j++;
        k++;
        if (j > m) {
            break;
        }
    }
}
