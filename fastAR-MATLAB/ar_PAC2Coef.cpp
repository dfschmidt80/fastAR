#include "math.h"
#include "mex.h"

// Call with: alpha = PACF_d_cpp(rho)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Check that the function is called properly
    if (nrhs != 1) {
        mexErrMsgTxt("One input argument required.");
    } 
    if (nlhs > 1) {
        mexErrMsgTxt("Too many output arguments.");
    }

    int i, j;
    
    const mwSize *dims_rho = mxGetDimensions(prhs[0]);
    int nc = dims_rho[0];
    int nr = dims_rho[1];
    int p;
    
    // Figure out what type of input vector has been passed
    if ((nc != 1 && nr != 1) && !(nr == 0 && nc == 0))
    {
        mexErrMsgTxt("Input argument must be either a row or column vector.");
    }
    if (nc == 1)
    {
        p = nr;
    }
    else
    {
        p = nc;
    }
    
    //mexPrintf("%d\n",p);
    
    double *rho = mxGetPr(prhs[0]);
    
    double *a = new double[p];
    plhs[0] = mxCreateDoubleMatrix(1, p+1, mxREAL);
    double *alpha = mxGetPr(plhs[0]);
    
    // If empty vector passed, just return an empty vector
    if (p == 0)
    {
        alpha[0] = 1;
        return;
    }
    
    // Initialise
    alpha[0] = 1;
    alpha[1] = -rho[0];
    
    // Recur ...
    for (j = 1; j<p; j++)
    {
        // "fliplr"
        for (i = 0; i<j; i++)
            a[i] = alpha[j-i];
        
        for (i = 0; i<j; i++)
        {
            alpha[i+1] = alpha[i+1] - a[i]*rho[j];
        }
        alpha[j+1] = -rho[j];
    }
    
    delete a;
    
    return;
}