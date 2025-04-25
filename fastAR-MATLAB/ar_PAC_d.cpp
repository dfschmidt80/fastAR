#include "math.h"
#include "mex.h"

// Call with: alpha = PACF_d_cpp(rho)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Check that the function is called properly
    if (nrhs != 2) {
        mexErrMsgTxt("Two input arguments required.");
    } 
    if (nlhs > 1) {
        mexErrMsgTxt("Too many output arguments.");
    }
    
/*
function alpha = PACF_d(rc, k)

p = length(rc);
alpha = zeros(1, p);

%% Initialise
if (k > 1)
    alpha(1) = -rc(1);
else
    alpha(1) = -1;
end

%% Recur ...
for j = 2:p,
    if (k > j)
        alpha(1:j) = [alpha(1:j-1) - fliplr(alpha(1:j-1))*rc(j), -rc(j)];
    elseif (k == j)
        alpha(1:j) = [-fliplr(alpha(1:j-1)), -1];
    else
        alpha(1:j) = [alpha(1:j-1) - fliplr(alpha(1:j-1))*rc(j), 0];
    end
end
*/    

    int i, j;
    
    const mwSize *dims_rho = mxGetDimensions(prhs[0]);
    int p = dims_rho[1];
    
    //mexPrintf("%d\n",p);
    
    double *rho = mxGetPr(prhs[0]);
    int k = (int)(*mxGetPr(prhs[1]));
    k--;
    //mexPrintf("k=%d\n",k);
    /*for (int i = 0; i<p; i++)
    {
        mexPrintf("%f, ", rho[i]);
    }*/
    
    double *a = new double[p];
    plhs[0] = mxCreateDoubleMatrix(1, p, mxREAL);
    double *alpha = mxGetPr(plhs[0]);
    
    // Initialise
    if (k > 0)
        alpha[0] = -rho[0];
    else
        alpha[0] = -1;
    
    // Recur ...
    for (j = 1; j<p; j++)
    {
        // "fliplr"
        for (i = 0; i<j; i++)
            a[i] = alpha[j-i-1];
        
        if (k > j)
        {
            for (i = 0; i<j; i++)
            {
                alpha[i] = alpha[i] - a[i]*rho[j];
            }
            alpha[j] = -rho[j];
        }
        else if (k == j)
        {
            for (i = 0; i<j; i++)
                alpha[i] = -a[i];
            alpha[j] = -1;
        }
        else 
        {
            for (i = 0; i<j; i++)
            {
                alpha[i] = alpha[i] - a[i]*rho[j];
            }
            alpha[j] = 0;
        }
    }
    
    delete a;
    
    
    //plhs[0] = mxCreateDoubleMatrix(p, 1, mxREAL);
    
    //double *alpha = mxGetPr(plhs[0]);
    
    // Check matrix dimensions; X~2D matrix and Y~1D matrix
/*    int numdims_X = mxGetNumberOfDimensions(prhs[0]);
    int numdims_y = mxGetNumberOfDimensions(prhs[1]);
    int numdims_delta = mxGetNumberOfDimensions(prhs[2]);
    int numdims_nsteps = mxGetNumberOfDimensions(prhs[3]);
    int numdims_SavePath = mxGetNumberOfDimensions(prhs[4]);
    
    const mwSize *dims_X = mxGetDimensions(prhs[0]);
    const mwSize *dims_y = mxGetDimensions(prhs[1]);
    const mwSize *dims_delta = mxGetDimensions(prhs[2]);
    const mwSize *dims_nsteps = mxGetDimensions(prhs[3]);
    const mwSize *dims_SavePath = mxGetDimensions(prhs[4]);
        
    if(numdims_X != 2) {
        mexErrMsgTxt("Argument 'X' must be a 2D matrix.");
    }
    if(numdims_y != 2) {
        mexErrMsgTxt("Argument 'y' must be a column vector.");
    }
    if(numdims_delta != 2) {
        mexErrMsgTxt("Argument 'delta' must be a scalar.");
    }  
    if(numdims_nsteps != 2) {
        mexErrMsgTxt("Argument 'nsteps' must be a scalar.");
    }
    if(numdims_SavePath != 2) {
        mexErrMsgTxt("Argument 'SavePath' must be a scalar.");
    }  
    
    if(dims_y[1] != 1){
        mexErrMsgTxt("Argument 'y' must be a column vector.");
    }
    if(dims_delta[0] != 1 || dims_delta[1] != 1){
        mexErrMsgTxt("Argument 'delta' must be a scalar.");
    }
    if(dims_nsteps[0] != 1 || dims_nsteps[1] != 1){
        mexErrMsgTxt("Argument 'nsteps' must be a scalar.");
    }
    if(dims_SavePath[0] != 1 || dims_SavePath[1] != 1){
        mexErrMsgTxt("Argument 'nsteps' must be a scalar.");
    }
    if(dims_X[0] != dims_y[0]) {
        mexErrMsgTxt("'X' and 'y' have incompatible dimensions.");
    }

    // get the data sizes
    mwSize n = dims_X[0];                   // number of samples
    mwSize d = dims_X[1];                   // number of regressors
    
    double delta = *mxGetPr(prhs[2]);       // step-size
    double nsteps = *mxGetPr(prhs[3]);      // number of steps
    
    // Create output matrices for beta, beta0 and likelihood paths
    bool SavePath = (*mxGetPr(prhs[4])>0);
    if (SavePath)
        plhs[0] = mxCreateDoubleMatrix(d, d+1, mxREAL);
    else
        plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
            
    plhs[1] = mxCreateDoubleMatrix(d+1, 1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(d+1, 1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(d+1, NUM_CRITERIA, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[5] = mxCreateDoubleMatrix(d, NUM_CRITERIA, mxREAL);
    
    double *bpath=mxGetPr(plhs[0]);
    double *b0path=mxGetPr(plhs[1]);
    double *Lpath=mxGetPr(plhs[2]);
    double *scores=mxGetPr(plhs[3]);
    double *PathSize=mxGetPr(plhs[4]);
    double *beta_best=mxGetPr(plhs[5]);
 
    // call GPS
    double *X=mxGetPr(prhs[0]);
    double *y=mxGetPr(prhs[1]);
    
    *PathSize = GPS(X,y,delta,nsteps,(long)n,(long)d,bpath,b0path,Lpath,scores,beta_best,SavePath);*/
    
    return;
}