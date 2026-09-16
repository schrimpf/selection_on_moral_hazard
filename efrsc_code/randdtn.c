/* 
   Mex-function to generate random truncated normal variables.
   
   to compile: use make.m
*/

#include <math.h>
#include <string.h>
#include "mex.h" /* header file for matlab API */
#include "alcoa.h"

static RANDSTATE rs;
static int rsInit=0;

void mexFunction
(int nlhs, /* number of left hand side arguments */
 mxArray *plhs[], /* pointer to lhs*/
 int nrhs, /* number of rhs arguments*/
 const mxArray *prhs[]) /* pointer to rhs arguments*/
{
  double *mu, *sig, *lo, *hi;
  double *x;
  int n,i;
  if (nrhs!=4) {
    mexPrintf("nrhs=%d\n",nrhs);
    mxAssert(nrhs==4, "Wrong number of arguments");
  }
  /* Get inputs */
  mu = mxGetPr(prhs[0]);
  sig = mxGetPr(prhs[1]);
  lo = mxGetPr(prhs[2]);
  hi = mxGetPr(prhs[3]);
  n = mxGetNumberOfElements(prhs[0]);
  mxAssert(n==mxGetNumberOfElements(prhs[1]),"numel(sig)!=numel(mu)");
  mxAssert(n==mxGetNumberOfElements(prhs[2]),"numel(lo)!=numel(mu)");
  mxAssert(n==mxGetNumberOfElements(prhs[3]),"numel(hi)!=numel(mu)");
  if (rsInit==0) {
    set_time_seed(&rs);
    rsInit=1;
  }  
  plhs[0] = mxCreateDoubleMatrix(1,n,mxREAL);
  x = mxGetPr(plhs[0]);
  for (i=0;i<n;i++) {
    //mexPrintf("i=%d mu=%g sig=%g lo=%g hi=%g\n",i,mu[i],sig[i],lo[i],hi[i]);
    x[i] = Sample_Double_Truncated_Normal(mu[i],sig[i]*sig[i],lo[i],hi[i],&rs);
    //mexPrintf("lo=%g <? x=%g <? hi=%g\n",lo[i],x[i],hi[i]);
  }
  return;
}

