/* 
   Mex-function that given latent variables and plan characteristics,
   calculates each person's choices, expected utility, and expected
   spending.  

   Input: 
    - logOmega = N x 1
    - loglambda = T x N
    - lamlo (kappa) = N x 1
    - sigma_omega (not used, only present for consistency with other mex functions)
    - mu_lambda = T x N
    - sigma_lambda = N x 1
    - psi = N x 1
    - xint = nint x 1 (integration points)
    - wint = nint x 1 (integration weights)
    - deduct = nplan x T x N 
    - maxoop = nplan x T x N 
    - premium = nplan x T x N
    - choice = T x N (not used, only present for consistency with other mex functions)
    - spend = T x N (not used, only present for consistency with other mex functions)
    - maxIter = scalar (not used, only present for consistency with other mex functions)
   Output: 
    - choiceOut = T x N, calculated choices
    - expectedValue = T x N , expected utility
    - expectedSpend = T x N, expected spending 
    - expectedOop = T x N, expected out of pocket spending
    
*/

#include <math.h>
#include <string.h>
#include <omp.h>
#include "mex.h" /* header file for matlab API */
#include "alcoa.h"

//#define copay 0.1

void mexFunction
(int nlhs, /* number of left hand side arguments */
 mxArray *plhs[], /* pointer to lhs*/
 int nrhs, /* number of rhs arguments*/
 const mxArray *prhs[]) /* pointer to rhs arguments*/
{
  /* rhs args */
  double *logomega;
  double *loglambda;
  double *lamlo;
  double sigO;
  double *muL;
  double *sigL;
  double *psi;

  double *xint;
  double *wint;
  double *deduct;
  double *maxoop;
  double *prem;
  double *avail;
  double *choice;
  double *spend;
  double copay;
  
  int nint; 
  int maxIter;
  int N; /* number of observations */
  int T; /* number of time periods */
  int t,i,np,nplan;
  /* output */
  double *choiceOut;
  double *expectedValue, *expectedSpend,*expectedOop;
  double *ev, *es, *eo;
  if (nrhs!=16 && nrhs!=17 && nrhs!=18) {
    mexPrintf("nrhs=%d\n",nrhs);
    mxAssert(nrhs==16 || nrhs==17, "Wrong number of arguments");
  }
  if (nrhs<17) copay = 0.1;
  else copay = mxGetScalar(prhs[16]);
  if (nrhs<18) setNegLambdaZeroUtil(0);
  else setNegLambdaZeroUtil((int) mxGetScalar(prhs[17]));

  /* Get inputs */
  logomega = mxGetPr(prhs[0]);
  loglambda = mxGetPr(prhs[1]);
  lamlo = mxGetPr(prhs[2]);
  sigO = mxGetScalar(prhs[3]);
  muL = mxGetPr(prhs[4]);
  sigL = mxGetPr(prhs[5]);
  psi = mxGetPr(prhs[6]);
                
  nint = mxGetNumberOfElements(prhs[7]);
  xint = mxGetPr(prhs[7]);
  wint = mxGetPr(prhs[8]);
  np = mxGetDimensions(prhs[9])[0];
  T = mxGetDimensions(prhs[9])[1];
  N = mxGetDimensions(prhs[9])[2];
  deduct = mxGetPr(prhs[9]);
  maxoop = mxGetPr(prhs[10]);
  prem = mxGetPr(prhs[11]);
  avail = mxGetPr(prhs[12]);
  choice = mxGetPr(prhs[13]);
  //mxAssert(mxGetNumberOfElements(prhs[13])==N*T,"deduct and/or choice wrong size");
  spend = mxGetPr(prhs[14]);
  //mxAssert(N==mxGetN(prhs[14]),"spend N!=N");
  //mxAssert(T==mxGetM(prhs[14]),"spend T!=N");
  maxIter = (int) mxGetScalar(prhs[15]);
  /* Create and get outputs */
  plhs[0] = mxCreateDoubleMatrix(T,N,mxREAL);
  if (nlhs>1) {
    int dims[3];
    dims[0] = np;;
    dims[1] = T;
    dims[2] = N;
    plhs[1] = mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
    expectedValue = mxGetPr(plhs[1]);
  } 
  if (nlhs>2) {
    int dims[3];
    dims[0] = np;;
    dims[1] = T;
    dims[2] = N;
    plhs[2] = mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
    expectedSpend = mxGetPr(plhs[2]);
  } 
  if (nlhs>3) {
    int dims[3];
    dims[0] = np;;
    dims[1] = T;
    dims[2] = N;
    plhs[3] = mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
    expectedOop = mxGetPr(plhs[3]);
  } 
  choiceOut = mxGetPr(plhs[0]);
  //omp_set_num_threads(20);
#pragma omp parallel private(t,ev,es,eo) 
  {
    if (nlhs<2) ev = myCalloc(np,sizeof(double));
    if (nlhs<3) es = myCalloc(np,sizeof(double));
    if (nlhs<4) eo = myCalloc(np,sizeof(double));
#pragma omp for 
    for(i=0;i<N;i++) {
      for(t=0;t<T;t++) {
        if (nlhs>1) ev = expectedValue+i*np*T+t*np;
        if (nlhs>2) es = expectedSpend+i*np*T+t*np;
        if (nlhs>3) eo = expectedOop+i*np*T+t*np;
        //mexPrintf("%d %d\n",i,t);
        choiceOut[i*T + t] = (double) choiceFnEx(exp(logomega[i]),psi[i],muL[i*T+t],sigL[i], 
                                                 xint, wint,nint,
                                                 deduct+i*np*T+t*np, maxoop+i*np*T+t*np, 
                                                 prem+i*np*T+t*np, np,copay, -lamlo[i], 0, 
                                                 ev, es, eo);
      }
    }
    if (nlhs<4) myFree(eo);
    if (nlhs<3) myFree(es);
    if (nlhs<2) myFree(ev);
  } // end parallel section
  return;
}


