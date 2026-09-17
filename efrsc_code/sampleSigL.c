/* 
   Mex-function to sample omega and lambda
   
   to compile: use make.m
*/

#include <math.h>
#include <string.h>
#include "mex.h" /* header file for matlab API */
#include "alcoa.h"
#include <omp.h>

#define copay 0.1

static RANDSTATE rs[MAXTHREADS];
static int rsInit=0;

struct sampleArgs {
  double *logomega;
  double *loglambda;
  double *shapei;
  double *scalei; 
  double *muL;
  double *sigL; // not needed but lazy...
  double *psi;

  double *xint;
  double *wint;
  double *deduct;
  double *maxoop;
  double *prem;
  double *avail;
  double *choice;
  double *lamlo;
  double varHi; 
  double varLo;

  int nint; 
  int maxIter;
  int T; // number of time periods 
  int np;
  int tottries, nbad;
  int multMH; 

  // output 
  double *sigLout;
  // thread management 
  int tid;
  int start;
  int end;
};

void mexFunction
(int nlhs, /* number of left hand side arguments */
 mxArray *plhs[], /* pointer to lhs*/
 int nrhs, /* number of rhs arguments*/
 const mxArray *prhs[]) /* pointer to rhs arguments*/
{
  /* rhs args */
  struct sampleArgs sa;
  int N; /* number of observations */
  int i;
  int tottries = 0, nbad = 0;
  int nthread;

  if (nrhs!=20) {
    mexPrintf("nrhs=%d\n",nrhs);
    mxAssert(nrhs==20, "Wrong number of arguments");
  }
  
  /* Get inputs */
  sa.logomega = mxGetPr(prhs[0]);
  sa.loglambda = mxGetPr(prhs[1]);
  sa.shapei = mxGetPr(prhs[2]);
  sa.scalei = mxGetPr(prhs[3]);
  sa.muL = mxGetPr(prhs[4]);
  sa.sigL = mxGetPr(prhs[5]);
  sa.psi = mxGetPr(prhs[6]);
                
  sa.nint = mxGetNumberOfElements(prhs[7]);
  sa.xint = mxGetPr(prhs[7]);
  sa.wint = mxGetPr(prhs[8]);
  sa.np = mxGetDimensions(prhs[9])[0];
  sa.T = mxGetDimensions(prhs[9])[1];
  N = mxGetDimensions(prhs[9])[2];
  sa.deduct = mxGetPr(prhs[9]);
  sa.maxoop = mxGetPr(prhs[10]);
  sa.prem = mxGetPr(prhs[11]);
  sa.avail = mxGetPr(prhs[12]);
  sa.choice = mxGetPr(prhs[13]);
  mxAssert(mxGetNumberOfElements(prhs[13])==N*sa.T,"deduct and/or choice wrong size");
  sa.lamlo = mxGetPr(prhs[14]);
  sa.maxIter = (int) mxGetScalar(prhs[15]);
  nthread = (int) mxGetScalar(prhs[16]);
  mxAssert(nthread<=MAXTHREADS,"too many threads\n");
  /* Create and get outputs */
  if (rsInit==0) {
    for (i=0;i<MAXTHREADS;i++) set_time_seed(rs+i);
    rsInit=1;
  }  
  // limits on variance
  sa.varLo = mxGetScalar(prhs[17]);
  sa.varHi = mxGetScalar(prhs[18]);
  sa.multMH = (int) mxGetScalar(prhs[19]);
  /* Create and get outputs */
  plhs[0] = mxCreateDoubleMatrix(1,N,mxREAL);
  sa.sigLout = mxGetPr(plhs[0]);
  memcpy(sa.sigLout,sa.sigL,sizeof(double)*N);
  
  if (nthread > 0) {
    omp_set_num_threads(nthread);
  }

  #pragma omp parallel for schedule(dynamic, 16) reduction(+:tottries) reduction(+:nbad)
  for (i = 0; i < N; i++) {
    int tid = omp_get_thread_num();
    int accept = 0;
    int tries = 0;
    double sl;
    while (!accept) {
      int t; 
      accept = 1;
      sl = sqrt(1.0/Sample_Truncated_Gamma(sa.shapei[i], 1.0/sa.scalei[i], 1/sa.varHi, 1/sa.varLo, &rs[tid]));
      for (t = 0; t < sa.T && accept; t++) {
        if (sa.choice[t+i*sa.T] >= 0) {
          double exv[5];
          int J = choiceFnEx(exp(sa.logomega[i]), sa.psi[i], sa.muL[i*sa.T+t], sl, sa.xint, sa.wint, sa.nint,
                             sa.deduct+t*sa.np+i*sa.T*sa.np, sa.maxoop+t*sa.np+i*sa.T*sa.np, 
                             sa.prem+t*sa.np+i*sa.T*sa.np, sa.np, copay, -sa.lamlo[i],
                             sa.multMH, 
                             exv, exv, exv);        
          accept = (((int) sa.choice[t+i*sa.T]) == J);
          if (!accept) break;
        }
      }
      tries++;
      if (tries > sa.maxIter) {
        nbad++;
        break;
      }
    }
    if (tries <= sa.maxIter) {
      sa.sigLout[i] = sl;
    }
    tottries += tries;
  }
  { 
    int verbosity = (int) mxGetScalar(mexGetVariable("global","verbosity"));
    if (verbosity>0) {
      
      mexPrintf("sampleSigL: average number of tries=%6.3g, reached maximum %d times\n",
                ((double) tottries)/((double) N),nbad);
    }
  }
  return;
}
