/* 
   Mex-function to sample omega and lambda
   
   to compile: use make.m
*/

#include <math.h>
#include <string.h>
#include "mex.h" /* header file for matlab API */
#include "alcoa.h"
#include <pthread.h>

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

void *sampleSigL(void *voidptr)
{
  struct sampleArgs *sa = (struct sampleArgs*) voidptr;
  int i;
  sa->tottries=0;
  sa->nbad=0;
  for(i=sa->start;i<sa->end;i++) {
    int accept=0;
    int tries=0;
    double sl;
    while (!accept) {
      int t; 
      accept=1;
      sl = sqrt(1.0/Sample_Truncated_Gamma(sa->shapei[i],1.0/sa->scalei[i],1/sa->varHi,1/sa->varLo,&rs[sa->tid]));
      //mexPrintf("%g %g\n",sl,sigL[i]);
      for(t=0;t<sa->T && accept;t++) {
        if (sa->choice[t+i*sa->T]>=0) {
          double exv[5];
          int J = choiceFnEx(exp(sa->logomega[i]),sa->psi[i],sa->muL[i*sa->T+t],sl,sa->xint, sa->wint,sa->nint,
                             sa->deduct+t*sa->np+i*sa->T*sa->np, sa->maxoop+t*sa->np+i*sa->T*sa->np, 
                             sa->prem+t*sa->np+i*sa->T*sa->np, sa->np, copay,-sa->lamlo[i],
                             sa->multMH, 
                             exv,exv,exv);        
          //mexPrintf("%d: %d %d\n",i,choice[t+i*T],J);
          accept = (((int) sa->choice[t+i*sa->T])== J);
          if (!accept) break;
        }
      }
      tries++;
      if (tries>sa->maxIter) {
        sa->nbad++;
        //mexPrintf("obs %d: too many rejections\n",i);
        break;
      }
      //accept = 1;
    }
    if (tries<=sa->maxIter) {
      sa->sigLout[i] = sl;
    }
    sa->tottries += tries;
  }    
  return NULL;
}

void mexFunction
(int nlhs, /* number of left hand side arguments */
 mxArray *plhs[], /* pointer to lhs*/
 int nrhs, /* number of rhs arguments*/
 const mxArray *prhs[]) /* pointer to rhs arguments*/
{
  /* rhs args */
  struct sampleArgs sa, *ta;
  int N; /* number of observations */
  int i;
  int tottries = 0, nbad=0;
  int nthread;
  pthread_t *thrptr;

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
  sa.multMH = (int) mxGetScalar(prhs[18]);
  /* Create and get outputs */
  plhs[0] = mxCreateDoubleMatrix(1,N,mxREAL);
  sa.sigLout = mxGetPr(plhs[0]);
  memcpy(sa.sigLout,sa.sigL,sizeof(double)*N);
  
  // spawn threads
  ta = myCalloc(nthread,sizeof(struct sampleArgs));
  thrptr = myCalloc(nthread,sizeof(pthread_t));
  for(i=0;i<nthread;i++) {
    memcpy(ta+i,&sa,sizeof(struct sampleArgs));
    ta[i].tid = i;
    ta[i].start = i*N/nthread;
    if (i<nthread-1) ta[i].end = (i+1)*N/nthread;
    else ta[i].end = N;       
    pthread_create(thrptr+i,NULL, sampleSigL, (void *) (ta+i)); 
  }
  // join threads
  for (i=0;i<nthread;i++) {
    void *status; 
    pthread_join(thrptr[i],&status);
    tottries += ta[i].tottries;
    nbad += ta[i].nbad;
  }
  //mexPrintf("joined\n");
  myFree(ta);
  //mexPrintf("freed args\n");
  myFree(thrptr);
  { 
    int verbosity = (int) mxGetScalar(mexGetVariable("global","verbosity"));
    if (verbosity>0) {
      
      mexPrintf("sampleSigL: average number of tries=%6.3g, reached maximum %d times\n",
                ((double) tottries)/((double) N),nbad);
    }
  }
  return;
}
