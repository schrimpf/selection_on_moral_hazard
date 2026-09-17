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
  double *muP;
  double *sigP; 
  int nsigP;
  double *muL;
  double *sigL;
  double *logpsi;

  double *xint;
  double *wint;
  double *deduct;
  double *maxoop;
  double *prem;
  double *avail;
  double *choice;
  double *meanMuL;
  double *lamlo;
  double *stdevi;

  int nint; 
  int maxIter;
  int T; /* number of time periods */
  int np;
  int multMH;
  
  int tries;
  int totover;
  /* output */
  double *logpout;
  double nbad;
  /* thread management */
  int tid;
  int start;
  int end;
};


void *sampleP(void *voidptr) 
{
  struct sampleArgs *sa = (struct sampleArgs*) voidptr;
  int i;
  sa->tries=0;
  sa->totover=0;
  sa->nbad=0;
  for(i=sa->start;i<sa->end;i++) {
    int accept = 0;
    int tries=0;
    double lpsi;
    double sp = sa->nsigP>1? sa->sigP[i]:sa->sigP[0];
    while (!accept) {
      int t; 
      accept = 1;
      lpsi = sa->muP[i] + Sample_Normal(0,1,&rs[sa->tid])*sp*sa->stdevi[i]; 
      for(t=0;t<sa->T && accept;t++) {
        if (sa->choice[t+i*sa->T]>=0) {
          double exv[5];
          int J = choiceFnEx(exp(sa->logomega[i]),exp(lpsi),sa->muL[i*sa->T+t],sa->sigL[i],sa->xint, sa->wint,sa->nint,
                             sa->deduct+t*sa->np+i*sa->T*sa->np, sa->maxoop+t*sa->np+i*sa->T*sa->np, 
                             sa->prem+t*sa->np+i*sa->T*sa->np, sa->np, copay,-sa->lamlo[i], 
                             sa->multMH, 
                             exv,exv,exv);        
          sa->totover += (J==0);
          accept = (((int) sa->choice[t+i*sa->T])== J);
        }
      }
      tries++;
      if (tries>sa->maxIter) {
        sa->nbad++;
        break;
      }
      //accept = 1;
    }
    sa->tries+= tries;
    if (tries<=sa->maxIter) {
      sa->logpout[i] = lpsi;
    }     
    //mexPrintf("%4d: %5.3g %5.3g\n      %5.3g %5.3g\n",i,rpm[0+i*4],rpm[1+i*4],rpm[2+i*4],rpm[3+i*4]);
    //mexPrintf("%4d: %d\n",
  }
  return NULL;
}

void mexFunction
(int nlhs, /* number of left hand side arguments */
 mxArray *plhs[], /* pointer to lhs*/
 int nrhs, /* number of rhs arguments*/
 const mxArray *prhs[]) /* pointer to rhs arguments*/
{
  struct sampleArgs sa, *ta;
  int tottries = 0,nbad=0,totover=0;
  int N;
  int i;
  int nthread;
  pthread_t *thrptr;
  if (nrhs!=19) {
    mexPrintf("nrhs=%d\n",nrhs);
    mxAssert(nrhs==19, "Wrong number of arguments");
  }
  /* Get inputs */
  sa.logomega = mxGetPr(prhs[0]);
  sa.loglambda = mxGetPr(prhs[1]);
  sa.muP = mxGetPr(prhs[2]);
  sa.sigP = mxGetPr(prhs[3]);
  sa.nsigP = mxGetNumberOfElements(prhs[3]);
  sa.muL = mxGetPr(prhs[4]);
  sa.sigL = mxGetPr(prhs[5]);
  sa.logpsi = mxGetPr(prhs[6]);
                
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
  mxAssert(N==mxGetNumberOfElements(prhs[14]),"lamlo N!=N");
  sa.stdevi = mxGetPr(prhs[15]);
  //mxAssert(sa.T==mxGetM(prhs[14]),"spend T!=N");
  sa.maxIter = (int) mxGetScalar(prhs[16]);
  nthread = (int) mxGetScalar(prhs[17]);
  sa.multMH = (int) mxGetScalar(prhs[18]);
  mxAssert(nthread<=MAXTHREADS,"too many threads\n");
  /* Create and get outputs */
  if (rsInit==0) {
    for (i=0;i<MAXTHREADS;i++) set_time_seed(rs+i);
    rsInit=1;
  }  
  plhs[0] = mxCreateDoubleMatrix(1,N,mxREAL);
  sa.logpout = mxGetPr(plhs[0]);
  memcpy(sa.logpout,sa.logpsi,sizeof(double)*N);
  ta = myCalloc(nthread,sizeof(struct sampleArgs));
  thrptr = myCalloc(nthread,sizeof(pthread_t));
  for(i=0;i<nthread;i++) {
    memcpy(ta+i,&sa,sizeof(struct sampleArgs));
    ta[i].tid = i;
    ta[i].start = i*N/nthread;
    if (i<nthread-1) ta[i].end = (i+1)*N/nthread;
    else ta[i].end = N;       
    pthread_create(thrptr+i,NULL, sampleP, (void *) (ta+i)); 
  }
  for (i=0;i<nthread;i++) {
    void *status; 
    pthread_join(thrptr[i],&status);
    tottries += ta[i].tries;
    totover += ta[i].totover;
    nbad += ta[i].nbad;
  }
  //mexPrintf("joined\n");
  myFree(ta);
  //mexPrintf("freed args\n");
  myFree(thrptr);
  { 
    const mxArray *v_arr = mexGetVariable("global","verbosity");
    int verbosity = v_arr ? (int) mxGetScalar(v_arr) : 0;
    if (verbosity>0) {
      mexPrintf("samplePsi: avg tries=%g, reached max %d times, overflow %d times (%.2lf)\n",
                ((double) tottries)/((double) N),nbad,totover,((double) totover)/((double) tottries));
    }
  } 
  return;
}

