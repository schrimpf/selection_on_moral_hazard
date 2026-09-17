/*
   Mex-function to sample omega and lambda

   to compile: use make.m
*/

#include <math.h>
#include <string.h>
#include "mex.h" /* header file for matlab API */
#include "alcoa.h"
//#include "lapack.h"
//#include "acml.h" // AMD core math library, includes lapack
#include <pthread.h>

// two lapack functions needed (can comment this out and include
// header file instead, but then need header in your include path)
//extern void dpotrf(char uplo, int n, double *a, int lda, int *info);
//extern void dpotri(char uplo, int n, double *a, int lda, int *info);

// above prototypes assume you have a C interface for lapack, if you
// want to use the fortran interface directly, use the following
// declarations, and change the calls below accordingly
extern void dpotrf_(char *uplo, int *n, double *a, int *lda, int *info, int uplo_len);
extern void dpotri_(char *uplo, int *n, double *a, int *lda, int *info, int uplo_len);

#define copay 0.1

static RANDSTATE rs[MAXTHREADS];
static int rsInit=0;

struct sampleArgs {
  double *logomega;
  double *loglambda;
  double *logpsi;
  double *muL;
  double *xbO;
  double *xbP;
  double *xbL;
  double *sigL;
  double *Sig;
  int nSig;

  double *xint;
  double *wint;
  double *deduct;
  double *maxoop;
  double *prem;
  double *avail;
  double *choice;
  double *lamlo;
  double *stdevi;
  double *meanMuL;

  int nint;
  int maxIter;
  int T; /* number of time periods */
  int np;
  int multMH;

  int tries;
  int totover;
  /* output */
  double *muLout;
  double nbad;
  /* thread management */
  int tid;
  int start;
  int end;
};


void *sampleMu(void *voidptr)
{
  double *Vi,*Ci;
  struct sampleArgs *sa = (struct sampleArgs*) voidptr;
  int i;
  Vi = myCalloc((2+sa->T)*(2+sa->T),sizeof(double));
  Ci = myCalloc((2+sa->T),sizeof(double));
  sa->tries=0;
  sa->totover=0;
  sa->nbad=0;
  for(i=sa->start;i<sa->end;i++) {
    int t;
    int tries=0;
    double ml;
    int Tp2 = sa->T+2;
    double *Sig = sa->nSig==1? sa->Sig:(sa->Sig+i*Tp2*Tp2);
    for(t=0;t<sa->T;t++) {
      int accept = 0;
      double sig2mul;
      double meanmul;
      int t1,t2,j1,j2;
      char uplo = 'U';
      int info,N=Tp2;
      // create V_{i,-t} -- first Sig_{-t}
      for(j1=0,t1=0;t1<(Tp2);t1++) {
        for(j2=0,t2=0;t2<(Tp2);t2++) {
          if (t1!=(2+t) && t2!=(t+2)) {
            Vi[j1*(Tp2) + j2] = Sig[t1*(Tp2)+t2];
            j2++;
          }
        }
        if (t1!=(t+2)) {
          Ci[j1] = Sig[t1*(Tp2)+(t+2)];
          j1++;
        }
      } // now the part for lambda
      for (t1=0;t1<Tp2;t1++) {
        Vi[(1+sa->T)*(Tp2)+t1] = Sig[(Tp2-1)*(Tp2)+t1];
        Vi[t1*(Tp2)+(1+sa->T)] = Sig[(Tp2-1)*(Tp2)+t1];
      }
      Vi[(1+sa->T)*(Tp2)+(1+sa->T)] = Sig[(t+2)*(Tp2)+(t+2)]+sa->sigL[i]*sa->sigL[i];
      Ci[1+sa->T] = Sig[2*(Tp2)+2];
      /*
        if (sa->tid==0 && i<=10) {
        mexPrintf("i=%d t=%d\nSig = \n",i,t);
        for(t1=0;t1<Tp2;t1++) {
          for(t2=0;t2<Tp2;t2++) {
            mexPrintf("  %6.3g",Sig[t1*Tp2+t2]);
          }
          mexPrintf("\n");
        }
        mexPrintf("Vi = \n");
        for(t1=0;t1<Tp2;t1++) {
          for(t2=0;t2<Tp2;t2++) {
            mexPrintf("  %6.3g",Vi[t1*Tp2+t2]);
          }
          mexPrintf("\n");
        }
        mexPrintf("Ci = \n");
        for(t1=0;t1<Tp2;t1++) {
          mexPrintf("  %6.3g",Ci[t1]);
        }
        mexPrintf("\n----------------------\n");
      }
      */
      // use lapack  to compute inv(V_i)
      //dpotrf(uplo,N,Vi,N,&info); // Vi = cholesky decomp(Vi)
      // or call fortran version directly
      dpotrf_(&uplo,&N,Vi,&N,&info,1);

      /*
      if (sa->tid==0 && i<=10) {
        mexPrintf("chol(Vi) = (info=%d N=%d)\n",info,(int) N);
        for(t1=0;t1<Tp2;t1++) {
          for(t2=0;t2<Tp2;t2++) {
            mexPrintf("  %6.3g",Vi[t1*Tp2+t2]);
          }
          mexPrintf("\n");
        }
      }
      */
      //dpotri(uplo,N,Vi,N,&info); // Vi = inv(Vi) (upper triangle)
      dpotri_(&uplo,&N,Vi,&N,&info,1);
      // fill in lower triangle
      for(j1=0;j1<Tp2;j1++) {
        for(j2=0;j2<j1;j2++) Vi[j1 + Tp2*j2] = Vi[j2 + Tp2*j1];
      }
      /*
        if (sa->tid==0 && i<=10) {
        mexPrintf("iVi = (info=%d N=%d)\n",info,(int) N);
        for(t1=0;t1<Tp2;t1++) {
          for(t2=0;t2<Tp2;t2++) {
            mexPrintf("  %6.3g",Vi[t1*Tp2+t2]);
          }
          mexPrintf("\n");
        }
      }
      */
      sig2mul = Sig[2*(Tp2)+2];
      //if (sa->tid==0) mexPrintf("sig2mul=%g, info=%d\n",sig2mul,info);
      meanmul  = sa->xbL[i*sa->T+t];
      for(j1=0;j1<Tp2;j1++) {
        meanmul += Ci[j1]*Vi[0 + Tp2*j1]*(sa->logomega[i]-sa->xbO[i]) +
          Ci[j1]*Vi[1 + Tp2*j1]*(sa->logpsi[i]-sa->xbP[i]);
        for(j2=2,t1=0;t1<sa->T;t1++) {
          if (t1!=t) meanmul += Ci[j1]*Vi[j2++ + Tp2*j1]*(sa->muLout[i*sa->T+t1]-sa->xbL[i*sa->T+t1]);
        }
        meanmul += Ci[j1]*Vi[Tp2-1+Tp2*j1]*(sa->loglambda[i*sa->T+t]-sa->xbL[i*sa->T+t]);
        for(j2=0;j2<Tp2;j2++) {
          sig2mul -= Ci[j1]*Vi[j1*Tp2+j2]*Ci[j2];
          //if (sa->tid==0) mexPrintf("sig2mul=%g j1=%d j2=%d\n",sig2mul,j1,j2);
        }
      }
      //if (sa->tid==0 && i<10) mexPrintf("%d: meanmul=%g sig2mul=%g\n",i,meanmul,sig2mul);
      tries = 0;
      while (!accept) {
        accept = 1;
        ml = Sample_Normal(meanmul,sig2mul*sa->stdevi[i]*sa->stdevi[i],
                           &rs[sa->tid]);
        if (sa->choice[t+i*sa->T]>=0) {
          double exv[5];
          int J = choiceFnEx(exp(sa->logomega[i]),exp(sa->logpsi[i]), ml, sa->sigL[i],sa->xint, sa->wint,sa->nint,
                             sa->deduct+t*sa->np+i*sa->T*sa->np, sa->maxoop+t*sa->np+i*sa->T*sa->np,
                             sa->prem+t*sa->np+i*sa->T*sa->np, sa->np, copay,-sa->lamlo[i], sa->multMH,
                             exv,exv,exv);
          sa->totover += (J==0);
          accept = (((int) sa->choice[t+i*sa->T])== J);
        }
        tries++;
        if (tries>sa->maxIter) {
          sa->nbad++;
          break;
        }
      } // end while
      //accept = 1;
      sa->tries+= tries;
      if (tries<=sa->maxIter) {
        sa->muLout[i*sa->T+t] = ml;
      }
    } // end t
    //mexPrintf("%4d: %5.3g %5.3g\n      %5.3g %5.3g\n",i,rpm[0+i*4],rpm[1+i*4],rpm[2+i*4],rpm[3+i*4]);
    //mexPrintf("%4d: %d\n",
  }
  myFree(Ci);
  myFree(Vi);
  return(NULL);
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
  if (nrhs!=21) {
    mexPrintf("nrhs=%d\n",nrhs);
    mxAssert(nrhs==21, "Wrong number of arguments");
  }
  /* Get inputs */
  sa.logomega = mxGetPr(prhs[0]);
  sa.loglambda = mxGetPr(prhs[1]);
  sa.logpsi = mxGetPr(prhs[2]);
  sa.muL = mxGetPr(prhs[3]);
  sa.xbO = mxGetPr(prhs[4]);
  sa.xbP = mxGetPr(prhs[5]);
  sa.xbL = mxGetPr(prhs[6]);
  sa.sigL = mxGetPr(prhs[7]);
  sa.Sig = mxGetPr(prhs[8]);
  sa.nSig = mxGetNumberOfElements(prhs[8]);

  sa.nint = mxGetNumberOfElements(prhs[9]);
  sa.xint = mxGetPr(prhs[9]);
  sa.wint = mxGetPr(prhs[10]);
  sa.np = mxGetDimensions(prhs[11])[0];
  sa.T = mxGetDimensions(prhs[11])[1];
  sa.nSig /= ((sa.T+2)*(sa.T+2));
  N = mxGetDimensions(prhs[11])[2];
  sa.deduct = mxGetPr(prhs[11]);
  sa.maxoop = mxGetPr(prhs[12]);
  sa.prem = mxGetPr(prhs[13]);
  sa.avail = mxGetPr(prhs[14]);
  sa.choice = mxGetPr(prhs[15]);
  mxAssert(mxGetNumberOfElements(prhs[15])==N*sa.T,"deduct and/or choice wrong size");
  sa.lamlo = mxGetPr(prhs[16]);
  sa.stdevi = mxGetPr(prhs[17]);
  //mxAssert(N==mxGetN(prhs[16]),"spend N!=N");
  //mxAssert(sa.T==mxGetM(prhs[16]),"spend T!=N");
  sa.maxIter = (int) mxGetScalar(prhs[18]);
  nthread = (int) mxGetScalar(prhs[19]);
  sa.multMH = (int) mxGetScalar(prhs[20]);
  mxAssert(nthread<=MAXTHREADS,"too many threads\n");
  /* Create and get outputs */
  if (rsInit==0) {
    for (i=0;i<MAXTHREADS;i++) set_time_seed(rs+i);
    rsInit=1;
  }
  plhs[0] = mxCreateDoubleMatrix(sa.T,N,mxREAL);
  sa.muLout = mxGetPr(plhs[0]);
  memcpy(sa.muLout,sa.muL,sizeof(double)*N*sa.T);
  ta = myCalloc(nthread,sizeof(struct sampleArgs));
  thrptr = myCalloc(nthread,sizeof(pthread_t));
  for(i=0;i<nthread;i++) {
    memcpy(ta+i,&sa,sizeof(struct sampleArgs));
    ta[i].tid = i;
    ta[i].start = i*N/nthread;
    if (i<nthread-1) ta[i].end = (i+1)*N/nthread;
    else ta[i].end = N;
    pthread_create(thrptr+i,NULL, sampleMu, (void *) (ta+i));
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
  int verbosity = (int) mxGetScalar(mexGetVariable("global","verbosity"));
  if (verbosity>0) {
    mexPrintf("sampleMu: avg tries=%g, reached max %d times, overflow %d times (%.2lf)\n",
              ((double) tottries)/((double) N),nbad,totover,((double) totover)/((double) tottries));
  }
  return;
}
