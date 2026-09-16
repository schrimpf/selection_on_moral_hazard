/* 
   Mex-function to sample kappa
   
   to compile: make.m
*/

#include <math.h>
#include <string.h>
#include <time.h>
#include "mex.h" /* header file for matlab API */
#include "alcoa.h"
#include <omp.h>
#include <signal.h>


//#define normcdf(x) 0.5*(1+erf(x/sqrt(2)));

double sampleLamLo(double lamlo0, double *lam, double *muL, double sigL, double mull, 
                   double sigll, int T, RANDSTATE *rs, int *try, int *accept,
                   const double *choice, const double omega, const double psi, 
                   const double *xint, const double *wint, const int nint,
                   const double *deduct, const double *maxoop, const double *prem, const int nplan,
                   const double copay, const double hi, const int multMH)
{
  const int niter = 5;
  int iter;
  double ll = lamlo0;
  double ld = 0; //-0.5*(ll-mull)*(ll-mull)/(sigll*sigll);  (because need P(accept) = P(new) P(draw old|new)/(P(old) P(draw new|old))
  double minll = -HUGE_VAL;
  for(int t=0;t<T;t++) {
    minll = minll < -lam[t]? (-lam[t]):minll;
    double x = (log(lam[t]+ll)-muL[t])/sigL;
    ld += -0.5*x*x - log(lam[t]+ll);    
  }
  for(iter=0;iter<niter;iter++) {
    double lltry = Sample_Double_Truncated_Normal(mull,sigll*sigll,minll, hi, rs);
    double ldtry = 0; //-0.5*(lltry-mull)*(lltry-mull)/(sigll*sigll);
    for(int t=0;t<T;t++) {
      double x = (log(lam[t]+lltry)-muL[t])/sigL;
      if (lltry+lam[t]<=0) ldtry += -HUGE_VAL;
      else ldtry += -0.5*x*x - log(lam[t]+lltry);
    }
    double u = Sample_Uniform(0,1,rs);
    (*try)++;
    //mexPrintf("u=%g ll=%g lltry=%g ld=%g ldtry=%g p=%g\n",
    //          u,ll,lltry,ld,ldtry,exp(ldtry-ld));
    if (u<exp(ldtry-ld)) {
      int good =1;
      for (int t=0;t<T;t++) {
        if (choice[t]>0) {
          double nothing1[10],nothing2[10];
          int c = choiceFnEx(omega,psi,muL[t],sigL, 
                             xint, wint,nint,
                             deduct+t*nplan, maxoop+t*nplan, 
                             prem+t*nplan, nplan,copay,-lltry,multMH,
                             nothing1,nothing2,nothing1);
          if (c!=(int) choice[t]) { good=0; break; }
        }
      }
      if (good) {
        ll = lltry;
        ld = ldtry;
        (*accept)++;
      }
    }
  }
  return(ll);
}

void mexFunction
(int nlhs, /* number of left hand side arguments */
 mxArray *plhs[], /* pointer to lhs*/
 int nrhs, /* number of rhs arguments*/
 const mxArray *prhs[]) /* pointer to rhs arguments*/
{
  /* rhs args */
  int j,N, T;
  double *lam;
  double *muL;
  double *sigL;
  double *mull;
  double sigll;
  double *lamlo0;
  double *choice,*omega, *psi, *xint, *wint, *deduct, *maxoop, *prem;
  double *stdevi;
  double hi;
  int nthread;
  int multMH;
  int nplan,nint;
  static RANDSTATE rstate[MAXTHREADS];
  static int rsInit=0;
  /* output */
  double *lamlo, *loglambda;
  /* Get inputs */
  signal(SIGINT,sigintHandler);
  j = 0;
  lam = mxGetPr(prhs[j++]);
  T = mxGetM(prhs[j-1]);
  N = mxGetN(prhs[j-1]);
  muL = mxGetPr(prhs[j++]);
  sigL = mxGetPr(prhs[j++]);
  mull = mxGetPr(prhs[j++]);
  sigll = mxGetScalar(prhs[j++]);
  lamlo0 = mxGetPr(prhs[j++]);
  choice = mxGetPr(prhs[j++]);
  omega = mxGetPr(prhs[j++]);
  psi = mxGetPr(prhs[j++]);
  xint = mxGetPr(prhs[j++]);
  wint = mxGetPr(prhs[j++]);
  nint = mxGetNumberOfElements(prhs[j-1]);
  deduct = mxGetPr(prhs[j++]);
  nplan = mxGetDimensions(prhs[j-1])[0];
  maxoop = mxGetPr(prhs[j++]);
  prem = mxGetPr(prhs[j++]);
  hi = mxGetScalar(prhs[j++]);
  stdevi = mxGetPr(prhs[j++]);
  nthread = (int) mxGetScalar(prhs[j++]);
  multMH = (int) mxGetScalar(prhs[j++]);
  if (j!=nrhs) {
    mexErrMsgTxt("sampleLamlo: wrong number of input arguments");
  }
  /* Create and get outputs */
  plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL);
  lamlo = mxGetPr(plhs[0]);
  plhs[1] = mxCreateDoubleMatrix(T,N,mxREAL);
  loglambda = mxGetPr(plhs[1]);
  memcpy(lamlo,lamlo0,sizeof(double)*N);
  if (nthread>MAXTHREADS) mexErrMsgTxt("Too many threads.\n");
  if (rsInit==0) {
    time_t timeval;      
    for (int k=0;k<nthread;k++) {
      set_seed((unsigned int) time(&timeval) + k,rstate+k);
    }
    rsInit=1;
  }
  //mexPrintf("sampleLamlo: %d %d\n",N,T);
  //omp_set_num_threads(nthread);
#pragma omp parallel 
  {
    int i;
    int tid;
    int try = 0, accept=0;
    time_t t;      
    tid = omp_get_thread_num();
#pragma omp for 
    for(i=0;i<N;i++) {
      //mexPrintf("sampleLamlo: %d start | ",i);
      int t;
      //mexPrintf("%d: c[t=2]=%g\n",i,choice[i*T+1]);
      lamlo[i] = sampleLamLo(lamlo0[i], lam+i*T, muL+i*T, sigL[i]*stdevi[i], mull[i], sigll*stdevi[i],T, rstate+tid,
                             &try, &accept,
                             choice+i*T, omega[i],  psi[i],  xint,  wint, nint,
                             deduct+i*T*nplan,  maxoop+i*T*nplan,  prem+i*T*nplan, nplan,
                             0.1,hi,multMH);
      for(t=0;t<T;t++) {
        if (lamlo[i]+lam[t+i*T]<1e-12) mexPrintf("WARNING: lambda[i=%d,t=%d] = %.3g, ll0=%.3g l0=%.3g sum=%.3g\n",
                                                 i,t,lamlo[i]+lam[t+i*T],lamlo0[i],lam[t+i*T],lamlo0[i]+lam[t+i*T]);
        loglambda[t+i*T] = log(lamlo[i]+lam[t+i*T]);
      }
      //mexPrintf(" %d done\n",i);
    }
    if (omp_get_thread_num()==0) {
      int verbosity = (int) mxGetScalar(mexGetVariable("global","verbosity"));
      if (verbosity>0) 
        mexPrintf("sampleLamlo: accept rate=%.3f\n",((double) accept)/((double) try));
    }
  }
}

