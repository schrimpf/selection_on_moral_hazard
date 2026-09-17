/* 
   Mex-function to find latent variables that rationalize choices.  
*/
void utilDiff(double *utilDiff, const double *usl, 
              const double *deduct,const double *maxoop,
              const double *prem,const double copay, const int np,
              const double *xint, const double *wint, const int nint,
              const double *choice, const int T);
void utilDiffHuge(double *diff, const double *usl,const double *deduct,const double *maxoop,
                  const double *prem,const double copay, const int np,
                  const double *xint, const double *wint, const int nint,
                  const double *choice, const int T);

#include <math.h>
#include "mex.h" /* header file for matlab API */
#include "alcoa.h"
#include <string.h>
#include <omp.h>
#include <mpfr.h>
#include <gmp.h>
#include <nlopt.h>
#include <signal.h>

static struct randState randStateVec[MAXTHREADS];
static int seeds_initialized =0;

#define NU (T+4)
typedef struct utilDiffArgs {
  int tid;
  int start;
  int end;
  double *muL;
  double *logomega;
  double *logpsi;
  double *logsigL;
  double *lamlo;
  double *deduct;
  double *maxoop;
  double *prem;
  int np;
  double *xint;
  double *wint;
  int nint;
  double *choice;
  int T;
  int nbad;
  double *uout;
  double *sigLout;
  double *lamloout;
  mxLogical *fail;
  int i;
  double *lb;
  double *ub;
  double *spend;
} UARGS;

typedef struct conArgs {
  UARGS *uargs;
  int t;
  int T;
} CONARGS;

double conFunc(int n, const double *x, double *grad, void *data) 
{
  UARGS *uargs = ((CONARGS*) data)->uargs;
  int t = ((CONARGS*) data)->t;
  int i = uargs->i;
  int np = uargs->np;
  int T = uargs->T;
  if (t<uargs->T) {
    double d, m, s,lamlo,omega;    
    int c=(int) uargs->choice[i*T+t]-1;
    if (c<0)  return(-1);
    d = uargs->deduct[i*np*T+t*np+c];
    m = uargs->maxoop[i*np*T+t*np+c];
    s = uargs->spend[i*T+t];
    omega = exp(x[0]);
    lamlo = x[3+T];
    if (s==0) return(-lamlo);
    else if (s>0 && s<=d) return(-s-lamlo);
    else if (s>d && d+0.1*(s-d)<m) return(-s-lamlo+0.9*omega);
    else if (d+0.1*(s-d)>=m) return(-s-lamlo+omega);
    mexErrMsgTxt("This should never happen");
    return(0);
  } else {
    int j = t-T;
    double diff[5*5];
    int i = uargs->i;
    int np = uargs->np;
    int T = uargs->T;
    mxAssert(n<=10 && uargs->T<=5 && uargs->np<=5,"something is too big");
    utilDiff(diff,x, uargs->deduct+i*np*T, uargs->maxoop+i*np*T, uargs->prem+i*np*T,
             0.1, np, uargs->xint, uargs->wint, uargs->nint, uargs->choice+i*T, T);
    return(diff[j]);
  }  
}

double objFunc(int n, const double *x, double *grad, void *data)
{
  CONARGS *ca = (CONARGS*) data;
  UARGS *uargs = ca[0].uargs;
  double val=0;
  double diff[5*5];
  int i;
  int np = uargs->np;
  int T = uargs->T;
  for(i=0;i<T;i++) {
    double cf = conFunc(n,x,grad,ca+i);
    if (cf>0) val += cf*cf;
  }
  i = uargs->i;
  utilDiff(diff,x, uargs->deduct+i*np*T, uargs->maxoop+i*np*T, uargs->prem+i*np*T,
           0.1, np, uargs->xint, uargs->wint, uargs->nint, uargs->choice+i*T, T);
  for(i=0;i<T*np;i++) {
    if (diff[i]>0) val+=diff[i]*diff[i];
  }  
  return(val);
}

/*
double objFunc(int n, const double *x, double *grad, void *x0) 
{
  int i;
  double val=0;
  for(i=0;i<n;i++) val+=(x[i]-x0[i])*(x[i]-x0[i]);
  return(sqrt(val));
}
*/

static int optimizeOneObs(int i, UARGS *uargs, CONARGS *cargs, double *u, double *u0)
{
  int T = uargs->T;
  int returnCode = 0;
  double diff;
  int t;
  uargs->i = i;
  u[0] = uargs->logomega[i];
  u[1] = uargs->logpsi[i];
  for (t = 0; t < T; t++) u[2+t] = uargs->muL[i*T + t];
  u[2+T] = uargs->logsigL[i];
  u[3+T] = uargs->lamlo[i];
  memcpy(u0, u, (4+T)*sizeof(double));
  for (t = 0; t < (4+T); t++) {
    if (u[t] < uargs->lb[i*NU + t]) u[t] = uargs->lb[i*NU + t] + 0.001;
    else if (u[t] > uargs->ub[i*NU + t]) u[t] = uargs->ub[i*NU + t] - 0.001;
  }
  diff = objFunc(4+T, u, NULL, cargs);
  if (diff > 0) {
    returnCode = nlopt_minimize_constrained(NLOPT_GN_DIRECT_L, 4+T, objFunc, cargs,
                                            0, NULL, NULL, 0,
                                            uargs->lb + i*NU, uargs->ub + i*NU, u, &diff,
                                            0, 0, 0, 1e-8, NULL, 0, 0.05);
  }
  if (diff <= 0) {
    memcpy(uargs->uout + i*(2+T), u, (2+T)*sizeof(double));
    uargs->sigLout[i] = exp(u[2+T]);
    uargs->lamloout[i] = u[3+T];
    uargs->fail[i] = 0;
    return 0;
  } else {
    uargs->fail[i] = 1;
    memcpy(uargs->uout + i*(2+T), u, (2+T)*sizeof(double));
    uargs->sigLout[i] = exp(u[2+T]);
    uargs->lamloout[i] = u[3+T];
    return 1;
  }
}

void mexFunction
(int nlhs, /* number of left hand side arguments */
 mxArray *plhs[], /* pointer to lhs*/
 int nrhs, /* number of rhs arguments*/
 const mxArray *prhs[]) /* pointer to rhs arguments*/
{
  int N, nbad = 0;
  UARGS uargs;
  int threads;
  signal(SIGINT, sigintHandler);
  uargs.nbad = 0;
  uargs.muL = mxGetPr(prhs[0]);
  uargs.logsigL = mxGetPr(prhs[1]);
  uargs.logomega = mxGetPr(prhs[2]);
  uargs.logpsi = mxGetPr(prhs[3]);
  uargs.deduct = mxGetPr(prhs[4]);
  uargs.np = mxGetDimensions(prhs[4])[0];
  uargs.T = mxGetDimensions(prhs[4])[1];
  N = mxGetDimensions(prhs[4])[2];
  uargs.maxoop = mxGetPr(prhs[5]);
  uargs.prem = mxGetPr(prhs[6]);
  uargs.xint = mxGetPr(prhs[7]);
  uargs.wint = mxGetPr(prhs[8]);
  uargs.nint = mxGetNumberOfElements(prhs[8]);
  uargs.choice = mxGetPr(prhs[9]);
  threads = (int) mxGetScalar(prhs[10]);
  uargs.lb = mxGetPr(prhs[11]);
  uargs.ub = mxGetPr(prhs[12]);
  uargs.lamlo = mxGetPr(prhs[13]); 
  uargs.spend = mxGetPr(prhs[14]);
  mxAssert(threads <= MAXTHREADS, "too many threads\n");

  plhs[0] = mxCreateDoubleMatrix(2+uargs.T, N, mxREAL); 
  plhs[1] = mxCreateDoubleMatrix(N, 1, mxREAL);
  plhs[2] = mxCreateDoubleMatrix(N, 1, mxREAL);
  plhs[3] = mxCreateLogicalMatrix(N, 1);
  uargs.uout = mxGetPr(plhs[0]);
  uargs.sigLout = mxGetPr(plhs[1]);
  uargs.lamloout = mxGetPr(plhs[2]);
  uargs.fail = mxGetLogicals(plhs[3]);

  /* Pre-populate all outputs with current values */
  for (int i = 0; i < N; i++) {
    uargs.uout[i*(2+uargs.T) + 0] = uargs.logomega[i];
    uargs.uout[i*(2+uargs.T) + 1] = uargs.logpsi[i];
    for (int t = 0; t < uargs.T; t++) {
      uargs.uout[i*(2+uargs.T) + 2 + t] = uargs.muL[i*uargs.T + t];
    }
    uargs.sigLout[i] = exp(uargs.logsigL[i]);
    uargs.lamloout[i] = uargs.lamlo[i];
    uargs.fail[i] = 0;
  }

  /* Optional badindex argument */
  int nbadobs = -1;
  double *badindices = NULL;
  if (nrhs > 15 && !mxIsEmpty(prhs[15])) {
    nbadobs = (int) mxGetNumberOfElements(prhs[15]);
    badindices = mxGetPr(prhs[15]);
  }

  int totalTasks = (nbadobs >= 0) ? nbadobs : N;
  if (threads > 0) {
    omp_set_num_threads(threads);
  }

  #pragma omp parallel reduction(+:nbad)
  {
    UARGS my_uargs = uargs;
    int T = my_uargs.T;
    int np = my_uargs.np;
    CONARGS cargs[1 + 5 + 25];
    double u[4 + 5];
    double u0[4 + 5];
    for (int j = 0; j < T + (np * T); j++) {
      cargs[j].uargs = &my_uargs;
      cargs[j].t = j;
      cargs[j].T = (T + np * T);
    }
    #pragma omp for schedule(dynamic, 1)
    for (int k = 0; k < totalTasks; k++) {
      int i = (badindices != NULL) ? ((int)(badindices[k] - 1)) : k;
      if (i >= 0 && i < N) {
        nbad += optimizeOneObs(i, &my_uargs, cargs, u, u0);
      }
    }
  }

  if (nbad > 0) mexPrintf("WARNING: failed to find valid latent for %d observations\n", nbad);
}
     

void utilDiff(double *utilDiff, const double *usl, 
              const double *deduct,const double *maxoop,
              const double *prem,const double copay, const int np,
              const double *xint, const double *wint, const int nint,
              const double *choice, const int T)
{
  double omega, psi, sigL, lamlo;
  long double util[5];
  double muL[5];
  int t;
  omega = exp(usl[0]);
  psi = exp(usl[1]);
  for(t=0;t<T;t++) muL[t] = usl[2+t];
  sigL = exp(usl[2+T]);
  lamlo = usl[3+T];
  for (t=0;t<np;t++) {
    int k; 
    for (k=0;k<np;k++) utilDiff[k+t*np]=0;
  }
  for(t=0;t<T;t++) {
    int k;
    if (choice[t]<=0) continue;
    for (k=0;k<np;k++) {
      int j;
      long double ev=0;      
      if (isnan(deduct[k+t*np]) || deduct[k + t*np]<0) {
        util[k] = -1e300;
        continue;
      }

      for(j=0;j<nint;j++) {
        double lambda = exp(muL[t]+xint[j]*sigL)-lamlo;
        double maxu = utilfn(0,lambda,omega,0.0,0);
        double spend = lambda;
        if (spend>0 && spend<deduct[k+t*np]) {
          double u = utilfn(spend,lambda,omega,spend,0);
          if (u>maxu) {
            maxu = u; 
          }
        }
        spend = lambda + (1-copay)*omega;
        if (spend>=deduct[k+t*np] && spend<deduct[k+t*np]+(maxoop[k+t*np]-deduct[k+t*np])/copay && spend>0) {
          double u = utilfn(spend,lambda,omega,deduct[k+t*np]+copay*(spend-deduct[k+t*np]),0);
          if (u>maxu) {
            maxu = u; 
          }
        }
        spend = lambda+omega;
        if (spend>=deduct[k+t*np]+(maxoop[k+t*np]-deduct[k+t*np])/copay) {
          double u=utilfn(spend,lambda,omega,maxoop[k+t*np],0);
          if (u>maxu) {
            maxu = u; 
          }
        }
        ev += expl(-psi*maxu)*wint[j];
        if (!isfinite(ev)) {
          mexPrintf("WARNING(findValidLatent): overflow, switching to utilDiffHuge\n");
          utilDiffHuge(utilDiff, usl, deduct,maxoop, prem, copay, np, xint, wint,nint, choice, T); 
          return;
        }
      }
      util[k] = -logl(ev) - psi*prem[k+t*np];
    } // for(k=0;k<np;k++) 
    for (k=0;k<np;k++) {
      int c = (int) choice[t]-1;
      utilDiff[k+t*np] = util[k]-util[c];
    }
  } // for(t)
}

#define PREC 53
void utilDiffHuge(double *utilDiff, const double *usl, 
                  const double *deduct,const double *maxoop,
                  const double *prem,const double copay, const int np,
                  const double *xint, const double *wint, const int nint,
                  const double *choice, const int T)
{
  double omega, psi, *muL, sigL, lamlo;
  //long double utilDiff=0;
  mpfr_t *util;
  int t;
  util = myCalloc(np,sizeof(mpfr_t));
  muL = myCalloc(T,sizeof(double));
  for(t=0;t<np;t++) {
    mpfr_init2(util[t],PREC);
  }
  omega = exp(usl[0]);
  psi = exp(usl[1]);
  for(t=0;t<T;t++) muL[t] = usl[2+t];
  sigL = exp(usl[2+T]);
  lamlo = usl[3+T];
  for (t=0;t<np;t++) {
    int k; 
    for (k=0;k<np;k++) utilDiff[k+t*np]=0;
  }
  for(t=0;t<T;t++) {
    int k;
    if (choice[t]<=0) continue;
    for (k=0;k<np;k++) {
      int j;
      mpfr_t ev;      
      mpfr_t euw;
      if (isnan(deduct[k+t*np]) || deduct[k + t*np]<0) {
        mpfr_set_d(util[k],-1e307,GMP_RNDN);
        continue;
      }
      mpfr_init2(ev,PREC);
      mpfr_init2(euw,PREC);
      mpfr_set_d(ev,0.0,GMP_RNDN);

      for(j=0;j<nint;j++) {
        double lambda = exp(muL[t]+xint[j]*sigL)-lamlo;
        double maxu = utilfn(0,lambda,omega,0.0,0);
        double spend = lambda;
        if (spend>0 && spend<deduct[k+t*np]) {
          double u = utilfn(spend,lambda,omega,spend,0);
          if (u>maxu) {
            maxu = u; 
          }
        }
        spend = lambda + (1-copay)*omega;
        if (spend>=deduct[k+t*np] && spend<deduct[k+t*np]+(maxoop[k+t*np]-deduct[k+t*np])/copay && spend>0) {
          double u = utilfn(spend,lambda,omega,deduct[k+t*np]+copay*(spend-deduct[k+t*np]),0);
          if (u>maxu) {
            maxu = u; 
          }
        }
        spend = lambda+omega;
        if (spend>=deduct[k+t*np]+(maxoop[k+t*np]-deduct[k+t*np])/copay) {
          double u=utilfn(spend,lambda,omega,maxoop[k+t*np],0);
          if (u>maxu) {
            maxu = u; 
          }
        }
        mpfr_set_d(euw,-psi*maxu,GMP_RNDN);
        mpfr_exp(euw,euw,GMP_RNDN);
        mpfr_mul_d(euw,euw,wint[j],GMP_RNDN);
        mpfr_add(ev,euw,ev,GMP_RNDN);
        //ev += expl(-psi*maxu)*wint[j];
      }
      //util[k] = -logl(ev) - psi*prem[k+t*np];
      mpfr_log(ev,ev,GMP_RNDN);
      mpfr_mul_si(ev,ev,-1,GMP_RNDN);
      mpfr_sub_d(util[k],ev,psi*prem[k+t*np],GMP_RNDN);
      //mexPrintf("%d %d: %g %g %g\n",t,k,deduct[k+t*np],maxoop[k+t*np],prem[k+t*np]);
      mpfr_clear(ev);
      mpfr_clear(euw);
    } // for(k=0;k<np;k++) 
    for (k=0;k<np;k++) {
      int c = (int) choice[t]-1;
      //if (mpfr_cmp(util[k],util[c])>0) {
        mpfr_t reld,dom;
        double diff;
        mpfr_init2(reld,PREC);        
        mpfr_init2(dom,PREC);        
        //double reld = (util[k]-util[c])/(fabs(util[c])+1);
        //mpfr_abs(dom,util[c],GMP_RNDN);
        //mpfr_add_d(dom,dom,1,GMP_RNDN);
        mpfr_sub(reld,util[k],util[c],GMP_RNDN);
        //mpfr_div(reld,reld,dom,GMP_RNDN);
        diff = mpfr_get_d(reld,GMP_RNDN);
        //if (diff>=0) {
        utilDiff[k+t*np] = diff; //mpfr_get_d(reld,GMP_RNDN);
        //}
        //mexPrintf("diff=%g\n",diff);
        //mexPrintf("ud=%g u=%g %g %g sigL=%g\n",utilDiff, 
        //          log(omega),log(psi),muL,sigL);  
        //if (utilDiff<=0) mexPrintf("utilDiff=%g loop\n",utilDiff);
        mpfr_clear(reld);
        mpfr_clear(dom);
        //}
      //if (k==3 || k==4) mexPrintf("u(%d)=%Lg\n",k,util[k]);
    }
  } // for(t)
  for(t=0;t<np;t++) mpfr_clear(util[t]);
  myFree(util);                  
  myFree(muL);
  //if (utilDiff<=0) mexPrintf("utilDiff=%g\n",utilDiff);
  //return(utilDiff);
}
