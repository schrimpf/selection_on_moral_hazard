/*
   Mex-function to sample lambda and omega

   to compile: make.m
*/

#include <math.h>
#include <string.h>
#include "mex.h" /* header file for matlab API */
#include "alcoa.h"
#include <pthread.h>
#include <signal.h>

/* Declaration of Function that computes indirect utility and spending */
int sampleLambdaOmega
(double *logomega, double *loglambda, const double muO, const double sigO,
 const double *muL, const double sigL,  const double psi,
 const double *xint, const double *wint, const int nint,
 const double *deduct, const double *maxoop, const double *prem, const int nplan,
 const double *choice, const double *spend, const int T, const int maxIter,const double logOmegaHi,
 int *tries, int *nover, RANDSTATE *rs, double lamlo);

struct u0parm {
  double *spend,c, *muL, sigL, muO, sigO, omega, lamlo;
  int T, *useS;
};
static RANDSTATE rstate[MAXTHREADS];
static int rsInit=0;

int utils(double *u, double lam, double om, double ded, double max, double co)
{
  double m, oop;
  double maxu;
  int seg;
  const int multMH = 0;
  maxu = u[0] = utilfn(0,lam,om,0,multMH);
  seg = 0;
  oop = m = lam;
  if (m<0 || m>ded) u[1] = -1e300;
  else   {
    u[1] = utilfn(m,lam,om,oop,multMH);
    if (u[1]>maxu) {
      maxu = u[1];
      seg = 1;
    }
  }
  m = lam+(1-co)*om;
  oop = ded + co*(m - ded);
  if (m<ded || oop>max)
    u[2] = -1e300;
  else   {
    u[2] = utilfn(m,lam,om,oop,multMH);
    if (u[2]>maxu) {
      maxu = u[2];
      seg = 2;
    }
  }
  m = lam+om;
  if (max>ded+co*(m-ded) || m<0) u[3] = -1e300;
  else {
    u[3] = utilfn(m,lam,om,max,multMH);
    if (u[3]>maxu) {
      maxu = u[3];
      seg = 3;
    }
  }
  return(seg);
}

double logDensity(double logo,void *udata) {
  struct u0parm *d=(struct u0parm*) udata;
  double logd= 0 ;//-.5*(logo - d->muO)*(logo-d->muO)/(d->sigO*d->sigO);
  int t;
  for(t=0;t<d->T;t++) {
    if (d->useS[t]==3) {
      double logl = d->spend[t]-exp(logo)+d->lamlo;
      if (logl<0) {
        return(-HUGE_VAL);
      } else {
        logl = log(logl);
      }
      logd += -.5*(logl-d->muL[t])*(logl-d->muL[t])/(d->sigL*d->sigL) - logl;
    }
    else if (d->useS[t]==2) {
      double logl = (-(1-d->c)*exp(logo)+d->spend[t]+d->lamlo);
      if (logl<0) logl=0;
      logl = log(logl);
      logd += -.5*(logl-d->muL[t])*(logl-d->muL[t])/(d->sigL*d->sigL) - logl;
    }
  }
  return(logd);
}

struct sampleArgs {
  double *logomega;
  double *loglambda;
  double *muO;
  double *sigO;
  int nsigO;
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
  double logOmegaHi;
  double *lamlo;

  double *stdevi;

  int nint;
  int maxIter;
  int T; /* number of time periods */
  int np;

  int tries;
  int totover;
  /* output */
  double *logo,*logl;
  double nbad;
  /* thread management */
  int tid;
  int start;
  int end;
};

void *sampleLO(void *voidptr)
{
  struct sampleArgs *sa = (struct sampleArgs*) voidptr;
  int i;
  sa->tries=0;
  sa->totover=0;
  for(i=sa->start;i<sa->end;i++) {
    int t;
    int nover = 0;
    double so = sa->nsigO>1? sa->sigO[i]:sa->sigO[0];
    sa->nbad += sampleLambdaOmega(sa->logo+i,sa->logl+i*sa->T,sa->muO[i],so*sa->stdevi[i],
                                  sa->muL+i*sa->T,
                                  sa->sigL[i]*sa->stdevi[i],
                                  sa->psi[i],
                                  sa->xint,sa->wint,sa->nint,
                                  sa->deduct+i*sa->np*sa->T,sa->maxoop+i*sa->np*sa->T,sa->prem+i*sa->np*sa->T,
                                  sa->np,sa->choice+i*sa->T, sa->spend+i*sa->T, sa->T,
                                  sa->maxIter, sa->logOmegaHi,&t, &nover, rstate+sa->tid, sa->lamlo[i]);
    sa->tries+=t;
    sa->totover+=nover;
  }
  pthread_exit(NULL);
}


void mexFunction
(int nlhs, /* number of left hand side arguments */
 mxArray *plhs[], /* pointer to lhs*/
 int nrhs, /* number of rhs arguments*/
 const mxArray *prhs[]) /* pointer to rhs arguments*/
{
  /* rhs args */
  struct sampleArgs sarg,*targ;
  int i,N;
  int tries = 0;
  int totover = 0;
  int nthread;
  pthread_t *thrptr;
  /* output */
  double *nbad;
  if (nrhs!=20) {
    mexPrintf("nrhs=%d\n",nrhs);
    mxAssert(nrhs==20, "Wrong number of arguments");
  }
  signal(SIGINT,sigintHandler);
  /* Get inputs */
  sarg.logomega = mxGetPr(prhs[0]);
  sarg.loglambda = mxGetPr(prhs[1]);
  sarg.muO = mxGetPr(prhs[2]);
  sarg.sigO = mxGetPr(prhs[3]);
  sarg.nsigO = mxGetNumberOfElements(prhs[3]);
  sarg.muL = mxGetPr(prhs[4]);
  sarg.sigL = mxGetPr(prhs[5]);
  sarg.psi = mxGetPr(prhs[6]);

  sarg.nint = mxGetNumberOfElements(prhs[7]);
  sarg.xint = mxGetPr(prhs[7]);
  sarg.wint = mxGetPr(prhs[8]);
  sarg.np = mxGetDimensions(prhs[9])[0];
  sarg.T = mxGetDimensions(prhs[9])[1];
  N = mxGetDimensions(prhs[9])[2];
  sarg.deduct = mxGetPr(prhs[9]);
  sarg.maxoop = mxGetPr(prhs[10]);
  sarg.prem = mxGetPr(prhs[11]);
  sarg.avail = mxGetPr(prhs[12]);
  sarg.choice = mxGetPr(prhs[13]);
  mxAssert(mxGetNumberOfElements(prhs[13])==N*sarg.T,"deduct and/or choice wrong size");
  sarg.spend = mxGetPr(prhs[14]);
  mxAssert(N==mxGetN(prhs[14]),"spend N!=N");
  mxAssert(sarg.T==mxGetM(prhs[14]),"spend T!=N");
  sarg.maxIter = (int) mxGetScalar(prhs[15]);
  nthread = (int) mxGetScalar(prhs[16]);
  sarg.logOmegaHi= mxGetScalar(prhs[17]);
  sarg.lamlo = mxGetPr(prhs[18]);
  sarg.stdevi = mxGetPr(prhs[19]);

  mxAssert(nthread<=MAXTHREADS,"too many threads\n");
  if (rsInit==0) {
    for (i=0;i<nthread;i++) set_time_seed(rstate+i);
    rsInit=1;
  }
  /* Create and get outputs */
  plhs[0] = mxCreateDoubleMatrix(1,N,mxREAL);
  sarg.logo = mxGetPr(plhs[0]);
  plhs[1] = mxCreateDoubleMatrix(sarg.T,N,mxREAL);
  sarg.logl = mxGetPr(plhs[1]);
  if (nlhs>2) {
    plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
    nbad = mxGetPr(plhs[2]);
    *(nbad) = 0;
  } else {
    if (NULL==(nbad = myCalloc(1,sizeof(double)))) mexErrMsgTxt("Allocation failure");
  }
  memcpy(sarg.logo,sarg.logomega,sizeof(double)*N);
  memcpy(sarg.logl,sarg.loglambda,sizeof(double)*N*sarg.T);
  targ = myCalloc(nthread,sizeof(struct sampleArgs));
  thrptr = myCalloc(nthread,sizeof(pthread_t));
  for(i=0;i<nthread;i++) {
    memcpy(targ+i,&sarg,sizeof(struct sampleArgs));
    targ[i].tid = i;
    targ[i].start = i*N/nthread;
    if (i<nthread-1) targ[i].end = (i+1)*N/nthread;
    else targ[i].end = N;
    pthread_create(thrptr+i,NULL, sampleLO, (void *) (targ+i));
  }
  for (i=0;i<nthread;i++) {
    void *status;
    pthread_join(thrptr[i],&status);
    tries += targ[i].tries;
    totover += targ[i].totover;
    nbad[0] += targ[i].nbad;
  }
  //mexPrintf("joined\n");
  myFree(targ);
  //mexPrintf("freed args\n");
  myFree(thrptr);
  int verbosity = (int) mxGetScalar(mexGetVariable("global","verbosity"));
  if (verbosity>0) {
    mexPrintf("sampleLambdaOmega: avg tries=%g, overflowed %d times (%.2lf of tries), reached max %d times\n", \
              ((double) tries)/((double) N),totover,
              ((double) totover)/(1e-300+(double) tries) , (int) nbad[0]);
  }
  if (nlhs<=2) myFree(nbad);
}


#define copay 0.1
int sampleLambdaOmega
(double *logomega, double *loglambda, const double muO, const double sigO,
 const double *muL, const double sigL,  const double psi,
 const double *xint, const double *wint, const int nint,
 const double *deduct, const double *maxoop, const double *prem, const int nplan,
 const double *choice, const double *spend, const int T, const int maxIter, const double logOmegaHi,
 int *tries, int *noverflow, RANDSTATE *rs, double lamlo)
{
  int accept=0;
  int bad=0;
  struct u0parm upm;
  int t;
  double *logl = myCalloc(T,sizeof(double));
  double logo = *logomega;
  int nsamp = 5;
  double usamp[5];
  int s=nsamp+1;
  int reason[5] = {0,0,0,0,0};
  double xinit[20];
  double lohi = logOmegaHi;
  int ninit = 20;
  int mhA = 0;
  int mhT = 0;

  upm.spend = spend;
  upm.T = T;
  upm.c = copay;
  upm.muL = muL;
  upm.sigL = sigL;
  upm.muO = muO;
  upm.sigO = sigO;
  upm.useS = myCalloc(T,sizeof(int));
  upm.lamlo = lamlo;
  for(t=0;t<T;t++) {
    if (choice[t]>0) {
      int k = ((int) choice[t])-1+t*nplan;
      if (spend[t]>0 && spend[t]<deduct[k]) upm.useS[t] = 1;
      else if (spend[t]>=deduct[k] &&
               spend[t]<deduct[k] + (maxoop[k] - deduct[k])/copay &&
               spend[t]>0)  {
        upm.useS[t] = 2;
        lohi = lohi<log((spend[t]+lamlo)/(1-copay))? lohi:log((spend[t]+lamlo)/(1-copay));
      }
      else if (spend[t]>0 && spend[t]>deduct[k] + (maxoop[k] - deduct[k])/copay) {
        upm.useS[t] = 3;
        lohi = lohi<log(spend[t]+lamlo)? lohi:log(spend[t]+lamlo);
      }
      else upm.useS[t] = 0;
    } else upm.useS[t] = 0;
  }

  if (logo>lohi) mexPrintf("crap %g %g\n",logo,lohi);

  (*tries)=0;
  //for (t=0;t<ninit;t++) xinit[t] = muO + sigO*6*((double) t-nint/2)/((double) nint);
  while (!accept) {
    double lo=muO-9*sigO;
    double hi=muO+9*sigO; // (inf)
    double convex = 5;
    double s2;
    double om;
    int neval;
    int out;
    accept=1;
    // do five M-H steps
    s2 = logDensity(logo,&upm);
    //if (!isfinite(s2)) mexPrintf("s2 is not finite\n");
    for(t=0;t<5;t++) {
      double lonew;
      double pnew;
      double u = Sample_Uniform(0,1,rs);
      if (isfinite(lohi)) lonew = Sample_Truncated_Normal(upm.muO,upm.sigO*upm.sigO,lohi,0,rs);
      else lonew = Sample_Normal(upm.muO,upm.sigO*upm.sigO,rs);
      pnew = logDensity(lonew,&upm);
      if (u<exp(pnew - s2)) {
        logo = lonew;
        s2 = pnew;
        mhA++;
      }
      mhT++;
    }

    //mexPrintf("old=%g sample=%g\n",*logomega,logo);
    //logo = Sample_Normal(muO,sigO*sigO);
    om = exp(logo);
    for(t=0;t<T;t++) {
      int k = ((int) choice[t])-1+t*nplan;
      double lam;
      double u[4];
      if (choice[t]<=0) {
        logl[t] = Sample_Normal(muL[t],sigL*sigL,rs);
        continue; // person missing this period
      }
      if (spend[t]<=0) {
        logl[t] = Sample_Normal(muL[t],sigL*sigL,rs);
        lam = exp(logl[t]) - lamlo;
        accept = accept && 0==utils(u, lam,om, deduct[k], maxoop[k], copay);
        if (!accept) reason[0]++;
      } else if (upm.useS[t]==1) {
        lam = spend[t];
        if (lam+lamlo<0) { accept=0;}
        logl[t] = log(lam+lamlo);
        accept = accept && 1==utils(u, lam,om, deduct[k], maxoop[k], copay);
        if (!accept) reason[1]++;
      } else if (upm.useS[t]==2) {
        lam = spend[t]-(1-copay)*om;
        if (lam+lamlo<0) { accept=0;}
        logl[t] = log(lam+lamlo);
        accept = accept && 2==utils(u, lam,om, deduct[k], maxoop[k], copay);
        if (!accept) reason[2]++;
      } else if (upm.useS[t]==3) {
        lam = spend[t] - om;
        if (lam+lamlo<0) { accept=0;}
        logl[t] = log(lam+lamlo);
        accept = accept && 3==utils(u, lam,om, deduct[k], maxoop[k], copay);
        if (!accept) reason[3]++;
      }
      if (accept) {
        double exv[5],exo[5];
        int newChoice=choiceFnEx(om,psi,muL[t],sigL,
                                 xint, wint,nint,
                                 deduct+t*nplan, maxoop+t*nplan,
                                 prem+t*nplan, nplan,copay,-lamlo,0,
                                 exv,exv,exo);
        if (newChoice==0) {
          accept = 0;
          (*noverflow)++;
        } else {
          accept = accept && ((int) choice[t])== newChoice;
        }
        //accept=1;
        if (!accept) reason[4]++;
      }
      if (!accept) break;
    }
    (*tries)++;
    if (*tries>maxIter) {
      int j;
      bad = 1;
      break;
    }
    //accept = 1;
  } /* while(!accept) */
  if (!bad) {
    *logomega = logo;
    memcpy(loglambda,logl,sizeof(double)*T);
  } //else {
    //mexPrintf("reached maxIter %d %d %d %d %d\n",reason[0],reason[1],reason[2],reason[3],reason[4]);
    //mexPrintf("Accepted %g of draws\n",((double) mhA)/((double) mhT));
  //}
  /*
  if (*tries>100) {
    int j;
    mexPrintf("%d tries, reason=", *tries);
    for(j=0;j<5;j++) mexPrintf(" %d ",reason[j]);
    mexPrintf("\n");
    if (reason[4]>maxIter/2) {
      int t, J;
      for(t=0;t<T;t++) {
        int J = choiceFn(exp(*logomega),psi,muL[t],sigL,
                         xint, wint,nint,
                         deduct+t*nplan, maxoop+t*nplan,
                         prem+t*nplan, nplan,copay);
          mexPrintf("t=%d choice=%d J=%d nplan=%d\n",t,((int) choice[t]),J,nplan);
      }
    }
  }
  */
  myFree(upm.useS);
  myFree(logl);
  return(bad);
}
