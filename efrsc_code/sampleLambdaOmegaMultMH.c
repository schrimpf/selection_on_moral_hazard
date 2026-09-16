/* 
   Mex-function to sample psi and muL
   
   to compile: make.m
*/

#include <math.h>
#include <string.h>
#include <time.h>
#include "mex.h" /* header file for matlab API */
#include "alcoa.h"
#include <signal.h>

#include <omp.h>

/* Conditional on spending, log omega can be bounded. This structure
   contains information about the bounds */
struct bounds {
  int n;      // number of observatoins
  int maxint; // maximum number of intervals
  double *ninterval; // number of intervals describing bounds for each person
  double **interval; // bound for ith observation is the union over k of 
                     // (interval[k][i], interval[k][i + bnd.n])
};

double normcdf(double x,double mu, double sig) {
  int xinf=isinf(x);
  if (xinf<0) return(0);
  else if (xinf>0) return(1);
  else return( 0.5*(1+erf( (x-mu)/(sig*sqrt(2)) )) );
}

/* Returns a draw from a normal distribution with mean mu and variance
   var, truncated to lie in the area described by bounds.*[i]
   intersected with (-inf,xhi) */ 
double sampleNormalFromBounds(double mu,double var,
			      const struct bounds *bnd, int i, 
                              double xhi,
			      RANDSTATE *rs)
{
  double *p, sump=0;
  double u;
  int nint = bnd->ninterval[i];
  int k;
  double sig = sqrt(var);
  double lo, hi;
  double x;
  p = (double*) myCalloc(nint,sizeof(double));
  if (nint>1) {
    for(k=0;k<nint;k++) {
      double up = bnd->interval[k][i+bnd->n];
      up = up < xhi? up:xhi;
      sump += (p[k] = normcdf(up,mu,sig) - 
               normcdf(bnd->interval[k][i],mu,sig) );
    }
    //mexPrintf("sampleNormalFromBounds: i=%d, sump=%.4f\n",i,sump);
    u = Sample_Uniform(0,sump,rs);
    k = -1;
    sump = 0;
    while(u>sump) {
      sump += p[++k]; 
      //mexPrintf("u = %g p[%d]=%g\n",u,k,p[k]);
    }
  } else k = 0;
  lo = bnd->interval[k][i];
  hi = bnd->interval[k][i+bnd->n];
  hi = hi<xhi? hi:xhi;
  if ( isfinite(lo) ) {
    if (isfinite(hi)) x = Sample_Double_Truncated_Normal(mu,var,lo,hi,rs) ;
    else x = Sample_Truncated_Normal(mu,var,lo,1,rs);
  } else {
    if (isfinite(hi)) x = Sample_Truncated_Normal(mu,var,hi,0,rs) ;
    else x = Sample_Normal(mu,var,rs);
  }
  //mexPrintf("about to free p\n");
  myFree(p);
  return(x);
}

static RANDSTATE rstate[MAXTHREADS];
static int rsInit=0;


/* 
   Given omega, spending, lamlo, and other information, returns
   loglambda and log(density(m|omega))
*/
void loglAndlogDens(double *logl, double *ldens, 
		    const double *m, double logo, double lamlo, double c, 
		    const int *seg,const double *muL,double sigL, int T,
		    RANDSTATE *rs)
{
  double varL = sigL*sigL;
  *ldens = 0;
  for (int t=0;t<T;t++) {
    switch(seg[t]) {
    case -1: // missing
      logl[t] = Sample_Normal(muL[t],sigL,rs); 
      break;
    case 0: // spend = 0
      logl[t] = Sample_Truncated_Normal(muL[t],sigL,
                                          log(lamlo), 0, rs);
      break;
    case 1: // 0 < spend < deduct
      logl[t] = log(m[t]+lamlo);
      *ldens += 0;
      break;
    case 2: // deduct < spend< max
      logl[t] = log(m[t]/(1+(1-c)*exp(logo)) + lamlo);
      *ldens += -0.5*(logl[t]-muL[t])*(logl[t]-muL[t])/varL 
	- log(m[t]+lamlo*(1+(1-c)*exp(logo)));
      break;
    case 3: // max < spend
      logl[t] = log(m[t]/(1+exp(logo)) + lamlo);
      *ldens += -0.5*(logl[t]-muL[t])*(logl[t]-muL[t])/varL 
	- log(m[t]+lamlo*(1+exp(logo)));
      break;
    }
  }
}


// returns number of acceptances
#define copay 0.1
int sampleLambdaOmega
(double *logomega, double *loglambda, const double muO, const double sigO, 
 const double *muL, const double sigL,  const double psi,
 const double *xint, const double *wint, const int nint,
 const double *deduct, const double *maxoop, 
 const double *prem, const int nplan, 
 const double *choice, const double *spend, const int T, 
 const int mhTries, const double logOmegaHi,
 int *noverflow, int *mha, RANDSTATE *rs, double lamlo,
 const struct bounds *bnd, const int i)
{
  int accept=0;
  int naccept = 0;
  int t;
  double *logl = myCalloc(T,sizeof(double));
  double logo = *logomega;  
  double pold;
  int *seg = myCalloc(T,sizeof(int));
  int tid = omp_get_thread_num();
  double xtraHi = INFINITY; 
  (*noverflow) = (*mha) = 0;
  // add lower bound need omega > -(s+lamlo)/(lamlo*(1+(1-c)) to take
  // log(s+lamlo*(1+(1-c)om))  
  for(t=0;t<T;t++) {
    int p;
    if (choice[t]<=0) seg[t] = -1;
    else {
      if (spend[t]==0) { seg[t] = 0; continue; }
      p = t*nplan+(int) (choice[t]-1);
      if (spend[t]>0 && spend[t] <= deduct[p]) seg[t]=1;
      else if (spend[t]>deduct[p] && 
	       deduct[p]+copay*(spend[t] - deduct[p]) <= maxoop[p]) {
        double lb =  -(spend[t]+lamlo)/(lamlo*(1-copay));
        seg[t]=2;
        if (lb>0) {
          lb = log(lb);
          if (lamlo<0) {
            for(int k=0;k<bnd->ninterval[i];k++) {
              if (bnd->interval[k][i+bnd->n] > lb) xtraHi = xtraHi<lb? xtraHi:lb;
              if (bnd->interval[k][i]>lb) {
                mexPrintf("bound is empty %d %d\n",i,k);
                mexPrintf(" [%g , %g ], s=%g, lamlo=%g d=%g max=%g\n",
                          bnd->interval[k][i], lb, spend[t], lamlo,
                          deduct[p], maxoop[p]);
                mexErrMsgTxt("Cannot handle this");
              }            
            } 
          }
          else mexErrMsgTxt("lb>0 & lamlo>0 should be impossible");
        }
      }
      else {
        double lb =  -(spend[t]+lamlo)/lamlo;
        seg[t] = 3;	     
        if (lb>0) {
          lb = log(lb);
          if (lamlo<0) {
            for(int k=0;k<bnd->ninterval[i];k++) {
              if (bnd->interval[k][i+bnd->n]>lb) xtraHi = xtraHi<lb? xtraHi:lb;
              if (bnd->interval[k][i]>lb) {
                mexPrintf("bound is empty %d %d\n",i,k);
                mexPrintf(" [%g , %g ], s=%g, lamlo=%g d=%g max=%g\n",
                          bnd->interval[k][i], lb, spend[t], lamlo,
                          deduct[p], maxoop[p]);
                mexErrMsgTxt("Cannot handle this");
              }            
            } 
          } 
          else mexErrMsgTxt("lb>0 & lamlo>0 should be impossible");
        } // end if lb>0        
      } // end if else (seg)
    } // end else not missing
  } // end loop over t
  //mexPrintf("%d: setup seg\n",tid);
  loglAndlogDens(logl, &pold, 
		 spend, logo, lamlo, copay, 
		 seg, muL, sigL, T, rs);  
  for (int rep = 0;rep<mhTries;rep++) { // repeated metropolis draw
    double pnew,lonew;
    double om;
    double u = Sample_Uniform(0,1,rs);
    //mexPrintf("%d: mh rep %d\n",tid,rep);
    lonew = sampleNormalFromBounds(muO,sigO*sigO,bnd,i,xtraHi,rs); 
    //if (lonew<-20) {
    //mexPrintf("%d: obs %d, drew logom = %g. mu=%g sig=%g\n",tid,i,lonew,muO,sigO);
    //mexPrintf("%d:  bounds = %g %g \n", bnd->interval[0][i], bnd->interval[0][i+bnd->n]);
    //mexErrMsgTxt("sending error so can debug");
    //}
    loglAndlogDens(logl, &pnew, 
		   spend, lonew, lamlo, copay, 
		   seg, muL, sigL, T, rs);
    if (u<exp(pnew - pold)) {
      (*mha)++;
      accept = 1;      
      om = exp(lonew);    
      for(t = 0;t<T;t++) {
	double exv[5],exo[5];
	int newChoice;
	if (choice[t]>0) { // i.e. not missing
	  newChoice=choiceFnEx(om,psi,muL[t],sigL, 
			       xint, wint,nint,
			       deduct+t*nplan, maxoop+t*nplan, 
			       prem+t*nplan, nplan,copay,-lamlo,1, exv,exv,exo);
	  if (newChoice==0) {
	    accept = 0; 
	    (*noverflow)++;
	  } else {
	    accept = accept && ((int) choice[t])== newChoice;
	  }
	  if (!accept)  break;
	} // choice[t]>0
      }      
      if (accept){ 
        if (isnan(lonew) || isinf(lonew)) mexPrintf("lonew = %g accepted\n",lonew);
        for(t=0;t<T;t++) if (isnan(logl[t]) || isinf(logl[t])) 
          mexPrintf("logl[%d,%d] = %g  accepted, lonew=%g lamlo=%g | loold=%g loglold=%g spend=%g\n",
                    i,t,logl[t],lonew,lamlo,logomega[0],loglambda[t],spend[t]);
      	*logomega = lonew;
        memcpy(loglambda,logl,sizeof(double)*T);
	(naccept)++;
	logo = lonew;
	pold = pnew;
      }      
    } // if accept mh
  } // end loop over m-h draws
  myFree(logl);
  myFree(seg);
  return(naccept);
}

void mexFunction
(int nlhs, /* number of left hand side arguments */
 mxArray *plhs[], /* pointer to lhs*/
 int nrhs, /* number of rhs arguments*/
 const mxArray *prhs[]) /* pointer to rhs arguments*/
{
  /* rhs args */
  double *logomega;
  double *loglambda;
  double *muO;
  double *sigO;
  int nsigO;
  double *muL;
  double *sigL;
  double *psi;
  double *xint, *wint; // integration points and weights
  int nint;
  int nplan; // number of plans
  int T; // time periods
  int N; // observations
  double *deduct;  // deductible for obs i, time t, plan k is
                   // deduct[i*nplan*T+t*nplan+k] 
  double *maxoop;  
  double *prem;
  double *avail;   // indicator for whether each plan available for
                   // each obs/time
  double *choice;  // choice[i*T+t] \in \{1,..,nplan}
  double *spend; // [i*T + t]
  int mhTries; // maximum number of accept-reject attemps
  int nthread; 
  double logOmegaHi; 
  double *lamlo;
  double *stdevi;
  struct bounds bnd;
  int j;
  int accept = 0; // number accepted completely
  int mha = 0; // number accepted by mh step 
  int totover = 0;
  /* output */
  double *logo,*logl;
  if (nrhs!=21) {
    mexPrintf("nrhs=%d\n",nrhs);
    mexErrMsgTxt("Wrong number of arguments");
  }
  signal(SIGINT,sigintHandler);
  /* Get inputs */
  logomega = mxGetPr(prhs[0]);
  loglambda = mxGetPr(prhs[1]);
  muO = mxGetPr(prhs[2]);
  sigO = mxGetPr(prhs[3]); 
  nsigO = mxGetNumberOfElements(prhs[3]);
  muL = mxGetPr(prhs[4]);
  sigL = mxGetPr(prhs[5]);
  psi = mxGetPr(prhs[6]);
  nint = mxGetNumberOfElements(prhs[7]);
  xint = mxGetPr(prhs[7]);
  wint = mxGetPr(prhs[8]);
  nplan = mxGetDimensions(prhs[9])[0];
  T = mxGetDimensions(prhs[9])[1];
  N = mxGetDimensions(prhs[9])[2];
  deduct = mxGetPr(prhs[9]);
  maxoop = mxGetPr(prhs[10]);
  prem = mxGetPr(prhs[11]);
  avail = mxGetPr(prhs[12]);
  choice = mxGetPr(prhs[13]);
  if (mxGetNumberOfElements(prhs[13])!=N*T)
    mexErrMsgTxt("deduct and/or choice wrong size");
  spend = mxGetPr(prhs[14]);
  mxAssert(N==mxGetN(prhs[14]),"spend N!=N");
  mxAssert(T==mxGetM(prhs[14]),"spend T!=N");
  mhTries = (int) mxGetScalar(prhs[15]);
  nthread = (int) mxGetScalar(prhs[16]);
  logOmegaHi= mxGetScalar(prhs[17]);
  lamlo = mxGetPr(prhs[18]);
  stdevi = mxGetPr(prhs[19]);
  j  = 20;
  bnd.ninterval = mxGetPr(mxGetField(prhs[j],0,"ninterval"));
  bnd.n = mxGetNumberOfElements(mxGetField(prhs[j],0,"ninterval"));
  bnd.maxint = mxGetNumberOfElements(mxGetField(prhs[j],0,"interval"));
  bnd.interval = (double **) myCalloc(bnd.maxint,sizeof(double*));
  if (nrhs!=j+1) {
    mexErrMsgTxt("wrong number of arguments");
  }
  for(int k=0;k<bnd.maxint;k++) {
    bnd.interval[k] = 
      mxGetPr(mxGetCell(mxGetField(prhs[j],0,"interval"),k));
  }    
  mxAssert(nthread<=MAXTHREADS,"too many threads\n");
  if (rsInit==0) {
    time_t timeval;      
    for (int k=0;k<nthread;k++) {
      set_seed((unsigned int) time(&timeval) + k,rstate+k);
    }
    rsInit=1;
  }
  /* Create and get outputs */
  plhs[0] = mxCreateDoubleMatrix(1,N,mxREAL);
  logo = mxGetPr(plhs[0]);
  plhs[1] = mxCreateDoubleMatrix(T,N,mxREAL);
  logl = mxGetPr(plhs[1]);
  memcpy(logo,logomega,sizeof(double)*N);
  memcpy(logl,loglambda,sizeof(double)*N*T);
  //omp_set_num_threads(nthread);
#pragma omp parallel default(shared) 
  {
    int i;
    int tid = omp_get_thread_num();
    //if (tid==0) mexPrintf("using %d threads (asked for %d)\n",
    //          omp_get_num_threads(),nthread);
#pragma omp for reduction(+:totover,accept,mha) 
    for (i=0;i<N;i++) {
      int nover=0, mymha=0;
      double sigo = (nsigO>1)? sigO[i]:sigO[0];
      //mexPrintf("%d: working on %d\n",tid,i);
      accept += 
	sampleLambdaOmega(logo+i,logl+i*T, 
			  muO[i],(sigo)*stdevi[i],muL+i*T,sigL[i],psi[i],
			  xint,wint, nint,
			  deduct+i*T*nplan,maxoop+i*T*nplan,prem+i*T*nplan,
			  nplan,choice+i*T,spend+i*T,T,mhTries,logOmegaHi,
			  &nover,&mymha,rstate+tid,lamlo[i],&bnd,i);
      //mexPrintf("%d: finished %d\n",tid,i);
      totover += nover;
      mha += mymha;
    }
  } // end parallel region
  
  int verbosity = (int) mxGetScalar(mexGetVariable("global","verbosity"));
  if (verbosity>0) {
    mexPrintf("sampleLambdaOmega: total accept rate=%.3f, ", 
              ((double) accept)/((double) N*mhTries));
    mexPrintf(" mh accept rate = %.3f\n",((double) mha)/((double) N*mhTries));
    mexPrintf(" overflowed %d times (%.2lf of tries)\n",
              totover, ((double) totover)/(1e-300+(double) N*mhTries) );
  }
  if (nlhs>2) {
    plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
    (*mxGetPr(plhs[2])) = ((double) accept)/((double) N*mhTries);
  } 
  myFree(bnd.interval);
}
  
