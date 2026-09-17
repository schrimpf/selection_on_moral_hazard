/* 
   Mex-function to sample psi and muL
   
   to compile: make.m
*/

#include <math.h>
#include <string.h>
#include "mex.h" /* header file for matlab API */
//#include "arms.h"
#include "alcoa.h"

int choiceFnNMHuge
(const double omega, const double psi, const double muL, const double sigL, 
 const double *xint, const double *wint, const int nint,
 const double *deduct, const double *maxoop, const double *prem, const int nplan,
 const double copay, const double dl, double *exVal) ;

int choiceFnNMH
(const double omega, const double psi, const double muL, const double sigL, 
 const double *xint, const double *wint, const int nint,
 const double *deduct, const double *maxoop, const double *prem, const int nplan,
 const double copay, const double dl, double *exVal)
{
  int J=0;
  long double mev=-HUGE_VAL;
  int k;
  for(k=0;k<nplan;k++) {
    if (deduct[k] < 0) continue;
    int i;
    long double ev=0;
    for(i=0;i<nint;i++) {
      double lambda, oop, orate,m,util;
      int c;
      lambda = exp(xint[i]*sigL + muL) + dl;
      util = utilfn(0,lambda,omega,0);
      m = spendFn(lambda,omega,1.0);
      if (m > 0 && m < deduct[k]) {
        double u0 = utilfn(m,lambda,omega,m);
        if (u0>util) util = u0;          
      } 
      if (m>=deduct[k] && m < (deduct[k] + (maxoop[k]-deduct[k])/copay) && m>0) {  
        oop = deduct[k] + copay*(m-deduct[k]);
        double u0 = utilfn(m,lambda,omega,oop);
        if (u0>util) util = u0;  
      } 
      if (m>=(deduct[k] + (maxoop[k]-deduct[k])/copay)) {
        double u0 = utilfn(m,lambda,omega,maxoop[k]);
        if (u0>util) util = u0;;
      }       
      //mexPrintf("util[%2d] = %6.3g\n",i,util);
      ev += wint[i]*expl(-psi*util);
      if (!isfinite(ev)) {
        mexPrintf("WARNING(choiceFnNMH): overflow, switching to choiceFnHuge\n");
        return(choiceFnNMHuge(omega,  psi,  muL,  sigL,  xint,  wint, nint,
                              deduct,  maxoop,  prem, nplan,  copay,dl,exVal));
      }      
    } // for(i<nint) 
    ev = -logl(ev)-psi*prem[k];
    exVal[k]=ev;
    if (k==0 || ev>mev) { 
      mev = ev; 
      J = k+1; 
    }
  } // for(k<nplan) 
  //mexPrintf("\n");
  //if (!finitel(mev)) {
  //  J = 0; // mev is -inf or nan
  //mexPrintf("WARNING(choiceFnHuge): max E[V] = %Lg\n",mev);
  //}
  return(J);
}


#include <mpfr.h>
#define MPFR_USE_FILE
// number of bits for mantissa of large numbers -- set to 53 to match double
#define PREC 53 
int choiceFnNMHuge
(const double omega, const double psi, const double muL, const double sigL, 
 const double *xint, const double *wint, const int nint,
 const double *deduct, const double *maxoop, const double *prem, const int nplan,
 const double copay, const double dl, double *exVal)
{
  int J=0;
  mpfr_t mev;
  int k;
  mpfr_init2(mev,PREC);
  for(k=0;k<nplan;k++) {
    if (deduct[k] < 0) continue;
    int i;
    mpfr_t ev;
    mpfr_init2(ev,PREC);
    mpfr_set_d(ev,0.0,GMP_RNDN);
    for(i=0;i<nint;i++) {
      double lambda, oop, orate,m,util;
      int c;
      lambda = exp(xint[i]*sigL + muL) + dl;
      util = utilfn(0,lambda,omega,0);
      m = spendFn(lambda,omega,1.0);
      if (m > 0 && m < deduct[k]) {
        double u0 = utilfn(m,lambda,omega, m);
        if (u0>util) util = u0;          
      } 
      if (m>=deduct[k] && m < deduct[k] + (maxoop[k]-deduct[k])/copay && m>0) {  
        oop = deduct[k] + copay*(m-deduct[k]);
        double u0 = utilfn(m,lambda,omega, oop);
        if (u0>util) util = u0;  
      } 
      if (m>=deduct[k] + (maxoop[k]-deduct[k])/copay) {
        if ((-maxoop[k])>util) util = utilfn(m,lambda,omega,maxoop[k]);
      }          
      { 
        mpfr_t euw;
        mpfr_init2(euw,PREC); // euw = exp(-psi*util)
        mpfr_set_d(euw,-psi*util,GMP_RNDN);
        mpfr_exp(euw,euw,GMP_RNDN);
        if (i==0) mpfr_mul_d(ev,euw,wint[i],GMP_RNDN); // ev=wint[i]*exp(-psi*util)
        else {
          // ev += wint[i]*exp(-psi*util)
          mpfr_mul_d(euw,euw,wint[i],GMP_RNDN);
          mpfr_add(ev,euw,ev,GMP_RNDN);
        }
        mpfr_clear(euw);
      }
    } // for(i<nint) 
    // ev = -log(ev)-psi*prem[k]
    mpfr_log(ev,ev,GMP_RNDN);
    mpfr_mul_si(ev,ev,-1,GMP_RNDN);
    mpfr_sub_d(ev,ev,psi*prem[k],GMP_RNDN);
    //mexPrintf("huge: u(%d)=%Lg \n",k,lev);
    //mexPrintf("huge: u(%d) = %Lg \n",k,mpfr_get_ld(ev,GMP_RNDN));
    exVal[k]=mpfr_get_d(ev,GMP_RNDN);
    if (k==0 || mpfr_cmp(ev,mev)>0) //(k==0 || ev>mev) 
      { 
        mpfr_set(mev,ev,GMP_RNDN); //mev = ev; 
        J = k+1; 
      }
    mpfr_clear(ev);
  } // for(k<nplan) 
  mpfr_clear(mev);
  //mexPrintf("\n");
  //if (!finitel(mev)) {
  //  J = 0; // mev is -inf or nan
  //mexPrintf("WARNING(choiceFnHuge): max E[V] = %Lg\n",mev);
  //}
  return(J);
}
#undef PREC


#define copay 0.1
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
  
  int nint; 
  int maxIter;
  int N; /* number of observations */
  int T; /* number of time periods */
  int t,i,np,nplan;
  /* output */
  double *choiceOut;
  double *expectedValue;
  double *ev;
  if (nrhs!=16 && nrhs!=17) {
    mexPrintf("nrhs=%d\n",nrhs);
    mxAssert(nrhs==16, "Wrong number of arguments");
  }
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
  mxAssert(mxGetNumberOfElements(prhs[13])==N*T,"deduct and/or choice wrong size");
  spend = mxGetPr(prhs[14]);
  //mxAssert(N==mxGetN(prhs[14]),"spend N!=N");
  //mxAssert(T==mxGetM(prhs[14]),"spend T!=N");
  maxIter = (int) mxGetScalar(prhs[15]);
  if (nrhs<17) setNegLambdaZeroUtil(0);
  else setNegLambdaZeroUtil((int) mxGetScalar(prhs[16]));
  /* Create and get outputs */
  plhs[0] = mxCreateDoubleMatrix(T,N,mxREAL);
  if (nlhs>1) {
    mwSize dims[3];
    dims[0] = np;
    dims[1] = T;
    dims[2] = N;
    plhs[1] = mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
    expectedValue = mxGetPr(plhs[1]);
  } else {
    ev = mxCalloc(np,sizeof(double));
  }
  choiceOut = mxGetPr(plhs[0]);
  for(i=0;i<N;i++) {
    for(t=0;t<T;t++) {
      if (nlhs>1) ev = expectedValue+i*np*T+t*np;
      //if (i%100==0) mexPrintf("%d %d\n",i,t);
      choiceOut[i*T + t] = (double) choiceFnNMH(exp(logomega[i]),psi[i],muL[i*T+t],sigL[i], 
                                                xint, wint,nint,
                                                deduct+i*np*T+t*np, maxoop+i*np*T+t*np, 
                                                prem+i*np*T+t*np, np,copay, -lamlo[i], ev);
    }
  }
  if (nlhs<2) mxFree(ev);
  return;
}


