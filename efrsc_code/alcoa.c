#include <math.h>
#include <time.h>
#include "mex.h"
#include "alcoa.h"

void sigintHandler(int sig) {
  mexErrMsgTxt("Interrupted by ctrl-c\n");
}

int choiceFnHuge
(const double omega, const double psi, const double muL, const double sigL, 
 const double *xint, const double *wint, const int nint,
 const double *deduct, const double *maxoop, const double *prem, const int nplan,
 const double copay, const double dl, const int multMH, double *exVal, double *exSpend, double *exOop);

int choiceFn
(const double omega, const double psi, const double muL, const double sigL, 
 const double *xint, const double *wint, const int nint,
 const double *deduct, const double *maxoop, const double *prem, const int nplan,
 const double copay, const int multMH, double *exVal, double *exSpend, double *exOop)
{
  return(choiceFnEx(omega,  psi,  muL,  sigL,  xint,  wint, nint,
                    deduct,  maxoop,  prem, nplan,  copay,0, multMH, exVal, exSpend, exOop));
}

static int negLambdaZeroUtil = 0;
void setNegLambdaZeroUtil(int val) { 
  negLambdaZeroUtil = val; 
}

double utilfn(const double spend,const double lambda,const double omega,
              const double cost, const int multMH) {
  if (multMH) { // multiplicative
    //return((1+lambda/omega)*spend - spend*spend/(2*omega));
    if (lambda<=0) return 0;
    else return((spend-lambda)-(spend-lambda)*(spend-lambda)/(2*omega*lambda) - cost);
  } else { // additive 
    if (negLambdaZeroUtil && lambda<0) return 0;
    else return((spend-lambda)-(spend-lambda)*(spend-lambda)/(2*omega) - cost);
  }
}
double spendFn(const double lambda, const double omega, const double copay, 
               const int multMH) { 
  if (multMH) {
    if (lambda<=0) return 0;
    else return(lambda*(1+omega*(1-copay)));
  } else {
    return(lambda + (1-copay)*omega);
  }
}

void solveSpend(double *spend, double *util, double *oop, 
                const double lambda, const double omega, const double copay,
                const double deduct, const double maxoop,
                const int multMH)
{
  double m;
  *util = utilfn(0,lambda,omega,0,multMH);
  *spend = 0;
  //un = util;
  m = spendFn(lambda,omega,1.0,multMH);
  //if (m>0) un = omega/2-lambda;
  if (m > 0 && m < deduct) {
    double u0 = utilfn(m,lambda,omega,m,multMH);
    if (u0>(*util)) { (*util) = u0; *spend = m; *oop=m;}
  } 
  m = spendFn(lambda,omega,copay,multMH);
  if (m>=deduct && m < (deduct + (maxoop-deduct)/copay) && m>0) {  
    double thisOop = deduct + copay*(m-deduct);
    double u0 = utilfn(m,lambda,omega,thisOop,multMH);
    if (u0>*util) { *util = u0;  *spend=m; *oop=thisOop; }
  } 
  m = spendFn(lambda,omega,0,multMH);
  if (m>=deduct + (maxoop-deduct)/copay) {
    double u0 = utilfn(m,lambda,omega,maxoop,multMH);
    if (u0>*util) { *util = u0; *spend=m; *oop=maxoop;}
  }                
  return;  
}


int choiceFnEx
(const double omega, const double psi, const double muL, const double sigL, 
 const double *xint, const double *wint, const int nint,
 const double *deduct, const double *maxoop, const double *prem, const int nplan,
 const double copay, const double dl, const int multMH, 
 double *exVal, double *exSpend, double *exOop)
{
  int J=0;
  long double mev=-HUGE_VAL;
  int k;
  for(k=0;k<nplan;k++) {   
    if (deduct[k] < 0) continue;
    int i;
    long double ev = 0;
    exSpend[k] =0;
    exOop[k] = 0;
    for(i=0;i<nint;i++) {
      double lambda, oop, orate,m,util,spend,coop;
      int c;
      lambda = exp(xint[i]*sigL + muL) + dl;
      solveSpend(&spend,&util,&coop, 
                 lambda, omega, copay, deduct[k], maxoop[k], multMH);      
      ev += expl(-psi*util)*wint[i];
      exSpend[k] += spend*wint[i];
      exOop[k] += coop*wint[i];
      //evn += expl(-psi*un)*wint[i];
      //mexPrintf("%i: u=%g un=%g\n",i,util,un);
      if (!isfinite(ev)) {
        //mexPrintf("WARNING(choiceFn): overflow, switching to choiceFnHuge\n");
        return(choiceFnHuge(omega,  psi,  muL,  sigL,  xint,  wint, nint,
                            deduct,  maxoop,  prem, nplan,  copay, dl, multMH, 
                            exVal, exSpend,exOop));
      }      
    } // for(i<nint) 
    ev = -logl(ev) - psi*prem[k];
    //mexPrintf("%d: v=%Lg vn=%Lg\n",k,ev,-logl(evn));
    exVal[k]= ev; // + logl(evn))/psi;
    if (!isfinite(ev)) {
      //mexPrintf("WARNING(choiceFn): overflow, switching to choiceFnHuge\n");
      return(choiceFnHuge(omega,  psi,  muL,  sigL,  xint,  wint, nint,
                          deduct,  maxoop,  prem, nplan,  copay, dl, multMH, 
                          exVal, exSpend, exOop));
    }      
    //mexPrintf(" u(%d)=%Lg \n",k,ev);
    if (k==0 || ev>mev) { mev = ev; J = k+1; }    
  } // for(k<nplan) 
  //mexPrintf("\n");
  if (!isfinite(mev)) {
    J = 0; // mev is -inf or nan
    //mexPrintf("WARNING(choiceFn): max E[V] = %Lg\n",mev);
  }
  
  // for testing only
  /*{ 
    int Jhuge =choiceFnHuge(omega,  psi,  muL,  sigL,  xint,  wint, nint,
                            deduct,  maxoop,  prem, nplan,  copay); 
    if (J!=Jhuge) {
      mexPrintf("warning: choiceFn=%d, choiceFnHuge=%d psi=%g\n",J,Jhuge,psi);
    }
    } */ 
  return(J);
}


#include <mpfr.h>
#define MPFR_USE_FILE
// number of bits for mantissa of large numbers -- set to 53 to match double
#define PREC 53 
int choiceFnHuge
(const double omega, const double psi, const double muL, const double sigL, 
 const double *xint, const double *wint, const int nint,
 const double *deduct, const double *maxoop, const double *prem, const int nplan,
 const double copay, const double dl, const int multMH, 
 double *exVal, double *exSpend, double *exOop)
{
  int J=0;
  mpfr_t mev;
  int k;
  if (mpfr_buildopt_tls_p()==0)
    mexErrMsgTxt("This version of MPFR is not thread safe.\n");
  mpfr_init2(mev,PREC);
  for(k=0;k<nplan;k++) {
    if (deduct[k] < 0) continue;
    int i;
    mpfr_t ev;
    mpfr_init2(ev,PREC);
    mpfr_set_d(ev,0.0,GMP_RNDN);
    exSpend[k] = 0;
    exOop[k] = 0;
    for(i=0;i<nint;i++) {
      double lambda, oop, orate,m,util, spend,coop;
      int c;
      lambda = exp(xint[i]*sigL + muL) + dl;
      //if (lambda<0) lambda = 0;
      solveSpend(&spend,&util,&coop, 
                 lambda, omega, copay, deduct[k], maxoop[k], multMH);      
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
      exSpend[k] += wint[i]*spend;
      exOop[k] += wint[i]*coop;
    } // for(i<nint) 
    // ev = -log(ev)-psi*prem[k]
    mpfr_log(ev,ev,GMP_RNDN);
    mpfr_mul_si(ev,ev,-1,GMP_RNDN);
    mpfr_sub_d(ev,ev,psi*prem[k],GMP_RNDN);
    //mexPrintf("huge: u(%d)=%Lg \n",k,lev);
    //mexPrintf("huge: u(%d) = %Lg \n",k,mpfr_get_ld(ev,GMP_RNDN));
    exVal[k]=mpfr_get_d(ev,GMP_RNDN);
    //if (!isfinite(exVal[k])) mexPrintf("warning: exVal[k] not finite\n");
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


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////// RANDOM NUMBER GENERATION ///////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*  Parameters for: 
    Marsaglia & Tsang generator for random normals & random exponentials.
    Marsaglia, G. & Tsang, W.W. (2000) `The ziggurat method for generating random variables', J. Statist. Software, 
    v5(8). This is an electronic journal which can be downloaded from:  http://www.jstatsoft.org/v05/i08 */

static inline unsigned int SHR3(struct randState *rs) 
{ return(rs->jz=rs->jsr, 
         rs->jsr^=(rs->jsr<<13), 
         rs->jsr^=(rs->jsr>>17), 
         rs->jsr^=(rs->jsr<<5),
         rs->jz+rs->jsr); }
static inline double UNI(struct randState *rs) 
{ return(.5 + (signed) SHR3(rs)*.2328306e-9); }
#define IUNI(rs) SHR3(rs)
#define EPSILON  1.110326e-16

static double nfix(struct randState *rs)
{ // nfix() generates variates from the residue when rejection in RNOR occurs.
  const double r = 3.442620f;     /* The start of the right tail */
  static double x, y;
  for(;;) {
    x=rs->hz*rs->wn[rs->iz];      /* iz==0, handles the base strip */
    if(rs->iz==0) {
      do {
	x=-log(UNI(rs))*0.2904764; /* .2904764 is 1/r */
        y=-log(UNI(rs));
      } while(y+y<x*x);
      return (rs->hz>0)? r+x : -r-x;
      }
    /* iz>0, handle the wedges of other strips */
    if ( rs->fn[rs->iz]+UNI(rs)*(rs->fn[rs->iz-1]-rs->fn[rs->iz]) < exp(-.5*x*x) ) return x;
    /* initiate, try to exit for(;;) for loop*/
    rs->hz=SHR3(rs);
    rs->iz=rs->hz&127;
    if(fabs(rs->hz)<rs->kn[rs->iz]) return (rs->hz*rs->wn[rs->iz]);
  }
}

static double efix(struct randState *rs)
{ // efix() generates variates from the residue when rejection in REXP occurs. 
  double x;
  for(;;) {
    if(rs->iz==0) return (7.69711-log(UNI(rs)));          /* iz==0 */
    x=rs->jz*rs->we[rs->iz];
    if( rs->fe[rs->iz]+UNI(rs)*(rs->fe[rs->iz-1]-rs->fe[rs->iz]) < exp(-x) ) return (x);
    /* initiate, try to exit for(;;) loop */
    rs->jz=SHR3(rs);
    rs->iz=(rs->jz&255);
    if(rs->jz<rs->ke[rs->iz]) return (rs->jz*rs->we[rs->iz]);
  }
}

static inline double RNOR(struct randState *rs) 
{ return(rs->hz=SHR3(rs), rs->iz=rs->hz&127, 
         (fabs(rs->hz)<rs->kn[rs->iz])? rs->hz*rs->wn[rs->iz] : nfix(rs)); }
static inline double REXP(struct randState *rs) 
{ return(rs->jz=SHR3(rs), rs->iz=rs->jz&255, 
         (rs->jz <rs->ke[rs->iz])? rs->jz*rs->we[rs->iz] : efix(rs)); }

unsigned int get_jsr(struct randState *rs)
{
  return(rs->jsr);
}

void set_jsr(unsigned int jsrin, struct randState *rs)
{
  rs->jsr = jsrin;
  return;
}

// This helps when you want to start the seuqence at the same place again
void reset_seed(unsigned int jsrseed, struct randState *rs)
{ // This procedure sets the seed and creates the tables.  
  //************************** THIS IS WHAT MAKES IT DIFFERENT FROM SET_SEED ****************************************//
  rs->jsr = 123456789; 
  //*****************************************************************************************************************//
  rs->jsr^=jsrseed;
}

void set_seed(unsigned int jsrseed, struct randState *rs)
{ // This procedure sets the seed and creates the tables.  
  const double m1 = 2147483648.0, m2 = 4294967296.;
  double dn=3.442619855899,tn=dn,vn=9.91256303526217e-3, q;
  double de=7.697117470131487, te=de, ve=3.949659822581572e-3;
  int i;
  rs->jsr = 123456789;
  rs->jsr^=jsrseed;
  /* Set up tables for RNOR */
  q=vn/exp(-.5*dn*dn);
  rs->kn[0]=(int)((dn/q)*m1);
  rs->kn[1]=0;
  rs->wn[0]=q/m1;
  rs->wn[127]=dn/m1;
  rs->fn[0]=1.;
  rs->fn[127]=exp(-.5*dn*dn);
  for(i=126;i>=1;i--) {
    dn=sqrt(-2.*log(vn/dn+exp(-.5*dn*dn)));
    rs->kn[i+1]=(int)((dn/tn)*m1);
    tn=dn;
    rs->fn[i]=exp(-.5*dn*dn);
    rs->wn[i]=dn/m1;
  }
  /* Set up tables for REXP */
  q = ve/exp(-de);
  rs->ke[0]=(int)((de/q)*m2);
  rs->ke[1]=0;
  rs->we[0]=q/m2;
  rs->we[255]=de/m2;
  rs->fe[0]=1.;
  rs->fe[255]=exp(-de);
  for(i=254;i>=1;i--) {
    de=-log(ve/de+exp(-de));
    rs->ke[i+1]= (int)((de/te)*m2);
    te=de;
    rs->fe[i]=exp(-de);
    rs->we[i]=de/m2;
  }
  rs->seed_initialized=1;
}

void set_time_seed(RANDSTATE *rs)
{
  time_t t;
  set_seed((unsigned int) time(&t),rs);
}

double Sample_Uniform(const double a,const double b, struct randState *rs)
{ 
  time_t t;
  if (rs->seed_initialized!=1) set_seed((unsigned int) time(&t),rs);
  return(UNI(rs)*(b-a) + a);
}

double Sample_Normal(const double mu,const double var, struct randState *rs)
{
  time_t t;
  if (rs->seed_initialized!=1) set_seed((unsigned int) time(&t),rs);
  return(RNOR(rs)*sqrt(var) + mu);
}


#define T4 0.45
double Sample_Truncated_Normal(const double mu,const double var,const double a,const int lb, struct randState *rs) 
{ // Returns one draw from a truncated normal with underlying mean=mu and VARIANCE=var with truncation
  //           (a,+infty)    IF lb=TRUE(not 0)
  //           (-infty,a)    IF lb=FALSE(0)
  double u,z,phi_z,c;
  time_t t;
  if (rs->seed_initialized!=1) set_seed((unsigned int) time(&t),rs);
  c=(a-mu)/(sqrt(var));
  if (!lb) c=-c;
  if (c < T4) {
    // normal rejection sampling
    do {
      u=RNOR(rs);
    } while(u<=c);
  }
  else {
    // exponential rejection sampling
    do {
      u = UNI(rs);
      z = REXP(rs)/c;
      phi_z=exp(-.5*(z*z));
    } while(u>=phi_z);
    u=c + z;
  }
  return !lb ? (mu-(u*sqrt(var))): (mu+(u*sqrt(var)));
}
#undef T4

#define T1 0.375
#define T2 2.18
#define T3 0.725
#define T4 0.45
#define F(x) exp(-.5*(x*x))
double Sample_Double_Truncated_Normal(const double mu,const double var,const double a,const double b, struct randState *rs)
{ // Generates one draw from truncated standard normal distribution (mu,sigma) on (a,b)
  double c,c1,c2,u[2],x,cdel,f1,f2,z,az,bz;
  int lflip,j;
  time_t t;
  if (rs->seed_initialized!=1) set_seed((unsigned int) time(&t),rs);
  if (a>=b) mexPrintf("Sample_DTN: !(a<b), a=%.3g b=%.3g\n",a,b);
  c1=az=(a-mu)/sqrt(var);
  c2=bz=(b-mu)/sqrt(var);
  lflip=0;
  if (c1*c2<0.0) {
    if ((F(c1)>T1) && (F(c2)>T1)) {
      cdel=c2-c1;
      do {
	for(j=0;j<=1;j++) u[j] = UNI(rs);
	x=c1+cdel*u[0];
      } while (!(u[1]<F(x)));
    }
    else
      do {
	x=RNOR(rs);
      } while (!((x>c1) && (x<c2)));
  }
  else {
    if (c1<0.0) {
      c=c1;
      c1=-c2;
      c2 = -c;
      lflip=1;
    }
    f1=F(c1);
    f2=F(c2);
    if ( (f2<EPSILON) || (f1/f2>T2)) {
      if (c1>T3){
	// exponential rejection sampling
	c=c2-c1;
	do {
	  u[0] = UNI(rs);
	  z = REXP(rs)/c1;
	} while (!((z<c) && (u[0]<F(z))));
	x=c1+z;
      }
      else {
	// half-normal rejection sampling
	do {
	  x=RNOR(rs);
	  x=fabs(x);
	} while (!((x>c1) && (x<c2)));
      }
    }
    else {
      // uniform rejection sampling
      cdel=c2-c1;
      do {
	for(j = 0;j<= 1;j++) u[j] = UNI(rs);
	x=c1+cdel*u[0];
      } while (!(u[1] < F(x)/f1));
    }
  }
  return lflip ? (mu-sqrt(var)*x):(mu+sqrt(var)*x);
}
#undef F
#undef T1
#undef T2
#undef T3
#undef T4

double Sample_Gamma(const double ain,const double b, struct randState *rs)
{ // Returns a number distributed as a gamma distribution with shape parameter a and inverse scale b.
  // such that mean=a/b and var=a/(b*b). That is P(x) = [(b^a)/Gamma(a)] * [x^(a-1)] * exp(-x*b) for x>0 and a>=1.
  // Uses the algorithm in Marsaglia, G. and Tsang, W.W. (2000) `A simple method for generating
  // gamma variables', Trans. om Math. Software (TOMS), vol.26(3), pp.363-372.
  double d,c,x,v,u,a=ain,correction=1.0;
  time_t t;
  if (rs->seed_initialized==0) set_seed((unsigned int) time(&t),rs);
  if (a<1.0) {
    correction=pow(UNI(rs),1.0/a);
    a+=1.0;
  }
  d = a-1.0/3.0;
  c = 1.0/sqrt(9.0*d);
  for(;;) {
    do {
      x=RNOR(rs);
      v=1.0+c*x;
    } while(v<=0.0);
    v=v*v*v;
    u=UNI(rs);
    if( u<1.0-.0331*(x*x)*(x*x) ) return (correction*d*v/b);
    if( log(u)<0.5*x*x + d*(1.0-v+log(v)) ) return (correction*d*v/b);
  }
}

double gamcdf(const double x,const double a,const double b)
{ // Gives me the gamma cdf such that mean=a/b and var=a/(b*b)
  // That is F(x) = integral(0,x) of [(b^a)/Gamma(a)] * [x^(a-1)] * exp(-x*b) for x>0 and a>=1.
  double p;
  double _gammp(double, double);
  if (a<0.0) mexPrintf("a has to be positive: gamcdf");
  if (b<0.0) mexPrintf("b has to be positive: gamcdf");
  p=_gammp(a,x*b);
  return p>1.0 ? 1.0:p;
}

////////////////////////////////////////// Gamma Auxiliary Functions ////////////////////////////////////////////
#define SMALL 1.0e-250
double _gammln(double xx)
{
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,
			  0.1208650973866179e-2,-0.5395239384953e-5};
    int j;
    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}

double _gammp(double a,double x)
{
    void _gcf(double *gammcf, double a, double x, double *gln);
    void _gser(double *gamser, double a, double x, double *gln);
    double gamser,gammcf,gln;
    if (x < 0.0 || a <= 0.0) mexPrintf("Invalid arguments in routine gammp");
    if (x < (a+1.0)) {
	_gser(&gamser,a,x,&gln);
	return gamser;
    } else {
	_gcf(&gammcf,a,x,&gln);
	return 1.0-gammcf;
    }
}

#define ITMAX 2000
void _gser(double *gamser,double a,double x,double *gln)
{
    double _gammln(double xx);
    int n;
    double sum,del,ap;
    *gln=_gammln(a);
    if (x <= 0.0) {
	if (x < 0.0) mexPrintf("x less than 0 in routine gser");
	*gamser=0.0;
	return;
    } else {
	ap=a;
	del=sum=1.0/a;
	for (n=1;n<=ITMAX;n++) {
	    ++ap;
	    del *= x/ap;
	    sum += del;
	    if (fabs(del) < fabs(sum)*EPSILON) {
		*gamser=sum*exp(-x+a*log(x)-(*gln));
		return;
	    }
	}
	mexPrintf("a too large, ITMAX too small in routine gser");
	return;
    }
}

#define FPMIN SMALL/EPSILON
void _gcf(double *gammcf, double a, double x, double *gln)
{
    double _gammln(double xx);
    int i;
    double an,b,c,d,del,h;
    *gln=_gammln(a);
    b=x+1.0-a;
    c=1.0/FPMIN;
    d=1.0/b;
    h=d;
    for (i=1;i<=ITMAX;i++) {
	an = -i*(i-a);
	b += 2.0;
	d=an*d+b;
	if (fabs(d) < FPMIN) d=FPMIN;
	c=b+an/c;
	if (fabs(c) < FPMIN) c=FPMIN;
	d=1.0/d;
	del=d*c;
	h *= del;
	if (fabs(del-1.0) < EPSILON) break;
    }
    if (i > ITMAX) mexPrintf("a=%g too large, ITMAX=%d too small in gcf\n",a,ITMAX);
    *gammcf=exp(-x+a*log(x)-(*gln))*h;
}
#undef FPMIN
#undef ITMAX

#undef SMALL


float genbet(float aa,float bb, struct randState *rs)
/*
**********************************************************************
     float genbet(float aa,float bb)
               GeNerate BETa random deviate
                              Function
     Returns a single random deviate from the beta distribution with
     parameters A and B.  The density of the beta is
               x^(a-1) * (1-x)^(b-1) / B(a,b) for 0 < x < 1
                              Arguments
     aa --> First parameter of the beta distribution
       
     bb --> Second parameter of the beta distribution
       
                              Method
     R. C. H. Cheng
     Generating Beta Variatew with Nonintegral Shape Parameters
     Communications of the ACM, 21:317-322  (1978)
     (Algorithms BB and BC)
**********************************************************************
Modified by Paul Schrimpf to use thread-safe Sample_Uniform(0,1,rs) as RNG
*/
{
#define expmax 89.0
#define infnty 1.0E38
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
static float olda = -1.0;
static float oldb = -1.0;
static float genbet,a,alpha,b,beta,delta,gamma,k1,k2,r,s,t,u1,u2,v,w,y,z;
static long qsame;

    qsame = olda == aa && oldb == bb;
    if(qsame) goto S20;
    if(!(aa <= 0.0 || bb <= 0.0)) goto S10;
    //fputs(" AA or BB <= 0 in GENBET - Abort!",stderr);
    mexPrintf(" AA: %16.6E BB %16.6E\n",aa,bb);
    exit(1);
S10:
    olda = aa;
    oldb = bb;
S20:
    if(!(min(aa,bb) > 1.0)) goto S100;
/*
     Alborithm BB
     Initialize
*/
    if(qsame) goto S30;
    a = min(aa,bb);
    b = max(aa,bb);
    alpha = a+b;
    beta = sqrt((alpha-2.0)/(2.0*a*b-alpha));
    gamma = a+1.0/beta;
S30:
S40:
    u1 = Sample_Uniform(0,1,rs);
/*
     Step 1
*/
    u2 = Sample_Uniform(0,1,rs);
    v = beta*log(u1/(1.0-u1));
    if(!(v > expmax)) goto S50;
    w = infnty;
    goto S60;
S50:
    w = a*exp(v);
S60:
    z = pow(u1,2.0)*u2;
    r = gamma*v-1.3862944;
    s = a+r-w;
/*
     Step 2
*/
    if(s+2.609438 >= 5.0*z) goto S70;
/*
     Step 3
*/
    t = log(z);
    if(s > t) goto S70;
/*
     Step 4
*/
    if(r+alpha*log(alpha/(b+w)) < t) goto S40;
S70:
/*
     Step 5
*/
    if(!(aa == a)) goto S80;
    genbet = w/(b+w);
    goto S90;
S80:
    genbet = b/(b+w);
S90:
    goto S230;
S100:
/*
     Algorithm BC
     Initialize
*/
    if(qsame) goto S110;
    a = max(aa,bb);
    b = min(aa,bb);
    alpha = a+b;
    beta = 1.0/b;
    delta = 1.0+a-b;
    k1 = delta*(1.38889E-2+4.16667E-2*b)/(a*beta-0.777778);
    k2 = 0.25+(0.5+0.25/delta)*b;
S110:
S120:
    u1 = Sample_Uniform(0,1,rs);
/*
     Step 1
*/
    u2 = Sample_Uniform(0,1,rs);
    if(u1 >= 0.5) goto S130;
/*
     Step 2
*/
    y = u1*u2;
    z = u1*y;
    if(0.25*u2+z-y >= k1) goto S120;
    goto S170;
S130:
/*
     Step 3
*/
    z = pow(u1,2.0)*u2;
    if(!(z <= 0.25)) goto S160;
    v = beta*log(u1/(1.0-u1));
    if(!(v > expmax)) goto S140;
    w = infnty;
    goto S150;
S140:
    w = a*exp(v);
S150:
    goto S200;
S160:
    if(z >= k2) goto S120;
S170:
/*
     Step 4
     Step 5
*/
    v = beta*log(u1/(1.0-u1));
    if(!(v > expmax)) goto S180;
    w = infnty;
    goto S190;
S180:
    w = a*exp(v);
S190:
    if(alpha*(log(alpha/(b+w))+v)-1.3862944 < log(z)) goto S120;
S200:
/*
     Step 6
*/
    if(!(a == aa)) goto S210;
    genbet = w/(b+w);
    goto S220;
S210:
    genbet = b/(b+w);
S230:
S220:
    return genbet;
#undef expmax
#undef infnty
#undef ABS
#undef min
#undef max
}

/* Simulation of right and left truncated Gamma distributions by mixtures     */
/* by  Anne PHILIPPE                                                          */
/* CNRS UPRES-A 60 85, Universit\'e de Rouen                                  */
/* 76821 Mont Saint Aignan Cedex, France                                      */
/* Anne.Philippe@univ-rouen.fr                                                */
/******************************************************************************/
/*     algorithm to generate right truncated gamma distribution               */
/******************************************************************************/

/* Modified by Paul Schrimpf to use thread-safe Sample_Uniform(0,1,rs) as RNG */
/*  and genbet in place of beta() and Sample_Gamma in place of gamma()

/* the following function gives the optimal number of components              */
/* for p=0.95 fixed                                                           */
int n_optimal(b)                                              
     double b;
{
  double nr,q;
  q=1.64;                                                  
  nr=floor(pow(q+sqrt(q*q+4*b),2.0)/4.0-1.0); 
  return((int) nr);
}

/* the following function returns a random number from TG^-(a,b,1)            */ 
double inter_ri(a,b,rs)
     double a,b;
     struct randState *rs;
{
  int n,i,j;
  double x,u,test=0,y,z,yy,zz;
  double *wl,*wlc;
  n=n_optimal(b);
  wl=(double *)myCalloc((int)(n+2),sizeof(double));
  wlc=(double *)myCalloc((int)(n+2),sizeof(double));
  wl[0]=1.0;
  wlc[0]=1.0;
  for (i=1; i<=n;i++)
    {
      wl[i]=wl[i-1]*b/(a+i);
      wlc[i]=wlc[i-1]+wl[i];
    };
  for(i=0; i<=n; i++) 
    { 
      wlc[i]=wlc[i]/wlc[n];
    };
  y=1.0; yy=1.0;
  for (i=1; i<=n; i++) 
    {
      yy=yy*b/i; y=y+yy;
    };

  while (test == 0.0)
    {
      u=Sample_Uniform(0,1,rs);
      j=0;
      while(u>wlc[j]){j=j+1;};
      x=(double) genbet((float) a, (float) j+1,rs);
      u=Sample_Uniform(0,1,rs);
      zz=1.0; z=1.0;
      for (i=1; i<=n; i++)
        {
          zz=zz*(1-x)*b/i;
          z=z+zz;
        }; 
      z=exp(-b*x)*y/z;
      if (u <= z ) {test=1.0;};
    };
  myFree(wl);
  myFree(wlc);
  return(x);
};

/* the following function returns a random number from TG^-(a,b,1)            */
/* algorithm [A_1N]  N fixed such that P(N)> .95                              */
/* N fixed such that P(N)> .95                                                */
double gamma_right(a,b,t,rs)                             
     double a,b,t;         
     struct randState *rs;                        
{ 
  return(inter_ri(a,b*t,rs)*t);
};

/****************************************************************************/
/*algorithm to generate left truncated gamma distribution  */
/****************************************************************************/
/* the following function returns a random number from TG^+(a,b,1)  */
/* when a is an interger  */
/* the function gamma(a) returns a random number from the gamma  */
/* distribution G(a,1).q  */
/* See Devroye, L. (1985) Non-Uniform Random Variate Generation  */
/* Springer-Verlag, New-York.  */

double integer(a,b,rs)                             
     double a,b;
     struct randState *rs;
{
  double u,x;
  double *wl,*wlc;
  int i;
  wl=(double *)myCalloc((int)(a+1),sizeof(double));
  wlc=(double *)myCalloc((int)(a+1),sizeof(double));
  wl[1]=1.0;
  wlc[1]=1.0;
  for(i=2; i<= (int)a ;i++ )
    {
      wl[i]=wl[i-1]*(a-i+1)/b;
      wlc[i]=wlc[i-1]+wl[i];
    };
  for(i=1; i<= (int) a; i++)
    {
      wlc[i]=wlc[i]/wlc[(int)a];
    };
  u=Sample_Uniform(0,1,rs);
  i=1;
  while(u>wlc[i]) {
    i=i+1;
    if (i>a) {
      mexPrintf("ERROR (integer): invalid read."
                " i = %d a=%d u=%.4f a=%g  b=%g \n",i,(int) a, u,a,b);
      for(int k=0;k<=(int) a;k++) {
        mexPrintf("    wlc[%2d] = %.4f \n",k,wlc[k]);
      }
      mexErrMsgTxt("Going to segfault.");
    }
  }
  x=Sample_Gamma( (double) i,1,rs)/b+1.0;        
  myFree(wl);
  myFree(wlc);               
  return(x);
};

/* the following function returns a random number from TG^+(a,b,1)            */   
double inter_le(a,b,rs)
     double a,b;
     struct randState *rs;
{
  double test=0,u,x,y,M;
  if (a<1.0)
    {
      M=1.0;
      while (test == 0)
        {
          x=1-(1/b)*log(1-Sample_Uniform(0,1,rs));
          y=1/pow(x,1-a);
          if (Sample_Uniform(0,1,rs)< y/M) test=1.0;
        };
    }
  else
    {
      if (a<b) 
        {
          M=exp(floor(a)-a);
          while (test == 0)
            {
              x=integer(floor(a), b*floor(a)/a, rs);
              y=pow(x, a-floor(a))*exp(-x*b*(1-floor(a)/a));
              if (Sample_Uniform(0,1,rs)< y/M) test=1.0;
            };
        }
      else
        {
          M=exp(floor(a)-a)*pow(a/b,a-floor(a));
          while (test == 0)
            {
              x=integer(floor(a), b+floor(a)-a, rs);
              y=pow(x, a-floor(a))*exp(-x*(-floor(a)+a));
              if (Sample_Uniform(0,1,rs)< y/M) test=1.0;
            };
        };
    };
  return(x);
};



/* the following function returns a random number from TG^+(a,b,t)            */
/* algorithm [A_5]                                                            */
double gamma_left(a,b,t,rs)                                 
     double a,b,t;
     struct randState *rs;
{
  return(inter_le(a, b*t,rs)*t);
}


#define MAXITER 100
double Sample_Truncated_Gamma(const double ain, const double b, 
                              const double lo, const double hi, struct randState *rs) 
{ /*
    Returns a number distributed as a truncated gamma distribution
   with shape parameter a and inverse scale b.  such that mean=a/b and
   var=a/(b*b). That is P(x) = [(b^a)/Gamma(a)] * [x^(a-1)] *
   exp(-x*b) for x>0 and a>=1.  The distribution is truncated to
   [lo,hi].  Set lo<=0 and hi=INFINITY or NaN for no left/right
   truncation
  */
  double x;
  int iter=0;
  time_t t;
  if (rs->seed_initialized!=1) set_seed((unsigned int) time(&t),rs);
  if (lo<=0 && !isfinite(hi)) return(Sample_Gamma(ain,b,rs));
  else if (lo<=0) { // only right truncation
    // use method of Phillipe above
    x=gamma_right(ain,b,hi,rs);
  } else if (!isfinite(hi)) { // only left truncation
    if (lo<(ain/b + ain/(b*b))) {
      x = lo-1;
      while (x<lo && ++iter<MAXITER) x = Sample_Gamma(ain,b,rs);
      if (iter>=MAXITER) x = gamma_left(ain,b,lo,rs);
      return(x);
    } else {
      // use method of Phillipe above
      x = gamma_left(ain,b,lo,rs);     
      return(x);
    }
  } else { // truncated on both sides
    // simple rejection sampling
    x = hi+1;
    while ((x>hi || x<lo) && ++iter<MAXITER) {
      if ((iter%2)==0) x = gamma_right(ain,b,hi,rs);
      else x = gamma_left(ain,b,lo,rs);
    }
    // NOTE: this may be very efficient, but it isn't needed in the current code anyway ...    
  }
  if (iter>=MAXITER) {
    mexPrintf("Sample_Truncated_Gamma: max iter reached a=%g b=%g, lo=%g hi=%g\n",ain,b,lo,hi);
    return(lo*1.01);
  } else return(x);
}
#undef MAXITER


