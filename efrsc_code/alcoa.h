#ifndef _ALCOA_H
#define _ALCOA_H

int choiceFn
(const double omega, const double psi, const double muL, const double sigL, 
 const double *xint, const double *wint, const int nint,
 const double *deduct, const double *maxoop, const double *prem, const int nplan,
 const double copay, const int multMH, double *exVal, double *exSpend, double *exOop);


int choiceFnEx
(const double omega, const double psi, const double muL, const double sigL, 
 const double *xint, const double *wint, const int nint,
 const double *deduct, const double *maxoop, const double *prem, const int nplan,
 const double copay, const double dl, const int multMH, 
 double *exVal, double *exSpend, double *exOop);

void sigintHandler(int sig);


typedef struct randState {
  unsigned int jz;
  unsigned int jsr; //=123456789;
  int seed_initialized;
  int hz;
  unsigned int iz, kn[128], ke[256];
  double wn[128],fn[128], we[256],fe[256];
} RANDSTATE;

double Sample_Normal(const double mu,const double var, RANDSTATE *rs);
double Sample_Truncated_Normal(const double mu,const double var,const double a,const int lb, RANDSTATE *rs);
double Sample_Double_Truncated_Normal(const double mu,const double var,const double a,const double b, RANDSTATE *rs);
double Sample_Uniform(const double a,const double b, RANDSTATE *rs);
double Sample_Gamma(const double ain,const double b, RANDSTATE *rs);
double Sample_Truncated_Gamma(const double ain, const double b, 
                              const double lo, const double hi, struct randState *rs);
double gamcdf(const double x,const double a,const double b);

void set_seed(unsigned int jsrseed, struct randState *rs);
void set_time_seed(RANDSTATE *rs);
// Returns a number distributed as a gamma distribution with shape parameter a and inverse scale b.
// such that mean=a/b and var=a/(b*b). That is P(x) = [(b^a)/Gamma(a)] * [x^(a-1)] * exp(-x*b) for x>0 and a>=1. 

#define MAXTHREADS 50
#define myCalloc calloc
#define myFree free

// utility ignoring cost of spending
double utilfn(const double spend,const double lambda,const double omega,
              const double cost, const int multMH);
double spendFn(const double lambda, const double omega, const double copay, 
               const int multMH);
void setNegLambdaZeroUtil(int val);
#endif // ifndef _ALCOA_H
