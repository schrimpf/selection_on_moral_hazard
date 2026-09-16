/* 
   Mex-function to sample shape parameter of gamma distribution
   
   to compile: make.m
*/

#include <math.h>
#include <string.h>
#include "mex.h" /* header file for matlab API */
#include "arms.h"
#include "alcoa.h"

struct gammaParm {
  double logScale;
  double scale;
  double sumlogx;
  double lo;
  int N;
};
double logDensity(double shape,void *data) {
  struct gammaParm *gp=(struct gammaParm*) data;
  return(shape * (gp->sumlogx  - gp->N*gp->logScale) - gp->sumlogx 
         - gp->N * lgamma(shape) 
         - gp->N*log(1-gamcdf(gp->lo,shape,1/exp(gp->scale))) );
}


void mexFunction
(int nlhs, /* number of left hand side arguments */
 mxArray *plhs[], /* pointer to lhs*/
 int nrhs, /* number of rhs arguments*/
 const mxArray *prhs[]) /* pointer to rhs arguments*/
{
  double shape0, *shape;
  struct gammaParm gp;
  double *x;
  int i;
  
  // stuff for arms
  double xinit[5] = {1, 2, 3, 4, 5};
  int ninit = 5;
  int neval = 0;
  double lo = 0.01, hi=100;
  double convex = 1.0;
  
  shape0 = mxGetScalar(prhs[0]);
  gp.logScale = log(mxGetScalar(prhs[1]));
  gp.scale=mxGetScalar(prhs[1]);
  x = mxGetPr(prhs[2]);
  gp.N = mxGetNumberOfElements(prhs[2]);
  gp.lo = mxGetScalar(prhs[3]);
  plhs[0] = mxCreateDoubleScalar(0);
  shape = mxGetPr(plhs[0]);
  for(i=0,gp.sumlogx=0;i<gp.N;i++) {
    gp.sumlogx += log(x[i]);
  }
  arms(xinit,ninit,&lo,&hi, logDensity, &gp, &convex,
       50, 1, &shape0, shape, 1, NULL, NULL, 0,
       &neval);
  return;
}
