#include <string.h>
#include "mex.h" /* header file for matlab API */

/* returns kronEye(x,N) = kron(x,eye(N)) 
   assumes memset(double *x,0, size_t bytes) sets x to array of 0.0's
    (this is true on x86 *nix and windows, but may not be true everywhere)
*/
void mexFunction
(int nlhs, /* number of left hand side arguments */
 mxArray *plhs[], /* pointer to lhs*/
 int nrhs, /* number of rhs arguments*/
 const mxArray *prhs[]) /* pointer to rhs arguments*/
{
  double *x, *kxi;
  int N;
  int ix,jx,mx,nx,n;
  mxAssert(nrhs==2,"2 rhs needed");
  mxAssert(nlhs==1,"1 lhs needed");
  x = mxGetPr(prhs[0]);
  mx = mxGetM(prhs[0]);
  nx = mxGetN(prhs[0]);
  N = (int) mxGetScalar(prhs[1]);
  plhs[0] = mxCreateDoubleMatrix(mx*N,nx*N,mxREAL);
  kxi = mxGetPr(plhs[0]);
  memset(kxi,0,sizeof(double)*mx*N*nx*N); 
  for(ix=0;ix<mx;ix++) {
    for(jx=0;jx<nx;jx++) {
      for(n=0;n<N;n++) {
        kxi[(ix*N+n) + N*mx*(jx*N+n)] = x[ix+mx*jx];
      }
    }
  }
  return;
}
