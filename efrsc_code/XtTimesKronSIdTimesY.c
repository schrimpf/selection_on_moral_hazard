#include <string.h>
#include "mex.h" /* header file for matlab API */

/* returns f(x,S,N,y) = x'*kron(S,eye(N))*y
*/
void mexFunction
(int nlhs, /* number of left hand side arguments */
 mxArray *plhs[], /* pointer to lhs*/
 int nrhs, /* number of rhs arguments*/
 const mxArray *prhs[]) /* pointer to rhs arguments*/
{
  double *x, *S, *prod,*y;
  int N;
  int jx,jy,Mx,Nx,n,Ns,is,js,Ny;
  mxAssert(nrhs==4,"4 rhs needed");
  mxAssert(nlhs==1,"1 lhs needed");
  x = mxGetPr(prhs[0]);
  Mx = mxGetM(prhs[0]);
  Nx = mxGetN(prhs[0]);
  S = mxGetPr(prhs[1]);
  Ns = mxGetN(prhs[1]);
  mxAssert(Ns==mxGetM(prhs[1]),"S must be square");
  N = (int) mxGetScalar(prhs[2]);
  mxAssert(N*Ns==Mx,"Bad dimensions");
  y = mxGetPr(prhs[3]);
  mxAssert(Mx==mxGetM(prhs[3]),"rows(y)!=rows(x)");
  Ny = mxGetN(prhs[3]);
  plhs[0] = mxCreateDoubleMatrix(Nx,Ny,mxREAL);
  prod = mxGetPr(plhs[0]);
  for(jx=0;jx<Nx;jx++) {
    for(jy=0;jy<Ny;jy++) {
      prod[jx + jy*Nx] = 0;
      for(is=0;is<Ns;is++) {
        for(js=0;js<Ns;js++) {
          for(n=0;n<N;n++) {
            int ix = (is*N+n);
            int iy = (js*N+n);
            prod[jx+jy*Nx] += x[ix + jx*Mx]*S[is+Ns*js]*y[iy+jy*Mx];
          }
        }
      }
    }
  }
  return;
}
