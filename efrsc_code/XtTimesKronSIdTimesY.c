#include <string.h>
#include "mex.h" /* header file for matlab API */

/* returns f(x,S,N,y) = x'*kron(S,eye(N))*y
*/
/* Fortran BLAS dgemm declaration */
extern void dgemm_(const char *transa, const char *transb,
                   const int *m, const int *n, const int *k,
                   const double *alpha, const double *a, const int *lda,
                   const double *b, const int *ldb,
                   const double *beta, double *c, const int *ldc);

/* returns f(x,S,N,y) = x'*kron(S,eye(N))*y
   Optimized using block linear combination and BLAS DGEMM:
   1. Z = kron(S, eye(N)) * y
      For each block a in 0..Ns-1:
        Za = sum_{b=0}^{Ns-1} S(a, b) * Yb
   2. prod = x' * Z  (computed via dgemm)
*/
void mexFunction
(int nlhs, /* number of left hand side arguments */
 mxArray *plhs[], /* pointer to lhs*/
 int nrhs, /* number of rhs arguments*/
 const mxArray *prhs[]) /* pointer to rhs arguments*/
{
  double *x, *S, *prod, *y;
  int N;
  int Mx, Nx, Ns, Ny;
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

  /* Allocate intermediate matrix Z of size Mx x Ny */
  double *Z = (double*) mxCalloc(Mx * Ny, sizeof(double));

  /* Compute Z = kron(S, eye(N)) * y block by block */
  for (int jy = 0; jy < Ny; jy++) {
    const double *y_col = y + jy * Mx;
    double *Z_col = Z + jy * Mx;
    for (int a = 0; a < Ns; a++) {
      double *Za = Z_col + a * N;
      for (int b = 0; b < Ns; b++) {
        double Sab = S[a + Ns * b];
        if (Sab == 0.0) continue;
        const double *Yb = y_col + b * N;
        for (int n = 0; n < N; n++) {
          Za[n] += Sab * Yb[n];
        }
      }
    }
  }

  /* Compute prod = x' * Z using BLAS DGEMM:
     x is Mx x Nx, so x' is Nx x Mx.
     Z is Mx x Ny.
     prod is Nx x Ny.
  */
  char transa = 'T';
  char transb = 'N';
  int m = Nx;
  int n = Ny;
  int k = Mx;
  double alpha = 1.0;
  double beta_val = 0.0;
  int lda = Mx;
  int ldb = Mx;
  int ldc = Nx;
  dgemm_(&transa, &transb, &m, &n, &k, &alpha, x, &lda, Z, &ldb, &beta_val, prod, &ldc);

  mxFree(Z);
  return;
}
