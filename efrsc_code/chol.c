// computes cholesky decomposition using lapack, created as a
// workaround for octave's chol, which hangs unexpectedly
// given symmetric, real A, returns upper triangular U such that U'U = A;
#include <string.h>
#include "mex.h" /* header file for matlab API */

void dpotrf_(char*,int*,double*,int*,int*);

void mexFunction
(int nlhs, /* number of left hand side arguments */
 mxArray *plhs[], /* pointer to lhs*/
 int nrhs, /* number of rhs arguments*/
 const mxArray *prhs[]) /* pointer to rhs arguments*/
{
  char uplo='U';
  int N;
  double *A; 
  int info;
  mxAssert(nrhs==1,"Exactly one input required.");
  mxAssert(nlhs==1 || nlhs==2,"One or two output required.");
  N = mxGetN(prhs[0]);
  mxAssert(N==mxGetM(prhs[0]),"Input matrix must be square");  
  plhs[0] = mxCreateDoubleMatrix(N,N,mxREAL);
  A = mxGetPr(plhs[0]);
  memcpy(A,mxGetPr(prhs[0]),N*N*sizeof(double)); // Output=Input
  dpotrf_(&uplo,&N,A,&N,&info); // A = chol(A)
  if (nlhs==2) {
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    *(mxGetPr(plhs[1])) = (double) info;
  } else mxAssert(info==0,"dpotrf_ returned info!=0");
  for(int i=1;i<N;i++)
    for(int j=0;j<i;j++) A[i+j*N] = 0;
}
