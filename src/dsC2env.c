// dsC2env-conversion;
//   at present limited to upper triangle storage (A@uplo = "U").
//  Note that Cholesky factors are of type 'dCholCMatrix' and stored in lower 
//   triangle form (uplo = "L").
//

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

void dsC2env(int *pnrows, double* Ax, int* Ai, int* Ap, double *diagvec, double *env, int *xenv) {
  
  int n = *pnrows;
  int i, j, k;

  //compute xenv --- obsolete as passed from R ----
//  int* xenv = calloc((n + 1), sizeof(int));
  int diff = 0;
//  if (uplo == "U") {
//    for (i = 0; i < n; i++) {  
//      diff += i - Ai[Ap[i]];
//      xenv[i + 1] = diff;
//    }
// }
//  if (uplo == "L") {  
//    for (i = n; i > 0; i--) {
//      diff += Ai[Ap[i] - 1] - (i - 1);
//      xenv[n - i + 1] = diff; 
//      //  printf("diff %d occurs at %d entry of xenv \n", diff, n-i+1);
//    }
//  }

  //compute env and diag
//  double* diagvec = calloc(n, sizeof(double));
//  double* env = calloc(xenv[n], sizeof(double));
  int pos = 0;
  int nenv;
  for (i = 0; i < n; i++) {  
//    if (uplo == "U") {
    diagvec[i] = Ax[Ap[i + 1] - 1];
//    } 
//    if (uplo == "L") {
//      diagvec[i] = Ax[Ap[i]];
//    }
    nenv = xenv[i + 1] - xenv[i];
    k = 0;
    for (j = 0; j < nenv; j++) {
      if ((Ai[Ap[i]] + j) == Ai[Ap[i] + k]) {
	env[pos + j] = Ax[Ap[i] + k];
	k += 1;
      } //endif
    } //endfor
    pos += nenv;
  }//endfor

//  free(diagvec); free(xenv); free(env);
return;

}
