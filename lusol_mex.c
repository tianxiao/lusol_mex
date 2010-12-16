#include "mex.h"
#include <string.h>
#include <stdlib.h>

#define LU1FAC_FUNC lu1fac_
#define LU6SOL_FUNC lu6sol_
#define LU6MUL_FUNC lu6mul_
#define LU8RPC_FUNC lu8rpc_
#define LU8ADC_FUNC lu8adc_
#define LU8ADR_FUNC lu8adr_
#define LU8DLC_FUNC lu8dlc_
#define LU8DLR_FUNC lu8dlr_
#define LU8MOD_FUNC lu8mod_
#define LU8RPR_FUNC lu8rpr_

void gateway_lu1fac(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void gateway_lu6sol(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void gateway_lu6mul(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void gateway_lu8rpc(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void gateway_lu8adc(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void gateway_lu8adr(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void gateway_lu8dlc(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void gateway_lu8dlr(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void gateway_lu8mod(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void gateway_lu8rpr(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

void LU1FAC_FUNC(int *m, int *n, int *nelem, int *lena, int *luparm, double *parmlu, double *a, int *indc, int *indr, int *ip, int *iq, int *lenc, int *lenr, int *locc, int *locr, int *iploc, int *iqloc, int *ipinv, int *iqinv, double *w, int *inform);
void LU6SOL_FUNC(int *mode, int *m, int *n, double *v, double *w, int *lena, int *luparm, double *parmlu, double *a, int *indc, int *indr, int *ip, int *iq, int *lenc, int *lenr, int *locc, int *locr, int *inform);
void LU6MUL_FUNC(int *mode, int *m, int *n, double *v, double *w, int *lena, int *luparm, double *parmlu, double *a, int *indc, int *indr, int *ip, int *iq, int *lenc, int *lenr, int *locc, int *locr);
void LU8RPC_FUNC(int *mode1, int *mode2, int *m, int *n, int *jrep, double *v, double *w, int *lena, int *luparm, double *parmlu, double *a, int *indc, int *indr, int *ip, int *iq, int *lenc, int *lenr, int *locc, int *locr, int *inform, double *diag, double *vnorm);
void LU8ADC_FUNC(int *mode, int *m, int *n, double *v, double *w, int *lena, int *luparm, double *parmlu, double *a, int *indc, int *indr, int *ip, int *iq, int *lenc, int *lenr, int *locc, int *locr, int *inform, double *diag, double *vnorm);
void LU8ADR_FUNC(int *m, int *n, double *w, int *lena, int *luparm, double *parmlu, double *a, int *indc, int *indr, int *ip, int *iq, int *lenc, int *lenr, int *locc, int *locr, int *inform, double *diag);
void LU8DLC_FUNC(int *m, int *n, int *jdel, int *lena, int *luparm, double *parmlu, double *a, int *indc, int *indr, int *ip, int *iq, int *lenc, int *lenr, int *locc, int *locr, int *inform);
void LU8DLR_FUNC(int *mode, int *m, int *n, int *idel, double *v, double *w, int *lena, int *luparm, double *parmlu, double *a, int *indc, int *indr, int *ip, int *iq, int *lenc, int *lenr, int *locc, int *locr, int *inform);
void LU8MOD_FUNC(int *mode, int *m, int *n, double *beta, double *v, double *w, int *lena, int *luparm, double *parmlu, double *a, int *indc, int *indr, int *ip, int *iq, int *lenc, int *lenr, int *locc, int *locr, int *inform);
void LU8RPR_FUNC(int *mode1, int *mode2, int *m, int *n, int *irep, double *v, double *w, double *wnew, int *lena, int *luparm, double *parmlu, double *a, int *indc, int *indr, int *ip, int *iq, int *lenc, int *lenr, int *locc, int *locr, int *inform);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { 
  if (nrhs < 1) {
    mexErrMsgIdAndTxt("lusol_mex:gateway","there are no input arguments.");
  }

  if (!mxIsChar(prhs[0])) {
    mexErrMsgIdAndTxt("lusol_mex:gateway","the first argument must be a string indicating the subroutine to call.");
  }

  mwSize strM = mxGetM(prhs[0]);
  mwSize strN = mxGetN(prhs[0]);
  if ( strM != 1 || strN < 1 ) {
    mexErrMsgIdAndTxt("lusol_mex:gateway","the input string has incorrect size.");
  }
  char *subroutine_name = (char*) mxMalloc(sizeof(char)*(strN+1));
  if (!subroutine_name) mexErrMsgIdAndTxt("lusol_mex:gateway","could not allocate memory for string.");
  int check;
  check = mxGetString(prhs[0],subroutine_name,strN+1);
  if (check) {
    mexErrMsgIdAndTxt("lusol_mex:gateway","could not read input string.");
  }

  if (strcmp(subroutine_name,"lu1fac") == 0) {
    gateway_lu1fac(nlhs, plhs, nrhs, prhs);
    mxFree(subroutine_name);
    return;
  }
  if (strcmp(subroutine_name,"lu6sol") == 0) {
    gateway_lu6sol(nlhs, plhs, nrhs, prhs);
    mxFree(subroutine_name);
    return;
  }
  if (strcmp(subroutine_name,"lu6mul") == 0) {
    gateway_lu6mul(nlhs, plhs, nrhs, prhs);
    mxFree(subroutine_name);
    return;
  }
  if (strcmp(subroutine_name,"lu8rpc") == 0) {
    gateway_lu8rpc(nlhs, plhs, nrhs, prhs);
    mxFree(subroutine_name);
    return;
  }
  if (strcmp(subroutine_name,"lu8adc") == 0) {
    gateway_lu8adc(nlhs, plhs, nrhs, prhs);
    mxFree(subroutine_name);
    return;
  }
  if (strcmp(subroutine_name,"lu8adr") == 0) {
    gateway_lu8adr(nlhs, plhs, nrhs, prhs);
    mxFree(subroutine_name);
    return;
  }
  if (strcmp(subroutine_name,"lu8dlc") == 0) {
    gateway_lu8dlc(nlhs, plhs, nrhs, prhs);
    mxFree(subroutine_name);
    return;
  }
  if (strcmp(subroutine_name,"lu8dlr") == 0) {
    gateway_lu8dlr(nlhs, plhs, nrhs, prhs);
    mxFree(subroutine_name);
    return;
  }
  if (strcmp(subroutine_name,"lu8mod") == 0) {
    gateway_lu8mod(nlhs, plhs, nrhs, prhs);
    mxFree(subroutine_name);
    return;
  }
  if (strcmp(subroutine_name,"lu8rpr") == 0) {
    gateway_lu8rpr(nlhs, plhs, nrhs, prhs);
    mxFree(subroutine_name);
    return;
  }
  mxFree(subroutine_name);
  mexErrMsgIdAndTxt("lusol_mex:gateway","the input string does not match any subroutine.");
}

void gateway_lu1fac(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  /* check nrhs */
  if (nrhs != 22) {
    mexErrMsgIdAndTxt("lusol_mex:lu1fac","incorrect number of inputs for lu1fac.");
  }

  /* check prhs[1], variable: m */
  if (!mxIsClass(prhs[1],"int32") || mxGetM(prhs[1]) < 1 || mxGetN(prhs[1]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu1fac","m must have type int32 and size at least (1,1)");
  }
  int *m = (int*) mxGetData(prhs[1]);

  /* check prhs[2], variable: n */
  if (!mxIsClass(prhs[2],"int32") || mxGetM(prhs[2]) < 1 || mxGetN(prhs[2]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu1fac","n must have type int32 and size at least (1,1)");
  }
  int *n = (int*) mxGetData(prhs[2]);

  /* check prhs[3], variable: nelem */
  if (!mxIsClass(prhs[3],"int32") || mxGetM(prhs[3]) < 1 || mxGetN(prhs[3]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu1fac","nelem must have type int32 and size at least (1,1)");
  }
  int *nelem = (int*) mxGetData(prhs[3]);

  /* check prhs[4], variable: lena */
  if (!mxIsClass(prhs[4],"int32") || mxGetM(prhs[4]) < 1 || mxGetN(prhs[4]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu1fac","lena must have type int32 and size at least (1,1)");
  }
  int *lena = (int*) mxGetData(prhs[4]);

  /* check prhs[5], variable: luparm */
  if (!mxIsClass(prhs[5],"int32") || mxGetM(prhs[5]) < 30 || mxGetN(prhs[5]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu1fac","luparm must have type int32 and size at least (30,1)");
  }
  int *luparm = (int*) mxGetData(prhs[5]);

  /* check prhs[6], variable: parmlu */
  if (!mxIsClass(prhs[6],"double") || mxGetM(prhs[6]) < 30 || mxGetN(prhs[6]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu1fac","parmlu must have type double and size at least (30,1)");
  }
  double *parmlu = (double*) mxGetData(prhs[6]);

  /* check prhs[7], variable: a */
  if (!mxIsClass(prhs[7],"double") || mxGetM(prhs[7]) < *lena || mxGetN(prhs[7]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu1fac","a must have type double and size at least (lena,1)");
  }
  double *a = (double*) mxGetData(prhs[7]);

  /* check prhs[8], variable: indc */
  if (!mxIsClass(prhs[8],"int32") || mxGetM(prhs[8]) < *lena || mxGetN(prhs[8]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu1fac","indc must have type int32 and size at least (lena,1)");
  }
  int *indc = (int*) mxGetData(prhs[8]);

  /* check prhs[9], variable: indr */
  if (!mxIsClass(prhs[9],"int32") || mxGetM(prhs[9]) < *lena || mxGetN(prhs[9]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu1fac","indr must have type int32 and size at least (lena,1)");
  }
  int *indr = (int*) mxGetData(prhs[9]);

  /* check prhs[10], variable: ip */
  if (!mxIsClass(prhs[10],"int32") || mxGetM(prhs[10]) < *m || mxGetN(prhs[10]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu1fac","ip must have type int32 and size at least (m,1)");
  }
  int *ip = (int*) mxGetData(prhs[10]);

  /* check prhs[11], variable: iq */
  if (!mxIsClass(prhs[11],"int32") || mxGetM(prhs[11]) < *n || mxGetN(prhs[11]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu1fac","iq must have type int32 and size at least (n,1)");
  }
  int *iq = (int*) mxGetData(prhs[11]);

  /* check prhs[12], variable: lenc */
  if (!mxIsClass(prhs[12],"int32") || mxGetM(prhs[12]) < *n || mxGetN(prhs[12]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu1fac","lenc must have type int32 and size at least (n,1)");
  }
  int *lenc = (int*) mxGetData(prhs[12]);

  /* check prhs[13], variable: lenr */
  if (!mxIsClass(prhs[13],"int32") || mxGetM(prhs[13]) < *m || mxGetN(prhs[13]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu1fac","lenr must have type int32 and size at least (m,1)");
  }
  int *lenr = (int*) mxGetData(prhs[13]);

  /* check prhs[14], variable: locc */
  if (!mxIsClass(prhs[14],"int32") || mxGetM(prhs[14]) < *n || mxGetN(prhs[14]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu1fac","locc must have type int32 and size at least (n,1)");
  }
  int *locc = (int*) mxGetData(prhs[14]);

  /* check prhs[15], variable: locr */
  if (!mxIsClass(prhs[15],"int32") || mxGetM(prhs[15]) < *m || mxGetN(prhs[15]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu1fac","locr must have type int32 and size at least (m,1)");
  }
  int *locr = (int*) mxGetData(prhs[15]);

  /* check prhs[16], variable: iploc */
  if (!mxIsClass(prhs[16],"int32") || mxGetM(prhs[16]) < *n || mxGetN(prhs[16]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu1fac","iploc must have type int32 and size at least (n,1)");
  }
  int *iploc = (int*) mxGetData(prhs[16]);

  /* check prhs[17], variable: iqloc */
  if (!mxIsClass(prhs[17],"int32") || mxGetM(prhs[17]) < *m || mxGetN(prhs[17]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu1fac","iqloc must have type int32 and size at least (m,1)");
  }
  int *iqloc = (int*) mxGetData(prhs[17]);

  /* check prhs[18], variable: ipinv */
  if (!mxIsClass(prhs[18],"int32") || mxGetM(prhs[18]) < *m || mxGetN(prhs[18]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu1fac","ipinv must have type int32 and size at least (m,1)");
  }
  int *ipinv = (int*) mxGetData(prhs[18]);

  /* check prhs[19], variable: iqinv */
  if (!mxIsClass(prhs[19],"int32") || mxGetM(prhs[19]) < *n || mxGetN(prhs[19]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu1fac","iqinv must have type int32 and size at least (n,1)");
  }
  int *iqinv = (int*) mxGetData(prhs[19]);

  /* check prhs[20], variable: w */
  if (!mxIsClass(prhs[20],"double") || mxGetM(prhs[20]) < *n || mxGetN(prhs[20]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu1fac","w must have type double and size at least (n,1)");
  }
  double *w = (double*) mxGetData(prhs[20]);

  /* check prhs[21], variable: inform */
  if (!mxIsClass(prhs[21],"int32") || mxGetM(prhs[21]) < 1 || mxGetN(prhs[21]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu1fac","inform must have type int32 and size at least (1,1)");
  }
  int *inform = (int*) mxGetData(prhs[21]);

  LU1FAC_FUNC(m, n, nelem, lena, luparm, parmlu, a, indc, indr, ip, iq, lenc, lenr, locc, locr, iploc, iqloc, ipinv, iqinv, w, inform);
}

void gateway_lu6sol(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  /* check nrhs */
  if (nrhs != 19) {
    mexErrMsgIdAndTxt("lusol_mex:lu6sol","incorrect number of inputs for lu6sol.");
  }

  /* check prhs[1], variable: mode */
  if (!mxIsClass(prhs[1],"int32") || mxGetM(prhs[1]) < 1 || mxGetN(prhs[1]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu6sol","mode must have type int32 and size at least (1,1)");
  }
  int *mode = (int*) mxGetData(prhs[1]);

  /* check prhs[2], variable: m */
  if (!mxIsClass(prhs[2],"int32") || mxGetM(prhs[2]) < 1 || mxGetN(prhs[2]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu6sol","m must have type int32 and size at least (1,1)");
  }
  int *m = (int*) mxGetData(prhs[2]);

  /* check prhs[3], variable: n */
  if (!mxIsClass(prhs[3],"int32") || mxGetM(prhs[3]) < 1 || mxGetN(prhs[3]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu6sol","n must have type int32 and size at least (1,1)");
  }
  int *n = (int*) mxGetData(prhs[3]);

  /* check prhs[4], variable: v */
  if (!mxIsClass(prhs[4],"double") || mxGetM(prhs[4]) < *m || mxGetN(prhs[4]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu6sol","v must have type double and size at least (m,1)");
  }
  double *v = (double*) mxGetData(prhs[4]);

  /* check prhs[5], variable: w */
  if (!mxIsClass(prhs[5],"double") || mxGetM(prhs[5]) < *n || mxGetN(prhs[5]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu6sol","w must have type double and size at least (n,1)");
  }
  double *w = (double*) mxGetData(prhs[5]);

  /* check prhs[6], variable: lena */
  if (!mxIsClass(prhs[6],"int32") || mxGetM(prhs[6]) < 1 || mxGetN(prhs[6]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu6sol","lena must have type int32 and size at least (1,1)");
  }
  int *lena = (int*) mxGetData(prhs[6]);

  /* check prhs[7], variable: luparm */
  if (!mxIsClass(prhs[7],"int32") || mxGetM(prhs[7]) < 30 || mxGetN(prhs[7]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu6sol","luparm must have type int32 and size at least (30,1)");
  }
  int *luparm = (int*) mxGetData(prhs[7]);

  /* check prhs[8], variable: parmlu */
  if (!mxIsClass(prhs[8],"double") || mxGetM(prhs[8]) < 30 || mxGetN(prhs[8]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu6sol","parmlu must have type double and size at least (30,1)");
  }
  double *parmlu = (double*) mxGetData(prhs[8]);

  /* check prhs[9], variable: a */
  if (!mxIsClass(prhs[9],"double") || mxGetM(prhs[9]) < *lena || mxGetN(prhs[9]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu6sol","a must have type double and size at least (lena,1)");
  }
  double *a = (double*) mxGetData(prhs[9]);

  /* check prhs[10], variable: indc */
  if (!mxIsClass(prhs[10],"int32") || mxGetM(prhs[10]) < *lena || mxGetN(prhs[10]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu6sol","indc must have type int32 and size at least (lena,1)");
  }
  int *indc = (int*) mxGetData(prhs[10]);

  /* check prhs[11], variable: indr */
  if (!mxIsClass(prhs[11],"int32") || mxGetM(prhs[11]) < *lena || mxGetN(prhs[11]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu6sol","indr must have type int32 and size at least (lena,1)");
  }
  int *indr = (int*) mxGetData(prhs[11]);

  /* check prhs[12], variable: ip */
  if (!mxIsClass(prhs[12],"int32") || mxGetM(prhs[12]) < *m || mxGetN(prhs[12]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu6sol","ip must have type int32 and size at least (m,1)");
  }
  int *ip = (int*) mxGetData(prhs[12]);

  /* check prhs[13], variable: iq */
  if (!mxIsClass(prhs[13],"int32") || mxGetM(prhs[13]) < *n || mxGetN(prhs[13]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu6sol","iq must have type int32 and size at least (n,1)");
  }
  int *iq = (int*) mxGetData(prhs[13]);

  /* check prhs[14], variable: lenc */
  if (!mxIsClass(prhs[14],"int32") || mxGetM(prhs[14]) < *n || mxGetN(prhs[14]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu6sol","lenc must have type int32 and size at least (n,1)");
  }
  int *lenc = (int*) mxGetData(prhs[14]);

  /* check prhs[15], variable: lenr */
  if (!mxIsClass(prhs[15],"int32") || mxGetM(prhs[15]) < *m || mxGetN(prhs[15]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu6sol","lenr must have type int32 and size at least (m,1)");
  }
  int *lenr = (int*) mxGetData(prhs[15]);

  /* check prhs[16], variable: locc */
  if (!mxIsClass(prhs[16],"int32") || mxGetM(prhs[16]) < *n || mxGetN(prhs[16]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu6sol","locc must have type int32 and size at least (n,1)");
  }
  int *locc = (int*) mxGetData(prhs[16]);

  /* check prhs[17], variable: locr */
  if (!mxIsClass(prhs[17],"int32") || mxGetM(prhs[17]) < *m || mxGetN(prhs[17]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu6sol","locr must have type int32 and size at least (m,1)");
  }
  int *locr = (int*) mxGetData(prhs[17]);

  /* check prhs[18], variable: inform */
  if (!mxIsClass(prhs[18],"int32") || mxGetM(prhs[18]) < 1 || mxGetN(prhs[18]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu6sol","inform must have type int32 and size at least (1,1)");
  }
  int *inform = (int*) mxGetData(prhs[18]);

  LU6SOL_FUNC(mode, m, n, v, w, lena, luparm, parmlu, a, indc, indr, ip, iq, lenc, lenr, locc, locr, inform);
}

void gateway_lu6mul(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  /* check nrhs */
  if (nrhs != 18) {
    mexErrMsgIdAndTxt("lusol_mex:lu6mul","incorrect number of inputs for lu6mul.");
  }

  /* check prhs[1], variable: mode */
  if (!mxIsClass(prhs[1],"int32") || mxGetM(prhs[1]) < 1 || mxGetN(prhs[1]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu6mul","mode must have type int32 and size at least (1,1)");
  }
  int *mode = (int*) mxGetData(prhs[1]);

  /* check prhs[2], variable: m */
  if (!mxIsClass(prhs[2],"int32") || mxGetM(prhs[2]) < 1 || mxGetN(prhs[2]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu6mul","m must have type int32 and size at least (1,1)");
  }
  int *m = (int*) mxGetData(prhs[2]);

  /* check prhs[3], variable: n */
  if (!mxIsClass(prhs[3],"int32") || mxGetM(prhs[3]) < 1 || mxGetN(prhs[3]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu6mul","n must have type int32 and size at least (1,1)");
  }
  int *n = (int*) mxGetData(prhs[3]);

  /* check prhs[4], variable: v */
  if (!mxIsClass(prhs[4],"double") || mxGetM(prhs[4]) < *m || mxGetN(prhs[4]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu6mul","v must have type double and size at least (m,1)");
  }
  double *v = (double*) mxGetData(prhs[4]);

  /* check prhs[5], variable: w */
  if (!mxIsClass(prhs[5],"double") || mxGetM(prhs[5]) < *n || mxGetN(prhs[5]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu6mul","w must have type double and size at least (n,1)");
  }
  double *w = (double*) mxGetData(prhs[5]);

  /* check prhs[6], variable: lena */
  if (!mxIsClass(prhs[6],"int32") || mxGetM(prhs[6]) < 1 || mxGetN(prhs[6]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu6mul","lena must have type int32 and size at least (1,1)");
  }
  int *lena = (int*) mxGetData(prhs[6]);

  /* check prhs[7], variable: luparm */
  if (!mxIsClass(prhs[7],"int32") || mxGetM(prhs[7]) < 30 || mxGetN(prhs[7]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu6mul","luparm must have type int32 and size at least (30,1)");
  }
  int *luparm = (int*) mxGetData(prhs[7]);

  /* check prhs[8], variable: parmlu */
  if (!mxIsClass(prhs[8],"double") || mxGetM(prhs[8]) < 30 || mxGetN(prhs[8]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu6mul","parmlu must have type double and size at least (30,1)");
  }
  double *parmlu = (double*) mxGetData(prhs[8]);

  /* check prhs[9], variable: a */
  if (!mxIsClass(prhs[9],"double") || mxGetM(prhs[9]) < *lena || mxGetN(prhs[9]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu6mul","a must have type double and size at least (lena,1)");
  }
  double *a = (double*) mxGetData(prhs[9]);

  /* check prhs[10], variable: indc */
  if (!mxIsClass(prhs[10],"int32") || mxGetM(prhs[10]) < *lena || mxGetN(prhs[10]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu6mul","indc must have type int32 and size at least (lena,1)");
  }
  int *indc = (int*) mxGetData(prhs[10]);

  /* check prhs[11], variable: indr */
  if (!mxIsClass(prhs[11],"int32") || mxGetM(prhs[11]) < *lena || mxGetN(prhs[11]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu6mul","indr must have type int32 and size at least (lena,1)");
  }
  int *indr = (int*) mxGetData(prhs[11]);

  /* check prhs[12], variable: ip */
  if (!mxIsClass(prhs[12],"int32") || mxGetM(prhs[12]) < *m || mxGetN(prhs[12]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu6mul","ip must have type int32 and size at least (m,1)");
  }
  int *ip = (int*) mxGetData(prhs[12]);

  /* check prhs[13], variable: iq */
  if (!mxIsClass(prhs[13],"int32") || mxGetM(prhs[13]) < *n || mxGetN(prhs[13]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu6mul","iq must have type int32 and size at least (n,1)");
  }
  int *iq = (int*) mxGetData(prhs[13]);

  /* check prhs[14], variable: lenc */
  if (!mxIsClass(prhs[14],"int32") || mxGetM(prhs[14]) < *n || mxGetN(prhs[14]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu6mul","lenc must have type int32 and size at least (n,1)");
  }
  int *lenc = (int*) mxGetData(prhs[14]);

  /* check prhs[15], variable: lenr */
  if (!mxIsClass(prhs[15],"int32") || mxGetM(prhs[15]) < *m || mxGetN(prhs[15]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu6mul","lenr must have type int32 and size at least (m,1)");
  }
  int *lenr = (int*) mxGetData(prhs[15]);

  /* check prhs[16], variable: locc */
  if (!mxIsClass(prhs[16],"int32") || mxGetM(prhs[16]) < *n || mxGetN(prhs[16]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu6mul","locc must have type int32 and size at least (n,1)");
  }
  int *locc = (int*) mxGetData(prhs[16]);

  /* check prhs[17], variable: locr */
  if (!mxIsClass(prhs[17],"int32") || mxGetM(prhs[17]) < *m || mxGetN(prhs[17]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu6mul","locr must have type int32 and size at least (m,1)");
  }
  int *locr = (int*) mxGetData(prhs[17]);

  LU6MUL_FUNC(mode, m, n, v, w, lena, luparm, parmlu, a, indc, indr, ip, iq, lenc, lenr, locc, locr);
}

void gateway_lu8rpc(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  /* check nrhs */
  if (nrhs != 23) {
    mexErrMsgIdAndTxt("lusol_mex:lu8rpc","incorrect number of inputs for lu8rpc.");
  }

  /* check prhs[1], variable: mode1 */
  if (!mxIsClass(prhs[1],"int32") || mxGetM(prhs[1]) < 1 || mxGetN(prhs[1]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpc","mode1 must have type int32 and size at least (1,1)");
  }
  int *mode1 = (int*) mxGetData(prhs[1]);

  /* check prhs[2], variable: mode2 */
  if (!mxIsClass(prhs[2],"int32") || mxGetM(prhs[2]) < 1 || mxGetN(prhs[2]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpc","mode2 must have type int32 and size at least (1,1)");
  }
  int *mode2 = (int*) mxGetData(prhs[2]);

  /* check prhs[3], variable: m */
  if (!mxIsClass(prhs[3],"int32") || mxGetM(prhs[3]) < 1 || mxGetN(prhs[3]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpc","m must have type int32 and size at least (1,1)");
  }
  int *m = (int*) mxGetData(prhs[3]);

  /* check prhs[4], variable: n */
  if (!mxIsClass(prhs[4],"int32") || mxGetM(prhs[4]) < 1 || mxGetN(prhs[4]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpc","n must have type int32 and size at least (1,1)");
  }
  int *n = (int*) mxGetData(prhs[4]);

  /* check prhs[5], variable: jrep */
  if (!mxIsClass(prhs[5],"int32") || mxGetM(prhs[5]) < 1 || mxGetN(prhs[5]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpc","jrep must have type int32 and size at least (1,1)");
  }
  int *jrep = (int*) mxGetData(prhs[5]);

  /* check prhs[6], variable: v */
  if (!mxIsClass(prhs[6],"double") || mxGetM(prhs[6]) < *m || mxGetN(prhs[6]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpc","v must have type double and size at least (m,1)");
  }
  double *v = (double*) mxGetData(prhs[6]);

  /* check prhs[7], variable: w */
  if (!mxIsClass(prhs[7],"double") || mxGetM(prhs[7]) < *n || mxGetN(prhs[7]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpc","w must have type double and size at least (n,1)");
  }
  double *w = (double*) mxGetData(prhs[7]);

  /* check prhs[8], variable: lena */
  if (!mxIsClass(prhs[8],"int32") || mxGetM(prhs[8]) < 1 || mxGetN(prhs[8]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpc","lena must have type int32 and size at least (1,1)");
  }
  int *lena = (int*) mxGetData(prhs[8]);

  /* check prhs[9], variable: luparm */
  if (!mxIsClass(prhs[9],"int32") || mxGetM(prhs[9]) < 30 || mxGetN(prhs[9]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpc","luparm must have type int32 and size at least (30,1)");
  }
  int *luparm = (int*) mxGetData(prhs[9]);

  /* check prhs[10], variable: parmlu */
  if (!mxIsClass(prhs[10],"double") || mxGetM(prhs[10]) < 30 || mxGetN(prhs[10]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpc","parmlu must have type double and size at least (30,1)");
  }
  double *parmlu = (double*) mxGetData(prhs[10]);

  /* check prhs[11], variable: a */
  if (!mxIsClass(prhs[11],"double") || mxGetM(prhs[11]) < *lena || mxGetN(prhs[11]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpc","a must have type double and size at least (lena,1)");
  }
  double *a = (double*) mxGetData(prhs[11]);

  /* check prhs[12], variable: indc */
  if (!mxIsClass(prhs[12],"int32") || mxGetM(prhs[12]) < *lena || mxGetN(prhs[12]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpc","indc must have type int32 and size at least (lena,1)");
  }
  int *indc = (int*) mxGetData(prhs[12]);

  /* check prhs[13], variable: indr */
  if (!mxIsClass(prhs[13],"int32") || mxGetM(prhs[13]) < *lena || mxGetN(prhs[13]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpc","indr must have type int32 and size at least (lena,1)");
  }
  int *indr = (int*) mxGetData(prhs[13]);

  /* check prhs[14], variable: ip */
  if (!mxIsClass(prhs[14],"int32") || mxGetM(prhs[14]) < *m || mxGetN(prhs[14]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpc","ip must have type int32 and size at least (m,1)");
  }
  int *ip = (int*) mxGetData(prhs[14]);

  /* check prhs[15], variable: iq */
  if (!mxIsClass(prhs[15],"int32") || mxGetM(prhs[15]) < *n || mxGetN(prhs[15]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpc","iq must have type int32 and size at least (n,1)");
  }
  int *iq = (int*) mxGetData(prhs[15]);

  /* check prhs[16], variable: lenc */
  if (!mxIsClass(prhs[16],"int32") || mxGetM(prhs[16]) < *n || mxGetN(prhs[16]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpc","lenc must have type int32 and size at least (n,1)");
  }
  int *lenc = (int*) mxGetData(prhs[16]);

  /* check prhs[17], variable: lenr */
  if (!mxIsClass(prhs[17],"int32") || mxGetM(prhs[17]) < *m || mxGetN(prhs[17]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpc","lenr must have type int32 and size at least (m,1)");
  }
  int *lenr = (int*) mxGetData(prhs[17]);

  /* check prhs[18], variable: locc */
  if (!mxIsClass(prhs[18],"int32") || mxGetM(prhs[18]) < *n || mxGetN(prhs[18]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpc","locc must have type int32 and size at least (n,1)");
  }
  int *locc = (int*) mxGetData(prhs[18]);

  /* check prhs[19], variable: locr */
  if (!mxIsClass(prhs[19],"int32") || mxGetM(prhs[19]) < *m || mxGetN(prhs[19]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpc","locr must have type int32 and size at least (m,1)");
  }
  int *locr = (int*) mxGetData(prhs[19]);

  /* check prhs[20], variable: inform */
  if (!mxIsClass(prhs[20],"int32") || mxGetM(prhs[20]) < 1 || mxGetN(prhs[20]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpc","inform must have type int32 and size at least (1,1)");
  }
  int *inform = (int*) mxGetData(prhs[20]);

  /* check prhs[21], variable: diag */
  if (!mxIsClass(prhs[21],"double") || mxGetM(prhs[21]) < 1 || mxGetN(prhs[21]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpc","diag must have type double and size at least (1,1)");
  }
  double *diag = (double*) mxGetData(prhs[21]);

  /* check prhs[22], variable: vnorm */
  if (!mxIsClass(prhs[22],"double") || mxGetM(prhs[22]) < 1 || mxGetN(prhs[22]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpc","vnorm must have type double and size at least (1,1)");
  }
  double *vnorm = (double*) mxGetData(prhs[22]);

  LU8RPC_FUNC(mode1, mode2, m, n, jrep, v, w, lena, luparm, parmlu, a, indc, indr, ip, iq, lenc, lenr, locc, locr, inform, diag, vnorm);
}

void gateway_lu8adc(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  /* check nrhs */
  if (nrhs != 21) {
    mexErrMsgIdAndTxt("lusol_mex:lu8adc","incorrect number of inputs for lu8adc.");
  }

  /* check prhs[1], variable: mode */
  if (!mxIsClass(prhs[1],"int32") || mxGetM(prhs[1]) < 1 || mxGetN(prhs[1]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8adc","mode must have type int32 and size at least (1,1)");
  }
  int *mode = (int*) mxGetData(prhs[1]);

  /* check prhs[2], variable: m */
  if (!mxIsClass(prhs[2],"int32") || mxGetM(prhs[2]) < 1 || mxGetN(prhs[2]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8adc","m must have type int32 and size at least (1,1)");
  }
  int *m = (int*) mxGetData(prhs[2]);

  /* check prhs[3], variable: n */
  if (!mxIsClass(prhs[3],"int32") || mxGetM(prhs[3]) < 1 || mxGetN(prhs[3]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8adc","n must have type int32 and size at least (1,1)");
  }
  int *n = (int*) mxGetData(prhs[3]);

  /* check prhs[4], variable: v */
  if (!mxIsClass(prhs[4],"double") || mxGetM(prhs[4]) < *m || mxGetN(prhs[4]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8adc","v must have type double and size at least (m,1)");
  }
  double *v = (double*) mxGetData(prhs[4]);

  /* check prhs[5], variable: w */
  if (!mxIsClass(prhs[5],"double") || mxGetM(prhs[5]) < *n || mxGetN(prhs[5]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8adc","w must have type double and size at least (n,1)");
  }
  double *w = (double*) mxGetData(prhs[5]);

  /* check prhs[6], variable: lena */
  if (!mxIsClass(prhs[6],"int32") || mxGetM(prhs[6]) < 1 || mxGetN(prhs[6]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8adc","lena must have type int32 and size at least (1,1)");
  }
  int *lena = (int*) mxGetData(prhs[6]);

  /* check prhs[7], variable: luparm */
  if (!mxIsClass(prhs[7],"int32") || mxGetM(prhs[7]) < 30 || mxGetN(prhs[7]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8adc","luparm must have type int32 and size at least (30,1)");
  }
  int *luparm = (int*) mxGetData(prhs[7]);

  /* check prhs[8], variable: parmlu */
  if (!mxIsClass(prhs[8],"double") || mxGetM(prhs[8]) < 30 || mxGetN(prhs[8]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8adc","parmlu must have type double and size at least (30,1)");
  }
  double *parmlu = (double*) mxGetData(prhs[8]);

  /* check prhs[9], variable: a */
  if (!mxIsClass(prhs[9],"double") || mxGetM(prhs[9]) < *lena || mxGetN(prhs[9]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8adc","a must have type double and size at least (lena,1)");
  }
  double *a = (double*) mxGetData(prhs[9]);

  /* check prhs[10], variable: indc */
  if (!mxIsClass(prhs[10],"int32") || mxGetM(prhs[10]) < *lena || mxGetN(prhs[10]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8adc","indc must have type int32 and size at least (lena,1)");
  }
  int *indc = (int*) mxGetData(prhs[10]);

  /* check prhs[11], variable: indr */
  if (!mxIsClass(prhs[11],"int32") || mxGetM(prhs[11]) < *lena || mxGetN(prhs[11]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8adc","indr must have type int32 and size at least (lena,1)");
  }
  int *indr = (int*) mxGetData(prhs[11]);

  /* check prhs[12], variable: ip */
  if (!mxIsClass(prhs[12],"int32") || mxGetM(prhs[12]) < *m || mxGetN(prhs[12]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8adc","ip must have type int32 and size at least (m,1)");
  }
  int *ip = (int*) mxGetData(prhs[12]);

  /* check prhs[13], variable: iq */
  if (!mxIsClass(prhs[13],"int32") || mxGetM(prhs[13]) < *n || mxGetN(prhs[13]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8adc","iq must have type int32 and size at least (n,1)");
  }
  int *iq = (int*) mxGetData(prhs[13]);

  /* check prhs[14], variable: lenc */
  if (!mxIsClass(prhs[14],"int32") || mxGetM(prhs[14]) < *n || mxGetN(prhs[14]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8adc","lenc must have type int32 and size at least (n,1)");
  }
  int *lenc = (int*) mxGetData(prhs[14]);

  /* check prhs[15], variable: lenr */
  if (!mxIsClass(prhs[15],"int32") || mxGetM(prhs[15]) < *m || mxGetN(prhs[15]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8adc","lenr must have type int32 and size at least (m,1)");
  }
  int *lenr = (int*) mxGetData(prhs[15]);

  /* check prhs[16], variable: locc */
  if (!mxIsClass(prhs[16],"int32") || mxGetM(prhs[16]) < *n || mxGetN(prhs[16]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8adc","locc must have type int32 and size at least (n,1)");
  }
  int *locc = (int*) mxGetData(prhs[16]);

  /* check prhs[17], variable: locr */
  if (!mxIsClass(prhs[17],"int32") || mxGetM(prhs[17]) < *m || mxGetN(prhs[17]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8adc","locr must have type int32 and size at least (m,1)");
  }
  int *locr = (int*) mxGetData(prhs[17]);

  /* check prhs[18], variable: inform */
  if (!mxIsClass(prhs[18],"int32") || mxGetM(prhs[18]) < 1 || mxGetN(prhs[18]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8adc","inform must have type int32 and size at least (1,1)");
  }
  int *inform = (int*) mxGetData(prhs[18]);

  /* check prhs[19], variable: diag */
  if (!mxIsClass(prhs[19],"double") || mxGetM(prhs[19]) < 1 || mxGetN(prhs[19]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8adc","diag must have type double and size at least (1,1)");
  }
  double *diag = (double*) mxGetData(prhs[19]);

  /* check prhs[20], variable: vnorm */
  if (!mxIsClass(prhs[20],"double") || mxGetM(prhs[20]) < 1 || mxGetN(prhs[20]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8adc","vnorm must have type double and size at least (1,1)");
  }
  double *vnorm = (double*) mxGetData(prhs[20]);

  LU8ADC_FUNC(mode, m, n, v, w, lena, luparm, parmlu, a, indc, indr, ip, iq, lenc, lenr, locc, locr, inform, diag, vnorm);
}

void gateway_lu8adr(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  /* check nrhs */
  if (nrhs != 18) {
    mexErrMsgIdAndTxt("lusol_mex:lu8adr","incorrect number of inputs for lu8adr.");
  }

  /* check prhs[1], variable: m */
  if (!mxIsClass(prhs[1],"int32") || mxGetM(prhs[1]) < 1 || mxGetN(prhs[1]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8adr","m must have type int32 and size at least (1,1)");
  }
  int *m = (int*) mxGetData(prhs[1]);

  /* check prhs[2], variable: n */
  if (!mxIsClass(prhs[2],"int32") || mxGetM(prhs[2]) < 1 || mxGetN(prhs[2]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8adr","n must have type int32 and size at least (1,1)");
  }
  int *n = (int*) mxGetData(prhs[2]);

  /* check prhs[3], variable: w */
  if (!mxIsClass(prhs[3],"double") || mxGetM(prhs[3]) < *n || mxGetN(prhs[3]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8adr","w must have type double and size at least (n,1)");
  }
  double *w = (double*) mxGetData(prhs[3]);

  /* check prhs[4], variable: lena */
  if (!mxIsClass(prhs[4],"int32") || mxGetM(prhs[4]) < 1 || mxGetN(prhs[4]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8adr","lena must have type int32 and size at least (1,1)");
  }
  int *lena = (int*) mxGetData(prhs[4]);

  /* check prhs[5], variable: luparm */
  if (!mxIsClass(prhs[5],"int32") || mxGetM(prhs[5]) < 30 || mxGetN(prhs[5]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8adr","luparm must have type int32 and size at least (30,1)");
  }
  int *luparm = (int*) mxGetData(prhs[5]);

  /* check prhs[6], variable: parmlu */
  if (!mxIsClass(prhs[6],"double") || mxGetM(prhs[6]) < 30 || mxGetN(prhs[6]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8adr","parmlu must have type double and size at least (30,1)");
  }
  double *parmlu = (double*) mxGetData(prhs[6]);

  /* check prhs[7], variable: a */
  if (!mxIsClass(prhs[7],"double") || mxGetM(prhs[7]) < *lena || mxGetN(prhs[7]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8adr","a must have type double and size at least (lena,1)");
  }
  double *a = (double*) mxGetData(prhs[7]);

  /* check prhs[8], variable: indc */
  if (!mxIsClass(prhs[8],"int32") || mxGetM(prhs[8]) < *lena || mxGetN(prhs[8]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8adr","indc must have type int32 and size at least (lena,1)");
  }
  int *indc = (int*) mxGetData(prhs[8]);

  /* check prhs[9], variable: indr */
  if (!mxIsClass(prhs[9],"int32") || mxGetM(prhs[9]) < *lena || mxGetN(prhs[9]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8adr","indr must have type int32 and size at least (lena,1)");
  }
  int *indr = (int*) mxGetData(prhs[9]);

  /* check prhs[10], variable: ip */
  if (!mxIsClass(prhs[10],"int32") || mxGetM(prhs[10]) < *m || mxGetN(prhs[10]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8adr","ip must have type int32 and size at least (m,1)");
  }
  int *ip = (int*) mxGetData(prhs[10]);

  /* check prhs[11], variable: iq */
  if (!mxIsClass(prhs[11],"int32") || mxGetM(prhs[11]) < *n || mxGetN(prhs[11]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8adr","iq must have type int32 and size at least (n,1)");
  }
  int *iq = (int*) mxGetData(prhs[11]);

  /* check prhs[12], variable: lenc */
  if (!mxIsClass(prhs[12],"int32") || mxGetM(prhs[12]) < *n || mxGetN(prhs[12]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8adr","lenc must have type int32 and size at least (n,1)");
  }
  int *lenc = (int*) mxGetData(prhs[12]);

  /* check prhs[13], variable: lenr */
  if (!mxIsClass(prhs[13],"int32") || mxGetM(prhs[13]) < *m || mxGetN(prhs[13]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8adr","lenr must have type int32 and size at least (m,1)");
  }
  int *lenr = (int*) mxGetData(prhs[13]);

  /* check prhs[14], variable: locc */
  if (!mxIsClass(prhs[14],"int32") || mxGetM(prhs[14]) < *n || mxGetN(prhs[14]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8adr","locc must have type int32 and size at least (n,1)");
  }
  int *locc = (int*) mxGetData(prhs[14]);

  /* check prhs[15], variable: locr */
  if (!mxIsClass(prhs[15],"int32") || mxGetM(prhs[15]) < *m || mxGetN(prhs[15]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8adr","locr must have type int32 and size at least (m,1)");
  }
  int *locr = (int*) mxGetData(prhs[15]);

  /* check prhs[16], variable: inform */
  if (!mxIsClass(prhs[16],"int32") || mxGetM(prhs[16]) < 1 || mxGetN(prhs[16]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8adr","inform must have type int32 and size at least (1,1)");
  }
  int *inform = (int*) mxGetData(prhs[16]);

  /* check prhs[17], variable: diag */
  if (!mxIsClass(prhs[17],"double") || mxGetM(prhs[17]) < 1 || mxGetN(prhs[17]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8adr","diag must have type double and size at least (1,1)");
  }
  double *diag = (double*) mxGetData(prhs[17]);

  LU8ADR_FUNC(m, n, w, lena, luparm, parmlu, a, indc, indr, ip, iq, lenc, lenr, locc, locr, inform, diag);
}

void gateway_lu8dlc(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  /* check nrhs */
  if (nrhs != 17) {
    mexErrMsgIdAndTxt("lusol_mex:lu8dlc","incorrect number of inputs for lu8dlc.");
  }

  /* check prhs[1], variable: m */
  if (!mxIsClass(prhs[1],"int32") || mxGetM(prhs[1]) < 1 || mxGetN(prhs[1]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8dlc","m must have type int32 and size at least (1,1)");
  }
  int *m = (int*) mxGetData(prhs[1]);

  /* check prhs[2], variable: n */
  if (!mxIsClass(prhs[2],"int32") || mxGetM(prhs[2]) < 1 || mxGetN(prhs[2]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8dlc","n must have type int32 and size at least (1,1)");
  }
  int *n = (int*) mxGetData(prhs[2]);

  /* check prhs[3], variable: jdel */
  if (!mxIsClass(prhs[3],"int32") || mxGetM(prhs[3]) < 1 || mxGetN(prhs[3]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8dlc","jdel must have type int32 and size at least (1,1)");
  }
  int *jdel = (int*) mxGetData(prhs[3]);

  /* check prhs[4], variable: lena */
  if (!mxIsClass(prhs[4],"int32") || mxGetM(prhs[4]) < 1 || mxGetN(prhs[4]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8dlc","lena must have type int32 and size at least (1,1)");
  }
  int *lena = (int*) mxGetData(prhs[4]);

  /* check prhs[5], variable: luparm */
  if (!mxIsClass(prhs[5],"int32") || mxGetM(prhs[5]) < 30 || mxGetN(prhs[5]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8dlc","luparm must have type int32 and size at least (30,1)");
  }
  int *luparm = (int*) mxGetData(prhs[5]);

  /* check prhs[6], variable: parmlu */
  if (!mxIsClass(prhs[6],"double") || mxGetM(prhs[6]) < 30 || mxGetN(prhs[6]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8dlc","parmlu must have type double and size at least (30,1)");
  }
  double *parmlu = (double*) mxGetData(prhs[6]);

  /* check prhs[7], variable: a */
  if (!mxIsClass(prhs[7],"double") || mxGetM(prhs[7]) < *lena || mxGetN(prhs[7]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8dlc","a must have type double and size at least (lena,1)");
  }
  double *a = (double*) mxGetData(prhs[7]);

  /* check prhs[8], variable: indc */
  if (!mxIsClass(prhs[8],"int32") || mxGetM(prhs[8]) < *lena || mxGetN(prhs[8]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8dlc","indc must have type int32 and size at least (lena,1)");
  }
  int *indc = (int*) mxGetData(prhs[8]);

  /* check prhs[9], variable: indr */
  if (!mxIsClass(prhs[9],"int32") || mxGetM(prhs[9]) < *lena || mxGetN(prhs[9]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8dlc","indr must have type int32 and size at least (lena,1)");
  }
  int *indr = (int*) mxGetData(prhs[9]);

  /* check prhs[10], variable: ip */
  if (!mxIsClass(prhs[10],"int32") || mxGetM(prhs[10]) < *m || mxGetN(prhs[10]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8dlc","ip must have type int32 and size at least (m,1)");
  }
  int *ip = (int*) mxGetData(prhs[10]);

  /* check prhs[11], variable: iq */
  if (!mxIsClass(prhs[11],"int32") || mxGetM(prhs[11]) < *n || mxGetN(prhs[11]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8dlc","iq must have type int32 and size at least (n,1)");
  }
  int *iq = (int*) mxGetData(prhs[11]);

  /* check prhs[12], variable: lenc */
  if (!mxIsClass(prhs[12],"int32") || mxGetM(prhs[12]) < *n || mxGetN(prhs[12]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8dlc","lenc must have type int32 and size at least (n,1)");
  }
  int *lenc = (int*) mxGetData(prhs[12]);

  /* check prhs[13], variable: lenr */
  if (!mxIsClass(prhs[13],"int32") || mxGetM(prhs[13]) < *m || mxGetN(prhs[13]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8dlc","lenr must have type int32 and size at least (m,1)");
  }
  int *lenr = (int*) mxGetData(prhs[13]);

  /* check prhs[14], variable: locc */
  if (!mxIsClass(prhs[14],"int32") || mxGetM(prhs[14]) < *n || mxGetN(prhs[14]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8dlc","locc must have type int32 and size at least (n,1)");
  }
  int *locc = (int*) mxGetData(prhs[14]);

  /* check prhs[15], variable: locr */
  if (!mxIsClass(prhs[15],"int32") || mxGetM(prhs[15]) < *m || mxGetN(prhs[15]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8dlc","locr must have type int32 and size at least (m,1)");
  }
  int *locr = (int*) mxGetData(prhs[15]);

  /* check prhs[16], variable: inform */
  if (!mxIsClass(prhs[16],"int32") || mxGetM(prhs[16]) < 1 || mxGetN(prhs[16]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8dlc","inform must have type int32 and size at least (1,1)");
  }
  int *inform = (int*) mxGetData(prhs[16]);

  LU8DLC_FUNC(m, n, jdel, lena, luparm, parmlu, a, indc, indr, ip, iq, lenc, lenr, locc, locr, inform);
}

void gateway_lu8dlr(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  /* check nrhs */
  if (nrhs != 20) {
    mexErrMsgIdAndTxt("lusol_mex:lu8dlr","incorrect number of inputs for lu8dlr.");
  }

  /* check prhs[1], variable: mode */
  if (!mxIsClass(prhs[1],"int32") || mxGetM(prhs[1]) < 1 || mxGetN(prhs[1]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8dlr","mode must have type int32 and size at least (1,1)");
  }
  int *mode = (int*) mxGetData(prhs[1]);

  /* check prhs[2], variable: m */
  if (!mxIsClass(prhs[2],"int32") || mxGetM(prhs[2]) < 1 || mxGetN(prhs[2]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8dlr","m must have type int32 and size at least (1,1)");
  }
  int *m = (int*) mxGetData(prhs[2]);

  /* check prhs[3], variable: n */
  if (!mxIsClass(prhs[3],"int32") || mxGetM(prhs[3]) < 1 || mxGetN(prhs[3]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8dlr","n must have type int32 and size at least (1,1)");
  }
  int *n = (int*) mxGetData(prhs[3]);

  /* check prhs[4], variable: idel */
  if (!mxIsClass(prhs[4],"int32") || mxGetM(prhs[4]) < 1 || mxGetN(prhs[4]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8dlr","idel must have type int32 and size at least (1,1)");
  }
  int *idel = (int*) mxGetData(prhs[4]);

  /* check prhs[5], variable: v */
  if (!mxIsClass(prhs[5],"double") || mxGetM(prhs[5]) < *m || mxGetN(prhs[5]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8dlr","v must have type double and size at least (m,1)");
  }
  double *v = (double*) mxGetData(prhs[5]);

  /* check prhs[6], variable: w */
  if (!mxIsClass(prhs[6],"double") || mxGetM(prhs[6]) < *n || mxGetN(prhs[6]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8dlr","w must have type double and size at least (n,1)");
  }
  double *w = (double*) mxGetData(prhs[6]);

  /* check prhs[7], variable: lena */
  if (!mxIsClass(prhs[7],"int32") || mxGetM(prhs[7]) < 1 || mxGetN(prhs[7]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8dlr","lena must have type int32 and size at least (1,1)");
  }
  int *lena = (int*) mxGetData(prhs[7]);

  /* check prhs[8], variable: luparm */
  if (!mxIsClass(prhs[8],"int32") || mxGetM(prhs[8]) < 30 || mxGetN(prhs[8]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8dlr","luparm must have type int32 and size at least (30,1)");
  }
  int *luparm = (int*) mxGetData(prhs[8]);

  /* check prhs[9], variable: parmlu */
  if (!mxIsClass(prhs[9],"double") || mxGetM(prhs[9]) < 30 || mxGetN(prhs[9]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8dlr","parmlu must have type double and size at least (30,1)");
  }
  double *parmlu = (double*) mxGetData(prhs[9]);

  /* check prhs[10], variable: a */
  if (!mxIsClass(prhs[10],"double") || mxGetM(prhs[10]) < *lena || mxGetN(prhs[10]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8dlr","a must have type double and size at least (lena,1)");
  }
  double *a = (double*) mxGetData(prhs[10]);

  /* check prhs[11], variable: indc */
  if (!mxIsClass(prhs[11],"int32") || mxGetM(prhs[11]) < *lena || mxGetN(prhs[11]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8dlr","indc must have type int32 and size at least (lena,1)");
  }
  int *indc = (int*) mxGetData(prhs[11]);

  /* check prhs[12], variable: indr */
  if (!mxIsClass(prhs[12],"int32") || mxGetM(prhs[12]) < *lena || mxGetN(prhs[12]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8dlr","indr must have type int32 and size at least (lena,1)");
  }
  int *indr = (int*) mxGetData(prhs[12]);

  /* check prhs[13], variable: ip */
  if (!mxIsClass(prhs[13],"int32") || mxGetM(prhs[13]) < *m || mxGetN(prhs[13]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8dlr","ip must have type int32 and size at least (m,1)");
  }
  int *ip = (int*) mxGetData(prhs[13]);

  /* check prhs[14], variable: iq */
  if (!mxIsClass(prhs[14],"int32") || mxGetM(prhs[14]) < *n || mxGetN(prhs[14]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8dlr","iq must have type int32 and size at least (n,1)");
  }
  int *iq = (int*) mxGetData(prhs[14]);

  /* check prhs[15], variable: lenc */
  if (!mxIsClass(prhs[15],"int32") || mxGetM(prhs[15]) < *n || mxGetN(prhs[15]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8dlr","lenc must have type int32 and size at least (n,1)");
  }
  int *lenc = (int*) mxGetData(prhs[15]);

  /* check prhs[16], variable: lenr */
  if (!mxIsClass(prhs[16],"int32") || mxGetM(prhs[16]) < *m || mxGetN(prhs[16]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8dlr","lenr must have type int32 and size at least (m,1)");
  }
  int *lenr = (int*) mxGetData(prhs[16]);

  /* check prhs[17], variable: locc */
  if (!mxIsClass(prhs[17],"int32") || mxGetM(prhs[17]) < *n || mxGetN(prhs[17]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8dlr","locc must have type int32 and size at least (n,1)");
  }
  int *locc = (int*) mxGetData(prhs[17]);

  /* check prhs[18], variable: locr */
  if (!mxIsClass(prhs[18],"int32") || mxGetM(prhs[18]) < *m || mxGetN(prhs[18]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8dlr","locr must have type int32 and size at least (m,1)");
  }
  int *locr = (int*) mxGetData(prhs[18]);

  /* check prhs[19], variable: inform */
  if (!mxIsClass(prhs[19],"int32") || mxGetM(prhs[19]) < 1 || mxGetN(prhs[19]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8dlr","inform must have type int32 and size at least (1,1)");
  }
  int *inform = (int*) mxGetData(prhs[19]);

  LU8DLR_FUNC(mode, m, n, idel, v, w, lena, luparm, parmlu, a, indc, indr, ip, iq, lenc, lenr, locc, locr, inform);
}

void gateway_lu8mod(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  /* check nrhs */
  if (nrhs != 20) {
    mexErrMsgIdAndTxt("lusol_mex:lu8mod","incorrect number of inputs for lu8mod.");
  }

  /* check prhs[1], variable: mode */
  if (!mxIsClass(prhs[1],"int32") || mxGetM(prhs[1]) < 1 || mxGetN(prhs[1]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8mod","mode must have type int32 and size at least (1,1)");
  }
  int *mode = (int*) mxGetData(prhs[1]);

  /* check prhs[2], variable: m */
  if (!mxIsClass(prhs[2],"int32") || mxGetM(prhs[2]) < 1 || mxGetN(prhs[2]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8mod","m must have type int32 and size at least (1,1)");
  }
  int *m = (int*) mxGetData(prhs[2]);

  /* check prhs[3], variable: n */
  if (!mxIsClass(prhs[3],"int32") || mxGetM(prhs[3]) < 1 || mxGetN(prhs[3]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8mod","n must have type int32 and size at least (1,1)");
  }
  int *n = (int*) mxGetData(prhs[3]);

  /* check prhs[4], variable: beta */
  if (!mxIsClass(prhs[4],"double") || mxGetM(prhs[4]) < 1 || mxGetN(prhs[4]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8mod","beta must have type double and size at least (1,1)");
  }
  double *beta = (double*) mxGetData(prhs[4]);

  /* check prhs[5], variable: v */
  if (!mxIsClass(prhs[5],"double") || mxGetM(prhs[5]) < *m || mxGetN(prhs[5]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8mod","v must have type double and size at least (m,1)");
  }
  double *v = (double*) mxGetData(prhs[5]);

  /* check prhs[6], variable: w */
  if (!mxIsClass(prhs[6],"double") || mxGetM(prhs[6]) < *n || mxGetN(prhs[6]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8mod","w must have type double and size at least (n,1)");
  }
  double *w = (double*) mxGetData(prhs[6]);

  /* check prhs[7], variable: lena */
  if (!mxIsClass(prhs[7],"int32") || mxGetM(prhs[7]) < 1 || mxGetN(prhs[7]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8mod","lena must have type int32 and size at least (1,1)");
  }
  int *lena = (int*) mxGetData(prhs[7]);

  /* check prhs[8], variable: luparm */
  if (!mxIsClass(prhs[8],"int32") || mxGetM(prhs[8]) < 30 || mxGetN(prhs[8]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8mod","luparm must have type int32 and size at least (30,1)");
  }
  int *luparm = (int*) mxGetData(prhs[8]);

  /* check prhs[9], variable: parmlu */
  if (!mxIsClass(prhs[9],"double") || mxGetM(prhs[9]) < 30 || mxGetN(prhs[9]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8mod","parmlu must have type double and size at least (30,1)");
  }
  double *parmlu = (double*) mxGetData(prhs[9]);

  /* check prhs[10], variable: a */
  if (!mxIsClass(prhs[10],"double") || mxGetM(prhs[10]) < *lena || mxGetN(prhs[10]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8mod","a must have type double and size at least (lena,1)");
  }
  double *a = (double*) mxGetData(prhs[10]);

  /* check prhs[11], variable: indc */
  if (!mxIsClass(prhs[11],"int32") || mxGetM(prhs[11]) < *lena || mxGetN(prhs[11]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8mod","indc must have type int32 and size at least (lena,1)");
  }
  int *indc = (int*) mxGetData(prhs[11]);

  /* check prhs[12], variable: indr */
  if (!mxIsClass(prhs[12],"int32") || mxGetM(prhs[12]) < *lena || mxGetN(prhs[12]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8mod","indr must have type int32 and size at least (lena,1)");
  }
  int *indr = (int*) mxGetData(prhs[12]);

  /* check prhs[13], variable: ip */
  if (!mxIsClass(prhs[13],"int32") || mxGetM(prhs[13]) < *m || mxGetN(prhs[13]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8mod","ip must have type int32 and size at least (m,1)");
  }
  int *ip = (int*) mxGetData(prhs[13]);

  /* check prhs[14], variable: iq */
  if (!mxIsClass(prhs[14],"int32") || mxGetM(prhs[14]) < *n || mxGetN(prhs[14]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8mod","iq must have type int32 and size at least (n,1)");
  }
  int *iq = (int*) mxGetData(prhs[14]);

  /* check prhs[15], variable: lenc */
  if (!mxIsClass(prhs[15],"int32") || mxGetM(prhs[15]) < *n || mxGetN(prhs[15]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8mod","lenc must have type int32 and size at least (n,1)");
  }
  int *lenc = (int*) mxGetData(prhs[15]);

  /* check prhs[16], variable: lenr */
  if (!mxIsClass(prhs[16],"int32") || mxGetM(prhs[16]) < *m || mxGetN(prhs[16]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8mod","lenr must have type int32 and size at least (m,1)");
  }
  int *lenr = (int*) mxGetData(prhs[16]);

  /* check prhs[17], variable: locc */
  if (!mxIsClass(prhs[17],"int32") || mxGetM(prhs[17]) < *n || mxGetN(prhs[17]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8mod","locc must have type int32 and size at least (n,1)");
  }
  int *locc = (int*) mxGetData(prhs[17]);

  /* check prhs[18], variable: locr */
  if (!mxIsClass(prhs[18],"int32") || mxGetM(prhs[18]) < *m || mxGetN(prhs[18]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8mod","locr must have type int32 and size at least (m,1)");
  }
  int *locr = (int*) mxGetData(prhs[18]);

  /* check prhs[19], variable: inform */
  if (!mxIsClass(prhs[19],"int32") || mxGetM(prhs[19]) < 1 || mxGetN(prhs[19]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8mod","inform must have type int32 and size at least (1,1)");
  }
  int *inform = (int*) mxGetData(prhs[19]);

  LU8MOD_FUNC(mode, m, n, beta, v, w, lena, luparm, parmlu, a, indc, indr, ip, iq, lenc, lenr, locc, locr, inform);
}

void gateway_lu8rpr(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  /* check nrhs */
  if (nrhs != 22) {
    mexErrMsgIdAndTxt("lusol_mex:lu8rpr","incorrect number of inputs for lu8rpr.");
  }

  /* check prhs[1], variable: mode1 */
  if (!mxIsClass(prhs[1],"int32") || mxGetM(prhs[1]) < 1 || mxGetN(prhs[1]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpr","mode1 must have type int32 and size at least (1,1)");
  }
  int *mode1 = (int*) mxGetData(prhs[1]);

  /* check prhs[2], variable: mode2 */
  if (!mxIsClass(prhs[2],"int32") || mxGetM(prhs[2]) < 1 || mxGetN(prhs[2]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpr","mode2 must have type int32 and size at least (1,1)");
  }
  int *mode2 = (int*) mxGetData(prhs[2]);

  /* check prhs[3], variable: m */
  if (!mxIsClass(prhs[3],"int32") || mxGetM(prhs[3]) < 1 || mxGetN(prhs[3]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpr","m must have type int32 and size at least (1,1)");
  }
  int *m = (int*) mxGetData(prhs[3]);

  /* check prhs[4], variable: n */
  if (!mxIsClass(prhs[4],"int32") || mxGetM(prhs[4]) < 1 || mxGetN(prhs[4]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpr","n must have type int32 and size at least (1,1)");
  }
  int *n = (int*) mxGetData(prhs[4]);

  /* check prhs[5], variable: irep */
  if (!mxIsClass(prhs[5],"int32") || mxGetM(prhs[5]) < 1 || mxGetN(prhs[5]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpr","irep must have type int32 and size at least (1,1)");
  }
  int *irep = (int*) mxGetData(prhs[5]);

  /* check prhs[6], variable: v */
  if (!mxIsClass(prhs[6],"double") || mxGetM(prhs[6]) < *m || mxGetN(prhs[6]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpr","v must have type double and size at least (m,1)");
  }
  double *v = (double*) mxGetData(prhs[6]);

  /* check prhs[7], variable: w */
  if (!mxIsClass(prhs[7],"double") || mxGetM(prhs[7]) < *n || mxGetN(prhs[7]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpr","w must have type double and size at least (n,1)");
  }
  double *w = (double*) mxGetData(prhs[7]);

  /* check prhs[8], variable: wnew */
  if (!mxIsClass(prhs[8],"double") || mxGetM(prhs[8]) < *n || mxGetN(prhs[8]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpr","wnew must have type double and size at least (n,1)");
  }
  double *wnew = (double*) mxGetData(prhs[8]);

  /* check prhs[9], variable: lena */
  if (!mxIsClass(prhs[9],"int32") || mxGetM(prhs[9]) < 1 || mxGetN(prhs[9]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpr","lena must have type int32 and size at least (1,1)");
  }
  int *lena = (int*) mxGetData(prhs[9]);

  /* check prhs[10], variable: luparm */
  if (!mxIsClass(prhs[10],"int32") || mxGetM(prhs[10]) < 30 || mxGetN(prhs[10]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpr","luparm must have type int32 and size at least (30,1)");
  }
  int *luparm = (int*) mxGetData(prhs[10]);

  /* check prhs[11], variable: parmlu */
  if (!mxIsClass(prhs[11],"double") || mxGetM(prhs[11]) < 30 || mxGetN(prhs[11]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpr","parmlu must have type double and size at least (30,1)");
  }
  double *parmlu = (double*) mxGetData(prhs[11]);

  /* check prhs[12], variable: a */
  if (!mxIsClass(prhs[12],"double") || mxGetM(prhs[12]) < *lena || mxGetN(prhs[12]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpr","a must have type double and size at least (lena,1)");
  }
  double *a = (double*) mxGetData(prhs[12]);

  /* check prhs[13], variable: indc */
  if (!mxIsClass(prhs[13],"int32") || mxGetM(prhs[13]) < *lena || mxGetN(prhs[13]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpr","indc must have type int32 and size at least (lena,1)");
  }
  int *indc = (int*) mxGetData(prhs[13]);

  /* check prhs[14], variable: indr */
  if (!mxIsClass(prhs[14],"int32") || mxGetM(prhs[14]) < *lena || mxGetN(prhs[14]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpr","indr must have type int32 and size at least (lena,1)");
  }
  int *indr = (int*) mxGetData(prhs[14]);

  /* check prhs[15], variable: ip */
  if (!mxIsClass(prhs[15],"int32") || mxGetM(prhs[15]) < *m || mxGetN(prhs[15]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpr","ip must have type int32 and size at least (m,1)");
  }
  int *ip = (int*) mxGetData(prhs[15]);

  /* check prhs[16], variable: iq */
  if (!mxIsClass(prhs[16],"int32") || mxGetM(prhs[16]) < *n || mxGetN(prhs[16]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpr","iq must have type int32 and size at least (n,1)");
  }
  int *iq = (int*) mxGetData(prhs[16]);

  /* check prhs[17], variable: lenc */
  if (!mxIsClass(prhs[17],"int32") || mxGetM(prhs[17]) < *n || mxGetN(prhs[17]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpr","lenc must have type int32 and size at least (n,1)");
  }
  int *lenc = (int*) mxGetData(prhs[17]);

  /* check prhs[18], variable: lenr */
  if (!mxIsClass(prhs[18],"int32") || mxGetM(prhs[18]) < *m || mxGetN(prhs[18]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpr","lenr must have type int32 and size at least (m,1)");
  }
  int *lenr = (int*) mxGetData(prhs[18]);

  /* check prhs[19], variable: locc */
  if (!mxIsClass(prhs[19],"int32") || mxGetM(prhs[19]) < *n || mxGetN(prhs[19]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpr","locc must have type int32 and size at least (n,1)");
  }
  int *locc = (int*) mxGetData(prhs[19]);

  /* check prhs[20], variable: locr */
  if (!mxIsClass(prhs[20],"int32") || mxGetM(prhs[20]) < *m || mxGetN(prhs[20]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpr","locr must have type int32 and size at least (m,1)");
  }
  int *locr = (int*) mxGetData(prhs[20]);

  /* check prhs[21], variable: inform */
  if (!mxIsClass(prhs[21],"int32") || mxGetM(prhs[21]) < 1 || mxGetN(prhs[21]) != 1) { 
    mexErrMsgIdAndTxt("lusol_mex:lu8rpr","inform must have type int32 and size at least (1,1)");
  }
  int *inform = (int*) mxGetData(prhs[21]);

  LU8RPR_FUNC(mode1, mode2, m, n, irep, v, w, wnew, lena, luparm, parmlu, a, indc, indr, ip, iq, lenc, lenr, locc, locr, inform);
}

