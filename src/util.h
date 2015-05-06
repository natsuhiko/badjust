#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <gsl_sf_gamma.h>
#include <gsl_sf_psi.h>


//#include <R.h>

double pchisq(double q, double k);
double Qchisq(double p0, double k);


double inverse(double* a, double* inv);
double getStepSize(double* H, double* g, int fixPi);


double pf1(double fval, int df1, int df2);
void pf(double* fvals, int N, int df1, int df2);

double sign1(double g);
void dcopy(long n, double* x, double* y);



long imax(double* x, long n);
//double max(double* x, long n);
long lomax(long* x, long n);
double doub(double x);

double esum(double *x, long n, double y);

double lgamma(double x);
double digamma(double x);
double trigamma(double x);
double psi_2(double x);

long ncol(FILE* fp, long L);
void scale(double* x, double* v, long L);

void progress(long j, long L);
void isna(double* v, long n, char tex[]);
void print(double* v, long n);
void printSep(double* v, long n, char sep);
void lprint(double* v, long n);
void print2(double* v, long n, long offs);
void printLong(long* v, long n);
void printInt(int* v, long n);
void printM(double* m, long n, long p);
void printL(double* m, long n, long p);
void clear(double* v, long N);
void clear1(double* v, long N);
void clearInt(int* v, int N);
void clearLong(long* v, long N);
void apply1(double* x, long n, long p, double* y);
void apply2(double* x, long n, long p, double* y);
void BackSolveOld(double* R, double* y, double* x, long p);
void QRold(double* X, long n, long p, double* R);
double ifelse(long a, long b, double* Beta);
double nrm2(double* x, long n, long offset, long inc);
double asum0(double* x, long n, long offset, long inc);
double asum(double* x, long n);
double sum(double *x, long n);
void wshift(double *x, double *k, long n);
double wsum(double *x, double *k, long n);
double iwsum(double *x, double *k, long n);
double iwlsum(long *x, double *k, long n);
long sumLong(long *x, long n);
double ddot(double* x, double* y, long n, long offset, long inc);

double ltrace(double* R, long n, long ldR);
double lsum(double* x, long n);

void WriteBin(double* x, long n, char* fname);

double logit(double x);
double expit(double x);


long minL(long i, long j);
//double min(double a, double b);

double mean(double* x, long n);
double lmean(double* x, long n);
double lsd(double* x, long n);
double sd(double* x, long n);
double wmean(double* x, double* w, long n);
double wvar(double* x, double* w, long n);

void sum2one(double* x, int n);
