double lgamma(double x);
double digamma(double x);
double trigamma(double x);
double logm1(double x);
double lbeta(double p, double q);
double logK(double x, double p, double q);
double hg2F1(double a1, double a2, double b, double x);
double logI(double x, double p, double q);
double logIn(double x, double p, double q, int n, double* res, double* work);
//void logIforR(double* v, double* work);
double getLkhd(double* x, int* ntests, int n, double* abc);
double nlm(double* x, int* ntests, int n, double* abc, double* grad, double* hess, double* work, int* fixParam, int verbose2);
void nlmForR(double* x, int* ntests, int* pn, double* abc, double* lkhd, double* grad, double* hess, double* work, int* fixParam, int* pverb);


double getI(double x, double a, double b);
double getI2(double x, double a, double b);
