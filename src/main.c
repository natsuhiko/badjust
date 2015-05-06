#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
//#include <gsl_cblas.h>
//#include <f2c.h>
//#include <clapack.h>
//#include <blaswrap.h>
//#include <cblas.h>
#include <gsl_sf_gamma.h>
#include <gsl_sf_hyperg.h>
#include <gsl_sf_psi.h>
#include <zlib.h>
#include "util.h"

#include "badjust.h"
int verbose=0;
int main(int argc, char** argv){
	verbose=1;
	double* v;
	double* res;
	double* w;
	v=(double*)calloc(3, sizeof(double));
	res=(double*)calloc(100, sizeof(double));
	w=(double*)calloc(24*101, sizeof(double));
	v[0]=atof(argv[1]);
	v[1]=atof(argv[2]);
	v[2]=atof(argv[3]);
	logIn(v[0], v[1], v[2], 15, res, w);
	printf("x=%lf, a=%lf, b=%lf\n\nlog I=%lf %lf %lf %lf %lf %lf\n", v[0], v[1], v[2], res[0], res[1], res[2] ,res[3], res[4], res[5]);
}
double lgamma(double x);
double digamma(double x);
double trigamma(double x);
double lbeta(double p, double q){
    return gsl_sf_lnbeta(p, q);
}
double logK(double x, double p, double q){
    //printf("beta=%lf %lf %lf %lf", x, log(x), log1p(-x), lbeta(p,q));
    return p*log(x) + (q-1)*log1p(-x) - log(p) - lbeta(p,q);
}
double hg2F1(double a1, double a2, double b, double x){
    return gsl_sf_hyperg_2F1(a1, a2, b, x);
}
double hg2F1rn(double a1, double a2, double b, double x){
    return gsl_sf_hyperg_2F1_renorm(a1, a2, b, x);
}
void getIforR(double* x, int* pn, double* p, double* q, double* res){
    int i;
    for(i=0; i<(*pn); i++){
        if(x[i]<=0.0){res[i]=0.0;}else if(x[i]>=1.0){res[i]=1.0;}else{
            res[i] = getI(x[i], *p, *q);
        }
    }
}
double getI(double x, double p, double q){
    if(x>p/(p+q)){
        return 1.0 - gsl_sf_beta_inc(q, p, 1.0-x);
    }else{
        return gsl_sf_beta_inc(p, q, x);
    }
}
double getI2(double x, double p, double q){
    if(fabs(-x/(1.0-x))>=1.0){
        return 1.0 - hg2F1rn(1.0-p, 1.0, q+1.0, -(1.0-x)/x) * exp(lgamma(q+1.0) + logK(1.0-x, q, p));
    }else{
        return hg2F1rn(1.0-q, 1.0, p+1.0, -x/(1.0-x)) * exp(lgamma(p+1.0) + logK(x, p, q));
    }
}
double logI(double x, double p, double q){
    if(x/(1.0-x)>=1.0){
        return log(1.0 - hg2F1(1.0-p, 1.0, q+1.0, -(1.0-x)/x)* exp(logK(1.0-x, q, p)));
    }else{
        return log(hg2F1(1.0-q, 1.0, p+1.0, -x/(1.0-x))) + logK(x, p, q);
    }
}
// work length: n * 24
double logIn(double x, double p, double q, int n, double* res, double* work){
    int flag=0;
    double tmp;
    if(x>p/(p+q)){
        flag=1;
        x = 1.0-x;
        tmp = p;
        p=q;
        q=tmp;
    }
    double dn; 
    double *a,    *b,    *A,    *B;
    double *dap,  *dbp,  *dAp,  *dBp;
    double *daq,  *dbq,  *dAq,  *dBq;
    double *dap2, *dbp2, *dAp2, *dBp2;
    double *dapq, *dbpq, *dApq, *dBpq;
    double *daq2, *dbq2, *dAq2, *dBq2;
    a = work;//(double*)calloc(n, sizeof(double));
    dap = work+n;//(double*)calloc(n, sizeof(double));
    daq = work+2*n;//(double*)calloc(n, sizeof(double));
    dap2 = work+3*n;//(double*)calloc(n, sizeof(double));
    daq2 = work+4*n;//(double*)calloc(n, sizeof(double));
    dapq = work+5*n;//(double*)calloc(n, sizeof(double));
    
    b = work+6*n;//(double*)calloc(n, sizeof(double));
    dbp = work+7*n;//(double*)calloc(n, sizeof(double));
    dbq = work+8*n;//(double*)calloc(n, sizeof(double));
    dbp2 = work+9*n;//(double*)calloc(n, sizeof(double));
    dbpq = work+10*n;//(double*)calloc(n, sizeof(double));
    dbq2 = work+11*n;//(double*)calloc(n, sizeof(double));
    
    
    A = work+12*n;//(double*)calloc(n, sizeof(double));
    dAp = work+13*n;//
    dAq = work+14*n;//
    dAp2 = work+15*n;//
    dApq = work+16*n;//
    dAq2 = work+17*n;//

    B = work+18*n;//(double*)calloc(n, sizeof(double));
    dBp = work+19*n;//
    dBq = work+20*n;//
    dBp2 = work+21*n;//
    dBpq = work+22*n;//
    dBq2 = work+23*n;//
    double lp = log(p);
    double lq = log(q);
    double ldn;
    
    //double lf = lq + log(x) - lp - log1p(-x);
    double xt = x/(1.0-x);
    //double lxt = log(x) - log1p(-x);
    
    verbose=0;
    
    int i;
    double tmp1, dltp, dltq, dltp2, dltpq, dltq2;
    double tmp2, tmp3, dlt1p, dlt1q, dlt1p2, dlt1pq, dlt1q2, dlt3p, dlt3p2, dta1, dta2, dt1p;
    for(i=0; i<n; i++){
        dn = (double)(i+1);
        ldn = log(dn);
        if(i==0){
            a[i]    = (q-1.0) * xt * exp(- log1p(p));
            dap[i]  = - a[i]/(1.0+p);
            daq[i]  = xt * exp(- log1p(p));
            dap2[i] = - dap[i]/(1.0+p) + a[i]/pow(1.0+p,2.0);
            dapq[i] = - daq[i]/(1.0+p);
            daq2[i] = 0.0;
        }else{
            // tmp1 can be dtq as well
            tmp1  = pow(xt,2.0)*exp( log(dn-1.0) + log(p+q+dn-2.0) + log(p+dn-1.0) - log(p+2.0*dn-3.0) - 2.0*log(p+2.0*dn-2.0) - log(p+2.0*dn-1.0) );
            dltp  = (                            1./(p+q+dn-2.0) + 1./(p+dn-1.0) - 1./(p+2.0*dn-3.0) - 2.0/(p+2.0*dn-2.0)    - 1./(p+2.0*dn-1.0) );
            dltq  = (                            1./(p+q+dn-2.0));
            dltp2 = ( - 1./pow(p+q+dn-2.0,2.0) - 1./pow(p+dn-1.0,2.0) + 1./pow(p+2.0*dn-3.0,2.0) + 2.0/pow(p+2.0*dn-2.0,2.0) + 1./pow(p+2.0*dn-1.0,2.0) );
            dltpq = - dltq * dltq;
            dltq2 = dltpq;
            // d tmp1 / dp = tmp1 * d log_tmp1 /dp = tmp1 * dltp
            // d tmp1 / dq = tmp1 * d log_tmp1 /dq = tmp1 * dltq
            a[i]    = (q-dn) * tmp1;
            dap[i]  = a[i] * dltp;        // (q-dn) * tmp1 * d log_tmp1/dp
            daq[i]  = a[i] * dltq + tmp1; // (q-dn) * tmp1 * d log_tmp1/dq + tmp1
            dap2[i] = dap[i]*dltp + a[i]*dltp2;
            dapq[i] = dap[i]*dltq + a[i]*dltpq + tmp1*dltp;
            daq2[i] = daq[i]*dltq + a[i]*dltq2 + tmp1*dltq;
            
        }
        tmp1   = exp( log(2.0) + ldn + log(dn+p-1.0) + log(2.0+xt) - log(p+2.0*dn-2.0)        - log(p+2.0*dn) );// f(p)
        dlt1p  =                       1./(dn+p-1.0)               - 1./(p+2.0*dn-2.0)        - 1./(p+2.0*dn);
        dlt1p2 =                     - 1./pow(dn+p-1.0,2.0)        + 1./pow(p+2.0*dn-2.0,2.0) + 1./pow(p+2.0*dn,2.0);
        
        tmp2   = (p-2.0-q*xt); // f(p,q)
        
        tmp3   = exp( lp                    - log(p+2.0*dn-2.0)        - log(p+2.0*dn) ); // f(p)
        dlt3p   =      1./p                 - 1./(p+2.0*dn-2.0)        - 1./(p+2.0*dn);
        dlt3p2  =    - 1./pow(dn+p-1.0,2.0) + 1./pow(p+2.0*dn-2.0,2.0) + 1./pow(p+2.0*dn,2.0);
	
        b[i]    = tmp1         + tmp2 * tmp3;
        dbp[i]  = tmp1 * dlt1p + tmp3 + tmp2 * tmp3 * dlt3p;
        dbq[i]  = (-xt) * tmp3;
        //dbp2[i] = tmp1 * pow(dlt1p,2.0) + tmp1 * dlt1p2     +     tmp3 * dlt3p      +       tmp3 * dlt3p + tmp2 * tmp3 * pow(dlt3p,2.0) + tmp2 * tmp3 * dlt3p2;
        dbp2[i] = tmp1 * (pow(dlt1p,2.0) + dlt1p2)     +     2.0 * tmp3 * dlt3p      +      tmp2 * tmp3 * (pow(dlt3p,2.0) + dlt3p2);
        dbpq[i] = dbq[i] * dlt3p;
        dbq2[i] = 0.0;
        if(verbose>1){printf("%d: a=(%lf %lf %lf %lf %lf %lf)  b=(%lf %lf %lf %lf %lf %lf)\n", i, a[i], dap[i], daq[i], dap2[i], dapq[i], daq2[i], b[i], dbp[i], dbq[i], dbp2[i], dbpq[i], dbq2[i]);}
    }
     A[0]    = a[0]    + b[0];
    dAp[0]  = dap[0]  + dbp[0];
    dAq[0]  = daq[0]  + dbq[0];
    dAp2[0] = dap2[0] + dbp2[0];
    dApq[0] = dapq[0] + dbpq[0];
    dAq2[0] = daq2[0] + dbq2[0];
    
     B[0]    = b[0];
    dBp[0]  = dbp[0];
    dBq[0]  = dbq[0];
    dBp2[0] = dbp2[0];
    dBpq[0] = dbpq[0];
    dBq2[0] = dbq2[0];
    
    A[1]    = a[1]    + b[1]*A[0];
    dAp[1]  = dap[1]  + dbp[1]*A[0] + b[1]*dAp[0];
    dAq[1]  = daq[1]  + dbq[1]*A[0] + b[1]*dAq[0];
    dAp2[1] = dap2[1] + dbp2[1]*A[0] + dbp[1]*dAp[0] + dbp[1]*dAp[0] + b[1]*dAp2[0];
    dApq[1] = dapq[1] + dbpq[1]*A[0] + dbq[1]*dAp[0] + dbp[1]*dAq[0] + b[1]*dApq[0];
    dAq2[1] = daq2[1] + dbq2[1]*A[0] + dbq[1]*dAq[0] + dbq[1]*dAq[0] + b[1]*dAq2[0];
    
    B[1]    = a[1]    + b[1]*B[0];
    dBp[1]  = dap[1]  + dbp[1]*B[0]  + b[1]*dBp[0];
    dBq[1]  = daq[1]  + dbq[1]*B[0]  + b[1]*dBq[0];
    dBp2[1] = dap2[1] + dbp2[1]*B[0] + dbp[1]*dBp[0] + dbp[1]*dBp[0] + b[1]*dBp2[0];
    dBpq[1] = dapq[1] + dbpq[1]*B[0] + dbq[1]*dBp[0] + dbp[1]*dBq[0] + b[1]*dBpq[0];
    dBq2[1] = daq2[1] + dbq2[1]*B[0] + dbq[1]*dBq[0] + dbq[1]*dBq[0] + b[1]*dBq2[0];
    
    for(i=2; i<n; i++){
         A[i]   =  a[i] * A[i-2] +  b[i]  *A[i-1];
        dAp[i]  = dap[i] *A[i-2] + dbp[i] *A[i-1] +  a[i] *dAp[i-2] +  b[i] *dAp[i-1];
        dAq[i]  = daq[i] *A[i-2] + dbq[i] *A[i-1] +  a[i] *dAq[i-2] +  b[i] *dAq[i-1];
        dAp2[i] = dap2[i]*A[i-2] + dbp2[i]*A[i-1] + dap[i]*dAp[i-2] + dbp[i]*dAp[i-1] + dap[i]*dAp[i-2] + dbp[i]*dAp[i-1] + a[i]*dAp2[i-2] + b[i]*dAp2[i-1];
        dApq[i] = dapq[i]*A[i-2] + dbpq[i]*A[i-1] + daq[i]*dAp[i-2] + dbq[i]*dAp[i-1] + dap[i]*dAq[i-2] + dbp[i]*dAq[i-1] + a[i]*dApq[i-2] + b[i]*dApq[i-1];
        dAq2[i] = daq2[i]*A[i-2] + dbq2[i]*A[i-1] + daq[i]*dAq[i-2] + dbq[i]*dAq[i-1] + daq[i]*dAq[i-2] + dbq[i]*dAq[i-1] + a[i]*dAq2[i-2] + b[i]*dAq2[i-1];
        
         B[i]   =  a[i]  *B[i-2] +  b[i]  *B[i-1];
        dBp[i]  = dap[i] *B[i-2] + dbp[i] *B[i-1]  + a[i] *dBp[i-2] +  b[i] *dBp[i-1];
        dBq[i]  = daq[i] *B[i-2] + dbq[i] *B[i-1]  + a[i] *dBq[i-2] +  b[i] *dBq[i-1];
        dBp2[i] = dap2[i]*B[i-2] + dbp2[i]*B[i-1] + dap[i]*dBp[i-2] + dbp[i]*dBp[i-1] + dap[i]*dBp[i-2] + dbp[i]*dBp[i-1] + a[i]*dBp2[i-2] + b[i]*dBp2[i-1];
        dBpq[i] = dapq[i]*B[i-2] + dbpq[i]*B[i-1] + daq[i]*dBp[i-2] + dbq[i]*dBp[i-1] + dap[i]*dBq[i-2] + dbp[i]*dBq[i-1] + a[i]*dBpq[i-2] + b[i]*dBpq[i-1];
        dBq2[i] = daq2[i]*B[i-2] + dbq2[i]*B[i-1] + daq[i]*dBq[i-2] + dbq[i]*dBq[i-1] + daq[i]*dBq[i-2] + dbq[i]*dBq[i-1] + a[i]*dBq2[i-2] + b[i]*dBq2[i-1];
        if(verbose>1)printf("%d: %lf %lf %lf %lf %lf %lf\n", i, A[i], B[i], dAp[i], dBp[i], dAq[i], dBq[i]);
    }
    if(verbose>0)printf("%d: %lf %lf %lf %lf %lf %lf\n", i, A[n-1], B[n-1], dAp[n-1], dBp[n-1], dAq[n-1], dBq[n-1]);
    double K, dKp, dKq, dKp2, dKpq, dKq2;
    
    // normal scale
    K = exp(logK(x, p, q));
    /*dKp  = K * (log(x)     - 1./p + digamma(p+q) - digamma(p));
    dKq  = K * (log1p(-x)         + digamma(p+q) - digamma(q));
    dKp2 = K * (pow(log(x) - 1./p + digamma(p+q) - digamma(p),2.0) + 1./p/p + trigamma(p+q) - trigamma(p));
    dKpq = K * (   (log(x) - 1./p + digamma(p+q) - digamma(p) )*( log1p(-x) + digamma(p+q) - digamma(q) ) + trigamma(p+q) );
    dKq2 = K * (pow(log1p(-x)     + digamma(p+q) - digamma(q),2.0) + trigamma(p+q) - trigamma(q));*/
    
    // log scale
    if(x<=0){
        dKp = dKq = dKp2 = dKpq = dKq2 = 0.0;
    }else{
        dKp  = (log(x)     - 1./p + digamma(p+q) - digamma(p));
        dKq  = (log1p(-x)         + digamma(p+q) - digamma(q));
        dKp2 = (pow(log(x) - 1./p + digamma(p+q) - digamma(p),2.0) + 1./p/p + trigamma(p+q) - trigamma(p));
        dKpq = (   (log(x) - 1./p + digamma(p+q) - digamma(p) )*( log1p(-x) + digamma(p+q) - digamma(q) ) + trigamma(p+q) );
        dKq2 = (pow(log1p(-x)     + digamma(p+q) - digamma(q),2.0) + trigamma(p+q) - trigamma(q));
    }
    // normal scale
    res[0] = exp(logK(x, p, q) + log(A[n-1]) - log(B[n-1]));
    //res[1] = K*dKp*A[n-1]/B[n-1] + K*(dAp[n-1]/B[n-1] - A[n-1]/B[n-1]*dBp[n-1]/B[n-1]);
    //res[2] = K*dKq*A[n-1]/B[n-1] + K*(dAq[n-1]/B[n-1] - A[n-1]/B[n-1]*dBq[n-1]/B[n-1]);
    res[1] = res[0] * (dKp + dAp[n-1]/A[n-1] - dBp[n-1]/B[n-1]);
    res[2] = res[0] * (dKq + dAq[n-1]/A[n-1] - dBq[n-1]/B[n-1]);
    
    res[3] = res[1] * (dKp + dAp[n-1]/A[n-1] - dBp[n-1]/B[n-1])
           + res[0] * (dKp2 + dAp2[n-1]/A[n-1] - dBp2[n-1]/B[n-1] - (pow(dKp,2.0) + pow(dAp[n-1]/A[n-1],2.0) - pow(dBp[n-1]/B[n-1],2.0)));
    
    res[4] = res[0] * (dKp + dAp[n-1]/A[n-1] - dBp[n-1]/B[n-1]) * (dKq + dAq[n-1]/A[n-1] - dBq[n-1]/B[n-1])
           + res[0] * (dKpq + dApq[n-1]/A[n-1] - dBpq[n-1]/B[n-1] - ((dKp)*(dKq)+ (dAp[n-1]/A[n-1])*(dAq[n-1]/A[n-1]) - (dBp[n-1]/B[n-1])*(dBq[n-1]/B[n-1])));
    
    res[5] = res[2] * (dKq + dAq[n-1]/A[n-1] - dBq[n-1]/B[n-1])
           + res[0] * (dKq2 + dAq2[n-1]/A[n-1] - dBq2[n-1]/B[n-1] - (pow(dKq,2.0) + pow(dAq[n-1]/A[n-1],2.0) - pow(dBq[n-1]/B[n-1],2.0)));
    
    // log scale
    /*res[0] = (logK(x, p, q) + log(A[n-1]) - log(B[n-1]));
    res[1] = dKp + dAp[n-1]/A[n-1] - dBp[n-1]/B[n-1];
    res[2] = dKq + dAq[n-1]/A[n-1] - dBq[n-1]/B[n-1];
    res[3] = dKp2 + dAp2[n-1]/A[n-1] - dBp2[n-1]/B[n-1] - (pow(dKp,2.0) + pow(dAp[n-1]/A[n-1],2.0) - pow(dBp[n-1]/B[n-1],2.0));
    res[4] = dKpq + dApq[n-1]/A[n-1] - dBpq[n-1]/B[n-1] - ((dKp)*(dKq)+ (dAp[n-1]/A[n-1])*(dAq[n-1]/A[n-1]) - (dBp[n-1]/B[n-1])*(dBq[n-1]/B[n-1]) );
    res[5] = dKq2 + dAq2[n-1]/A[n-1] - dBq2[n-1]/B[n-1] - (pow(dKq,2.0) + pow(dAq[n-1]/A[n-1],2.0) - pow(dBq[n-1]/B[n-1],2.0));*/
    
    if(flag){
        //if(verbose>0){printf("flipped! %lf %lf %lf %lf\n", logK(x, p, q), res[0], res[1], res[2]);}
        
        res[0] = 1.0-res[0];
        tmp    = res[1];
        res[1] = -res[2];
        res[2] = -tmp ;
        
        tmp    = res[3];
        res[3] = -res[5];
        res[4] = -res[4];
        res[5] = -tmp;
        
        //tmp1   = res[0];
        //res[0] = log1p(-exp(res[0]));
        //res[1] = - exp(tmp1-res[0])*res[2];
        //res[2] = - exp(tmp1-res[0])*tmp;
    } 
    return 1;
}
void logIforR(double* v, double* work){
    //v[3] = logI(v[0], v[1], v[2]);
    logIn(v[0], v[1], v[2], 100, v+3, work+1);
}

double kappa=11.0;
double omega=10.0;

double getLkhd(double* x, int* ntests, int n, double* abc){
    double lkhd=0.0;
    int i;
    double dn, dm;
    for(i=0; i<n; i++){
        dn = (double)(ntests[i]); dm = (pow(dn, abc[2])-1.0);
        //logIn(x[i], abc[0], abc[1], 100, lghI, work+27);
        //lkhd += abc[2]*log(dn) + dm*log(1.0-lghI[0]) + (abc[0]-1.0)*log(x[i]) + (abc[1]-1.0)*log1p(-x[i]) - lbeta(abc[0], abc[1]);
        lkhd += abc[2]*log(dn) + dm*log(getI(1.0-x[i], abc[1], abc[0])) + (abc[0]-1.0)*log(x[i]) + (abc[1]-1.0)*log1p(-x[i]) - lbeta(abc[0], abc[1]);
    }
    lkhd += 3.*kappa*log(omega) + (kappa-1.0)*(log(abc[0])+log(abc[1])+log(abc[2])) - omega*(abc[0]+abc[1]+abc[2]) - 3.*lgamma(kappa);
    return lkhd;
}

void nlmForR(double* x, int* ntests, int* pn, double* abc, double* lkhd, double* grad, double* hess, double* work, int* fixParam, int* pverb){
    (*lkhd) = nlm(x, ntests, (*pn), abc, grad, hess, work, fixParam, *pverb);
}


void gridSearch(double* x, int* ntests, int n, double* abc, double* abc1, int* fixParam, int verbose2){
    int itr, itr_step, j;
    double ssize=1.0;
    double lkhd=0.0, lkhd0=-1e20;
    for(itr=0; itr<10; itr++){
        for(j=2; j>=0; j--){if(fixParam[j]==0){
            ssize=0.1/pow((double)(itr+1),0.5);
            for(itr_step=0; itr_step<10; itr_step++){
                abc1[0]=abc[0]; abc1[1]=abc[1]; abc1[2]=abc[2];
                abc1[j] = exp( log(abc[j]) - ssize);
                if(abc1[j]<1e-2){abc1[j]=1e-2;}else if(abc1[j]>10.){abc1[j]=10.;}
                lkhd = getLkhd(x, ntests, n, abc1);
                if(verbose2>1){fprintf(stderr, " l=%lf l1=%lf a=%lf b=%lf c=%lf ss=%lf\n", lkhd0, lkhd, abc1[0], abc1[1], abc1[2], ssize);}
                if(lkhd > lkhd0){
                    lkhd0 = lkhd;
                    abc[j] = abc1[j];
                    break;
                }else if(ssize>0.0){ssize *= -1.0;}else{ssize/=-4.0;}
            }
        }}
    }
}

double nlm(double* x, int* ntests, int n, double* abc, double* grad, double* hess, double* work, int* fixParam, int verbose2){
    int i, itr, itr_step;
    double dn, dm;
    double lkhd=0.0, lkhd0=-1.0e20, lkhd1, xi=0.0;
    double* g;
    double* h;
    double* hinv;
    double* step;
    double* lghI;
    double* abc1;
    g    = work;   // 3
    h    = work+3; // 6
    hinv = work+9; // 6
    step = work+15;// 3
    lghI = work+18;// 6
    abc1 = work+24;// 3
    if(abc[0]<1e-2){abc[0]=1e-2;}else if(abc[0]>1e2){abc[0]=1e2;}
    if(abc[1]<1e-2){abc[1]=1e-2;}else if(abc[1]>1e2){abc[1]=1e2;}
    if(abc[2]<1e-2){abc[2]=1e-2;}else if(abc[2]>10.0){abc[2]=10.0;}
    gridSearch(x, ntests, n, abc, abc1, fixParam, verbose2);
    double ssize=1.0;
    int stuck=0;
    for(itr=0; itr<20; itr++){
        clear1(g,3);
        clear1(h,6);
        lkhd=0;
        for(i=0; i<n; i++){
            dn = (double)(ntests[i]);
            dm = (pow(dn, abc[2])-1.0);
            logIn(x[i], abc[0], abc[1], 10, lghI, work+27);
            lghI[0] = getI(x[i], abc[0], abc[1]);
            lkhd += abc[2]*log(dn) + dm*log(1.0-lghI[0]) + (abc[0]-1.0)*log(x[i]) + (abc[1]-1.0)*log1p(-x[i]) - lbeta(abc[0], abc[1]);
            g[0] += dm*(-lghI[1])/(1.0-lghI[0]) + log(x[i])    - digamma(abc[0]) + digamma(abc[0]+abc[1]);
            g[1] += dm*(-lghI[2])/(1.0-lghI[0]) + log1p(-x[i]) - digamma(abc[1]) + digamma(abc[0]+abc[1]);
            g[2] += log(dn) + pow(dn, abc[2])*log(dn)*log(1.0-lghI[0]);
            h[0] += dm*(-lghI[3])/(1.0-lghI[0]) - dm*pow(lghI[1],2.0)/pow(1.0-lghI[0],2.0) - trigamma(abc[0]) + trigamma(abc[0]+abc[1]);
            h[1] += dm*(-lghI[4])/(1.0-lghI[0]) - dm*lghI[1]*lghI[2] /pow(1.0-lghI[0],2.0)                    + trigamma(abc[0]+abc[1]);
            h[3] += dm*(-lghI[5])/(1.0-lghI[0]) - dm*pow(lghI[2],2.0)/pow(1.0-lghI[0],2.0) - trigamma(abc[1]) + trigamma(abc[0]+abc[1]);
            h[2] += pow(dn, abc[2])*log(dn)*(-lghI[1])/(1.0-lghI[0]);
            h[4] += pow(dn, abc[2])*log(dn)*(-lghI[2])/(1.0-lghI[0]);
            h[5] += pow(dn, abc[2])*log(dn)*log(dn)*log(1.0-lghI[0]);
        }
        // converting in log scale
        g[0] *= abc[0];
        g[1] *= abc[1];
        g[2] *= abc[2];
        h[0] *= abc[0]*abc[0];
        h[1] *= abc[0]*abc[1];
        h[2] *= abc[0]*abc[2];
        h[3] *= abc[1]*abc[1];
        h[4] *= abc[1]*abc[2];
        h[5] *= abc[2]*abc[2];
        h[0] += g[0];
        h[1] += g[1];
        h[2] += g[2];
        if(fixParam[0]==1){
            g[0]=0.0;
            h[0]=h[1]=h[2]=0.0;
        }
        if(fixParam[1]==1){
            g[1]=0.0;
            h[1]=h[3]=h[4]=0.0;
        }
        if(fixParam[2]==1){
            g[2]=0.0;
            h[2]=h[4]=h[5]=0.0;
        }
        
        //prior in log scale
        lkhd += 3.*kappa*log(omega) + (kappa-1.0)*(log(abc[0])+log(abc[1])+log(abc[2])) - omega*(abc[0]+abc[1]+abc[2]) - 3.*lgamma(kappa);
        g[0] += (kappa-1.0) - omega*abc[0];
        g[1] += (kappa-1.0) - omega*abc[1];
        g[2] += (kappa-1.0) - omega*abc[2];
        h[0] += - omega*abc[0];
        h[3] += - omega*abc[1];
        h[5] += - omega*abc[2];
        if(fabs(lkhd-lkhd0)<1e-5 && fabs(g[0])<1e-5 && fabs(g[1])<1e-5 && fabs(g[2])<1e-5){
            if(verbose2>0){fprintf(stderr, "l0=%lf l1=%lf a=%lf b=%lf c=%lf g1=%lf g2=%lf g3=%lf\n", lkhd0, lkhd, abc[0], abc[1], abc[2], g[0], g[1], g[2]);}
            return lkhd;
        }else{lkhd0=lkhd;}
        if(verbose2>0){fprintf(stderr, "l0=%lf l1=%lf a=%lf b=%lf c=%lf g1=%lf g2=%lf g3=%lf\n", lkhd0, lkhd, abc[0], abc[1], abc[2], g[0], g[1], g[2]);}
        //solve
        getStepSize(h, g, 0);
        ssize = 1.0;
        stuck++;
        for(itr_step=0; itr_step<30; itr_step++){
            abc1[0] = exp( log(abc[0]) - ssize*g[0] );
            abc1[1] = exp( log(abc[1]) - ssize*g[1] );
            abc1[2] = exp( log(abc[2]) - ssize*g[2] );
            if(abc1[0]<1e-2){abc1[0]=1e-2;}else if(abc1[0]>1e2){abc1[0]=1e2;}
            if(abc1[1]<1e-2){abc1[0]=1e-2;}else if(abc1[1]>1e2){abc1[1]=1e2;}
            if(abc1[2]<1e-2){abc1[0]=1e-2;}else if(abc1[2]>10.){abc1[2]=10.;}
            lkhd1 = getLkhd(x, ntests, n, abc1);
            if(verbose2>1){fprintf(stderr, "%d l=%lf l1=%lf a=%lf b=%lf c=%lf g1=%lf g2=%lf g3=%lf ss=%lf xi=%lf stuck=%d\n", itr_step, lkhd0, lkhd1, abc1[0], abc1[1], abc1[2], g[0], g[1], g[2], ssize, xi, stuck);}
            if(lkhd1 + xi > lkhd){
                stuck = 0;
                lkhd0 = lkhd1;
                abc[0] = abc1[0];
                abc[1] = abc1[1];
                abc[2] = abc1[2];
                break;
            }else if(ssize>0.0){ssize *= -1.0;}else{ssize/=-2.0;}
	}
        if(stuck>0){xi=(double)stuck;}else{xi=0.0;}
    }
    return -999999.9;
}






