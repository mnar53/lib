#ifndef ELEMENTARY
#define ELEMENTARY



#ifdef __cplusplus
extern "C" {
#endif

inline double rcp(const double x) { return 1 / x; }
inline double neg(const double x) { return -x; }
inline double sqr(const double x) { return x * x; }


//double poly(vector<double> cof, const double x);
double i0(const double x);
double i1(const double x);
double k0(const double x);
double k1(const double x);
double gammln(const double xx);
double factrl(const double x);
double factln(const double x);
double bico(const int n, const int k);
double beta(const double z, const double w);
double gser(const double a, const double x);
double gcf(const double a, const double x);
double gammpapprox(double a, double x, int psig);
double gammp(const double a, const double x);
double gammq(const double a, const double x);
double invgammp(double p, double a);
double betaiapprox(double a, double b, double x);
double betacf(const double a, const double b, const double x);
double betai(const double a, const double b, const double x);
double invbetai(double p, double a, double b);
double gammap(double alph, double bet, double x);
double gammacdf(double alph, double bet, double x);
double gammainvcdf(double alph, double bet, double p);
double betap(double alph, double bet, double x);
double betacdf(double alph, double bet, double x);
double betainvcdf(double alph, double bet, double p);
double studentp(double nu, double mu, double sig, double t);
double studentcdf(double nu, double mu, double sig, double t);
double studentinvcdf(double nu, double mu, double sig, double p);
double studentaa(double nu, double t);
double studentinvaa(double nu, double p);
double poissonp(double lam, int n);
double poissoncdf(double lam, int n);
int poissoninvcdf(double lam, double p);
double binomialp(int n, double pe, int k);
double binomialcdf(int n, double pe, int k);
int binomialinvcdf(int n, double pe, double p);
double chisqrp(double nu, double x2);
double chisqrcdf(double nu, double x2);
double chisqrinvcdf(double nu, double p);
double fdistp(double nu1, double nu2, double f);
double fdistcdf(double nu1, double nu2, double f);
double fdistinvcdf(double nu1, double nu2, double p);
double erfccheb(double z);
double inverfc(double p);
double erfcc(double x);
double normalp(double mu, double sig, double x);
double normalcdf(double mu, double sig, double x);
double normalinvcdf(double mu, double sig, double p);
double lognormalp(double mu, double sig, double x);
double lognormalcdf(double mu, double sig, double x);
double lognormalinvcdf(double mu, double sig, double p);
double expint(const int n, const double x);
double ei(double x);

inline double fact(double n) { return factrl((int) n); }
inline double factLn(double n) { return factln((int) n); }




#ifdef __cplusplus
}
#endif

#endif