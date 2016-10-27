
#include "specfunc.h"
#include <vector>
using namespace std;

template<class T> inline const T &MAX(const T &a, const T &b) {	return b > a ? (b) : (a); }
template<class T> inline const T &MIN(const T &a, const T &b) {	return b < a ? (b) : (a); }
template<class T> inline T SQR(const T a) { return a*a; }



#ifdef __cplusplus
extern "C" {
#endif



double poly (vector<double> cof, const double x) {
	int n = cof.size(); double ans = cof[n];
	for (int i = n-1; i >= 0; i--) 
		ans = ans*x + cof[i];
	return ans;
}



double i0(const double x) {
	vector<double> i0p = { 9.999999999999997e-1, 2.466405579426905e-1, 1.478980363444585e-2, 3.826993559940360e-4, 5.395676869878828e-6,
		4.700912200921704e-8, 2.733894920915608e-10, 1.115830108455192e-12, 3.301093025084127e-15, 7.209167098020555e-18, 1.166898488777214e-20,
		1.378948246502109e-23, 1.124884061857506e-26, 5.498556929587117e-30 },
	i0q = { 4.463598170691436e-1, 1.702205745042606e-3, 2.792125684538934e-6, 2.369902034785866e-9, 8.965900179621208e-13 },
	i0pp= { 1.192273748120670e-1, 1.947452015979746e-1, 7.629241821600588e-2, 8.474903580801549e-3, 2.023821945835647e-4 },
	i0qq = { 2.962898424533095e-1, 4.866115913196384e-1,1.938352806477617e-1, 2.261671093400046e-2, 6.450448095075585e-4, 1.529835782400450e-6 };
	double ax, z, y;
	if ((ax = abs(x)) < 15.0) {
		y = x*x;
		return poly(i0p, y) / poly(i0q, 225. - y);
	}
	else {
		z = 1.0 - 15.0 / ax;
		return exp(ax)*poly(i0pp, z) / (poly(i0qq,  z)*sqrt(ax));
	}
}



double i1(const double x) {
	vector<double> i1p = { 5.000000000000000e-1, 6.090824836578078e-2,2.407288574545340e-3, 4.622311145544158e-5, 5.161743818147913e-7,
		3.712362374847555e-9, 1.833983433811517e-11, 6.493125133990706e-14,	1.693074927497696e-16, 3.299609473102338e-19, 4.813071975603122e-22,
		5.164275442089090e-25, 3.846870021788629e-28, 1.712948291408736e-31 },
	i1q = { 4.665973211630446e-1, 1.677754477613006e-3,	2.583049634689725e-6, 2.045930934253556e-9, 7.166133240195285e-13 },
	i1pp = { 1.286515211317124e-1, 1.930915272916783e-1, 6.965689298161343e-2, 7.345978783504595e-3, 1.963602129240502e-4 },
	i1qq = { 3.309385098860755e-1, 4.878218424097628e-1, 1.663088501568696e-1, 1.473541892809522e-2, 1.964131438571051e-4, -1.034524660214173e-6 };
	double ax, z, y, ans;
	if ((ax = abs(x)) < 15.0) {
		y = x*x;
		return x*poly(i1p, y) / poly(i1q, 225. - y);
	}
	else {
		z = 1.0 - 15.0 / ax;
		ans = exp(ax)*poly(i1pp, z) / (poly(i1qq, z)*sqrt(ax));
		return x > 0.0 ? ans : -ans;
	}
}



double k0(const double x) {
	vector<double> k0pi = { 1.0, 2.346487949187396e-1, 1.187082088663404e-2, 2.150707366040937e-4, 1.425433617130587e-6 },
	k0qi = { 9.847324170755358e-1, 1.518396076767770e-2,8.362215678646257e-5 },
	k0p = { 1.159315156584126e-1, 2.770731240515333e-1,	2.066458134619875e-2, 4.574734709978264e-4, 3.454715527986737e-6 },
	k0q = { 9.836249671709183e-1, 1.627693622304549e-2,	9.809660603621949e-5 },
	k0pp = { 1.253314137315499, 1.475731032429900e1, 6.123767403223466e1, 1.121012633939949e2, 9.285288485892228e1,
		3.198289277679660e1, 3.595376024148513, 6.160228690102976e-2 },
	k0qq = { 1.0, 1.189963006673403e1, 5.027773590829784e1,	9.496513373427093e1, 8.318077493230258e1, 3.181399777449301e1,	4.443672926432041, 1.408295601966600e-1 };
	double z, term;
	if (x <= 1.0) {
		z = x*x;
		term = poly(k0pi, z)*log(x) / poly(k0qi, 1. - z);
		return poly(k0p, z) / poly(k0q, 1.- z) - term;
	}
	else {
		z = 1.0 / x;
		return exp(-x)*poly(k0pp, z) / (poly(k0qq, z)*sqrt(x));
	}
}


double k1(const double x) {
	vector<double> k1pi = { 0.5, 5.598072040178741e-2, 1.818666382168295e-3, 2.397509908859959e-5, 1.239567816344855e-7 },
	k1qi = { 9.870202601341150e-1, 1.292092053534579e-2, 5.881933053917096e-5 },
	k1p = { -3.079657578292062e-1, -8.109417631822442e-2, -3.477550948593604e-3, -5.385594871975406e-5, -3.110372465429008e-7 },
	k1q = { 9.861813171751389e-1, 1.375094061153160e-2,	6.774221332947002e-5 },
	k1pp = { 1.253314137315502, 1.457171340220454e1,6.063161173098803e1, 1.147386690867892e2, 1.040442011439181e2,
		4.356596656837691e1, 7.265230396353690, 3.144418558991021e-1 },
	k1qq = { 1.0, 1.125154514806458e1, 4.427488496597630e1,	7.616113213117645e1, 5.863377227890893e1, 1.850303673841586e1,
		1.857244676566022, 2.538540887654872e-2 };	
	double z, term;
	if (x <= 1.0) {
		z = x*x;
		term = poly(k1pi, z)*log(x) / poly(k1qi, 1. - z);
		return x*(poly(k1p, z) / poly(k1q, 1. - z) + term) + 1. / x;
	}
	else {
		z = 1.0 / x;
		return exp(-x)*poly(k1pp, z) / (poly(k1qq, z)*sqrt(x));
	}
}


double gammln(const double xx) {
	int j;	double x, tmp, y, ser;
	vector<double> cof = { 57.1562356658629235, -59.5979603554754912, 14.1360979747417471, -0.491913816097620199, .339946499848118887e-4,
		.465236289270485756e-4, -.983744753048795646e-4, .158088703224912494e-3,-.210264441724104883e-3, .217439618115212643e-3, -.164318106536763890e-3,
		.844182239838527433e-4, -.261908384015814087e-4, .368991826595316234e-5 };
	if (xx <= 0) throw("bad arg in gammln");
	y = x = xx;
	tmp = x + 5.24218750000000000;
	tmp = (x + 0.5)*log(tmp) - tmp;
	ser = 0.999999999999997092;
	for (j = 0; j<14; j++) 
		ser += cof[j] / ++y;
	return tmp + log(2.5066282746310005*ser / x);
}


double factrl(const double x) {
	static vector<double> a(171); static bool init = true;
	if (init) {
		init = false;
		a[0] = 1.;
		for (int i = 1; i<171; i++) a[i] = i*a[i - 1];
	}
	int n = floor(x);
	if (n < 0 || n > 170) throw("factrl out of range");
	return a[n];
}



double factln(const double x) {
	static const int NTOP = 2000;
	static vector<double> a(NTOP); static bool init = true;
	if (init) {
		init = false;
		for (int i = 0; i<NTOP; i++) 
			a[i] = gammln(i + 1.);
	}
	int n = floor(x);
	if (n < 0) throw("negative arg in factln");
	if (n < NTOP) return a[n];
	return gammln(n + 1.);
}


double bico(const int n, const int k) {
	if (n<0 || k<0 || k>n) throw("bad args in bico");
	if (n<171) return floor(0.5 + factrl(n) / (factrl(k)*factrl(n - k)));
	return floor(0.5 + exp(factln(n) - factln(k) - factln(n - k)));
}


double beta(const double z, const double w) {
	return exp(gammln(z) + gammln(w) - gammln(z + w));
}





#define ASWITCH 100
#define EPS  DBL_EPSILON
#define FPMIN  DBL_MIN / EPS
#define EULER  0.577215664901533
#define BIG  DBL_MAX*EPS

double gser(const double a, const double x) {
	double sum, del, ap, gln;
	gln = gammln(a);
	ap = a;
	del = sum = 1.0 / a;
	for (;;) {
		++ap;
		del *= x / ap;
		sum += del;
		if (fabs(del) < fabs(sum)*EPS) {
			return sum*exp(-x + a*log(x) - gln);
		}
	}
}

double gcf(const double a, const double x) {
	int i;	double an, b, c, d, del, h, gln;
	gln = gammln(a);
	b = x + 1.0 - a;
	c = 1.0 / FPMIN;
	d = 1.0 / b;
	h = d;
	for (i = 1;; i++) {
		an = -i*(i - a);
		b += 2.0;
		d = an*d + b;
		if (fabs(d) < FPMIN) d = FPMIN;
		c = b + an / c;
		if (fabs(c) < FPMIN) c = FPMIN;
		d = 1.0 / d;
		del = d*c;
		h *= del;
		if (fabs(del - 1.0) <= EPS) break;
	}
	return exp(-x + a*log(x) - gln)*h;
}


static vector<double> gamy = { 0.0021695375159141994, 0.011413521097787704, 0.027972308950302116, 0.051727015600492421,
0.082502225484340941, 0.12007019910960293, 0.16415283300752470, 0.21442376986779355, 0.27051082840644336, 0.33199876341447887,
0.39843234186401943, 0.46931971407375483, 0.54413605556657973, 0.62232745288031077, 0.70331500465597174, 0.78649910768313447,
0.87126389619061517, 0.95698180152629142 },
gamw = { 0.0055657196642445571, 0.012915947284065419, 0.020181515297735382, 0.027298621498568734, 0.034213810770299537,
0.040875750923643261, 0.047235083490265582, 0.053244713977759692, 0.058860144245324798, 0.064039797355015485, 0.068745323835736408,
0.072941885005653087, 0.076598410645870640, 0.079687828912071670, 0.082187266704339706, 0.084078218979661945, 0.085346685739338721, 0.085983275670394821 };


double gammpapprox(double a, double x, int psig) {
	int j; double xu, t, sum, ans, gln, a1 = a - 1.0, lna1 = log(a1), sqrta1 = sqrt(a1);
	gln = gammln(a);
	if (x > a1) xu = MAX(a1 + 11.5*sqrta1, x + 6.0*sqrta1);
	else xu = MAX(0., MIN(a1 - 7.5*sqrta1, x - 5.0*sqrta1));
	sum = 0;
	for (j = 0; j < gamy.size(); j++) {
		t = x + (xu - x) * gamy[j];
		sum += gamw[j] * exp(-(t - a1) + a1*(log(t) - lna1));
	}
	ans = sum*(xu - x)*exp(a1*(lna1 - 1.) - gln);
	return (psig ? (ans>0.0 ? 1.0 - ans : -ans) : (ans >= 0.0 ? ans : 1.0 + ans));
}

double gammp(const double a, const double x) {
	if (x < 0.0 || a <= 0.0) throw("bad args in gammp");
	if (x == 0.0) return 0.0;
	else if ((int)a >= ASWITCH) return gammpapprox(a, x, 1);
	else if (x < a + 1.0) return gser(a, x);
	else return 1.0 - gcf(a, x);
}

double gammq(const double a, const double x) {
	if (x < 0.0 || a <= 0.0) throw("bad args in gammq");
	if (x == 0.0) return 1.0;
	else if ((int)a >= ASWITCH) return gammpapprox(a, x, 0);
	else if (x < a + 1.0) return 1.0 - gser(a, x);
	else return gcf(a, x);
}


double invgammp(double p, double a) {
	int j;	double x, err, t, u, pp, lna1, afac, gln, a1 = a - 1;
	const double EPSLN = 1.e-8;
	gln = gammln(a);
	if (a <= 0.) throw("a must be pos in invgammap");
	if (p >= 1.) return MAX(100., a + 100.*sqrt(a));
	if (p <= 0.) return 0.0;
	if (a > 1.) {
		lna1 = log(a1);
		afac = exp(a1*(lna1 - 1.) - gln);
		pp = (p < 0.5) ? p : 1. - p;
		t = sqrt(-2.*log(pp));
		x = (2.30753 + t*0.27061) / (1. + t*(0.99229 + t*0.04481)) - t;
		if (p < 0.5) x = -x;
		x = MAX(1.e-3, a*pow(1. - 1. / (9.*a) - x / (3.*sqrt(a)), 3));
	}
	else {
		t = 1.0 - a*(0.253 + a*0.12);
		if (p < t) x = pow(p / t, 1. / a);
		else x = 1. - log(1. - (p - t) / (1. - t));
	}
	for (j = 0; j<12; j++) {
		if (x <= 0.0) return 0.0;
		err = gammp(a, x) - p;
		if (a > 1.) t = afac*exp(-(x - a1) + a1*(log(x) - lna1));
		else t = exp(-x + a1*log(x) - gln);
		u = err / t;
		x -= (t = u / (1. - 0.5*MIN(1., u*((a - 1.) / x - 1))));
		if (x <= 0.) x = 0.5*(x + t);
		if (fabs(t) < EPSLN*x) break;
	}
	return x;
}


#define SWITCH 3000


double betaiapprox(double a, double b, double x) {
	int j;
	double xu, t, sum, ans;
	double a1 = a - 1.0, b1 = b - 1.0, mu = a / (a + b);
	double lnmu = log(mu), lnmuc = log(1. - mu);
	t = sqrt(a*b / (SQR(a + b)*(a + b + 1.0)));
	if (x > a / (a + b)) {
		if (x >= 1.0) return 1.0;
		xu = MIN(1., MAX(mu + 10.*t, x + 5.0*t));
	}
	else {
		if (x <= 0.0) return 0.0;
		xu = MAX(0., MIN(mu - 10.*t, x - 5.0*t));
	}
	sum = 0;
	for (j = 0; j<18; j++) {
		t = x + (xu - x) * gamy[j];
		sum += gamw[j] * exp(a1*(log(t) - lnmu) + b1*(log(1 - t) - lnmuc));
	}
	ans = sum*(xu - x)*exp(a1*lnmu - gammln(a) + b1*lnmuc - gammln(b) + gammln(a + b));
	return ans>0.0 ? 1.0 - ans : -ans;
}


double betacf(const double a, const double b, const double x) {
	int m, m2;
	double aa, c, d, del, h, qab, qam, qap;
	qab = a + b;
	qap = a + 1.0;
	qam = a - 1.0;
	c = 1.0;
	d = 1.0 - qab*x / qap;
	if (fabs(d) < FPMIN) d = FPMIN;
	d = 1.0 / d;
	h = d;
	for (m = 1; m<10000; m++) {
		m2 = 2 * m;
		aa = m*(b - m)*x / ((qam + m2)*(a + m2));
		d = 1.0 + aa*d;
		if (fabs(d) < FPMIN) d = FPMIN;
		c = 1.0 + aa / c;
		if (fabs(c) < FPMIN) c = FPMIN;
		d = 1.0 / d;
		h *= d*c;
		aa = -(a + m)*(qab + m)*x / ((a + m2)*(qap + m2));
		d = 1.0 + aa*d;
		if (fabs(d) < FPMIN) d = FPMIN;
		c = 1.0 + aa / c;
		if (fabs(c) < FPMIN) c = FPMIN;
		d = 1.0 / d;
		del = d*c;
		h *= del;
		if (fabs(del - 1.0) <= EPS) break;
	}
	return h;
}


double betai(const double a, const double b, const double x) {
	double bt;
	if (a <= 0.0 || b <= 0.0) throw("Bad a or b in routine betai");
	if (x < 0.0 || x > 1.0) throw("Bad x in routine betai");
	if (x == 0.0 || x == 1.0) return x;
	if (a > SWITCH && b > SWITCH) return betaiapprox(a, b, x);
	bt = exp(gammln(a + b) - gammln(a) - gammln(b) + a*log(x) + b*log(1.0 - x));
	if (x < (a + 1.0) / (a + b + 2.0)) return bt*betacf(a, b, x) / a;
	else return 1.0 - bt*betacf(b, a, 1.0 - x) / b;
}





double invbetai(double p, double a, double b) {
	const double EPSLN = 1.e-8;
	double pp, t, u, err, x, al, h, w, afac, a1 = a - 1., b1 = b - 1.;
	int j;
	if (p <= 0.) return 0.;
	else if (p >= 1.) return 1.;
	else if (a >= 1. && b >= 1.) {
		pp = (p < 0.5) ? p : 1. - p;
		t = sqrt(-2.*log(pp));
		x = (2.30753 + t*0.27061) / (1. + t*(0.99229 + t*0.04481)) - t;
		if (p < 0.5) x = -x;
		al = (SQR(x) - 3.) / 6.;
		h = 2. / (1. / (2.*a - 1.) + 1. / (2.*b - 1.));
		w = (x*sqrt(al + h) / h) - (1. / (2.*b - 1) - 1. / (2.*a - 1.))*(al + 5. / 6. - 2. / (3.*h));
		x = a / (a + b*exp(2.*w));
	}
	else {
		double lna = log(a / (a + b)), lnb = log(b / (a + b));
		t = exp(a*lna) / a;
		u = exp(b*lnb) / b;
		w = t + u;
		if (p < t / w) x = pow(a*w*p, 1. / a);
		else x = 1. - pow(b*w*(1. - p), 1. / b);
	}
	afac = -gammln(a) - gammln(b) + gammln(a + b);
	for (j = 0; j<10; j++) {
		if (x == 0. || x == 1.) return x;
		err = betai(a, b, x) - p;
		t = exp(a1*log(x) + b1*log(1. - x) + afac);
		u = err / t;
		x -= (t = u / (1. - 0.5*MIN(1., u*(a1 / x - b1 / (1. - x)))));
		if (x <= 0.) x = 0.5*(x + t);
		if (x >= 1.) x = 0.5*(x + t + 1.);
		if (fabs(t) < EPSLN*x && j > 0) break;
	}
	return x;
}



double gammap(double alph, double bet, double x) {
	if (x <= 0.) throw("bad x in Gammadist");
	return exp(-bet*x + (alph - 1.)*log(x) + alph*log(bet) - gammln(alph));
}

double gammacdf(double alph, double bet, double x) {
	if (x < 0.) throw("bad x in Gammadist");
	return gammp(alph, bet*x);
}

double gammainvcdf(double alph, double bet, double p) {
	if (p < 0. || p >= 1.) throw("bad p in Gammadist");
	return invgammp(p, alph) / bet;
}




double betap(double alph, double bet, double x) {
	if (x <= 0. || x >= 1.) throw("bad x in Betadist");
	return exp((alph - 1.)*log(x) + (bet - 1.)*log(1. - x) + gammln(alph + bet) - gammln(alph) - gammln(bet));
}
double betacdf(double alph, double bet, double x) {
	if (x < 0. || x > 1.) throw("bad x in Betadist");
	return betai(alph, bet, x);
}
double betainvcdf(double alph, double bet, double p) {
	if (p < 0. || p > 1.) throw("bad p in Betadist");
	return invbetai(p, alph, bet);
}


double studentp(double nu, double mu, double sig,  double t) {
	double np = 0.5*(nu + 1.);
	return exp(-np*log(1. + SQR((t - mu) / sig) / nu) + gammln(np) - gammln(0.5*nu))
		/ (sqrt(3.14159265358979324*nu)*sig);
}
double studentcdf(double nu, double mu, double sig, double t) {
	double p = 0.5*betai(0.5*nu, 0.5, nu / (nu + SQR((t - mu) / sig)));
	if (t >= mu) return 1. - p;
	else return p;
}
double studentinvcdf(double nu, double mu, double sig, double p) {
	if (p <= 0. || p >= 1.) throw("bad p in Studentdist");
	double x = invbetai(2.*MIN(p, 1. - p), 0.5*nu, 0.5);
	x = sig*sqrt(nu*(1. - x) / x);
	return (p >= 0.5 ? mu + x : mu - x);
}
double studentaa(double nu, double t) {
	if (t < 0.) throw("bad t in Studentdist");
	return 1. - betai(0.5*nu, 0.5, nu / (nu + SQR(t)));
}
double studentinvaa(double nu, double p) {
	if (p < 0. || p >= 1.) throw("bad p in Studentdist");
	double x = invbetai(1. - p, 0.5*nu, 0.5);
	return sqrt(nu*(1. - x) / x);
}



double poissonp(double lam, int n) {
	if (n < 0) throw("bad n in Poissondist");
	return exp(-lam + n*log(lam) - gammln(n + 1.));
}
double poissoncdf(double lam, int n) {
	if (n < 0) throw("bad n in Poissondist");
	if (n == 0) return 0.;
	return gammq((double)n, lam);
}
int poissoninvcdf(double lam, double p) {
	int n, nl, nu, inc = 1;
	if (p <= 0. || p >= 1.) throw("bad p in Poissondist");
	if (p < exp(-lam)) return 0;
	n = (int)MAX(sqrt(lam), 5.);
	if (p < poissoncdf(lam, n)) {
		do {
			n = MAX(n - inc, 0);
			inc *= 2;
		} while (p < poissoncdf(lam, n));
		nl = n; nu = n + inc / 2;
	}
	else {
		do {
			n += inc;
			inc *= 2;
		} while (p > poissoncdf(lam, n));
		nu = n; nl = n - inc / 2;
	}
	while (nu - nl>1) {
		n = (nl + nu) / 2;
		if (p < poissoncdf(lam, n)) nu = n;
		else nl = n;
	}
	return nl;
}


double binomialp(int n, double pe, int k) {
	if (k < 0) throw("bad k in Binomialdist");
	if (k > n) return 0.;
	return exp(k*log(pe) + (n - k)*log(1. - pe)
		+ gammln(n + 1.) - gammln(k + 1.) - gammln(n - k + 1.));
}
double binomialcdf(int n, double pe, int k) {
	if (k < 0) throw("bad k in Binomialdist");
	if (k == 0) return 0.;
	if (k > n) return 1.;
	return 1. - betai((double)k, n - k + 1., pe);
}
int binomialinvcdf(int n, double pe, double p) {
	int k, kl, ku, inc = 1;
	if (p <= 0. || p >= 1.) throw("bad p in Binomialdist");
	k = MAX(0, MIN(n, (int)(n*pe)));
	if (p < binomialcdf(n,pe,k)) {
		do {
			k = MAX(k - inc, 0);
			inc *= 2;
		} while (p < binomialcdf(n,pe,k));
		kl = k; ku = k + inc / 2;
	}
	else {
		do {
			k = MIN(k + inc, n + 1);
			inc *= 2;
		} while (p > binomialcdf(n,pe,k));
		ku = k; kl = k - inc / 2;
	}
	while (ku - kl>1) {
		k = (kl + ku) / 2;
		if (p < binomialcdf(n,pe,k)) ku = k;
		else kl = k;
	}
	return kl;
}


double chisqrp(double nu, double x2) {
	if (x2 <= 0.) throw("bad x2 in Chisqdist");
	return exp(-0.5*(x2 - (nu - 2.)*log(x2)) - (0.693147180559945309*(0.5*nu) + gammln(0.5*nu)));
}
double chisqrcdf(double nu, double x2) {
	if (x2 < 0.) throw("bad x2 in Chisqdist");
	return gammp(0.5*nu, 0.5*x2);
}
double chisqrinvcdf(double nu, double p) {
	if (p < 0. || p >= 1.) throw("bad p in Chisqdist");
	return 2.*invgammp(p, 0.5*nu);
}

double fdistp(double nu1, double nu2, double f) {
	if (f <= 0.) throw("bad f in Fdist");
	return exp((0.5*nu1 - 1.)*log(f) - 0.5*(nu1 + nu2)*log(nu2 + nu1*f) + 
		0.5*(nu1*log(nu1) + nu2*log(nu2)) + gammln(0.5*(nu1 + nu2))	- gammln(0.5*nu1) - gammln(0.5*nu2));
}
double fdistcdf(double nu1, double nu2, double f) {
	if (f < 0.) throw("bad f in Fdist");
	return betai(0.5*nu1, 0.5*nu2, nu1*f / (nu2 + nu1*f));
}
double fdistinvcdf(double nu1, double nu2, double p) {
	if (p <= 0. || p >= 1.) throw("bad p in Fdist");
	double x = invbetai(p, 0.5*nu1, 0.5*nu2);
	return nu2*x / (nu1*(1. - x));
}





double erfccheb(double z){
	vector<double> cof = { -1.3026537197817094, 6.4196979235649026e-1, 1.9476473204185836e-2, -9.561514786808631e-3, -9.46595344482036e-4,
		3.66839497852761e-4, 4.2523324806907e-5, -2.0278578112534e-5, -1.624290004647e-6, 1.303655835580e-6, 1.5626441722e-8, -8.5238095915e-8,
		6.529054439e-9, 5.059343495e-9, -9.91364156e-10, -2.27365122e-10, 9.6467911e-11, 2.394038e-12, -6.886027e-12, 8.94487e-13, 3.13092e-13,
		-1.12708e-13, 3.81e-16, 7.106e-15, -1.523e-15, -9.4e-17, 1.21e-16, -2.8e-17 };
	int j;	double t, ty, tmp, d = 0., dd = 0.;
	if (z < 0.) throw("erfccheb requires nonnegative argument");
	t = 2. / (2. + z);
	ty = 4.*t - 2.;
	for (j = cof.size() - 1; j>0; j--) {
		tmp = d;
		d = ty*d - dd + cof[j];
		dd = tmp;
	}
	return t*exp(-z*z + 0.5*(cof[0] + ty*d) - dd);
}

double inverfc(double p) {
	double x, err, t, pp;
	if (p >= 2.0) return -100.;
	if (p <= 0.0) return 100.;
	pp = (p < 1.0) ? p : 2. - p;
	t = sqrt(-2.*log(pp / 2.));
	x = -0.70711*((2.30753 + t*0.27061) / (1. + t*(0.99229 + t*0.04481)) - t);
	for (int j = 0; j<2; j++) {
		err = erfc(x) - pp;
		x += err / (1.12837916709551257*exp(-SQR(x)) - x*err);
	}
	return (p < 1.0 ? x : -x);
}






double erfcc(double x) {
	double t, z = fabs(x), ans;
	t = 2. / (2. + z);
	ans = t*exp(-z*z - 1.26551223 + t*(1.00002368 + t*(0.37409196 + t*(0.09678418 +
		t*(-0.18628806 + t*(0.27886807 + t*(-1.13520398 + t*(1.48851587 +
		t*(-0.82215223 + t*0.17087277)))))))));
	return (x >= 0.0 ? ans : 2.0 - ans);
}


double normalp(double mu, double sig, double x) {
	return (0.398942280401432678 / sig)*exp(-0.5*SQR((x - mu) / sig));
}
double normalcdf(double mu, double sig, double x) {
	return 0.5*erfc(-0.707106781186547524*(x - mu) / sig);
}
double normalinvcdf(double mu, double sig, double p) {
	if (p <= 0. || p >= 1.) throw("bad p in Normaldist");
	return -1.41421356237309505*sig*inverfc(2.*p) + mu;
}



double lognormalp(double mu, double sig, double x) {
	if (x < 0.) throw("bad x in Lognormaldist");
	if (x == 0.) return 0.;
	return (0.398942280401432678 / (sig*x))*exp(-0.5*SQR((log(x) - mu) / sig));
}
double lognormalcdf(double mu, double sig, double x) {
	if (x < 0.) throw("bad x in Lognormaldist");
	if (x == 0.) return 0.;
	return 0.5*erfc(-0.707106781186547524*(log(x) - mu) / sig);
}
double lognormalinvcdf(double mu, double sig, double p) {
	if (p <= 0. || p >= 1.) throw("bad p in Lognormaldist");
	return exp(-1.41421356237309505*sig*inverfc(2.*p) + mu);
}


double expint(const int n, const double x)
{	static const int MAXIT = 100; int i, ii, nm1 = n - 1;
	double a, b, c, d, del, fact, h, psi, ans;
	if (n < 0 || x < 0.0 || (x == 0.0 && (n == 0 || n == 1)))
		throw("bad arguments in expint");
	if (n == 0) ans = exp(-x) / x;
	else {
		if (x == 0.0) ans = 1.0 / nm1;
		else {
			if (x > 1.0) {
				b = x + n;
				c = BIG;
				d = 1.0 / b;
				h = d;
				for (i = 1; i <= MAXIT; i++) {
					a = -i*(nm1 + i);
					b += 2.0;
					d = 1.0 / (a*d + b);
					c = b + a / c;
					del = c*d;
					h *= del;
					if (abs(del - 1.0) <= EPS) {
						ans = h*exp(-x);
						return ans;
					}
				}
				throw("continued fraction failed in expint");
			}
			else {
				ans = (nm1 != 0 ? 1.0 / nm1 : -log(x) - EULER);
				fact = 1.0;
				for (i = 1; i <= MAXIT; i++) {
					fact *= -x / i;
					if (i != nm1) del = -fact / (i - nm1);
					else {
						psi = -EULER;
						for (ii = 1; ii <= nm1; ii++) psi += 1.0 / ii;
						del = fact*(-log(x) + psi);
					}
					ans += del;
					if (abs(del) < abs(ans)*EPS) return ans;
				}
				throw("series failed in expint");
			}
		}
	}
	return ans;
}

double ei(double x) {
	static const int MAXIT = 100; int k;
	double fact, prev, sum, term;
	if (x <= 0.0) throw("Bad argument in ei");
	if (x < FPMIN) return log(x) + EULER;
	if (x <= -log(EPS)) {
		sum = 0.0;
		fact = 1.0;
		for (k = 1; k <= MAXIT; k++) {
			fact *= x / k;
			term = fact / k;
			sum += term;
			if (term < EPS*sum) break;
		}
		if (k > MAXIT) throw("Series failed in ei");
		return sum + log(x) + EULER;
	}
	else {
		sum = 0.0;
		term = 1.0;
		for (k = 1; k <= MAXIT; k++) {
			prev = term;
			term *= k / x;
			if (term < EPS) break;
			if (term < prev) sum += term;
			else {
				sum -= prev;
				break;
			}
		}
		return exp(x)*(1.0 + sum) / x;
	}
}

#ifdef __cplusplus
}
#endif