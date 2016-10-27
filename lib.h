#ifndef FP_INCLUDE
#define FP_INCLUDE

#include <Windows.h>
#include <atlsafe.h>
#include <math.h>
#include "xlcall.h"
#include <vector>
#include <string>
#include <iostream>
//#include <Eigen/Dense>

#include "specfunc.h"

using namespace std;
//using namespace Eigen;


LPXLOPER Arg(int i, const char* string);
LPXLOPER ArgNum(int i, double x);
void setCategory(const char* category);
int RegFcn(const char* procedure, const char* typeText, const char* argumentText, const char* functionHelp, const char* argumentHelp1 = NULL,
	const char* argumentHelp2 = NULL, const char* argumentHelp3 = NULL, const char* argumentHelp4 = NULL, const char* argumentHelp5 = NULL);
void RegFcnV(std::vector<char*> x);



inline double abs_diff(double x, double y) { return abs(x - y); }
inline double sqr_diff(double x, double y) { return (x - y)*(x - y); }


#define TOXL  'X'
#define INTERNAL  'I'


//#define UNASSIGNED 0
#define WHOLE 0
#define BYROW 1
#define BYCOL 2

#define xll_void		__declspec(dllexport) void          WINAPI 
#define xll_double		__declspec(dllexport) double        WINAPI 
#define xll_FP12		__declspec(dllexport) FP12 *        WINAPI 
#define xll_LPSAFEARRAY __declspec(dllexport) LPSAFEARRAY   WINAPI 

#define Q(x) #x
#define QUOTE(x) Q(x)


extern "C" xll_LPSAFEARRAY mat(VARIANT & src);
extern "C" xll_LPSAFEARRAY matcopy(LPSAFEARRAY * a);
extern "C" xll_LPSAFEARRAY size(int m, int n);

extern "C++" FP12 * Dim(int m, int n, double val=0.0, char assign=TOXL);
extern "C" xll_FP12 mCopy   (FP12 & A, char assign=TOXL);						// copies A
extern "C" xll_FP12 mConform(FP12 & A, int idir = WHOLE, char assign = TOXL);
extern "C" xll_LPSAFEARRAY vConform(LPSAFEARRAY * a, int idir);


extern "C++"
template<typename T>
class vbMatrix : public CComSafeArray<T> {
	void setvals() {
		rows = m_psa->rgsabound[1].cElements;
		cols = m_psa->rgsabound[0].cElements;
		data = static_cast<T *> (m_psa->pvData);
	}
public:
	T * data;
	UINT rows, cols;
	vbMatrix() : CComSafeArray<T>() {}
	vbMatrix(UINT a_rows, UINT a_cols) {
		size(a_rows, a_cols);
	}
	vbMatrix(LPSAFEARRAY * x)  {
		CComSafeArray<T>::Attach(*x);
		setvals();
	}
	void size(UINT a_rows, UINT a_cols) {
		LPSAFEARRAY x = ::size(a_rows, a_cols);
		CComSafeArray<T>::Attach(x);
		setvals();
	}
	inline T & operator[] (int i)  { return data[i - 1]; }											// 1 based
	~vbMatrix()  { Detach(); }
	inline T & operator()(UINT a_row, UINT a_col)	{ return data[rows * (a_col - 1) + (a_row - 1)]; }		// 1..rows, 1..cols
	inline int Els()                { return rows*cols; }
};

struct MATRIX;
struct VECTOR;
typedef MATRIX & (MATRIX::*MonFcn) ();
typedef MATRIX & (MATRIX::*DyadFcn) (MATRIX & B);
typedef double (VECTOR::*VecReduceFcn) ();


#ifdef __cplusplus
extern "C" {
#endif

enum UNARY_FCNS {			// define symbol values, for use in vba
#define CODE_DEF(x,y,z,t) x,
#include "codes.def"
MAX_UNI
};


enum BINARY_FCNS {			// define symbol values, for use in vba
#define CODE_DEF(x,y,z,t) x,
#include "bincodes.def"
MAX_BIN
};

enum REDUCE_FCNS {			// define symbol values, for use in vba
#define CODE_DEF(x,y,z,t) x,
#include "red_codes.def"
MAX_RED
};



struct lookupFcn {
	string name;
	union {
		MonFcn  fcn;
		DyadFcn fcn2;
		VecReduceFcn fcnred;
	};
	lookupFcn() {}
	lookupFcn(string nam, MonFcn f)	      { name = nam; fcn = f; }
	lookupFcn(string nam, DyadFcn f)      { name = nam; fcn2 = f; }
	lookupFcn(string nam, VecReduceFcn f) { name = nam; fcnred = f; }
};

extern lookupFcn uni_array[MAX_UNI];
extern lookupFcn bi_array[MAX_BIN];
extern lookupFcn red_array[MAX_RED];

int look(int n, lookupFcn fcnarr[], char * c);
inline int look_uni(char * c) { return look(MAX_UNI, uni_array, c); }
inline int look_bi (char * c) { return look(MAX_BIN, bi_array, c); }
inline int look_red(char * c) { return look(MAX_RED, red_array, c); }


inline int Els(FP12 & A) { return A.rows* A.cols;  }


#define SUM __sec_reduce_add

//#define MapFP12(A) Ref<MatrixXd> (Map<MatrixXd,RowMajor> (&(A).array[0],(A).cols,(A).rows))				// A = FP12 &
//#define MapSA(A)   Ref<MatrixXd> (Map<MatrixXd,ColMajor> ((double*) (A)->pvData, ROWS(A),COLS(A)))		// A = LPSAFEARRAY &

typedef double(*MonadicScalarFcn) (double x);
#define VEC(x)  (x).data[0:(x).n:(x).stride] 
#define SELF    VEC(*this)

struct VECTOR {
	int n;
	int stride;
	double * data;
	VECTOR & refMonadic(MonadicScalarFcn f) { data[0:n : stride] = f(data[0:n : stride]); return *this; }
	#define CODE_DEF(symb,fcn,op,desc) 	VECTOR & fcn() { return refMonadic(op); }
	#include "codes.def"
	#define CODE_DEF(symb,fcn,op,desc) VECTOR & fcn(VECTOR & B) BIN_OP(op)
	#include "bincodes.def"
	double Amin()     { return __sec_reduce_min(fabs(SELF)); }
	double Sum()      { return SUM(SELF); }
	double Amax()     { return __sec_reduce_max(fabs(SELF)); }
	double Sumsqr()   { return SUM(SELF * SELF); }
	double Asum()     { return SUM(fabs(SELF)); }
	double Nrm2()     { return sqrt(SUM(SELF*SELF)); }
	double Max()      { return __sec_reduce_max(SELF); }
	double Min()      { return __sec_reduce_min(SELF); }
	double iAmin()    { return (__sec_reduce_min_ind(fabs(SELF)) / stride + 1); }
	double iAmax()    { return (__sec_reduce_max_ind(fabs(SELF)) / stride + 1); }
	double iMin()     { return (__sec_reduce_min_ind(SELF) / stride + 1); }
	double iMax()     { return (__sec_reduce_max_ind(SELF) / stride + 1); }
	double Mean()     { return Sum() / n; }
	double Prod()     { return __sec_reduce_mul(SELF); }
	double AbsProd()  { return __sec_reduce_mul(fabs(SELF)); }
	double All()      { return __sec_reduce_all_nonzero(SELF); }
	double Any()      { return __sec_reduce_any_nonzero(SELF); }
	double Count()    { return SUM(fabs(SELF)>0); }
	double AbsDev()   { double mean = Mean(); return SUM(fabs(SELF - mean)) / n; }
	double Variance() { double mean = Mean(); return SUM(sqr(SELF - mean)) / n; }
};

double dot(VECTOR & x, VECTOR & y) { return SUM(VEC(x) * VEC(y)); }

#define COL(A,i) (A).data[(i-1)*(A).colstride : (A).rows : (A).rowstride]				// COL(A,1)
#define ROW(A,i) (A).data[(i-1)*(A).rowstride : (A).cols : (A).colstride]				// ROW(A,1)
#define ALL(A)   (A).data[0 : (A).rows * (A).cols]										// WHOLE MATRIX

#define ROWS(SA) SA->rgsabound[1].cElements
#define COLS(SA) SA->rgsabound[0].cElements
#define DATA(SA) (double*) SA->pvData

struct MATRIX {
	int rows;
	int cols;
	int rowstride;
	int colstride;
	double * data;
	inline double & operator () (int i, int j)	{ return data[i*rowstride + j*colstride]; }		// A(0,0) to speed up
	MATRIX(FP12 & A)        { rows = A.rows;  cols = A.cols;  data = A.array; rowstride = cols;	colstride = 1;	}
	MATRIX(LPSAFEARRAY & x) { rows = ROWS(x); cols = COLS(x); data = DATA(x); rowstride = 1;	colstride = rows; }
	int els() { return rows * cols; }
	bool isRowMajor() { return rowstride > 1; }
	VECTOR col(int i) { return {rows, rowstride, &data[(i - 1)*colstride] }; }		// Col(1..)
	VECTOR row(int i) { return {cols, colstride, &data[(i - 1)*rowstride] }; }		// Row(1..)
	VECTOR all()      { return {els(), 1, data}; }	
	void reduce (MATRIX & ret, VecReduceFcn func, int idir=WHOLE);
	#define CODE_DEF(symb,fcn,op,desc) 	MATRIX & fcn() { all().fcn(); return *this; }
	#include "codes.def"
	#define CODE_DEF(symb,fcn,op,desc) MATRIX & fcn(MATRIX & B) { all().fcn(B.all()); return *this; }
	#include "bincodes.def"
	#define CODE_DEF(symb,fcn,op,desc) double fcn(MATRIX & B) { return all().fcn(); }
	#include "red_codes.def"
	MATRIX & Transp(MATRIX & B) { for (int i = 1; i <= rows; i++) ROW(*this,i) = COL(B,i); return *this; }
};

#define CODE_DEF(symb,fcn,op,desc) 	inline xll_void m##fcn(FP12 & A) { MATRIX(A).fcn(); } 
#include "codes.def"
#define CODE_DEF(symb,fcn,op,desc) inline xll_void m##fcn(FP12 & A, FP12 & B) { MATRIX(A).fcn(MATRIX(B)); }
#include "bincodes.def"

inline xll_void mMonadic(FP12 & A, char * fcn) {int icode = look_uni(fcn); if (icode >= 0) (MATRIX(A).*(uni_array[icode].fcn)) (); }
inline xll_void mDyadic (FP12 & A, FP12 & B, char * fcn) { 
	int icode = look_bi(fcn);  if (icode >= 0) {
		DyadFcn f = (DyadFcn)bi_array[icode].fcn2;
		(MATRIX(A).*f) (MATRIX(B));
	}
}
inline xll_FP12 doReduce (FP12 & A, FP12 * ret, int icode, int idir = WHOLE) {
	if (icode >= 0)
		MATRIX(A).reduce(MATRIX(*ret), red_array[icode].fcnred, idir);
	return ret;
}
#define CODE_DEF(symb,fcn,op,desc) inline xll_FP12 m##fcn(FP12 & A, int idir=WHOLE) { return doReduce(A,mConform(A, idir),symb,idir); }
#include "red_codes.def"
inline xll_FP12 mReduce(FP12 & A, char * fcn, int idir = WHOLE) { return doReduce(A, mConform(A, idir), look_red(fcn), idir); }
inline xll_LPSAFEARRAY vReduce(LPSAFEARRAY * A, REDUCE_FCNS fcn, int idir = WHOLE) { 
	LPSAFEARRAY ret = vConform(A, idir);
	MATRIX(*A).reduce(MATRIX(ret), red_array[fcn].fcnred, idir);
	return ret;
}


inline xll_void vrefMonadic(LPSAFEARRAY & A, UNARY_FCNS icode) { (MATRIX(A).*(uni_array[icode].fcn)) (); }
inline xll_LPSAFEARRAY vMonadic(LPSAFEARRAY & A, UNARY_FCNS icode) { LPSAFEARRAY x = matcopy(&A); vrefMonadic(x, icode); return x; }
inline xll_void vrefDyadic(LPSAFEARRAY & A, LPSAFEARRAY & B, BINARY_FCNS icode)   { DyadFcn f = (DyadFcn)bi_array[icode].fcn2; (MATRIX(A).*f) (MATRIX(B)); }
inline xll_LPSAFEARRAY vDyadic(LPSAFEARRAY & A, LPSAFEARRAY & B, BINARY_FCNS icode) { LPSAFEARRAY x = matcopy(&A); vrefDyadic(x, B, icode); return x; }


inline xll_FP12 mTransp(FP12 & A) { FP12 * B = Dim(A.cols, A.rows); MATRIX(*B).Transp(MATRIX(A)); return B; }
inline xll_LPSAFEARRAY vTransp(LPSAFEARRAY & A) { LPSAFEARRAY B = size(COLS(A), ROWS(A)); MATRIX(B).Transp(MATRIX(A)); return B; }
inline xll_FP12 mCol(FP12 & A, int i) { FP12 *ret = mConform(A, BYCOL); ALL(MATRIX(*ret)) = COL(MATRIX(A), i); return ret; }
inline xll_FP12 mRow(FP12 & A, int i) { FP12 *ret = mConform(A, BYROW); ALL(MATRIX(*ret)) = ROW(MATRIX(A), i); return ret; }


inline void doAB(MATRIX & A, MATRIX & B, MATRIX & C) {		// A(m,in)   B(in,n)		rows by cols
	for (int i = 1; i <= C.rows; i++)
		for (int j = 1; j <= C.cols; j++)
			C(i - 1, j - 1) = SUM(ROW(A, i) * COL(B, j));
}
inline xll_FP12 mAB(FP12 & A, FP12 & B) { FP12 * C = Dim(A.rows, B.cols); doAB(MATRIX(A), MATRIX(B), MATRIX(*C)); return C; }
inline xll_LPSAFEARRAY vAB(LPSAFEARRAY & A, LPSAFEARRAY & B) { LPSAFEARRAY C = size(ROWS(A), COLS(B)); doAB(MATRIX(A), MATRIX(B), MATRIX(C)); return C; }


inline void registerMatrixFunc() {
	setCategory("Matrix Functions");
	std::vector<std::vector<char*>> x = {
		{ "mMonadic", "1K%C", "A,fcn", "generic monadic", "matrix A", "function " },
		{ "mDyadic", "1K%K%C", "A,B,fcn", "generic dyadic", "matrix A", "matrix B", "function " },
		{ "mReduce", "K%K%CJ", "A,fcn,idir", "generic reduction", "matrix A", "function", "direction " },
		{ "mTransp", "K%K%", "A", "transpose", "matrix A" },
		{ "mCopy", "K%K%", "A", "copy", "matrix A" },
		{ "mCol", "K%K%J", "A,idir", "copy", "matrix A" },
		{ "mRow", "K%K%J", "A,idir", "copy", "matrix A" },
		{ "mAB", "K%K%K%", "A,B", "copy", "matrix A", "matrix B " },
	};
	for (auto row : x) RegFcnV(row);

	setCategory("Monadic Functions");
	#define CODE_DEF(symb,f,addr,z)  uni_array[symb] = lookupFcn(#symb, (MonFcn) &MATRIX::f);
	#include "codes.def"
	x = {
		#define CODE_DEF(x,FN,addr,DESC)  {QUOTE(m##FN), "1K%", "A", #DESC, "matrix A" },
		#include "codes.def"
	};
	for (auto row : x)
		RegFcnV(row);

	setCategory("Matrix Binary Functions");
	#define CODE_DEF(symb,name,op,desc)  bi_array[symb] = lookupFcn(#symb, (DyadFcn) &MATRIX::name);
	#include "bincodes.def"
	x = {
		#define CODE_DEF(symb,name,op,desc)  { QUOTE(m##name), "1K%K%", "A,B", #desc, "matrix A", "matrix B" },
		#include "bincodes.def"
		{ "mLo", "1K%K%", "A,B", "clip low", "matrix A", "matrix B" },
		{ "mHi", "1K%K%", "A,B", "clip high", "matrix A", "matrix B" },
		{ "mSqrDiff", "1K%K%", "A,B", "squared difference", "matrix A", "matrix B" },
	};
	for (auto row : x)	RegFcnV(row);

	setCategory("Matrix Reduction Functions");
	#define CODE_DEF(symb,name,op,desc)  red_array[symb] = lookupFcn(#symb, (VecReduceFcn) &VECTOR::name);
	#include "red_codes.def"
	x = {
		#define CODE_DEF(symbol,name,func,desc) { QUOTE(m##name), "K%K%J", "A,idir", #desc, "matrix A", "direction "},
		#include "red_codes.def"
	};
	for (auto row : x)	RegFcnV(row);
}


#ifdef __cplusplus
}
#endif


#endif	//FP_INCLUDE

