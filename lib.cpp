
#include "lib.h"
//#include "specfunc.h"

/*
#include "linear.h"
#include "roots.h"
#include "polys.h"
#include "specfun.h"
#include "reduce.h"
#include "integrate.h"

//#include "comutil.h"
#include <comdef.h>
*/


HINSTANCE hDLL;               // Handle to DLL
lookupFcn uni_array[MAX_UNI];
lookupFcn bi_array[MAX_BIN];
lookupFcn red_array[MAX_RED];

BOOL APIENTRY DllMain(HMODULE hModule,	DWORD  ul_reason_for_call,	LPVOID lpReserved) {
	hDLL = NULL;
	switch (ul_reason_for_call)	{
	case DLL_PROCESS_ATTACH:
	case DLL_THREAD_ATTACH:
	case DLL_THREAD_DETACH:
		break;
	case DLL_PROCESS_DETACH:
		if (hDLL)
			FreeLibrary(hDLL);
	}
	hDLL = LoadLibrary(L"Lib.xll");
	return TRUE;
}


#define MAX_XLOPER_STR_LEN 1024 /* max len to use when converting xloper to string */
#define MAX_XL12_STR_LEN 32767u

static XLOPER xDLL;
static XLOPER args[14];
static XLOPER MISSING;

LPXLOPER ArgNum(int i, double x) {
	args[i].val.num = x;
	args[i].xltype = xltypeNum;
	return &args[i];
}


XLOPER ArgStr(const char* string) {
	XLOPER res;
	size_t len = strlen(string);
	if (len > 255) len = 255; // Excel strings are limited to 255 chars
	char * temp = (char *)malloc(len + 2);
	memcpy(temp + 1, string, len);
	temp[0] = (BYTE)len;
	temp[len + 1] = 0;
	res.val.str = temp;
	res.xltype = xltypeStr;
	return res;
}


wchar_t *new_xl12string(const char *text)
{
	size_t len;
	if (!text || !(len = strlen(text)))
		return NULL;
	if (len > MAX_XL12_STR_LEN)
		len = MAX_XL12_STR_LEN; // truncate
	wchar_t *p = (wchar_t *)malloc((len + 2) * sizeof(wchar_t));
	if (!p) return NULL;
	mbstowcs(p + 1, text, len);
	p[0] = len; // string p[1] is NOT null terminated
	p[len + 1] = 0; // now it is
	return p;
}


XLOPER12 ArgStr12(char * string) {
	XLOPER12 res;
	res.val.str = new_xl12string(string);
	res.xltype = xltypeStr;
	return res;
}


LPXLOPER Arg(int i, const char* string) {
	if (!args[i].val.str == NULL) free(args[i].val.str);
	if (string == NULL)
		args[i] = MISSING;
	else 
		args[i] = ArgStr(string);
	return &args[i];
}

void setDLL() {
	Excel4(xlGetName, &xDLL, 0);
	MISSING.val.str = NULL;
	MISSING.xltype = xltypeMissing;
	for (int i = 0; i < 14; i++)
		args[i] = MISSING;
	Arg(4, "1");
	Arg(6, NULL);
	Arg(7, NULL);
}

void setCategory(const char* category) {
	Arg(5, category);
}

int RegFcn(const char* procedure, const char* typeText, const char* argumentText, const char* functionHelp, const char* argumentHelp1,
	const char* argumentHelp2, const char* argumentHelp3, const char* argumentHelp4, const char* argumentHelp5) {
	Arg(0, procedure);
	int res = Excel4(xlfRegister, 0, 15, &xDLL, &args[0],
		Arg(1, typeText), &args[0], Arg(3, argumentText),
		&args[4], &args[5], &args[6], &args[7], Arg(8, functionHelp), Arg(9, argumentHelp1),
		Arg(10, argumentHelp2), Arg(11, argumentHelp3), Arg(12, argumentHelp4), Arg(13, argumentHelp5));
	if (res != 0) {
		//sprintf("Failed to register %s\n", procedure);
	}
	return res;
}

void RegFcnV(std::vector<char*> x){
	x.resize(9, NULL);
	RegFcn(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8]);
}


#ifdef __cplusplus
extern "C" {
#endif 

__declspec(dllexport) int WINAPI xlAutoClose(void)	{ return 1; }

__declspec(dllexport) LPXLOPER WINAPI xlAutoRegister(LPXLOPER pxName) {
	static XLOPER xDLL, xRegId;
	xRegId.xltype = xltypeErr;
	xRegId.val.err = xlerrValue;
	return (LPXLOPER)&xRegId;
}

__declspec(dllexport) int WINAPI xlAutoAdd(void) { return 1; }
__declspec(dllexport) int WINAPI xlAutoRemove(void) { return 1; }
__declspec(dllexport) void WINAPI xlAutoFree(LPXLOPER px) {}

__declspec(dllexport) LPXLOPER WINAPI xlAddInManagerInfo(LPXLOPER xAction) {
	static XLOPER xInfo, xIntAction, xIntType;
	xIntType.xltype = xltypeInt;
	xIntType.val.w = xltypeInt;
	xInfo.xltype = xltypeErr;
	xInfo.val.err = xlerrValue;
	return (LPXLOPER)&xInfo;
}



__declspec(dllexport) int WINAPI xlAutoOpen(void) {
	XLCallVer();
	setDLL();
	static XLOPER xDLL;
	Excel4(xlGetName, &xDLL, 0);
	registerMatrixFunc();
	/*
	registerLinear();
	registerRoots();
	registerPolys();
	registerReduce();
	registerMatrices();
	registerXlmatrix();
	registerIntegration();*/
	return 1;
}


// **************************************************************************************************************************

int look(int n, lookupFcn fcnarr[], char * c) {
	int x = -1;
	for (int i = 0; i < n; i++)
		if (c == fcnarr[i].name) {
			x = i; break;
		}
	return x;
}

FP12 * allocateFP(int m, int n, double val = 0)  {  // for internal use. Pointer after use must be FREEd
	FP12 * x; int size = m * n;
	if (size <= 0)
		return NULL;
	if ((x = (FP12 *)calloc(size + 1, sizeof(double)))) {
		x->rows = m;
		x->cols = n;
	}
	if (val != 0.0)
		std::fill_n((double *)&(x->array), m*n, val);
	return x;
}



extern "C++" FP12 * Dim(int m, int n, double val, char assign) {
	static FP12 * p_array = NULL; // Not thread-safe to use static
	if (assign == TOXL) {
		if (p_array) // free memory allocated on last call
		{
			free(p_array);
			p_array = NULL;
		}
		p_array = allocateFP(m, n, val);
		return p_array;
	} else {				// you must free it
		FP12 * p = allocateFP(m, n, val);
		return p;
	}
}


extern "C" xll_FP12 mCopy(FP12 & A, char assign) {
	FP12 * dest = Dim(A.rows, A.cols, 0.0, assign);
	memcpy(dest->array, A.array, A.rows * A.cols * sizeof(double));
	return dest;
}





xll_LPSAFEARRAY size(int m, int n) {
	SAFEARRAYBOUND bounds[2];
	bounds[0].cElements = m;
	bounds[0].lLbound = 1;
	bounds[1].cElements = n;
	bounds[1].lLbound = 1;
	LPSAFEARRAY x = SafeArrayCreate(VT_R8, 2, bounds);
	return x;
}


VARIANT getdata(VARIANT & ExcelArray) {
	VARIANT dvout;	EXCEPINFO excep; DISPPARAMS dispparams;	unsigned int uiArgErr; DISPID dispidValue;
	LPOLESTR XName = L"Value2";
	ExcelArray.pdispVal->GetIDsOfNames(IID_NULL, &XName, 1, LOCALE_SYSTEM_DEFAULT, &dispidValue);
	dispparams.cArgs = 0;
	dispparams.cNamedArgs = 0;
	ExcelArray.pdispVal->Invoke(dispidValue, IID_NULL, LOCALE_SYSTEM_DEFAULT, DISPATCH_PROPERTYGET, &dispparams, &dvout, &excep, &uiArgErr);
	ExcelArray.pdispVal->Release();
	return dvout;
}

inline LPSAFEARRAY getarray (VARIANT const &var) {
	// Extract the array (differs if it came in by reference)
	if (var.vt & VT_BYREF)
		return * var.pparray;
	else
		return var.parray;
}

xll_LPSAFEARRAY mat(VARIANT & src) {
	if (src.vt == VT_DISPATCH)		// if range, convert it
		src = getdata(src);
	if (src.vt == VT_R8) {		// a single value
		vbMatrix<double>  A(1, 1);
		A[1] = src.dblVal;
		return A;
	}
	else if ((src.vt & VT_VARIANT) != 0) {
		vbMatrix<VARIANT> v(&src.parray);
		vbMatrix<double>  A(v.rows, v.cols);
		for (int j = 1; j <= v.cols; j++)
			for (int i = 1; i <= v.rows; i++)
				A(i, j) = v(i, j).dblVal;
		return A;
	}
	else
		return getarray(src); 
}


xll_LPSAFEARRAY matcopy(LPSAFEARRAY * a) {
	vbMatrix<double>  A(a), dest(A.rows, A.cols);
	memcpy(dest.data, A.data, A.rows * A.cols * sizeof(double));
	return dest;
}

extern "C" xll_FP12 mConform(FP12 & A, int idir, char assign) {
	FP12 * dest;
	if (idir == WHOLE)
		dest = Dim(1, 1, assign);
	else if (idir == BYROW)
		dest = Dim(1, A.cols, assign);
	else
		dest = Dim(A.rows, 1, assign);
	return dest;
}

extern "C" xll_LPSAFEARRAY vConform(LPSAFEARRAY * a, int idir) {
	vbMatrix<double>  A(a);
	LPSAFEARRAY dest;
	if (idir == WHOLE)
		dest = size(1, 1);
	else if (idir == BYROW)
		dest = size(1, A.cols);
	else
		dest = size(A.rows, 1);
	return dest;
}

void MATRIX::reduce(MATRIX & ret, VecReduceFcn func, int idir) {
	if (idir == WHOLE)
		ret(0, 0) =	(all().*func)();
	else if (idir == BYROW)
		for (int i = 1; i <= cols; i++)
			ret(0, i - 1) = (col(i).*func)();
	else
		for (int i = 1; i <= rows; i++)
			ret(i - 1,0) = (row(i).*func)();
} 


#ifdef __cplusplus
}
#endif


