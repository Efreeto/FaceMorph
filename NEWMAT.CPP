//$$ newmat1.cpp   Matrix type functions

// Copyright (C) 1991,2,3,4: R B Davies

//#define WANT_STREAM

#include "newmat.h"

#ifdef use_namespace
namespace NEWMAT {
#endif

#ifdef DO_REPORT
#define REPORT { static ExeCounter ExeCount(__LINE__,1); ++ExeCount; }
#else
#define REPORT {}
#endif

#include <stdio.h>
#include <iostream>


	//$$ newmat8.cpp         Advanced LU transform, scalar functions

	// Copyright (C) 1991,2,3,4,8: R B Davies

#define WANT_MATH

//#include "include.h"
#include <string.h>
#include "newmat.h"

	//$$ newmatrc.h              definition file for row/column classes

	// Copyright (C) 1991,2,3,4,7: R B Davies

#ifndef NEWMATRC_LIB
#define NEWMATRC_LIB 0

#ifdef use_namespace
	namespace NEWMAT {
#endif

#include "controlw.h"


		/************** classes MatrixRowCol, MatrixRow, MatrixCol *****************/

		// Used for accessing the rows and columns of matrices
		// All matrix classes must provide routines for calculating matrix rows and
		// columns. Assume rows can be found very efficiently.

		enum LSF { LoadOnEntry=1,StoreOnExit=2,DirectPart=4,StoreHere=8,HaveStore=16 };


		class LoadAndStoreFlag : public ControlWord
		{
		public:
			LoadAndStoreFlag() {}
			LoadAndStoreFlag(int i) : ControlWord(i) {}
			LoadAndStoreFlag(LSF lsf) : ControlWord(lsf) {}
			LoadAndStoreFlag(const ControlWord& cwx) : ControlWord(cwx) {}
		};

		class MatrixRowCol
			// the row or column of a matrix
		{
		public:                                        // these are public to avoid
			// numerous friend statements
			int length;                                 // row or column length
			int skip;                                   // initial number of zeros
			int storage;                                // number of stored elements
			int rowcol;                                 // row or column number
			GeneralMatrix* gm;                          // pointer to parent matrix
			Real* data;                                 // pointer to local storage
			LoadAndStoreFlag cw;                        // Load? Store? Is a Copy?
			void IncrMat() { rowcol++; data += storage; }   // used by NextRow
			void IncrDiag() { rowcol++; skip++; data++; }
			void IncrId() { rowcol++; skip++; }
			void IncrUT() { rowcol++; data += storage; storage--; skip++; }
			void IncrLT() { rowcol++; data += storage; storage++; }

		public:
			void Zero();                                // set elements to zero
			void Add(const MatrixRowCol&);              // add a row/col
			void AddScaled(const MatrixRowCol&, Real);  // add a multiple of a row/col
			void Add(const MatrixRowCol&, const MatrixRowCol&);
			// add two rows/cols
			void Add(const MatrixRowCol&, Real);        // add a row/col
			void NegAdd(const MatrixRowCol&, Real);     // Real - a row/col
			void Sub(const MatrixRowCol&);              // subtract a row/col
			void Sub(const MatrixRowCol&, const MatrixRowCol&);
			// sub a row/col from another
			void RevSub(const MatrixRowCol&);           // subtract from a row/col
			void ConCat(const MatrixRowCol&, const MatrixRowCol&);
			// concatenate two row/cols
			void Multiply(const MatrixRowCol&);         // multiply a row/col
			void Multiply(const MatrixRowCol&, const MatrixRowCol&);
			// multiply two row/cols
			void KP(const MatrixRowCol&, const MatrixRowCol&);
			// Kronecker Product two row/cols
			void Copy(const MatrixRowCol&);             // copy a row/col
			void CopyCheck(const MatrixRowCol&);        // ... check for data loss
			void Check(const MatrixRowCol&);            // just check for data loss
			void Check();                               // check full row/col present
			void Copy(const double*&);                  // copy from an array
			void Copy(const float*&);                   // copy from an array
			void Copy(const int*&);                     // copy from an array
			void Copy(Real);                            // copy from constant
			void Add(Real);                             // add a constant
			void Multiply(Real);                        // multiply by constant
			Real SumAbsoluteValue();                    // sum of absolute values
			Real MaximumAbsoluteValue1(Real r, int& i); // maximum of absolute values
			Real MinimumAbsoluteValue1(Real r, int& i); // minimum of absolute values
			Real Maximum1(Real r, int& i);              // maximum
			Real Minimum1(Real r, int& i);              // minimum
			Real Sum();                                 // sum of values
			void Inject(const MatrixRowCol&);           // copy stored els of a row/col
			void Negate(const MatrixRowCol&);           // change sign of a row/col
			void Multiply(const MatrixRowCol&, Real);   // scale a row/col
			friend Real DotProd(const MatrixRowCol&, const MatrixRowCol&);
			// sum of pairwise product
			Real* Data() { return data; }
			int Skip() { return skip; }                 // number of elements skipped
			int Storage() { return storage; }           // number of elements stored
			int Length() { return length; }             // length of row or column
			void Skip(int i) { skip=i; }
			void Storage(int i) { storage=i; }
			void Length(int i) { length=i; }
			void SubRowCol(MatrixRowCol&, int, int) const;
			// get part of a row or column
			MatrixRowCol() {}                           // to stop warning messages
			~MatrixRowCol();
			FREE_CHECK(MatrixRowCol)
		};

		class MatrixRow : public MatrixRowCol
		{
		public:
			// bodies for these are inline at the end of this .h file
			MatrixRow(GeneralMatrix*, LoadAndStoreFlag, int=0);
			// extract a row
			~MatrixRow();
			void Next();                                // get next row
			FREE_CHECK(MatrixRow)
		};

		class MatrixCol : public MatrixRowCol
		{
		public:
			// bodies for these are inline at the end of this .h file
			MatrixCol(GeneralMatrix*, LoadAndStoreFlag, int=0);
			// extract a col
			MatrixCol(GeneralMatrix*, Real*, LoadAndStoreFlag, int=0);
			// store/retrieve a col
			~MatrixCol();
			void Next();                                // get next row
			FREE_CHECK(MatrixCol)
		};

		// MatrixColX is an alternative to MatrixCol where the complete
		// column is stored externally

		class MatrixColX : public MatrixRowCol
		{
		public:
			// bodies for these are inline at the end of this .h file
			MatrixColX(GeneralMatrix*, Real*, LoadAndStoreFlag, int=0);
			// store/retrieve a col
			~MatrixColX();
			void Next();                                // get next row
			Real* store;                                // pointer to local storage
			//    less skip
			FREE_CHECK(MatrixColX)
		};

		/**************************** inline bodies ****************************/

		inline MatrixRow::MatrixRow(GeneralMatrix* gmx, LoadAndStoreFlag cwx, int row)
		{ gm=gmx; cw=cwx; rowcol=row; gm->GetRow(*this); }

		inline void MatrixRow::Next() { gm->NextRow(*this); }

		inline MatrixCol::MatrixCol(GeneralMatrix* gmx, LoadAndStoreFlag cwx, int col)
		{ gm=gmx; cw=cwx; rowcol=col; gm->GetCol(*this); }

		inline MatrixCol::MatrixCol(GeneralMatrix* gmx, Real* r,
			LoadAndStoreFlag cwx, int col)
		{ gm=gmx; data=r; cw=cwx+StoreHere; rowcol=col; gm->GetCol(*this); }

		inline MatrixColX::MatrixColX(GeneralMatrix* gmx, Real* r,
			LoadAndStoreFlag cwx, int col)
		{ gm=gmx; store=data=r; cw=cwx+StoreHere; rowcol=col; gm->GetCol(*this); }


		inline void MatrixCol::Next() { gm->NextCol(*this); }

		inline void MatrixColX::Next() { gm->NextCol(*this); }

#ifdef use_namespace
	}
#endif

#endif

	//$$ precisio.h                          floating point constants

#ifndef PRECISION_LIB
#define PRECISION_LIB 0

#define WANT_MATH
//#include "include.h"              // in case being used as stand alone

#ifdef _STANDARD_                 // standard library available
#include <limits>
#endif

#ifdef use_namespace
	namespace NEWMAT {
#endif

#ifdef _STANDARD_                 // standard library available

#ifdef OPT_COMPATIBLE
#include <cfloat>                 // for FLT_MAX
#endif

		using namespace std;

		class FloatingPointPrecision
		{
		public:
			static int Dig()              // number of decimal digits or precision
			{ return numeric_limits<Real>::digits10 ; }

			static Real Epsilon()         // smallest number such that 1+Eps!=Eps
			{ return numeric_limits<Real>::epsilon(); }

			static int Mantissa()         // bits in mantisa
			{ return numeric_limits<Real>::digits; }

			static Real Maximum()         // maximum value
			{ return numeric_limits<Real>::max(); }

			static int MaximumDecimalExponent()  // maximum decimal exponent
			{ return numeric_limits<Real>::max_exponent10; }

			static int MaximumExponent()  // maximum binary exponent
			{ return numeric_limits<Real>::max_exponent; }

			static Real LnMaximum()       // natural log of maximum
			{ return (Real)log(Maximum()); }

			static Real Minimum()         // minimum positive value
			{ return numeric_limits<Real>::min(); } 

			static int MinimumDecimalExponent() // minimum decimal exponent
			{ return numeric_limits<Real>::min_exponent10; }

			static int MinimumExponent()  // minimum binary exponent
			{ return numeric_limits<Real>::min_exponent; }

			static Real LnMinimum()       // natural log of minimum
			{ return (Real)log(Minimum()); }

			static int Radix()            // exponent radix
			{ return numeric_limits<Real>::radix; }

			static int Rounds()           // addition rounding (1 = does round)
			{
				return numeric_limits<Real>::round_style ==
					round_to_nearest ? 1 : 0;
			}

		};


#else                              // _STANDARD_ not defined

#ifndef SystemV                    // if there is float.h

#ifdef USING_FLOAT

		class FloatingPointPrecision
		{
		public:
			static int Dig()
			{ return FLT_DIG; }        // number of decimal digits or precision

			static Real Epsilon()
			{ return FLT_EPSILON; }    // smallest number such that 1+Eps!=Eps

			static int Mantissa()
			{ return FLT_MANT_DIG; }   // bits in mantisa

			static Real Maximum()
			{ return FLT_MAX; }        // maximum value

			static int MaximumDecimalExponent()
			{ return FLT_MAX_10_EXP; } // maximum decimal exponent

			static int MaximumExponent()
			{ return FLT_MAX_EXP; }    // maximum binary exponent

			static Real LnMaximum()
			{ return (Real)log(Maximum()); } // natural log of maximum

			static Real Minimum()
			{ return FLT_MIN; }        // minimum positive value

			static int MinimumDecimalExponent()
			{ return FLT_MIN_10_EXP; } // minimum decimal exponent

			static int MinimumExponent()
			{ return FLT_MIN_EXP; }    // minimum binary exponent

			static Real LnMinimum()
			{ return (Real)log(Minimum()); } // natural log of minimum

			static int Radix()
			{ return FLT_RADIX; }      // exponent radix

			static int Rounds()
			{ return FLT_ROUNDS; }     // addition rounding (1 = does round)

		};

#endif                           // USING_FLOAT


#ifdef USING_DOUBLE

		class FloatingPointPrecision
		{
		public:

			static int Dig()
			{ return DBL_DIG; }        // number of decimal digits or precision

			static Real Epsilon()
			{ return DBL_EPSILON; }    // smallest number such that 1+Eps!=Eps

			static int Mantissa()
			{ return DBL_MANT_DIG; }   // bits in mantisa

			static Real Maximum()
			{ return DBL_MAX; }        // maximum value

			static int MaximumDecimalExponent()
			{ return DBL_MAX_10_EXP; } // maximum decimal exponent

			static int MaximumExponent()
			{ return DBL_MAX_EXP; }    // maximum binary exponent

			static Real LnMaximum()
			{ return (Real)log(Maximum()); } // natural log of maximum

			static Real Minimum()
			{
				//#ifdef __BCPLUSPLUS__
				//       return 2.225074e-308;     // minimum positive value
				//#else
				return DBL_MIN;
				//#endif
			}

			static int MinimumDecimalExponent()
			{ return DBL_MIN_10_EXP; } // minimum decimal exponent

			static int MinimumExponent()
			{ return DBL_MIN_EXP; }    // minimum binary exponent

			static Real LnMinimum()
			{ return (Real)log(Minimum()); } // natural log of minimum


			static int Radix()
			{ return FLT_RADIX; }      // exponent radix

			static int Rounds()
			{ return FLT_ROUNDS; }     // addition rounding (1 = does round)

		};

#endif                             // USING_DOUBLE

#else                              // if there is no float.h

#ifdef OPT_COMPATIBLE
#define FLT_MAX MAXFLOAT
#endif


#ifdef USING_FLOAT

		class FloatingPointPrecision
		{
		public:

			static Real Epsilon()
			{ return pow(2.0,(int)(1-FSIGNIF)); }
			// smallest number such that 1+Eps!=Eps

			static Real Maximum()
			{ return MAXFLOAT; }            // maximum value

			static Real LnMaximum()
			{ return (Real)log(Maximum()); }  // natural log of maximum

			static Real Minimum()
			{ return MINFLOAT; }             // minimum positive value

			static Real LnMinimum()
			{ return (Real)log(Minimum()); }  // natural log of minimum

		};

#endif                                  // USING_FLOAT


#ifdef USING_DOUBLE

		class FloatingPointPrecision
		{
		public:

			static Real Epsilon()
			{ return pow(2.0,(int)(1-DSIGNIF)); }
			// smallest number such that 1+Eps!=Eps

			static Real Maximum()
			{ return MAXDOUBLE; }           // maximum value

			static Real LnMaximum()
			{ return LN_MAXDOUBLE; }        // natural log of maximum

			static Real Minimum()
			{ return MINDOUBLE; }

			static Real LnMinimum()
			{ return LN_MINDOUBLE; }        // natural log of minimum
		};

#endif                                // USING_DOUBLE

#endif                                // SystemV

#endif                                // _STANDARD_




#ifdef use_namespace
	}
#endif                                // use_namespace



#endif                                // PRECISION_LIB

	static int tristore(int n)                    // elements in triangular matrix
	{ return (n*(n+1))/2; }

#ifdef use_namespace
	namespace NEWMAT {
#endif


#ifdef DO_REPORT
#define REPORT { static ExeCounter ExeCount(__LINE__,8); ++ExeCount; }
#else
#define REPORT {}
#endif


		// ************************* general utilities *************************/




		/************************** LU transformation ****************************/

		void CroutMatrix::ludcmp()
			// LU decomposition from Golub & Van Loan, algorithm 3.4.1, (the "outer
			// product" version).
			// This replaces the code derived from Numerical Recipes in C in previous
			// versions of newmat and being row oriented runs much faster with large
			// matrices.
		{
			REPORT
				Tracer tr( "Crout(ludcmp)" ); sing = false;
			Real* akk = store;                    // runs down diagonal

			Real big = fabs(*akk); int mu = 0; Real* ai = akk; int k;

			for (k = 1; k < nrows_val; k++)
			{
				ai += nrows_val; const Real trybig = fabs(*ai);
				if (big < trybig) { big = trybig; mu = k; }
			}


			if (nrows_val) for (k = 0;;)
			{
				/*
				int mu1;
				{
				Real big = fabs(*akk); mu1 = k; Real* ai = akk; int i;

				for (i = k+1; i < nrows_val; i++)
				{
				ai += nrows_val; const Real trybig = fabs(*ai);
				if (big < trybig) { big = trybig; mu1 = i; }
				}
				}
				if (mu1 != mu) cout << k << " " << mu << " " << mu1 << endl;
				*/

				indx[k] = mu;

				if (mu != k)                       //row swap
				{
					Real* a1 = store + nrows_val * k;
					Real* a2 = store + nrows_val * mu; d = !d;
					int j = nrows_val;
					while (j--) { const Real temp = *a1; *a1++ = *a2; *a2++ = temp; }
				}

				Real diag = *akk; big = 0; mu = k + 1;
				if (diag != 0)
				{
					ai = akk; int i = nrows_val - k - 1;
					while (i--)
					{
						ai += nrows_val; Real* al = ai;
						Real mult = *al / diag; *al = mult;
						int l = nrows_val - k - 1; Real* aj = akk;
						// work out the next pivot as part of this loop
						// this saves a column operation
						if (l-- != 0)
						{
							*(++al) -= (mult * *(++aj));
							const Real trybig = fabs(*al);
							if (big < trybig) { big = trybig; mu = nrows_val - i - 1; }
							while (l--) *(++al) -= (mult * *(++aj));
						}
					}
				}
				else sing = true;
				if (++k == nrows_val) break;          // so next line won't overflow
				akk += nrows_val + 1;
			}
		}

		void CroutMatrix::lubksb(Real* B, int mini)
		{
			REPORT
				// this has been adapted from Numerical Recipes in C. The code has been
				// substantially streamlined, so I do not think much of the original
				// copyright remains. However there is not much opportunity for
				// variation in the code, so it is still similar to the NR code.
				// I follow the NR code in skipping over initial zeros in the B vector.

				Tracer tr("Crout(lubksb)");
			if (sing) Throw(SingularException(*this));
			int i, j, ii = nrows_val;       // ii initialised : B might be all zeros


			// scan for first non-zero in B
			for (i = 0; i < nrows_val; i++)
			{
				int ip = indx[i]; Real temp = B[ip]; B[ip] = B[i]; B[i] = temp;
				if (temp != 0.0) { ii = i; break; }
			}

			Real* bi; Real* ai;
			i = ii + 1;

			if (i < nrows_val)
			{
				bi = B + ii; ai = store + ii + i * nrows_val;
				for (;;)
				{
					int ip = indx[i]; Real sum = B[ip]; B[ip] = B[i];
					Real* aij = ai; Real* bj = bi; j = i - ii;
					while (j--) sum -= *aij++ * *bj++;
					B[i] = sum;
					if (++i == nrows_val) break;
					ai += nrows_val;
				}
			}

			ai = store + nrows_val * nrows_val;

			for (i = nrows_val - 1; i >= mini; i--)
			{
				Real* bj = B+i; ai -= nrows_val; Real* ajx = ai+i;
				Real sum = *bj; Real diag = *ajx;
				j = nrows_val - i; while(--j) sum -= *(++ajx) * *(++bj);
				B[i] = sum / diag;
			}
		}

		/****************************** scalar functions ****************************/

		inline Real square(Real x) { return x*x; }

		Real GeneralMatrix::sum_square() const
		{
			REPORT
				Real sum = 0.0; int i = storage; Real* s = store;
			while (i--) sum += square(*s++);
			((GeneralMatrix&)*this).tDelete(); return sum;
		}

		Real GeneralMatrix::sum_absolute_value() const
		{
			REPORT
				Real sum = 0.0; int i = storage; Real* s = store;
			while (i--) sum += fabs(*s++);
			((GeneralMatrix&)*this).tDelete(); return sum;
		}

		Real GeneralMatrix::sum() const
		{
			REPORT
				Real sm = 0.0; int i = storage; Real* s = store;
			while (i--) sm += *s++;
			((GeneralMatrix&)*this).tDelete(); return sm;
		}

		// maxima and minima

		// There are three sets of routines
		// maximum_absolute_value, minimum_absolute_value, maximum, minimum
		// ... these find just the maxima and minima
		// maximum_absolute_value1, minimum_absolute_value1, maximum1, minimum1
		// ... these find the maxima and minima and their locations in a
		//     one dimensional object
		// maximum_absolute_value2, minimum_absolute_value2, maximum2, minimum2
		// ... these find the maxima and minima and their locations in a
		//     two dimensional object

		// If the matrix has no values throw an exception

		// If we do not want the location find the maximum or minimum on the
		// array stored by GeneralMatrix
		// This won't work for BandMatrices. We call ClearCorner for
		// maximum_absolute_value but for the others use the absolute_minimum_value2
		// version and discard the location.

		// For one dimensional objects, when we want the location of the
		// maximum or minimum, work with the array stored by GeneralMatrix

		// For two dimensional objects where we want the location of the maximum or
		// minimum proceed as follows:

		// For rectangular matrices use the array stored by GeneralMatrix and
		// deduce the location from the location in the GeneralMatrix

		// For other two dimensional matrices use the Matrix Row routine to find the
		// maximum or minimum for each row.

		static void NullMatrixError(const GeneralMatrix* gm)
		{
			((GeneralMatrix&)*gm).tDelete();
			Throw(ProgramException("Maximum or minimum of null matrix"));
		}

		Real GeneralMatrix::maximum_absolute_value() const
		{
			REPORT
				if (storage == 0) NullMatrixError(this);
			Real maxval = 0.0; int l = storage; Real* s = store;
			while (l--) { Real a = fabs(*s++); if (maxval < a) maxval = a; }
			((GeneralMatrix&)*this).tDelete(); return maxval;
		}

		Real GeneralMatrix::maximum_absolute_value1(int& i) const
		{
			REPORT
				if (storage == 0) NullMatrixError(this);
			Real maxval = 0.0; int l = storage; Real* s = store; int li = storage;
			while (l--)
			{ Real a = fabs(*s++); if (maxval <= a) { maxval = a; li = l; }  }
			i = storage - li;
			((GeneralMatrix&)*this).tDelete(); return maxval;
		}

		Real GeneralMatrix::minimum_absolute_value() const
		{
			REPORT
				if (storage == 0) NullMatrixError(this);
			int l = storage - 1; Real* s = store; Real minval = fabs(*s++);
			while (l--) { Real a = fabs(*s++); if (minval > a) minval = a; }
			((GeneralMatrix&)*this).tDelete(); return minval;
		}

		Real GeneralMatrix::minimum_absolute_value1(int& i) const
		{
			REPORT
				if (storage == 0) NullMatrixError(this);
			int l = storage - 1; Real* s = store; Real minval = fabs(*s++); int li = l;
			while (l--)
			{ Real a = fabs(*s++); if (minval >= a) { minval = a; li = l; }  }
			i = storage - li;
			((GeneralMatrix&)*this).tDelete(); return minval;
		}

		Real GeneralMatrix::maximum() const
		{
			REPORT
				if (storage == 0) NullMatrixError(this);
			int l = storage - 1; Real* s = store; Real maxval = *s++;
			while (l--) { Real a = *s++; if (maxval < a) maxval = a; }
			((GeneralMatrix&)*this).tDelete(); return maxval;
		}

		Real GeneralMatrix::maximum1(int& i) const
		{
			REPORT
				if (storage == 0) NullMatrixError(this);
			int l = storage - 1; Real* s = store; Real maxval = *s++; int li = l;
			while (l--) { Real a = *s++; if (maxval <= a) { maxval = a; li = l; } }
			i = storage - li;
			((GeneralMatrix&)*this).tDelete(); return maxval;
		}

		Real GeneralMatrix::minimum() const
		{
			REPORT
				if (storage == 0) NullMatrixError(this);
			int l = storage - 1; Real* s = store; Real minval = *s++;
			while (l--) { Real a = *s++; if (minval > a) minval = a; }
			((GeneralMatrix&)*this).tDelete(); return minval;
		}

		Real GeneralMatrix::minimum1(int& i) const
		{
			REPORT
				if (storage == 0) NullMatrixError(this);
			int l = storage - 1; Real* s = store; Real minval = *s++; int li = l;
			while (l--) { Real a = *s++; if (minval >= a) { minval = a; li = l; } }
			i = storage - li;
			((GeneralMatrix&)*this).tDelete(); return minval;
		}

		Real GeneralMatrix::maximum_absolute_value2(int& i, int& j) const
		{
			REPORT
				if (storage == 0) NullMatrixError(this);
			Real maxval = 0.0; int nr = Nrows();
			MatrixRow mr((GeneralMatrix*)this, LoadOnEntry+DirectPart);
			for (int r = 1; r <= nr; r++)
			{
				int c; maxval = mr.MaximumAbsoluteValue1(maxval, c);
				if (c > 0) { i = r; j = c; }
				mr.Next();
			}
			((GeneralMatrix&)*this).tDelete(); return maxval;
		}

		Real GeneralMatrix::minimum_absolute_value2(int& i, int& j) const
		{
			REPORT
				if (storage == 0)  NullMatrixError(this);
			Real minval = FloatingPointPrecision::Maximum(); int nr = Nrows();
			MatrixRow mr((GeneralMatrix*)this, LoadOnEntry+DirectPart);
			for (int r = 1; r <= nr; r++)
			{
				int c; minval = mr.MinimumAbsoluteValue1(minval, c);
				if (c > 0) { i = r; j = c; }
				mr.Next();
			}
			((GeneralMatrix&)*this).tDelete(); return minval;
		}

		Real GeneralMatrix::maximum2(int& i, int& j) const
		{
			REPORT
				if (storage == 0) NullMatrixError(this);
			Real maxval = -FloatingPointPrecision::Maximum(); int nr = Nrows();
			MatrixRow mr((GeneralMatrix*)this, LoadOnEntry+DirectPart);
			for (int r = 1; r <= nr; r++)
			{
				int c; maxval = mr.Maximum1(maxval, c);
				if (c > 0) { i = r; j = c; }
				mr.Next();
			}
			((GeneralMatrix&)*this).tDelete(); return maxval;
		}

		Real GeneralMatrix::minimum2(int& i, int& j) const
		{
			REPORT
				if (storage == 0) NullMatrixError(this);
			Real minval = FloatingPointPrecision::Maximum(); int nr = Nrows();
			MatrixRow mr((GeneralMatrix*)this, LoadOnEntry+DirectPart);
			for (int r = 1; r <= nr; r++)
			{
				int c; minval = mr.Minimum1(minval, c);
				if (c > 0) { i = r; j = c; }
				mr.Next();
			}
			((GeneralMatrix&)*this).tDelete(); return minval;
		}

		Real Matrix::maximum_absolute_value2(int& i, int& j) const
		{
			REPORT
				int k; Real m = GeneralMatrix::maximum_absolute_value1(k); k--;
			i = k / Ncols(); j = k - i * Ncols(); i++; j++;
			return m;
		}

		Real Matrix::minimum_absolute_value2(int& i, int& j) const
		{
			REPORT
				int k; Real m = GeneralMatrix::minimum_absolute_value1(k); k--;
			i = k / Ncols(); j = k - i * Ncols(); i++; j++;
			return m;
		}

		Real Matrix::maximum2(int& i, int& j) const
		{
			REPORT
				int k; Real m = GeneralMatrix::maximum1(k); k--;
			i = k / Ncols(); j = k - i * Ncols(); i++; j++;
			return m;
		}

		Real Matrix::minimum2(int& i, int& j) const
		{
			REPORT
				int k; Real m = GeneralMatrix::minimum1(k); k--;
			i = k / Ncols(); j = k - i * Ncols(); i++; j++;
			return m;
		}

		Real SymmetricMatrix::sum_square() const
		{
			REPORT
				Real sum1 = 0.0; Real sum2 = 0.0; Real* s = store; int nr = nrows_val;
			for (int i = 0; i<nr; i++)
			{
				int j = i;
				while (j--) sum2 += square(*s++);
				sum1 += square(*s++);
			}
			((GeneralMatrix&)*this).tDelete(); return sum1 + 2.0 * sum2;
		}

		Real SymmetricMatrix::sum_absolute_value() const
		{
			REPORT
				Real sum1 = 0.0; Real sum2 = 0.0; Real* s = store; int nr = nrows_val;
			for (int i = 0; i<nr; i++)
			{
				int j = i;
				while (j--) sum2 += fabs(*s++);
				sum1 += fabs(*s++);
			}
			((GeneralMatrix&)*this).tDelete(); return sum1 + 2.0 * sum2;
		}

		Real IdentityMatrix::sum_absolute_value() const
		{ REPORT  return fabs(trace()); }    // no need to do tDelete?

		Real SymmetricMatrix::sum() const
		{
			REPORT
				Real sum1 = 0.0; Real sum2 = 0.0; Real* s = store; int nr = nrows_val;
			for (int i = 0; i<nr; i++)
			{
				int j = i;
				while (j--) sum2 += *s++;
				sum1 += *s++;
			}
			((GeneralMatrix&)*this).tDelete(); return sum1 + 2.0 * sum2;
		}

		Real IdentityMatrix::sum_square() const
		{
			Real sum = *store * *store * nrows_val;
			((GeneralMatrix&)*this).tDelete(); return sum;
		}


		Real BaseMatrix::sum_square() const
		{
			REPORT GeneralMatrix* gm = ((BaseMatrix&)*this).Evaluate();
			Real s = gm->sum_square(); return s;
		}

		Real BaseMatrix::norm_Frobenius() const
		{ REPORT  return sqrt(sum_square()); }

		Real BaseMatrix::sum_absolute_value() const
		{
			REPORT GeneralMatrix* gm = ((BaseMatrix&)*this).Evaluate();
			Real s = gm->sum_absolute_value(); return s;
		}

		Real BaseMatrix::sum() const
		{
			REPORT GeneralMatrix* gm = ((BaseMatrix&)*this).Evaluate();
			Real s = gm->sum(); return s;
		}

		Real BaseMatrix::maximum_absolute_value() const
		{
			REPORT GeneralMatrix* gm = ((BaseMatrix&)*this).Evaluate();
			Real s = gm->maximum_absolute_value(); return s;
		}

		Real BaseMatrix::maximum_absolute_value1(int& i) const
		{
			REPORT GeneralMatrix* gm = ((BaseMatrix&)*this).Evaluate();
			Real s = gm->maximum_absolute_value1(i); return s;
		}

		Real BaseMatrix::maximum_absolute_value2(int& i, int& j) const
		{
			REPORT GeneralMatrix* gm = ((BaseMatrix&)*this).Evaluate();
			Real s = gm->maximum_absolute_value2(i, j); return s;
		}

		Real BaseMatrix::minimum_absolute_value() const
		{
			REPORT GeneralMatrix* gm = ((BaseMatrix&)*this).Evaluate();
			Real s = gm->minimum_absolute_value(); return s;
		}

		Real BaseMatrix::minimum_absolute_value1(int& i) const
		{
			REPORT GeneralMatrix* gm = ((BaseMatrix&)*this).Evaluate();
			Real s = gm->minimum_absolute_value1(i); return s;
		}

		Real BaseMatrix::minimum_absolute_value2(int& i, int& j) const
		{
			REPORT GeneralMatrix* gm = ((BaseMatrix&)*this).Evaluate();
			Real s = gm->minimum_absolute_value2(i, j); return s;
		}

		Real BaseMatrix::maximum() const
		{
			REPORT GeneralMatrix* gm = ((BaseMatrix&)*this).Evaluate();
			Real s = gm->maximum(); return s;
		}

		Real BaseMatrix::maximum1(int& i) const
		{
			REPORT GeneralMatrix* gm = ((BaseMatrix&)*this).Evaluate();
			Real s = gm->maximum1(i); return s;
		}

		Real BaseMatrix::maximum2(int& i, int& j) const
		{
			REPORT GeneralMatrix* gm = ((BaseMatrix&)*this).Evaluate();
			Real s = gm->maximum2(i, j); return s;
		}

		Real BaseMatrix::minimum() const
		{
			REPORT GeneralMatrix* gm = ((BaseMatrix&)*this).Evaluate();
			Real s = gm->minimum(); return s;
		}

		Real BaseMatrix::minimum1(int& i) const
		{
			REPORT GeneralMatrix* gm = ((BaseMatrix&)*this).Evaluate();
			Real s = gm->minimum1(i); return s;
		}

		Real BaseMatrix::minimum2(int& i, int& j) const
		{
			REPORT GeneralMatrix* gm = ((BaseMatrix&)*this).Evaluate();
			Real s = gm->minimum2(i, j); return s;
		}

		Real dotproduct(const Matrix& A, const Matrix& B)
		{
			REPORT
				int n = A.storage;
			if (n != B.storage)
			{
				Tracer tr("dotproduct");
				Throw(IncompatibleDimensionsException(A,B));
			}
			Real sum = 0.0; Real* a = A.store; Real* b = B.store;
			while (n--) sum += *a++ * *b++;
			return sum;
		}

		Real Matrix::trace() const
		{
			REPORT
				Tracer tr("trace");
			int i = nrows_val; int d = i+1;
			if (i != ncols_val) Throw(NotSquareException(*this));
			Real sum = 0.0; Real* s = store;
			//   while (i--) { sum += *s; s += d; }
			if (i) for (;;) { sum += *s; if (!(--i)) break; s += d; }
			((GeneralMatrix&)*this).tDelete(); return sum;
		}

		Real DiagonalMatrix::trace() const
		{
			REPORT
				int i = nrows_val; Real sum = 0.0; Real* s = store;
			while (i--) sum += *s++;
			((GeneralMatrix&)*this).tDelete(); return sum;
		}

		Real SymmetricMatrix::trace() const
		{
			REPORT
				int i = nrows_val; Real sum = 0.0; Real* s = store; int j = 2;
			// while (i--) { sum += *s; s += j++; }
			if (i) for (;;) { sum += *s; if (!(--i)) break; s += j++; }
			((GeneralMatrix&)*this).tDelete(); return sum;
		}

		Real LowerTriangularMatrix::trace() const
		{
			REPORT
				int i = nrows_val; Real sum = 0.0; Real* s = store; int j = 2;
			// while (i--) { sum += *s; s += j++; }
			if (i) for (;;) { sum += *s; if (!(--i)) break; s += j++; }
			((GeneralMatrix&)*this).tDelete(); return sum;
		}

		Real UpperTriangularMatrix::trace() const
		{
			REPORT
				int i = nrows_val; Real sum = 0.0; Real* s = store;
			while (i) { sum += *s; s += i--; }             // won t cause a problem
			((GeneralMatrix&)*this).tDelete(); return sum;
		}

		Real BandMatrix::trace() const
		{
			REPORT
				int i = nrows_val; int w = lower_val+upper_val+1;
			Real sum = 0.0; Real* s = store+lower_val;
			// while (i--) { sum += *s; s += w; }
			if (i) for (;;) { sum += *s; if (!(--i)) break; s += w; }
			((GeneralMatrix&)*this).tDelete(); return sum;
		}

		Real SymmetricBandMatrix::trace() const
		{
			REPORT
				int i = nrows_val; int w = lower_val+1;
			Real sum = 0.0; Real* s = store+lower_val;
			// while (i--) { sum += *s; s += w; }
			if (i) for (;;) { sum += *s; if (!(--i)) break; s += w; }
			((GeneralMatrix&)*this).tDelete(); return sum;
		}

		Real IdentityMatrix::trace() const
		{
			Real sum = *store * nrows_val;
			((GeneralMatrix&)*this).tDelete(); return sum;
		}


		Real BaseMatrix::trace() const
		{
			REPORT
				MatrixType Diag = MatrixType::Dg; Diag.SetDataLossOK();
			GeneralMatrix* gm = ((BaseMatrix&)*this).Evaluate(Diag);
			Real sum = gm->trace(); return sum;
		}

		void LogAndSign::operator*=(Real x)
		{
			if (x > 0.0) { log_val += log(x); }
			else if (x < 0.0) { log_val += log(-x); sign_val = -sign_val; }
			else sign_val = 0;
		}

		void LogAndSign::pow_eq(int k)
		{
			if (sign_val)
			{
				log_val *= k;
				if ( (k & 1) == 0 ) sign_val = 1;
			}
		}

		Real LogAndSign::value() const
		{
			Tracer et("LogAndSign::value");
			if (log_val >= FloatingPointPrecision::LnMaximum())
				Throw(OverflowException("Overflow in exponential"));
			return sign_val * exp(log_val);
		}

		LogAndSign::LogAndSign(Real f)
		{
			if (f == 0.0) { log_val = 0.0; sign_val = 0; return; }
			else if (f < 0.0) { sign_val = -1; f = -f; }
			else sign_val = 1;
			log_val = log(f);
		}

		LogAndSign DiagonalMatrix::log_determinant() const
		{
			REPORT
				int i = nrows_val; LogAndSign sum; Real* s = store;
			while (i--) sum *= *s++;
			((GeneralMatrix&)*this).tDelete(); return sum;
		}

		LogAndSign LowerTriangularMatrix::log_determinant() const
		{
			REPORT
				int i = nrows_val; LogAndSign sum; Real* s = store; int j = 2;
			// while (i--) { sum *= *s; s += j++; }
			if (i) for(;;) { sum *= *s; if (!(--i)) break; s += j++; }
			((GeneralMatrix&)*this).tDelete(); return sum;
		}

		LogAndSign UpperTriangularMatrix::log_determinant() const
		{
			REPORT
				int i = nrows_val; LogAndSign sum; Real* s = store;
			while (i) { sum *= *s; s += i--; }
			((GeneralMatrix&)*this).tDelete(); return sum;
		}

		LogAndSign IdentityMatrix::log_determinant() const
		{
			REPORT
				int i = nrows_val; LogAndSign sum;
			if (i > 0) { sum = *store; sum.PowEq(i); }
			((GeneralMatrix&)*this).tDelete(); return sum;
		}

		LogAndSign BaseMatrix::log_determinant() const
		{
			REPORT GeneralMatrix* gm = ((BaseMatrix&)*this).Evaluate();
			LogAndSign sum = gm->log_determinant(); return sum;
		}

		LogAndSign GeneralMatrix::log_determinant() const
		{
			REPORT
				Tracer tr("log_determinant");
			if (nrows_val != ncols_val) Throw(NotSquareException(*this));
			CroutMatrix C(*this); return C.log_determinant();
		}

		LogAndSign CroutMatrix::log_determinant() const
		{
			REPORT
				if (sing) return 0.0;
			int i = nrows_val; int dd = i+1; LogAndSign sum; Real* s = store;
			if (i) for(;;)
			{
				sum *= *s;
				if (!(--i)) break;
				s += dd;
			}
			if (!d) sum.ChangeSign(); return sum;

		}

		Real BaseMatrix::determinant() const
		{
			REPORT
				Tracer tr("determinant");
			REPORT GeneralMatrix* gm = ((BaseMatrix&)*this).Evaluate();
			LogAndSign ld = gm->log_determinant();
			return ld.Value();
		}

		LinearEquationSolver::LinearEquationSolver(const BaseMatrix& bm)
		{
			gm = ( ((BaseMatrix&)bm).Evaluate() )->MakeSolver();
			if (gm==&bm) { REPORT  gm = gm->Image(); }
			// want a copy if  *gm is actually bm
			else { REPORT  gm->Protect(); }
		}

		ReturnMatrix BaseMatrix::sum_square_rows() const
		{
			REPORT
				GeneralMatrix* gm = ((BaseMatrix&)*this).Evaluate();
			int nr = gm->nrows();
			ColumnVector ssq(nr);
			if (gm->size() == 0) { REPORT ssq = 0.0; }
			else
			{
				MatrixRow mr(gm, LoadOnEntry);
				for (int i = 1; i <= nr; ++i)
				{
					Real sum = 0.0;
					int s = mr.Storage();
					Real* in = mr.Data();
					while (s--) sum += square(*in++);
					ssq(i) = sum;   
					mr.Next();
				}
			}
			gm->tDelete();
			ssq.release(); return ssq.for_return();
		}

		ReturnMatrix BaseMatrix::sum_square_columns() const
		{
			REPORT
				GeneralMatrix* gm = ((BaseMatrix&)*this).Evaluate();
			int nr = gm->nrows(); int nc = gm->ncols();
			RowVector ssq(nc); ssq = 0.0;
			if (gm->size() != 0)
			{
				MatrixRow mr(gm, LoadOnEntry);
				for (int i = 1; i <= nr; ++i)
				{
					int s = mr.Storage();
					Real* in = mr.Data(); Real* out = ssq.data() + mr.Skip();
					while (s--) *out++ += square(*in++);
					mr.Next();
				}
			}
			gm->tDelete();
			ssq.release(); return ssq.for_return();
		}

		ReturnMatrix BaseMatrix::sum_rows() const
		{
			REPORT
				GeneralMatrix* gm = ((BaseMatrix&)*this).Evaluate();
			int nr = gm->nrows();
			ColumnVector sum_vec(nr);
			if (gm->size() == 0) { REPORT sum_vec = 0.0; }
			else
			{
				MatrixRow mr(gm, LoadOnEntry);
				for (int i = 1; i <= nr; ++i)
				{
					Real sum = 0.0;
					int s = mr.Storage();
					Real* in = mr.Data();
					while (s--) sum += *in++;
					sum_vec(i) = sum;   
					mr.Next();
				}
			}
			gm->tDelete();
			sum_vec.release(); return sum_vec.for_return();
		}

		ReturnMatrix BaseMatrix::sum_columns() const
		{
			REPORT
				GeneralMatrix* gm = ((BaseMatrix&)*this).Evaluate();
			int nr = gm->nrows(); int nc = gm->ncols();
			RowVector sum_vec(nc); sum_vec = 0.0;
			if (gm->size() != 0)
			{
				MatrixRow mr(gm, LoadOnEntry);
				for (int i = 1; i <= nr; ++i)
				{
					int s = mr.Storage();
					Real* in = mr.Data(); Real* out = sum_vec.data() + mr.Skip();
					while (s--) *out++ += *in++;
					mr.Next();
				}
			}
			gm->tDelete();
			sum_vec.release(); return sum_vec.for_return();
		}


#ifdef use_namespace
	}
#endif

//#include "include.h"

#include "newmat.h"
//#include "newmatrc.h"

#ifdef use_namespace
	namespace NEWMAT {
#endif


#ifdef DO_REPORT
#define REPORT { static ExeCounter ExeCount(__LINE__,7); ++ExeCount; }
#else
#define REPORT {}
#endif


		//***************************** solve routines ******************************/

		GeneralMatrix* GeneralMatrix::MakeSolver()
		{
			REPORT
				GeneralMatrix* gm = new CroutMatrix(*this);
			MatrixErrorNoSpace(gm); gm->ReleaseAndDelete(); return gm;
		}

		GeneralMatrix* Matrix::MakeSolver()
		{
			REPORT
				GeneralMatrix* gm = new CroutMatrix(*this);
			MatrixErrorNoSpace(gm); gm->ReleaseAndDelete(); return gm;
		}

		void CroutMatrix::Solver(MatrixColX& mcout, const MatrixColX& mcin)
		{
			REPORT
				int i = mcin.skip; Real* el = mcin.data-i; Real* el1 = el;
			while (i--) *el++ = 0.0;
			el += mcin.storage; i = nrows_val - mcin.skip - mcin.storage;
			while (i--) *el++ = 0.0;
			lubksb(el1, mcout.skip);
		}


		// Do we need check for entirely zero output?

		void UpperTriangularMatrix::Solver(MatrixColX& mcout,
			const MatrixColX& mcin)
		{
			REPORT
				int i = mcin.skip-mcout.skip; Real* elx = mcin.data-i;
			while (i-- > 0) *elx++ = 0.0;
			int nr = mcin.skip+mcin.storage;
			elx = mcin.data+mcin.storage; Real* el = elx;
			int j = mcout.skip+mcout.storage-nr;
			int nc = ncols_val-nr; i = nr-mcout.skip;
			while (j-- > 0) *elx++ = 0.0;
			Real* Ael = store + (nr*(2*ncols_val-nr+1))/2; j = 0;
			while (i-- > 0)
			{
				elx = el; Real sum = 0.0; int jx = j++; Ael -= nc;
				while (jx--) sum += *(--Ael) * *(--elx);
				elx--; *elx = (*elx - sum) / *(--Ael);
			}
		}

		void LowerTriangularMatrix::Solver(MatrixColX& mcout,
			const MatrixColX& mcin)
		{
			REPORT
				int i = mcin.skip-mcout.skip; Real* elx = mcin.data-i;
			while (i-- > 0) *elx++ = 0.0;
			int nc = mcin.skip; i = nc+mcin.storage; elx = mcin.data+mcin.storage;
			int nr = mcout.skip+mcout.storage; int j = nr-i; i = nr-nc;
			while (j-- > 0) *elx++ = 0.0;
			Real* el = mcin.data; Real* Ael = store + (nc*(nc+1))/2; j = 0;
			while (i-- > 0)
			{
				elx = el; Real sum = 0.0; int jx = j++; Ael += nc;
				while (jx--) sum += *Ael++ * *elx++;
				*elx = (*elx - sum) / *Ael++;
			}
		}

		//******************* carry out binary operations *************************/

		static GeneralMatrix*
			GeneralMult(GeneralMatrix*,GeneralMatrix*,MultipliedMatrix*,MatrixType);
		static GeneralMatrix*
			GeneralSolv(GeneralMatrix*,GeneralMatrix*,BaseMatrix*,MatrixType);
		static GeneralMatrix*
			GeneralSolvI(GeneralMatrix*,BaseMatrix*,MatrixType);
		static GeneralMatrix*
			GeneralKP(GeneralMatrix*,GeneralMatrix*,KPMatrix*,MatrixType);

		GeneralMatrix* MultipliedMatrix::Evaluate(MatrixType mt)
		{
			REPORT
				gm2 = ((BaseMatrix*&)bm2)->Evaluate();
			gm2 = gm2->Evaluate(gm2->type().MultRHS());     // no symmetric on RHS
			gm1=((BaseMatrix*&)bm1)->Evaluate();
			return GeneralMult(gm1, gm2, this, mt);
		}

		GeneralMatrix* SolvedMatrix::Evaluate(MatrixType mt)
		{
			REPORT
				gm1=((BaseMatrix*&)bm1)->Evaluate();
			gm2=((BaseMatrix*&)bm2)->Evaluate();
			return GeneralSolv(gm1,gm2,this,mt);
		}

		GeneralMatrix* KPMatrix::Evaluate(MatrixType mt)
		{
			REPORT
				gm1=((BaseMatrix*&)bm1)->Evaluate();
			gm2=((BaseMatrix*&)bm2)->Evaluate();
			return GeneralKP(gm1,gm2,this,mt);
		}

		// routines for adding or subtracting matrices of identical storage structure

		static void Add(GeneralMatrix* gm, GeneralMatrix* gm1, GeneralMatrix* gm2)
		{
			REPORT
				Real* s1=gm1->Store(); Real* s2=gm2->Store();
			Real* s=gm->Store(); int i=gm->Storage() >> 2;
			while (i--)
			{
				*s++ = *s1++ + *s2++; *s++ = *s1++ + *s2++;
				*s++ = *s1++ + *s2++; *s++ = *s1++ + *s2++;
			}
			i=gm->Storage() & 3; while (i--) *s++ = *s1++ + *s2++;
		}

		static void AddTo(GeneralMatrix* gm, const GeneralMatrix* gm2)
		{
			REPORT
				const Real* s2=gm2->Store(); Real* s=gm->Store(); int i=gm->Storage() >> 2;
			while (i--)
			{ *s++ += *s2++; *s++ += *s2++; *s++ += *s2++; *s++ += *s2++; }
			i=gm->Storage() & 3; while (i--) *s++ += *s2++;
		}

		void GeneralMatrix::PlusEqual(const GeneralMatrix& gm)
		{
			REPORT
				if (nrows_val != gm.nrows_val || ncols_val != gm.ncols_val)
					Throw(IncompatibleDimensionsException(*this, gm));
			AddTo(this, &gm);
		}

		static void Subtract(GeneralMatrix* gm, GeneralMatrix* gm1, GeneralMatrix* gm2)
		{
			REPORT
				Real* s1=gm1->Store(); Real* s2=gm2->Store();
			Real* s=gm->Store(); int i=gm->Storage() >> 2;
			while (i--)
			{
				*s++ = *s1++ - *s2++; *s++ = *s1++ - *s2++;
				*s++ = *s1++ - *s2++; *s++ = *s1++ - *s2++;
			}
			i=gm->Storage() & 3; while (i--) *s++ = *s1++ - *s2++;
		}

		static void SubtractFrom(GeneralMatrix* gm, const GeneralMatrix* gm2)
		{
			REPORT
				const Real* s2=gm2->Store(); Real* s=gm->Store(); int i=gm->Storage() >> 2;
			while (i--)
			{ *s++ -= *s2++; *s++ -= *s2++; *s++ -= *s2++; *s++ -= *s2++; }
			i=gm->Storage() & 3; while (i--) *s++ -= *s2++;
		}

		void GeneralMatrix::MinusEqual(const GeneralMatrix& gm)
		{
			REPORT
				if (nrows_val != gm.nrows_val || ncols_val != gm.ncols_val)
					Throw(IncompatibleDimensionsException(*this, gm));
			SubtractFrom(this, &gm);
		}

		static void ReverseSubtract(GeneralMatrix* gm, const GeneralMatrix* gm2)
		{
			REPORT
				const Real* s2=gm2->Store(); Real* s=gm->Store(); int i=gm->Storage() >> 2;
			while (i--)
			{
				*s = *s2++ - *s; s++; *s = *s2++ - *s; s++;
				*s = *s2++ - *s; s++; *s = *s2++ - *s; s++;
			}
			i=gm->Storage() & 3; while (i--) { *s = *s2++ - *s; s++; }
		}

		static void SP(GeneralMatrix* gm, GeneralMatrix* gm1, GeneralMatrix* gm2)
		{
			REPORT
				Real* s1=gm1->Store(); Real* s2=gm2->Store();
			Real* s=gm->Store(); int i=gm->Storage() >> 2;
			while (i--)
			{
				*s++ = *s1++ * *s2++; *s++ = *s1++ * *s2++;
				*s++ = *s1++ * *s2++; *s++ = *s1++ * *s2++;
			}
			i=gm->Storage() & 3; while (i--) *s++ = *s1++ * *s2++;
		}

		static void SP(GeneralMatrix* gm, GeneralMatrix* gm2)
		{
			REPORT
				Real* s2=gm2->Store(); Real* s=gm->Store(); int i=gm->Storage() >> 2;
			while (i--)
			{ *s++ *= *s2++; *s++ *= *s2++; *s++ *= *s2++; *s++ *= *s2++; }
			i=gm->Storage() & 3; while (i--) *s++ *= *s2++;
		}

		// routines for adding or subtracting matrices of different storage structure

		static void AddDS(GeneralMatrix* gm, GeneralMatrix* gm1, GeneralMatrix* gm2)
		{
			REPORT
				int nr = gm->Nrows();
			MatrixRow mr1(gm1, LoadOnEntry); MatrixRow mr2(gm2, LoadOnEntry);
			MatrixRow mr(gm, StoreOnExit+DirectPart);
			while (nr--) { mr.Add(mr1,mr2); mr1.Next(); mr2.Next(); mr.Next(); }
		}

		static void AddDS(GeneralMatrix* gm, GeneralMatrix* gm2)
			// Add into first argument
		{
			REPORT
				int nr = gm->Nrows();
			MatrixRow mr(gm, StoreOnExit+LoadOnEntry+DirectPart);
			MatrixRow mr2(gm2, LoadOnEntry);
			while (nr--) { mr.Add(mr2); mr.Next(); mr2.Next(); }
		}

		static void SubtractDS
			(GeneralMatrix* gm, GeneralMatrix* gm1, GeneralMatrix* gm2)
		{
			REPORT
				int nr = gm->Nrows();
			MatrixRow mr1(gm1, LoadOnEntry); MatrixRow mr2(gm2, LoadOnEntry);
			MatrixRow mr(gm, StoreOnExit+DirectPart);
			while (nr--) { mr.Sub(mr1,mr2); mr1.Next(); mr2.Next(); mr.Next(); }
		}

		static void SubtractDS(GeneralMatrix* gm, GeneralMatrix* gm2)
		{
			REPORT
				int nr = gm->Nrows();
			MatrixRow mr(gm, LoadOnEntry+StoreOnExit+DirectPart);
			MatrixRow mr2(gm2, LoadOnEntry);
			while (nr--) { mr.Sub(mr2); mr.Next(); mr2.Next(); }
		}

		static void ReverseSubtractDS(GeneralMatrix* gm, GeneralMatrix* gm2)
		{
			REPORT
				int nr = gm->Nrows();
			MatrixRow mr(gm, LoadOnEntry+StoreOnExit+DirectPart);
			MatrixRow mr2(gm2, LoadOnEntry);
			while (nr--) { mr.RevSub(mr2); mr2.Next(); mr.Next(); }
		}

		static void SPDS(GeneralMatrix* gm, GeneralMatrix* gm1, GeneralMatrix* gm2)
		{
			REPORT
				int nr = gm->Nrows();
			MatrixRow mr1(gm1, LoadOnEntry); MatrixRow mr2(gm2, LoadOnEntry);
			MatrixRow mr(gm, StoreOnExit+DirectPart);
			while (nr--) { mr.Multiply(mr1,mr2); mr1.Next(); mr2.Next(); mr.Next(); }
		}

		static void SPDS(GeneralMatrix* gm, GeneralMatrix* gm2)
			// SP into first argument
		{
			REPORT
				int nr = gm->Nrows();
			MatrixRow mr(gm, StoreOnExit+LoadOnEntry+DirectPart);
			MatrixRow mr2(gm2, LoadOnEntry);
			while (nr--) { mr.Multiply(mr2); mr.Next(); mr2.Next(); }
		}

		static GeneralMatrix* GeneralMult1(GeneralMatrix* gm1, GeneralMatrix* gm2,
			MultipliedMatrix* mm, MatrixType mtx)
		{
			REPORT
				Tracer tr("GeneralMult1");
			int nr=gm1->Nrows(); int nc=gm2->Ncols();
			if (gm1->Ncols() !=gm2->Nrows())
				Throw(IncompatibleDimensionsException(*gm1, *gm2));
			GeneralMatrix* gmx = mtx.New(nr,nc,mm);

			MatrixCol mcx(gmx, StoreOnExit+DirectPart);
			MatrixCol mc2(gm2, LoadOnEntry);
			while (nc--)
			{
				MatrixRow mr1(gm1, LoadOnEntry, mcx.Skip());
				Real* el = mcx.Data();                         // pointer to an element
				int n = mcx.Storage();
				while (n--) { *(el++) = DotProd(mr1,mc2); mr1.Next(); }
				mc2.Next(); mcx.Next();
			}
			gmx->ReleaseAndDelete(); gm1->tDelete(); gm2->tDelete(); return gmx;
		}

		static GeneralMatrix* GeneralMult2(GeneralMatrix* gm1, GeneralMatrix* gm2,
			MultipliedMatrix* mm, MatrixType mtx)
		{
			// version that accesses by row only - not good for thin matrices
			// or column vectors in right hand term.
			REPORT
				Tracer tr("GeneralMult2");
			int nr=gm1->Nrows(); int nc=gm2->Ncols();
			if (gm1->Ncols() !=gm2->Nrows())
				Throw(IncompatibleDimensionsException(*gm1, *gm2));
			GeneralMatrix* gmx = mtx.New(nr,nc,mm);

			MatrixRow mrx(gmx, LoadOnEntry+StoreOnExit+DirectPart);
			MatrixRow mr1(gm1, LoadOnEntry);
			while (nr--)
			{
				MatrixRow mr2(gm2, LoadOnEntry, mr1.Skip());
				Real* el = mr1.Data();                         // pointer to an element
				int n = mr1.Storage();
				mrx.Zero();
				while (n--) { mrx.AddScaled(mr2, *el++); mr2.Next(); }
				mr1.Next(); mrx.Next();
			}
			gmx->ReleaseAndDelete(); gm1->tDelete(); gm2->tDelete(); return gmx;
		}

		static GeneralMatrix* mmMult(GeneralMatrix* gm1, GeneralMatrix* gm2)
		{
			// matrix multiplication for type Matrix only
			REPORT
				Tracer tr("MatrixMult");

			int nr=gm1->Nrows(); int ncr=gm1->Ncols(); int nc=gm2->Ncols();
			if (ncr != gm2->Nrows()) Throw(IncompatibleDimensionsException(*gm1,*gm2));

			Matrix* gm = new Matrix(nr,nc); MatrixErrorNoSpace(gm);

			Real* s1=gm1->Store(); Real* s2=gm2->Store(); Real* s=gm->Store();

			if (ncr)
			{
				while (nr--)
				{
					Real* s2x = s2; int j = ncr;
					Real* sx = s; Real f = *s1++; int k = nc;
					while (k--) *sx++ = f * *s2x++;
					while (--j)
					{ sx = s; f = *s1++; k = nc; while (k--) *sx++ += f * *s2x++; }
					s = sx;
				}
			}
			else *gm = 0.0;

			gm->ReleaseAndDelete(); gm1->tDelete(); gm2->tDelete(); return gm;
		}

		static GeneralMatrix* GeneralMult(GeneralMatrix* gm1, GeneralMatrix* gm2,
			MultipliedMatrix* mm, MatrixType mtx)
		{
			if ( Rectangular(gm1->type(), gm2->type(), mtx))
			{
				REPORT
					return mmMult(gm1, gm2);
			}
			else
			{
				REPORT
					Compare(gm1->type() * gm2->type(),mtx);
				int nr = gm2->Nrows(); int nc = gm2->Ncols();
				if (nc <= 5 && nr > nc) { REPORT return GeneralMult1(gm1, gm2, mm, mtx); }
				else { REPORT return GeneralMult2(gm1, gm2, mm, mtx); }
			}
		}

		static GeneralMatrix* GeneralKP(GeneralMatrix* gm1, GeneralMatrix* gm2,
			KPMatrix* kp, MatrixType mtx)
		{
			REPORT
				Tracer tr("GeneralKP");
			int nr1 = gm1->Nrows(); int nc1 = gm1->Ncols();
			int nr2 = gm2->Nrows(); int nc2 = gm2->Ncols();
			Compare((gm1->type()).KP(gm2->type()),mtx);
			GeneralMatrix* gmx = mtx.New(nr1*nr2, nc1*nc2, kp);
			MatrixRow mrx(gmx, LoadOnEntry+StoreOnExit+DirectPart);
			MatrixRow mr1(gm1, LoadOnEntry);
			for (int i = 1; i <= nr1; ++i)
			{
				MatrixRow mr2(gm2, LoadOnEntry);
				for (int j = 1; j <= nr2; ++j)
				{ mrx.KP(mr1,mr2); mr2.Next(); mrx.Next(); }
				mr1.Next();
			}
			gmx->ReleaseAndDelete(); gm1->tDelete(); gm2->tDelete(); return gmx;
		}

		static GeneralMatrix* GeneralSolv(GeneralMatrix* gm1, GeneralMatrix* gm2,
			BaseMatrix* sm, MatrixType mtx)
		{
			REPORT
				Tracer tr("GeneralSolv");
			Compare(gm1->type().i() * gm2->type(),mtx);
			int nr = gm1->Nrows();
			if (nr != gm1->Ncols()) Throw(NotSquareException(*gm1));
			int nc = gm2->Ncols();
			if (gm1->Ncols() != gm2->Nrows())
				Throw(IncompatibleDimensionsException(*gm1, *gm2));
			GeneralMatrix* gmx = mtx.New(nr,nc,sm); MatrixErrorNoSpace(gmx);
			Real* r = new Real [nr]; MatrixErrorNoSpace(r);
			MONITOR_REAL_NEW("Make   (GenSolv)",nr,r)
				GeneralMatrix* gms = gm1->MakeSolver();
			Try
			{

				MatrixColX mcx(gmx, r, StoreOnExit+DirectPart);   // copy to and from r
				// this must be inside Try so mcx is destroyed before gmx
				MatrixColX mc2(gm2, r, LoadOnEntry);
				while (nc--) { gms->Solver(mcx, mc2); mcx.Next(); mc2.Next(); }
			}
			CatchAll
			{
				if (gms) gms->tDelete();
				delete gmx;                   // <--------------------
				gm2->tDelete();
				MONITOR_REAL_DELETE("Delete (GenSolv)",nr,r)
					// AT&T version 2.1 gives an internal error
					delete [] r;
				ReThrow;
			}
			gms->tDelete(); gmx->ReleaseAndDelete(); gm2->tDelete();
			MONITOR_REAL_DELETE("Delete (GenSolv)",nr,r)
				// AT&T version 2.1 gives an internal error
				delete [] r;
			return gmx;
		}

		// version for inverses - gm2 is identity
		static GeneralMatrix* GeneralSolvI(GeneralMatrix* gm1, BaseMatrix* sm,
			MatrixType mtx)
		{
			REPORT
				Tracer tr("GeneralSolvI");
			Compare(gm1->type().i(),mtx);
			int nr = gm1->Nrows();
			if (nr != gm1->Ncols()) Throw(NotSquareException(*gm1));
			int nc = nr;
			// DiagonalMatrix I(nr); I = 1;
			IdentityMatrix I(nr);
			GeneralMatrix* gmx = mtx.New(nr,nc,sm); MatrixErrorNoSpace(gmx);
			Real* r = new Real [nr]; MatrixErrorNoSpace(r);
			MONITOR_REAL_NEW("Make   (GenSolvI)",nr,r)
				GeneralMatrix* gms = gm1->MakeSolver();
			Try
			{

				MatrixColX mcx(gmx, r, StoreOnExit+DirectPart);   // copy to and from r
				// this must be inside Try so mcx is destroyed before gmx
				MatrixColX mc2(&I, r, LoadOnEntry);
				while (nc--) { gms->Solver(mcx, mc2); mcx.Next(); mc2.Next(); }
			}
			CatchAll
			{
				if (gms) gms->tDelete();
				delete gmx;
				MONITOR_REAL_DELETE("Delete (GenSolvI)",nr,r)
					// AT&T version 2.1 gives an internal error
					delete [] r;
				ReThrow;
			}
			gms->tDelete(); gmx->ReleaseAndDelete();
			MONITOR_REAL_DELETE("Delete (GenSolvI)",nr,r)
				// AT&T version 2.1 gives an internal error
				delete [] r;
			return gmx;
		}

		GeneralMatrix* InvertedMatrix::Evaluate(MatrixType mtx)
		{
			// Matrix Inversion - use solve routines
			Tracer tr("InvertedMatrix::Evaluate");
			REPORT
				gm=((BaseMatrix*&)bm)->Evaluate();
			return GeneralSolvI(gm,this,mtx);
		}

		//*************************** New versions ************************

		GeneralMatrix* AddedMatrix::Evaluate(MatrixType mtd)
		{
			REPORT
				Tracer tr("AddedMatrix::Evaluate");
			gm1=((BaseMatrix*&)bm1)->Evaluate(); gm2=((BaseMatrix*&)bm2)->Evaluate();
			int nr=gm1->Nrows(); int nc=gm1->Ncols();
			if (nr!=gm2->Nrows() || nc!=gm2->Ncols())
			{
				Try { Throw(IncompatibleDimensionsException(*gm1, *gm2)); }
				CatchAll
				{
					gm1->tDelete(); gm2->tDelete();
					ReThrow;
				}
			}
			MatrixType mt1 = gm1->type(), mt2 = gm2->type(); MatrixType mts = mt1 + mt2;
			if (!mtd) { REPORT mtd = mts; }
			else if (!(mtd.DataLossOK || mtd >= mts))
			{
				REPORT
					gm1->tDelete(); gm2->tDelete();
				Throw(ProgramException("Illegal Conversion", mts, mtd));
			}
			GeneralMatrix* gmx;
			bool c1 = (mtd == mt1), c2 = (mtd == mt2);
			if ( c1 && c2 && (gm1->SimpleAddOK(gm2) == 0) )
			{
				if (gm1->reuse()) { REPORT AddTo(gm1,gm2); gm2->tDelete(); gmx = gm1; }
				else if (gm2->reuse()) { REPORT AddTo(gm2,gm1); gmx = gm2; }
				else
				{
					REPORT
						// what if new throws an exception
						Try { gmx = mt1.New(nr,nc,this); }
					CatchAll
					{
						ReThrow;
					}
					gmx->ReleaseAndDelete(); Add(gmx,gm1,gm2);
				}
			}
			else
			{
				if (c1 && c2)
				{
					short SAO = gm1->SimpleAddOK(gm2);
					if (SAO & 1) { REPORT c1 = false; }
					if (SAO & 2) { REPORT c2 = false; }
				}
				if (c1 && gm1->reuse() )               // must have type test first
				{ REPORT AddDS(gm1,gm2); gm2->tDelete(); gmx = gm1; }
				else if (c2 && gm2->reuse() )
				{ REPORT AddDS(gm2,gm1); if (!c1) gm1->tDelete(); gmx = gm2; }
				else
				{
					REPORT
						Try { gmx = mtd.New(nr,nc,this); }
					CatchAll
					{
						if (!c1) gm1->tDelete(); if (!c2) gm2->tDelete();
						ReThrow;
					}
					AddDS(gmx,gm1,gm2);
					if (!c1) gm1->tDelete(); if (!c2) gm2->tDelete();
					gmx->ReleaseAndDelete();
				}
			}
			return gmx;
		}

		GeneralMatrix* SubtractedMatrix::Evaluate(MatrixType mtd)
		{
			REPORT
				Tracer tr("SubtractedMatrix::Evaluate");
			gm1=((BaseMatrix*&)bm1)->Evaluate(); gm2=((BaseMatrix*&)bm2)->Evaluate();
			int nr=gm1->Nrows(); int nc=gm1->Ncols();
			if (nr!=gm2->Nrows() || nc!=gm2->Ncols())
			{
				Try { Throw(IncompatibleDimensionsException(*gm1, *gm2)); }
				CatchAll
				{
					gm1->tDelete(); gm2->tDelete();
					ReThrow;
				}
			}
			MatrixType mt1 = gm1->type(), mt2 = gm2->type(); MatrixType mts = mt1 + mt2;
			if (!mtd) { REPORT mtd = mts; }
			else if (!(mtd.DataLossOK || mtd >= mts))
			{
				gm1->tDelete(); gm2->tDelete();
				Throw(ProgramException("Illegal Conversion", mts, mtd));
			}
			GeneralMatrix* gmx;
			bool c1 = (mtd == mt1), c2 = (mtd == mt2);
			if ( c1 && c2 && (gm1->SimpleAddOK(gm2) == 0) )
			{
				if (gm1->reuse())
				{ REPORT SubtractFrom(gm1,gm2); gm2->tDelete(); gmx = gm1; }
				else if (gm2->reuse()) { REPORT ReverseSubtract(gm2,gm1); gmx = gm2; }
				else
				{
					REPORT
						Try { gmx = mt1.New(nr,nc,this); }
					CatchAll
					{
						ReThrow;
					}
					gmx->ReleaseAndDelete(); Subtract(gmx,gm1,gm2);
				}
			}
			else
			{
				if (c1 && c2)
				{
					short SAO = gm1->SimpleAddOK(gm2);
					if (SAO & 1) { REPORT c1 = false; }
					if (SAO & 2) { REPORT c2 = false; }
				}
				if (c1 && gm1->reuse() )               // must have type test first
				{ REPORT SubtractDS(gm1,gm2); gm2->tDelete(); gmx = gm1; }
				else if (c2 && gm2->reuse() )
				{
					REPORT ReverseSubtractDS(gm2,gm1);
					if (!c1) gm1->tDelete(); gmx = gm2;
				}
				else
				{
					REPORT
						// what if New throws and exception
						Try { gmx = mtd.New(nr,nc,this); }
					CatchAll
					{
						if (!c1) gm1->tDelete(); if (!c2) gm2->tDelete();
						ReThrow;
					}
					SubtractDS(gmx,gm1,gm2);
					if (!c1) gm1->tDelete(); if (!c2) gm2->tDelete();
					gmx->ReleaseAndDelete();
				}
			}
			return gmx;
		}

		GeneralMatrix* SPMatrix::Evaluate(MatrixType mtd)
		{
			REPORT
				Tracer tr("SPMatrix::Evaluate");
			gm1=((BaseMatrix*&)bm1)->Evaluate(); gm2=((BaseMatrix*&)bm2)->Evaluate();
			int nr=gm1->Nrows(); int nc=gm1->Ncols();
			if (nr!=gm2->Nrows() || nc!=gm2->Ncols())
			{
				Try { Throw(IncompatibleDimensionsException(*gm1, *gm2)); }
				CatchAll
				{
					gm1->tDelete(); gm2->tDelete();
					ReThrow;
				}
			}
			MatrixType mt1 = gm1->type(), mt2 = gm2->type();
			MatrixType mts = mt1.SP(mt2);
			if (!mtd) { REPORT mtd = mts; }
			else if (!(mtd.DataLossOK || mtd >= mts))
			{
				gm1->tDelete(); gm2->tDelete();
				Throw(ProgramException("Illegal Conversion", mts, mtd));
			}
			GeneralMatrix* gmx;
			bool c1 = (mtd == mt1), c2 = (mtd == mt2);
			if ( c1 && c2 && (gm1->SimpleAddOK(gm2) == 0) )
			{
				if (gm1->reuse()) { REPORT SP(gm1,gm2); gm2->tDelete(); gmx = gm1; }
				else if (gm2->reuse()) { REPORT SP(gm2,gm1); gmx = gm2; }
				else
				{
					REPORT
						Try { gmx = mt1.New(nr,nc,this); }
					CatchAll
					{
						ReThrow;
					}
					gmx->ReleaseAndDelete(); SP(gmx,gm1,gm2);
				}
			}
			else
			{
				if (c1 && c2)
				{
					short SAO = gm1->SimpleAddOK(gm2);
					if (SAO & 1) { REPORT c2 = false; }    // c1 and c2 swapped
					if (SAO & 2) { REPORT c1 = false; }
				}
				if (c1 && gm1->reuse() )               // must have type test first
				{ REPORT SPDS(gm1,gm2); gm2->tDelete(); gmx = gm1; }
				else if (c2 && gm2->reuse() )
				{ REPORT SPDS(gm2,gm1); if (!c1) gm1->tDelete(); gmx = gm2; }
				else
				{
					REPORT
						// what if New throws and exception
						Try { gmx = mtd.New(nr,nc,this); }
					CatchAll
					{
						if (!c1) gm1->tDelete(); if (!c2) gm2->tDelete();
						ReThrow;
					}
					SPDS(gmx,gm1,gm2);
					if (!c1) gm1->tDelete(); if (!c2) gm2->tDelete();
					gmx->ReleaseAndDelete();
				}
			}
			return gmx;
		}



		//*************************** norm functions ****************************/

		Real BaseMatrix::norm1() const
		{
			// maximum of sum of absolute values of a column
			REPORT
				GeneralMatrix* gm = ((BaseMatrix&)*this).Evaluate();
			int nc = gm->Ncols(); Real value = 0.0;
			MatrixCol mc(gm, LoadOnEntry);
			while (nc--)
			{ Real v = mc.SumAbsoluteValue(); if (value < v) value = v; mc.Next(); }
			gm->tDelete(); return value;
		}

		Real BaseMatrix::norm_infinity() const
		{
			// maximum of sum of absolute values of a row
			REPORT
				GeneralMatrix* gm = ((BaseMatrix&)*this).Evaluate();
			int nr = gm->Nrows(); Real value = 0.0;
			MatrixRow mr(gm, LoadOnEntry);
			while (nr--)
			{ Real v = mr.SumAbsoluteValue(); if (value < v) value = v; mr.Next(); }
			gm->tDelete(); return value;
		}

		//********************** Concatenation and stacking *************************/

		GeneralMatrix* ConcatenatedMatrix::Evaluate(MatrixType mtx)
		{
			REPORT
				Tracer tr("Concatenate");
			gm2 = ((BaseMatrix*&)bm2)->Evaluate();
			gm1 = ((BaseMatrix*&)bm1)->Evaluate();
			Compare(gm1->type() | gm2->type(),mtx);
			int nr=gm1->Nrows(); int nc = gm1->Ncols() + gm2->Ncols();
			if (nr != gm2->Nrows())
				Throw(IncompatibleDimensionsException(*gm1, *gm2));
			GeneralMatrix* gmx = mtx.New(nr,nc,this);
			MatrixRow mr1(gm1, LoadOnEntry); MatrixRow mr2(gm2, LoadOnEntry);
			MatrixRow mr(gmx, StoreOnExit+DirectPart);
			while (nr--) { mr.ConCat(mr1,mr2); mr1.Next(); mr2.Next(); mr.Next(); }
			gmx->ReleaseAndDelete(); gm1->tDelete(); gm2->tDelete(); return gmx;
		}

		GeneralMatrix* StackedMatrix::Evaluate(MatrixType mtx)
		{
			REPORT
				Tracer tr("Stack");
			gm2 = ((BaseMatrix*&)bm2)->Evaluate();
			gm1 = ((BaseMatrix*&)bm1)->Evaluate();
			Compare(gm1->type() & gm2->type(),mtx);
			int nc=gm1->Ncols();
			int nr1 = gm1->Nrows(); int nr2 = gm2->Nrows();
			if (nc != gm2->Ncols())
				Throw(IncompatibleDimensionsException(*gm1, *gm2));
			GeneralMatrix* gmx = mtx.New(nr1+nr2,nc,this);
			MatrixRow mr1(gm1, LoadOnEntry); MatrixRow mr2(gm2, LoadOnEntry);
			MatrixRow mr(gmx, StoreOnExit+DirectPart);
			while (nr1--) { mr.Copy(mr1); mr1.Next(); mr.Next(); }
			while (nr2--) { mr.Copy(mr2); mr2.Next(); mr.Next(); }
			gmx->ReleaseAndDelete(); gm1->tDelete(); gm2->tDelete(); return gmx;
		}

		// ************************* equality of matrices ******************** //

		static bool RealEqual(Real* s1, Real* s2, int n)
		{
			int i = n >> 2;
			while (i--)
			{
				if (*s1++ != *s2++) return false; if (*s1++ != *s2++) return false;
				if (*s1++ != *s2++) return false; if (*s1++ != *s2++) return false;
			}
			i = n & 3; while (i--) if (*s1++ != *s2++) return false;
			return true;
		}

		static bool intEqual(int* s1, int* s2, int n)
		{
			int i = n >> 2;
			while (i--)
			{
				if (*s1++ != *s2++) return false; if (*s1++ != *s2++) return false;
				if (*s1++ != *s2++) return false; if (*s1++ != *s2++) return false;
			}
			i = n & 3; while (i--) if (*s1++ != *s2++) return false;
			return true;
		}


		bool operator==(const BaseMatrix& A, const BaseMatrix& B)
		{
			Tracer tr("BaseMatrix ==");
			REPORT
				GeneralMatrix* gmA = ((BaseMatrix&)A).Evaluate();
			GeneralMatrix* gmB = ((BaseMatrix&)B).Evaluate();

			if (gmA == gmB)                            // same matrix
			{ REPORT gmA->tDelete(); return true; }

			if ( gmA->Nrows() != gmB->Nrows() || gmA->Ncols() != gmB->Ncols() )
				// different dimensions
			{ REPORT gmA->tDelete(); gmB->tDelete(); return false; }

			// check for CroutMatrix or BandLUMatrix
			MatrixType AType = gmA->type(); MatrixType BType = gmB->type();
			if (AType.CannotConvert() || BType.CannotConvert() )
			{
				REPORT
					bool bx = gmA->IsEqual(*gmB);
				gmA->tDelete(); gmB->tDelete();
				return bx;
			}

			// is matrix storage the same
			// will need to modify if further matrix structures are introduced
			if (AType == BType && gmA->bandwidth() == gmB->bandwidth())
			{                                          // compare store
				REPORT
					bool bx = RealEqual(gmA->Store(),gmB->Store(),gmA->Storage());
				gmA->tDelete(); gmB->tDelete();
				return bx;
			}

			// matrix storage different - just subtract
			REPORT  return is_zero(*gmA-*gmB);
		}

		bool operator==(const GeneralMatrix& A, const GeneralMatrix& B)
		{
			Tracer tr("GeneralMatrix ==");
			// May or may not call tDeletes
			REPORT

				if (&A == &B)                              // same matrix
				{ REPORT return true; }

				if ( A.Nrows() != B.Nrows() || A.Ncols() != B.Ncols() )
				{ REPORT return false; }                // different dimensions

				// check for CroutMatrix or BandLUMatrix
				MatrixType AType = A.Type(); MatrixType BType = B.Type();
				if (AType.CannotConvert() || BType.CannotConvert() )
				{ REPORT  return A.IsEqual(B); }

				// is matrix storage the same
				// will need to modify if further matrix structures are introduced
				if (AType == BType && A.bandwidth() == B.bandwidth())
				{ REPORT return RealEqual(A.Store(),B.Store(),A.Storage()); }

				// matrix storage different - just subtract
				REPORT  return is_zero(A-B);
		}

		bool GeneralMatrix::is_zero() const
		{
			REPORT
				Real* s=store; int i = storage >> 2;
			while (i--)
			{
				if (*s++) return false; if (*s++) return false;
				if (*s++) return false; if (*s++) return false;
			}
			i = storage & 3; while (i--) if (*s++) return false;
			return true;
		}

		bool is_zero(const BaseMatrix& A)
		{
			Tracer tr("BaseMatrix::is_zero");
			REPORT
				GeneralMatrix* gm1 = 0; bool bx;
			Try { gm1=((BaseMatrix&)A).Evaluate(); bx = gm1->is_zero(); }
			CatchAll { if (gm1) gm1->tDelete(); ReThrow; }
			gm1->tDelete();
			return bx;
		}

		// IsEqual functions - insist matrices are of same type
		// as well as equal values to be equal

		bool GeneralMatrix::IsEqual(const GeneralMatrix& A) const
		{
			Tracer tr("GeneralMatrix IsEqual");
			if (A.type() != type())                       // not same types
			{ REPORT return false; }
			if (&A == this)                               // same matrix
			{ REPORT  return true; }
			if (A.nrows_val != nrows_val || A.ncols_val != ncols_val)
				// different dimensions
			{ REPORT return false; }
			// is matrix storage the same - compare store
			REPORT
				return RealEqual(A.store,store,storage);
		}

		bool CroutMatrix::IsEqual(const GeneralMatrix& A) const
		{
			Tracer tr("CroutMatrix IsEqual");
			if (A.type() != type())                       // not same types
			{ REPORT return false; }
			if (&A == this)                               // same matrix
			{ REPORT  return true; }
			if (A.nrows_val != nrows_val || A.ncols_val != ncols_val)
				// different dimensions
			{ REPORT return false; }
			// is matrix storage the same - compare store
			REPORT
				return RealEqual(A.store,store,storage)
				&& intEqual(((CroutMatrix&)A).indx, indx, nrows_val);
		}


		bool BandLUMatrix::IsEqual(const GeneralMatrix& A) const
		{
			Tracer tr("BandLUMatrix IsEqual");
			if (A.type() != type())                       // not same types
			{ REPORT  return false; }
			if (&A == this)                               // same matrix
			{ REPORT  return true; }
			if ( A.Nrows() != nrows_val || A.Ncols() != ncols_val
				|| ((BandLUMatrix&)A).m1 != m1 || ((BandLUMatrix&)A).m2 != m2 )
				// different dimensions
			{ REPORT  return false; }

			// matrix storage the same - compare store
			REPORT
				return RealEqual(A.Store(),store,storage)
				&& RealEqual(((BandLUMatrix&)A).store2,store2,storage2)
				&& intEqual(((BandLUMatrix&)A).indx, indx, nrows_val);
		}


		// ************************* cross products ******************** //

		inline void crossproduct_body(Real* a, Real* b, Real* c)
		{
			c[0] = a[1] * b[2] - a[2] * b[1];
			c[1] = a[2] * b[0] - a[0] * b[2];
			c[2] = a[0] * b[1] - a[1] * b[0];
		}

		Matrix crossproduct(const Matrix& A, const Matrix& B)
		{
			REPORT
				int ac = A.Ncols(); int ar = A.Nrows();
			int bc = B.Ncols(); int br = B.Nrows();
			Real* a = A.Store(); Real* b = B.Store();
			if (ac == 3)
			{
				if (bc != 3 || ar != 1 || br != 1)
				{ Tracer et("crossproduct"); IncompatibleDimensionsException(A, B); }
				REPORT
					RowVector C(3);  Real* c = C.Store(); crossproduct_body(a, b, c);
				return (Matrix&)C;
			}
			else
			{
				if (ac != 1 || bc != 1 || ar != 3 || br != 3)
				{ Tracer et("crossproduct"); IncompatibleDimensionsException(A, B); }
				REPORT
					ColumnVector C(3);  Real* c = C.Store(); crossproduct_body(a, b, c);
				return (Matrix&)C;
			}
		}

		ReturnMatrix crossproduct_rows(const Matrix& A, const Matrix& B)
		{
			REPORT
				int n = A.Nrows();
			if (A.Ncols() != 3 || B.Ncols() != 3 || n != B.Nrows())
			{
				Tracer et("crossproduct_rows"); IncompatibleDimensionsException(A, B);
			}
			Matrix C(n, 3);
			Real* a = A.Store(); Real* b = B.Store(); Real* c = C.Store();
			if (n--)
			{
				for (;;)
				{
					crossproduct_body(a, b, c);
					if (!(n--)) break;
					a += 3; b += 3; c += 3;
				}
			}

			return C.ForReturn();
		}

		ReturnMatrix crossproduct_columns(const Matrix& A, const Matrix& B)
		{
			REPORT
				int n = A.Ncols();
			if (A.Nrows() != 3 || B.Nrows() != 3 || n != B.Ncols())
			{
				Tracer et("crossproduct_columns");
				IncompatibleDimensionsException(A, B);
			}
			Matrix C(3, n);
			Real* a = A.Store(); Real* b = B.Store(); Real* c = C.Store();
			Real* an = a+n; Real* an2 = an+n;
			Real* bn = b+n; Real* bn2 = bn+n;
			Real* cn = c+n; Real* cn2 = cn+n;

			int i = n; 
			while (i--)
			{
				*c++   = *an    * *bn2   - *an2   * *bn;
				*cn++  = *an2++ * *b     - *a     * *bn2++;
				*cn2++ = *a++   * *bn++  - *an++  * *b++;
			}

			return C.ForReturn();
		}


#ifdef use_namespace
	}
#endif


	//$$ newmat6.cpp            Operators, element access, submatrices

	// Copyright (C) 1991,2,3,4: R B Davies

//#include "include.h"

#include "newmat.h"
//#include "newmatrc.h"

#ifdef use_namespace
	namespace NEWMAT {
#endif



#ifdef DO_REPORT
#define REPORT { static ExeCounter ExeCount(__LINE__,6); ++ExeCount; }
#else
#define REPORT {}
#endif

		/*************************** general utilities *************************/
/*
		static int tristore(int n)                      // els in triangular matrix
		{ return (n*(n+1))/2; }

*/
		/****************************** operators *******************************/

		Real& Matrix::operator()(int m, int n)
		{
			REPORT
				if (m<=0 || m>nrows_val || n<=0 || n>ncols_val)
					Throw(IndexException(m,n,*this));
			return store[(m-1)*ncols_val+n-1];
		}

		Real& SymmetricMatrix::operator()(int m, int n)
		{
			REPORT
				if (m<=0 || n<=0 || m>nrows_val || n>ncols_val)
					Throw(IndexException(m,n,*this));
			if (m>=n) return store[tristore(m-1)+n-1];
			else return store[tristore(n-1)+m-1];
		}

		Real& UpperTriangularMatrix::operator()(int m, int n)
		{
			REPORT
				if (m<=0 || n<m || n>ncols_val)
					Throw(IndexException(m,n,*this));
			return store[(m-1)*ncols_val+n-1-tristore(m-1)];
		}

		Real& LowerTriangularMatrix::operator()(int m, int n)
		{
			REPORT
				if (n<=0 || m<n || m>nrows_val)
					Throw(IndexException(m,n,*this));
			return store[tristore(m-1)+n-1];
		}

		Real& DiagonalMatrix::operator()(int m, int n)
		{
			REPORT
				if (n<=0 || m!=n || m>nrows_val || n>ncols_val)
					Throw(IndexException(m,n,*this));
			return store[n-1];
		}

		Real& DiagonalMatrix::operator()(int m)
		{
			REPORT
				if (m<=0 || m>nrows_val) Throw(IndexException(m,*this));
			return store[m-1];
		}

		Real& ColumnVector::operator()(int m)
		{
			REPORT
				if (m<=0 || m> nrows_val) Throw(IndexException(m,*this));
			return store[m-1];
		}

		Real& RowVector::operator()(int n)
		{
			REPORT
				if (n<=0 || n> ncols_val) Throw(IndexException(n,*this));
			return store[n-1];
		}

		Real& BandMatrix::operator()(int m, int n)
		{
			REPORT
				int w = upper_val+lower_val+1; int i = lower_val+n-m;
			if (m<=0 || m>nrows_val || n<=0 || n>ncols_val || i<0 || i>=w)
				Throw(IndexException(m,n,*this));
			return store[w*(m-1)+i];
		}

		Real& UpperBandMatrix::operator()(int m, int n)
		{
			REPORT
				int w = upper_val+1; int i = n-m;
			if (m<=0 || m>nrows_val || n<=0 || n>ncols_val || i<0 || i>=w)
				Throw(IndexException(m,n,*this));
			return store[w*(m-1)+i];
		}

		Real& LowerBandMatrix::operator()(int m, int n)
		{
			REPORT
				int w = lower_val+1; int i = lower_val+n-m;
			if (m<=0 || m>nrows_val || n<=0 || n>ncols_val || i<0 || i>=w)
				Throw(IndexException(m,n,*this));
			return store[w*(m-1)+i];
		}

		Real& SymmetricBandMatrix::operator()(int m, int n)
		{
			REPORT
				int w = lower_val+1;
			if (m>=n)
			{
				REPORT
					int i = lower_val+n-m;
				if ( m>nrows_val || n<=0 || i<0 )
					Throw(IndexException(m,n,*this));
				return store[w*(m-1)+i];
			}
			else
			{
				REPORT
					int i = lower_val+m-n;
				if ( n>nrows_val || m<=0 || i<0 )
					Throw(IndexException(m,n,*this));
				return store[w*(n-1)+i];
			}
		}


		Real Matrix::operator()(int m, int n) const
		{
			REPORT
				if (m<=0 || m>nrows_val || n<=0 || n>ncols_val)
					Throw(IndexException(m,n,*this));
			return store[(m-1)*ncols_val+n-1];
		}

		Real SymmetricMatrix::operator()(int m, int n) const
		{
			REPORT
				if (m<=0 || n<=0 || m>nrows_val || n>ncols_val)
					Throw(IndexException(m,n,*this));
			if (m>=n) return store[tristore(m-1)+n-1];
			else return store[tristore(n-1)+m-1];
		}

		Real UpperTriangularMatrix::operator()(int m, int n) const
		{
			REPORT
				if (m<=0 || n<m || n>ncols_val)
					Throw(IndexException(m,n,*this));
			return store[(m-1)*ncols_val+n-1-tristore(m-1)];
		}

		Real LowerTriangularMatrix::operator()(int m, int n) const
		{
			REPORT
				if (n<=0 || m<n || m>nrows_val)
					Throw(IndexException(m,n,*this));
			return store[tristore(m-1)+n-1];
		}

		Real DiagonalMatrix::operator()(int m, int n) const
		{
			REPORT
				if (n<=0 || m!=n || m>nrows_val || n>ncols_val)
					Throw(IndexException(m,n,*this));
			return store[n-1];
		}

		Real DiagonalMatrix::operator()(int m) const
		{
			REPORT
				if (m<=0 || m>nrows_val) Throw(IndexException(m,*this));
			return store[m-1];
		}

		Real ColumnVector::operator()(int m) const
		{
			REPORT
				if (m<=0 || m> nrows_val) Throw(IndexException(m,*this));
			return store[m-1];
		}

		Real RowVector::operator()(int n) const
		{
			REPORT
				if (n<=0 || n> ncols_val) Throw(IndexException(n,*this));
			return store[n-1];
		}

		Real BandMatrix::operator()(int m, int n) const
		{
			REPORT
				int w = upper_val+lower_val+1; int i = lower_val+n-m;
			if (m<=0 || m>nrows_val || n<=0 || n>ncols_val || i<0 || i>=w)
				Throw(IndexException(m,n,*this));
			return store[w*(m-1)+i];
		}

		Real UpperBandMatrix::operator()(int m, int n) const
		{
			REPORT
				int w = upper_val+1; int i = n-m;
			if (m<=0 || m>nrows_val || n<=0 || n>ncols_val || i<0 || i>=w)
				Throw(IndexException(m,n,*this));
			return store[w*(m-1)+i];
		}

		Real LowerBandMatrix::operator()(int m, int n) const
		{
			REPORT
				int w = lower_val+1; int i = lower_val+n-m;
			if (m<=0 || m>nrows_val || n<=0 || n>ncols_val || i<0 || i>=w)
				Throw(IndexException(m,n,*this));
			return store[w*(m-1)+i];
		}

		Real SymmetricBandMatrix::operator()(int m, int n) const
		{
			REPORT
				int w = lower_val+1;
			if (m>=n)
			{
				REPORT
					int i = lower_val+n-m;
				if ( m>nrows_val || n<=0 || i<0 )
					Throw(IndexException(m,n,*this));
				return store[w*(m-1)+i];
			}
			else
			{
				REPORT
					int i = lower_val+m-n;
				if ( n>nrows_val || m<=0 || i<0 )
					Throw(IndexException(m,n,*this));
				return store[w*(n-1)+i];
			}
		}


		Real BaseMatrix::as_scalar() const
		{
			REPORT
				GeneralMatrix* gm = ((BaseMatrix&)*this).Evaluate();

			if (gm->nrows_val!=1 || gm->ncols_val!=1)
			{
				Tracer tr("as_scalar");
				Try
				{ Throw(ProgramException("Cannot convert to scalar", *gm)); }
				CatchAll { gm->tDelete(); ReThrow; }
			}

			Real x = *(gm->store); gm->tDelete(); return x;
		}


		AddedMatrix BaseMatrix::operator+(const BaseMatrix& bm) const
		{ REPORT return AddedMatrix(this, &bm); }

		SPMatrix SP(const BaseMatrix& bm1,const BaseMatrix& bm2)
		{ REPORT return SPMatrix(&bm1, &bm2); }

		KPMatrix KP(const BaseMatrix& bm1,const BaseMatrix& bm2)
		{ REPORT return KPMatrix(&bm1, &bm2); }

		MultipliedMatrix BaseMatrix::operator*(const BaseMatrix& bm) const
		{ REPORT return MultipliedMatrix(this, &bm); }

		ConcatenatedMatrix BaseMatrix::operator|(const BaseMatrix& bm) const
		{ REPORT return ConcatenatedMatrix(this, &bm); }

		StackedMatrix BaseMatrix::operator&(const BaseMatrix& bm) const
		{ REPORT return StackedMatrix(this, &bm); }

		SolvedMatrix InvertedMatrix::operator*(const BaseMatrix& bmx) const
		{ REPORT return SolvedMatrix(bm, &bmx); }

		SubtractedMatrix BaseMatrix::operator-(const BaseMatrix& bm) const
		{ REPORT return SubtractedMatrix(this, &bm); }

		ShiftedMatrix BaseMatrix::operator+(Real f) const
		{ REPORT return ShiftedMatrix(this, f); }

		ShiftedMatrix operator+(Real f, const BaseMatrix& BM)
		{ REPORT return ShiftedMatrix(&BM, f); }

		NegShiftedMatrix operator-(Real f, const BaseMatrix& bm)
		{ REPORT return NegShiftedMatrix(f, &bm); }

		ScaledMatrix BaseMatrix::operator*(Real f) const
		{ REPORT return ScaledMatrix(this, f); }

		ScaledMatrix BaseMatrix::operator/(Real f) const
		{ REPORT return ScaledMatrix(this, 1.0/f); }

		ScaledMatrix operator*(Real f, const BaseMatrix& BM)
		{ REPORT return ScaledMatrix(&BM, f); }

		ShiftedMatrix BaseMatrix::operator-(Real f) const
		{ REPORT return ShiftedMatrix(this, -f); }

		TransposedMatrix BaseMatrix::t() const
		{ REPORT return TransposedMatrix(this); }

		NegatedMatrix BaseMatrix::operator-() const
		{ REPORT return NegatedMatrix(this); }

		ReversedMatrix BaseMatrix::reverse() const
		{ REPORT return ReversedMatrix(this); }

		InvertedMatrix BaseMatrix::i() const
		{ REPORT return InvertedMatrix(this); }


		RowedMatrix BaseMatrix::as_row() const
		{ REPORT return RowedMatrix(this); }

		ColedMatrix BaseMatrix::as_column() const
		{ REPORT return ColedMatrix(this); }

		DiagedMatrix BaseMatrix::as_diagonal() const
		{ REPORT return DiagedMatrix(this); }

		MatedMatrix BaseMatrix::as_matrix(int nrx, int ncx) const
		{ REPORT return MatedMatrix(this,nrx,ncx); }


		void GeneralMatrix::operator=(Real f)
		{ REPORT int i=storage; Real* s=store; while (i--) { *s++ = f; } }

		void Matrix::operator=(const BaseMatrix& X)
		{
			REPORT //CheckConversion(X);
				// MatrixConversionCheck mcc;
				Eq(X,MatrixType::Rt);
		} 

		void SquareMatrix::operator=(const BaseMatrix& X)
		{
			REPORT //CheckConversion(X);
				// MatrixConversionCheck mcc;
				Eq(X,MatrixType::Rt);
			if (nrows_val != ncols_val)
			{ Tracer tr("SquareMatrix(=)"); Throw(NotSquareException(*this)); }
		}

		void SquareMatrix::operator=(const Matrix& m)
		{
			REPORT
				if (m.nrows_val != m.ncols_val)
				{ Tracer tr("SquareMatrix(=Matrix)"); Throw(NotSquareException(*this)); }
				Eq(m);
		}

		void RowVector::operator=(const BaseMatrix& X)
		{
			REPORT  // CheckConversion(X);
				// MatrixConversionCheck mcc;
				Eq(X,MatrixType::RV);
			if (nrows_val!=1)
			{ Tracer tr("RowVector(=)"); Throw(VectorException(*this)); }
		}

		void ColumnVector::operator=(const BaseMatrix& X)
		{
			REPORT //CheckConversion(X);
				// MatrixConversionCheck mcc;
				Eq(X,MatrixType::CV);
			if (ncols_val!=1)
			{ Tracer tr("ColumnVector(=)"); Throw(VectorException(*this)); }
		}

		void SymmetricMatrix::operator=(const BaseMatrix& X)
		{
			REPORT // CheckConversion(X);
				// MatrixConversionCheck mcc;
				Eq(X,MatrixType::Sm);
		}

		void UpperTriangularMatrix::operator=(const BaseMatrix& X)
		{
			REPORT //CheckConversion(X);
				// MatrixConversionCheck mcc;
				Eq(X,MatrixType::UT);
		}

		void LowerTriangularMatrix::operator=(const BaseMatrix& X)
		{
			REPORT //CheckConversion(X);
				// MatrixConversionCheck mcc;
				Eq(X,MatrixType::LT);
		}

		void DiagonalMatrix::operator=(const BaseMatrix& X)
		{
			REPORT // CheckConversion(X);
				// MatrixConversionCheck mcc;
				Eq(X,MatrixType::Dg);
		}

		void IdentityMatrix::operator=(const BaseMatrix& X)
		{
			REPORT // CheckConversion(X);
				// MatrixConversionCheck mcc;
				Eq(X,MatrixType::Id);
		}


		void CroutMatrix::operator=(const CroutMatrix& gm)
		{
			if (&gm == this) { REPORT tag_val = -1; return; }
			REPORT
				if (indx > 0) { delete [] indx; indx = 0; }
				((CroutMatrix&)gm).get_aux(*this);
				Eq(gm);
		}





		void GeneralMatrix::operator<<(const double* r)
		{
			REPORT
				int i = storage; Real* s=store;
			while(i--) *s++ = (Real)*r++;
		}


		void GeneralMatrix::operator<<(const float* r)
		{
			REPORT
				int i = storage; Real* s=store;
			while(i--) *s++ = (Real)*r++;
		}


		void GeneralMatrix::operator<<(const int* r)
		{
			REPORT
				int i = storage; Real* s=store;
			while(i--) *s++ = (Real)*r++;
		}


		void GenericMatrix::operator=(const GenericMatrix& bmx)
		{
			if (&bmx != this) { REPORT if (gm) delete gm; gm = bmx.gm->Image();}
			else { REPORT }
			gm->Protect();
		}

		void GenericMatrix::operator=(const BaseMatrix& bmx)
		{
			if (gm)
			{
				int counter=bmx.search(gm);
				if (counter==0) { REPORT delete gm; gm=0; }
				else { REPORT gm->Release(counter); }
			}
			else { REPORT }
			GeneralMatrix* gmx = ((BaseMatrix&)bmx).Evaluate();
			if (gmx != gm) { REPORT if (gm) delete gm; gm = gmx->Image(); }
			else { REPORT }
			gm->Protect();
		}


		/*************************** += etc ***************************************/

		// will also need versions for SubMatrix



		// GeneralMatrix operators

		void GeneralMatrix::operator+=(const BaseMatrix& X)
		{
			REPORT
				Tracer tr("GeneralMatrix::operator+=");
			// MatrixConversionCheck mcc;
			Protect();                                   // so it cannot get deleted
			// during Evaluate
			GeneralMatrix* gm = ((BaseMatrix&)X).Evaluate();
			AddedMatrix am(this,gm);
			if (gm==this) Release(2); else Release();
			Eq2(am,type());
		}

		void GeneralMatrix::operator-=(const BaseMatrix& X)
		{
			REPORT
				Tracer tr("GeneralMatrix::operator-=");
			// MatrixConversionCheck mcc;
			Protect();                                   // so it cannot get deleted
			// during Evaluate
			GeneralMatrix* gm = ((BaseMatrix&)X).Evaluate();
			SubtractedMatrix am(this,gm);
			if (gm==this) Release(2); else Release();
			Eq2(am,type());
		}

		void GeneralMatrix::operator*=(const BaseMatrix& X)
		{
			REPORT
				Tracer tr("GeneralMatrix::operator*=");
			// MatrixConversionCheck mcc;
			Protect();                                   // so it cannot get deleted
			// during Evaluate
			GeneralMatrix* gm = ((BaseMatrix&)X).Evaluate();
			MultipliedMatrix am(this,gm);
			if (gm==this) Release(2); else Release();
			Eq2(am,type());
		}

		void GeneralMatrix::operator|=(const BaseMatrix& X)
		{
			REPORT
				Tracer tr("GeneralMatrix::operator|=");
			// MatrixConversionCheck mcc;
			Protect();                                   // so it cannot get deleted
			// during Evaluate
			GeneralMatrix* gm = ((BaseMatrix&)X).Evaluate();
			ConcatenatedMatrix am(this,gm);
			if (gm==this) Release(2); else Release();
			Eq2(am,type());
		}

		void GeneralMatrix::operator&=(const BaseMatrix& X)
		{
			REPORT
				Tracer tr("GeneralMatrix::operator&=");
			// MatrixConversionCheck mcc;
			Protect();                                   // so it cannot get deleted
			// during Evaluate
			GeneralMatrix* gm = ((BaseMatrix&)X).Evaluate();
			StackedMatrix am(this,gm);
			if (gm==this) Release(2); else Release();
			Eq2(am,type());
		}

		void GeneralMatrix::operator+=(Real r)
		{
			REPORT
				Tracer tr("GeneralMatrix::operator+=(Real)");
			// MatrixConversionCheck mcc;
			ShiftedMatrix am(this,r);
			Release(); Eq2(am,type());
		}

		void GeneralMatrix::operator*=(Real r)
		{
			REPORT
				Tracer tr("GeneralMatrix::operator*=(Real)");
			// MatrixConversionCheck mcc;
			ScaledMatrix am(this,r);
			Release(); Eq2(am,type());
		}


		// Generic matrix operators

		void GenericMatrix::operator+=(const BaseMatrix& X)
		{
			REPORT
				Tracer tr("GenericMatrix::operator+=");
			if (!gm) Throw(ProgramException("GenericMatrix is null"));
			gm->Protect();            // so it cannot get deleted during Evaluate
			GeneralMatrix* gmx = ((BaseMatrix&)X).Evaluate();
			AddedMatrix am(gm,gmx);
			if (gmx==gm) gm->Release(2); else gm->Release();
			GeneralMatrix* gmy = am.Evaluate();
			if (gmy != gm) { REPORT delete gm; gm = gmy->Image(); }
			else { REPORT }
			gm->Protect();
		}

		void GenericMatrix::operator-=(const BaseMatrix& X)
		{
			REPORT
				Tracer tr("GenericMatrix::operator-=");
			if (!gm) Throw(ProgramException("GenericMatrix is null"));
			gm->Protect();            // so it cannot get deleted during Evaluate
			GeneralMatrix* gmx = ((BaseMatrix&)X).Evaluate();
			SubtractedMatrix am(gm,gmx);
			if (gmx==gm) gm->Release(2); else gm->Release();
			GeneralMatrix* gmy = am.Evaluate();
			if (gmy != gm) { REPORT delete gm; gm = gmy->Image(); }
			else { REPORT }
			gm->Protect();
		}

		void GenericMatrix::operator*=(const BaseMatrix& X)
		{
			REPORT
				Tracer tr("GenericMatrix::operator*=");
			if (!gm) Throw(ProgramException("GenericMatrix is null"));
			gm->Protect();            // so it cannot get deleted during Evaluate
			GeneralMatrix* gmx = ((BaseMatrix&)X).Evaluate();
			MultipliedMatrix am(gm,gmx);
			if (gmx==gm) gm->Release(2); else gm->Release();
			GeneralMatrix* gmy = am.Evaluate();
			if (gmy != gm) { REPORT delete gm; gm = gmy->Image(); }
			else { REPORT }
			gm->Protect();
		}

		void GenericMatrix::operator|=(const BaseMatrix& X)
		{
			REPORT
				Tracer tr("GenericMatrix::operator|=");
			if (!gm) Throw(ProgramException("GenericMatrix is null"));
			gm->Protect();            // so it cannot get deleted during Evaluate
			GeneralMatrix* gmx = ((BaseMatrix&)X).Evaluate();
			ConcatenatedMatrix am(gm,gmx);
			if (gmx==gm) gm->Release(2); else gm->Release();
			GeneralMatrix* gmy = am.Evaluate();
			if (gmy != gm) { REPORT delete gm; gm = gmy->Image(); }
			else { REPORT }
			gm->Protect();
		}

		void GenericMatrix::operator&=(const BaseMatrix& X)
		{
			REPORT
				Tracer tr("GenericMatrix::operator&=");
			if (!gm) Throw(ProgramException("GenericMatrix is null"));
			gm->Protect();            // so it cannot get deleted during Evaluate
			GeneralMatrix* gmx = ((BaseMatrix&)X).Evaluate();
			StackedMatrix am(gm,gmx);
			if (gmx==gm) gm->Release(2); else gm->Release();
			GeneralMatrix* gmy = am.Evaluate();
			if (gmy != gm) { REPORT delete gm; gm = gmy->Image(); }
			else { REPORT }
			gm->Protect();
		}

		void GenericMatrix::operator+=(Real r)
		{
			REPORT
				Tracer tr("GenericMatrix::operator+= (Real)");
			if (!gm) Throw(ProgramException("GenericMatrix is null"));
			ShiftedMatrix am(gm,r);
			gm->Release();
			GeneralMatrix* gmy = am.Evaluate();
			if (gmy != gm) { REPORT delete gm; gm = gmy->Image(); }
			else { REPORT }
			gm->Protect();
		}

		void GenericMatrix::operator*=(Real r)
		{
			REPORT
				Tracer tr("GenericMatrix::operator*= (Real)");
			if (!gm) Throw(ProgramException("GenericMatrix is null"));
			ScaledMatrix am(gm,r);
			gm->Release();
			GeneralMatrix* gmy = am.Evaluate();
			if (gmy != gm) { REPORT delete gm; gm = gmy->Image(); }
			else { REPORT }
			gm->Protect();
		}


		/************************* element access *********************************/

		Real& Matrix::element(int m, int n)
		{
			REPORT
				if (m<0 || m>= nrows_val || n<0 || n>= ncols_val)
					Throw(IndexException(m,n,*this,true));
			return store[m*ncols_val+n];
		}

		Real Matrix::element(int m, int n) const
		{
			REPORT
				if (m<0 || m>= nrows_val || n<0 || n>= ncols_val)
					Throw(IndexException(m,n,*this,true));
			return store[m*ncols_val+n];
		}

		Real& SymmetricMatrix::element(int m, int n)
		{
			REPORT
				if (m<0 || n<0 || m >= nrows_val || n>=ncols_val)
					Throw(IndexException(m,n,*this,true));
			if (m>=n) return store[tristore(m)+n];
			else return store[tristore(n)+m];
		}

		Real SymmetricMatrix::element(int m, int n) const
		{
			REPORT
				if (m<0 || n<0 || m >= nrows_val || n>=ncols_val)
					Throw(IndexException(m,n,*this,true));
			if (m>=n) return store[tristore(m)+n];
			else return store[tristore(n)+m];
		}

		Real& UpperTriangularMatrix::element(int m, int n)
		{
			REPORT
				if (m<0 || n<m || n>=ncols_val)
					Throw(IndexException(m,n,*this,true));
			return store[m*ncols_val+n-tristore(m)];
		}

		Real UpperTriangularMatrix::element(int m, int n) const
		{
			REPORT
				if (m<0 || n<m || n>=ncols_val)
					Throw(IndexException(m,n,*this,true));
			return store[m*ncols_val+n-tristore(m)];
		}

		Real& LowerTriangularMatrix::element(int m, int n)
		{
			REPORT
				if (n<0 || m<n || m>=nrows_val)
					Throw(IndexException(m,n,*this,true));
			return store[tristore(m)+n];
		}

		Real LowerTriangularMatrix::element(int m, int n) const
		{
			REPORT
				if (n<0 || m<n || m>=nrows_val)
					Throw(IndexException(m,n,*this,true));
			return store[tristore(m)+n];
		}

		Real& DiagonalMatrix::element(int m, int n)
		{
			REPORT
				if (n<0 || m!=n || m>=nrows_val || n>=ncols_val)
					Throw(IndexException(m,n,*this,true));
			return store[n];
		}

		Real DiagonalMatrix::element(int m, int n) const
		{
			REPORT
				if (n<0 || m!=n || m>=nrows_val || n>=ncols_val)
					Throw(IndexException(m,n,*this,true));
			return store[n];
		}

		Real& DiagonalMatrix::element(int m)
		{
			REPORT
				if (m<0 || m>=nrows_val) Throw(IndexException(m,*this,true));
			return store[m];
		}

		Real DiagonalMatrix::element(int m) const
		{
			REPORT
				if (m<0 || m>=nrows_val) Throw(IndexException(m,*this,true));
			return store[m];
		}

		Real& ColumnVector::element(int m)
		{
			REPORT
				if (m<0 || m>= nrows_val) Throw(IndexException(m,*this,true));
			return store[m];
		}

		Real ColumnVector::element(int m) const
		{
			REPORT
				if (m<0 || m>= nrows_val) Throw(IndexException(m,*this,true));
			return store[m];
		}

		Real& RowVector::element(int n)
		{
			REPORT
				if (n<0 || n>= ncols_val)  Throw(IndexException(n,*this,true));
			return store[n];
		}

		Real RowVector::element(int n) const
		{
			REPORT
				if (n<0 || n>= ncols_val)  Throw(IndexException(n,*this,true));
			return store[n];
		}

		Real& BandMatrix::element(int m, int n)
		{
			REPORT
				int w = upper_val+lower_val+1; int i = lower_val+n-m;
			if (m<0 || m>= nrows_val || n<0 || n>= ncols_val || i<0 || i>=w)
				Throw(IndexException(m,n,*this,true));
			return store[w*m+i];
		}

		Real BandMatrix::element(int m, int n) const
		{
			REPORT
				int w = upper_val+lower_val+1; int i = lower_val+n-m;
			if (m<0 || m>= nrows_val || n<0 || n>= ncols_val || i<0 || i>=w)
				Throw(IndexException(m,n,*this,true));
			return store[w*m+i];
		}

		Real& UpperBandMatrix::element(int m, int n)
		{
			REPORT
				int w = upper_val+1; int i = n-m;
			if (m<0 || m>= nrows_val || n<0 || n>= ncols_val || i<0 || i>=w)
				Throw(IndexException(m,n,*this,true));
			return store[w*m+i];
		}

		Real UpperBandMatrix::element(int m, int n) const
		{
			REPORT
				int w = upper_val+1; int i = n-m;
			if (m<0 || m>= nrows_val || n<0 || n>= ncols_val || i<0 || i>=w)
				Throw(IndexException(m,n,*this,true));
			return store[w*m+i];
		}

		Real& LowerBandMatrix::element(int m, int n)
		{
			REPORT
				int w = lower_val+1; int i = lower_val+n-m;
			if (m<0 || m>= nrows_val || n<0 || n>= ncols_val || i<0 || i>=w)
				Throw(IndexException(m,n,*this,true));
			return store[w*m+i];
		}

		Real LowerBandMatrix::element(int m, int n) const
		{
			REPORT
				int w = lower_val+1; int i = lower_val+n-m;
			if (m<0 || m>= nrows_val || n<0 || n>= ncols_val || i<0 || i>=w)
				Throw(IndexException(m,n,*this,true));
			return store[w*m+i];
		}

		Real& SymmetricBandMatrix::element(int m, int n)
		{
			REPORT
				int w = lower_val+1;
			if (m>=n)
			{
				REPORT
					int i = lower_val+n-m;
				if ( m>=nrows_val || n<0 || i<0 )
					Throw(IndexException(m,n,*this,true));
				return store[w*m+i];
			}
			else
			{
				REPORT
					int i = lower_val+m-n;
				if ( n>=nrows_val || m<0 || i<0 )
					Throw(IndexException(m,n,*this,true));
				return store[w*n+i];
			}
		}

		Real SymmetricBandMatrix::element(int m, int n) const
		{
			REPORT
				int w = lower_val+1;
			if (m>=n)
			{
				REPORT
					int i = lower_val+n-m;
				if ( m>=nrows_val || n<0 || i<0 )
					Throw(IndexException(m,n,*this,true));
				return store[w*m+i];
			}
			else
			{
				REPORT
					int i = lower_val+m-n;
				if ( n>=nrows_val || m<0 || i<0 )
					Throw(IndexException(m,n,*this,true));
				return store[w*n+i];
			}
		}

#ifdef use_namespace
	}
#endif

	//$$ newmat5.cpp         Transpose, evaluate etc

	// Copyright (C) 1991,2,3,4: R B Davies

	//#define WANT_STREAM

//#include "include.h"

#include "newmat.h"
//#include "newmatrc.h"

#ifdef use_namespace
	namespace NEWMAT {
#endif


#ifdef DO_REPORT
#define REPORT { static ExeCounter ExeCount(__LINE__,5); ++ExeCount; }
#else
#define REPORT {}
#endif


		/************************ carry out operations ******************************/


		GeneralMatrix* GeneralMatrix::Transpose(TransposedMatrix* tm, MatrixType mt)
		{
			GeneralMatrix* gm1;

			if (Compare(Type().t(),mt))
			{
				REPORT
					gm1 = mt.New(ncols_val,nrows_val,tm);
				for (int i=0; i<ncols_val; i++)
				{
					MatrixRow mr(gm1, StoreOnExit+DirectPart, i);
					MatrixCol mc(this, mr.Data(), LoadOnEntry, i);
				}
			}
			else
			{
				REPORT
					gm1 = mt.New(ncols_val,nrows_val,tm);
				MatrixRow mr(this, LoadOnEntry);
				MatrixCol mc(gm1, StoreOnExit+DirectPart);
				int i = nrows_val;
				while (i--) { mc.Copy(mr); mr.Next(); mc.Next(); }
			}
			tDelete(); gm1->ReleaseAndDelete(); return gm1;
		}

		GeneralMatrix* SymmetricMatrix::Transpose(TransposedMatrix*, MatrixType mt)
		{ REPORT  return Evaluate(mt); }


		GeneralMatrix* DiagonalMatrix::Transpose(TransposedMatrix*, MatrixType mt)
		{ REPORT return Evaluate(mt); }

		GeneralMatrix* ColumnVector::Transpose(TransposedMatrix*, MatrixType mt)
		{
			REPORT
				GeneralMatrix* gmx = new RowVector; MatrixErrorNoSpace(gmx);
			gmx->nrows_val = 1; gmx->ncols_val = gmx->storage = storage;
			return BorrowStore(gmx,mt);
		}

		GeneralMatrix* RowVector::Transpose(TransposedMatrix*, MatrixType mt)
		{
			REPORT
				GeneralMatrix* gmx = new ColumnVector; MatrixErrorNoSpace(gmx);
			gmx->ncols_val = 1; gmx->nrows_val = gmx->storage = storage;
			return BorrowStore(gmx,mt);
		}

		GeneralMatrix* IdentityMatrix::Transpose(TransposedMatrix*, MatrixType mt)
		{ REPORT return Evaluate(mt); }

		GeneralMatrix* GeneralMatrix::Evaluate(MatrixType mt)
		{
			if (Compare(this->Type(),mt)) { REPORT return this; }
			REPORT
				GeneralMatrix* gmx = mt.New(nrows_val,ncols_val,this);
			MatrixRow mr(this, LoadOnEntry);
			MatrixRow mrx(gmx, StoreOnExit+DirectPart);
			int i=nrows_val;
			while (i--) { mrx.Copy(mr); mrx.Next(); mr.Next(); }
			tDelete(); gmx->ReleaseAndDelete(); return gmx;
		}

		GeneralMatrix* CroutMatrix::Evaluate(MatrixType mt)
		{
			if (Compare(this->Type(),mt)) { REPORT return this; }
			REPORT
				Tracer et("CroutMatrix::Evaluate");
			bool dummy = true;
			if (dummy) Throw(ProgramException("Illegal use of CroutMatrix", *this));
			return this;
		}

		GeneralMatrix* GenericMatrix::Evaluate(MatrixType mt)
		{ REPORT  return gm->Evaluate(mt); }

		GeneralMatrix* ShiftedMatrix::Evaluate(MatrixType mt)
		{
			gm=((BaseMatrix*&)bm)->Evaluate();
			int nr=gm->Nrows(); int nc=gm->Ncols();
			Compare(gm->Type().AddEqualEl(),mt);
			if (!(mt==gm->Type()))
			{
				REPORT
					GeneralMatrix* gmx = mt.New(nr,nc,this);
				MatrixRow mr(gm, LoadOnEntry);
				MatrixRow mrx(gmx, StoreOnExit+DirectPart);
				while (nr--) { mrx.Add(mr,f); mrx.Next(); mr.Next(); }
				gmx->ReleaseAndDelete(); gm->tDelete();
				return gmx;
			}
			else if (gm->reuse())
			{
				REPORT gm->Add(f);
				return gm;
			}
			else
			{
				REPORT GeneralMatrix* gmy = gm->Type().New(nr,nc,this);
				gmy->ReleaseAndDelete(); gmy->Add(gm,f);
				return gmy;
			}
		}

		GeneralMatrix* NegShiftedMatrix::Evaluate(MatrixType mt)
		{
			gm=((BaseMatrix*&)bm)->Evaluate();
			int nr=gm->Nrows(); int nc=gm->Ncols();
			Compare(gm->Type().AddEqualEl(),mt);
			if (!(mt==gm->Type()))
			{
				REPORT
					GeneralMatrix* gmx = mt.New(nr,nc,this);
				MatrixRow mr(gm, LoadOnEntry);
				MatrixRow mrx(gmx, StoreOnExit+DirectPart);
				while (nr--) { mrx.NegAdd(mr,f); mrx.Next(); mr.Next(); }
				gmx->ReleaseAndDelete(); gm->tDelete();
				return gmx;
			}
			else if (gm->reuse())
			{
				REPORT gm->NegAdd(f);
				return gm;
			}
			else
			{
				REPORT GeneralMatrix* gmy = gm->Type().New(nr,nc,this);
				gmy->ReleaseAndDelete(); gmy->NegAdd(gm,f);
				return gmy;
			}
		}

		GeneralMatrix* ScaledMatrix::Evaluate(MatrixType mt)
		{
			gm=((BaseMatrix*&)bm)->Evaluate();
			int nr=gm->Nrows(); int nc=gm->Ncols();
			if (Compare(gm->Type(),mt))
			{
				if (gm->reuse())
				{
					REPORT gm->Multiply(f);
					return gm;
				}
				else
				{
					REPORT GeneralMatrix* gmx = gm->Type().New(nr,nc,this);
					gmx->ReleaseAndDelete(); gmx->Multiply(gm,f);
					return gmx;
				}
			}
			else
			{
				REPORT
					GeneralMatrix* gmx = mt.New(nr,nc,this);
				MatrixRow mr(gm, LoadOnEntry);
				MatrixRow mrx(gmx, StoreOnExit+DirectPart);
				while (nr--) { mrx.Multiply(mr,f); mrx.Next(); mr.Next(); }
				gmx->ReleaseAndDelete(); gm->tDelete();
				return gmx;
			}
		}

		GeneralMatrix* NegatedMatrix::Evaluate(MatrixType mt)
		{
			gm=((BaseMatrix*&)bm)->Evaluate();
			int nr=gm->Nrows(); int nc=gm->Ncols();
			if (Compare(gm->Type(),mt))
			{
				if (gm->reuse())
				{
					REPORT gm->Negate();
					return gm;
				}
				else
				{
					REPORT
						GeneralMatrix* gmx = gm->Type().New(nr,nc,this);
					gmx->ReleaseAndDelete(); gmx->Negate(gm);
					return gmx;
				}
			}
			else
			{
				REPORT
					GeneralMatrix* gmx = mt.New(nr,nc,this);
				MatrixRow mr(gm, LoadOnEntry);
				MatrixRow mrx(gmx, StoreOnExit+DirectPart);
				while (nr--) { mrx.Negate(mr); mrx.Next(); mr.Next(); }
				gmx->ReleaseAndDelete(); gm->tDelete();
				return gmx;
			}
		}

		GeneralMatrix* ReversedMatrix::Evaluate(MatrixType mt)
		{
			gm=((BaseMatrix*&)bm)->Evaluate(); GeneralMatrix* gmx;

			if ((gm->Type()).is_band() && ! (gm->Type()).is_diagonal())
			{
				gm->tDelete();
				Throw(NotDefinedException("Reverse", "band matrices"));
			}

			if (gm->reuse()) { REPORT gm->ReverseElements(); gmx = gm; }
			else
			{
				REPORT
					gmx = gm->Type().New(gm->Nrows(), gm->Ncols(), this);
				gmx->ReverseElements(gm); gmx->ReleaseAndDelete();
			}
			return gmx->Evaluate(mt); // target matrix is different type?

		}

		GeneralMatrix* TransposedMatrix::Evaluate(MatrixType mt)
		{
			REPORT
				gm=((BaseMatrix*&)bm)->Evaluate();
			Compare(gm->Type().t(),mt);
			GeneralMatrix* gmx=gm->Transpose(this, mt);
			return gmx;
		}

		GeneralMatrix* RowedMatrix::Evaluate(MatrixType mt)
		{
			gm = ((BaseMatrix*&)bm)->Evaluate();
			GeneralMatrix* gmx = new RowVector; MatrixErrorNoSpace(gmx);
			gmx->nrows_val = 1; gmx->ncols_val = gmx->storage = gm->storage;
			return gm->BorrowStore(gmx,mt);
		}

		GeneralMatrix* ColedMatrix::Evaluate(MatrixType mt)
		{
			gm = ((BaseMatrix*&)bm)->Evaluate();
			GeneralMatrix* gmx = new ColumnVector; MatrixErrorNoSpace(gmx);
			gmx->ncols_val = 1; gmx->nrows_val = gmx->storage = gm->storage;
			return gm->BorrowStore(gmx,mt);
		}

		GeneralMatrix* DiagedMatrix::Evaluate(MatrixType mt)
		{
			gm = ((BaseMatrix*&)bm)->Evaluate();
			GeneralMatrix* gmx = new DiagonalMatrix; MatrixErrorNoSpace(gmx);
			gmx->nrows_val = gmx->ncols_val = gmx->storage = gm->storage;
			return gm->BorrowStore(gmx,mt);
		}

		GeneralMatrix* MatedMatrix::Evaluate(MatrixType mt)
		{
			Tracer tr("MatedMatrix::Evaluate");
			gm = ((BaseMatrix*&)bm)->Evaluate();
			GeneralMatrix* gmx = new Matrix; MatrixErrorNoSpace(gmx);
			gmx->nrows_val = nr; gmx->ncols_val = nc; gmx->storage = gm->storage;
			if (nr*nc != gmx->storage)
				Throw(IncompatibleDimensionsException());
			return gm->BorrowStore(gmx,mt);
		}

		GeneralMatrix* GetSubMatrix::Evaluate(MatrixType mt)
		{
			REPORT
				Tracer tr("SubMatrix(evaluate)");
			gm = ((BaseMatrix*&)bm)->Evaluate();
			if (row_number < 0) row_number = gm->Nrows();
			if (col_number < 0) col_number = gm->Ncols();
			if (row_skip+row_number > gm->Nrows() || col_skip+col_number > gm->Ncols())
			{
				gm->tDelete();
				Throw(SubMatrixDimensionException());
			}
			if (IsSym) Compare(gm->Type().ssub(), mt);
			else Compare(gm->Type().sub(), mt);
			GeneralMatrix* gmx = mt.New(row_number, col_number, this);
			int i = row_number;
			MatrixRow mr(gm, LoadOnEntry, row_skip); 
			MatrixRow mrx(gmx, StoreOnExit+DirectPart);
			MatrixRowCol sub;
			while (i--)
			{
				mr.SubRowCol(sub, col_skip, col_number);   // put values in sub
				mrx.Copy(sub); mrx.Next(); mr.Next();
			}
			gmx->ReleaseAndDelete(); gm->tDelete();
			return gmx;
		}


		GeneralMatrix* ReturnMatrix::Evaluate(MatrixType mt)
		{
			return gm->Evaluate(mt);
		}


		void GeneralMatrix::Add(GeneralMatrix* gm1, Real f)
		{
			REPORT
				Real* s1=gm1->store; Real* s=store; int i=(storage >> 2);
			while (i--)
			{ *s++ = *s1++ + f; *s++ = *s1++ + f; *s++ = *s1++ + f; *s++ = *s1++ + f; }
			i = storage & 3; while (i--) *s++ = *s1++ + f;
		}

		void GeneralMatrix::Add(Real f)
		{
			REPORT
				Real* s=store; int i=(storage >> 2);
			while (i--) { *s++ += f; *s++ += f; *s++ += f; *s++ += f; }
			i = storage & 3; while (i--) *s++ += f;
		}

		void GeneralMatrix::NegAdd(GeneralMatrix* gm1, Real f)
		{
			REPORT
				Real* s1=gm1->store; Real* s=store; int i=(storage >> 2);
			while (i--)
			{ *s++ = f - *s1++; *s++ = f - *s1++; *s++ = f - *s1++; *s++ = f - *s1++; }
			i = storage & 3; while (i--) *s++ = f - *s1++;
		}

		void GeneralMatrix::NegAdd(Real f)
		{
			REPORT
				Real* s=store; int i=(storage >> 2);
			while (i--)
			{
				*s = f - *s; s++; *s = f - *s; s++;
				*s = f - *s; s++; *s = f - *s; s++;
			}
			i = storage & 3; while (i--)  { *s = f - *s; s++; }
		}

		void GeneralMatrix::Negate(GeneralMatrix* gm1)
		{
			// change sign of elements
			REPORT
				Real* s1=gm1->store; Real* s=store; int i=(storage >> 2);
			while (i--)
			{ *s++ = -(*s1++); *s++ = -(*s1++); *s++ = -(*s1++); *s++ = -(*s1++); }
			i = storage & 3; while(i--) *s++ = -(*s1++);
		}

		void GeneralMatrix::Negate()
		{
			REPORT
				Real* s=store; int i=(storage >> 2);
			while (i--)
			{ *s = -(*s); s++; *s = -(*s); s++; *s = -(*s); s++; *s = -(*s); s++; }
			i = storage & 3; while(i--) { *s = -(*s); s++; }
		}

		void GeneralMatrix::Multiply(GeneralMatrix* gm1, Real f)
		{
			REPORT
				Real* s1=gm1->store; Real* s=store;  int i=(storage >> 2);
			while (i--)
			{ *s++ = *s1++ * f; *s++ = *s1++ * f; *s++ = *s1++ * f; *s++ = *s1++ * f; }
			i = storage & 3; while (i--) *s++ = *s1++ * f;
		}

		void GeneralMatrix::Multiply(Real f)
		{
			REPORT
				Real* s=store; int i=(storage >> 2);
			while (i--) { *s++ *= f; *s++ *= f; *s++ *= f; *s++ *= f; }
			i = storage & 3; while (i--) *s++ *= f;
		}


		/************************ MatrixInput routines ****************************/

		// int MatrixInput::n;          // number values still to be read
		// Real* MatrixInput::r;        // pointer to next location to be read to

		MatrixInput MatrixInput::operator<<(double f)
		{
			REPORT
				Tracer et("MatrixInput");
			if (n<=0) Throw(ProgramException("List of values too long"));
			*r = (Real)f; int n1 = n-1; n=0;   // n=0 so we won't trigger exception
			return MatrixInput(n1, r+1);
		}


		MatrixInput GeneralMatrix::operator<<(double f)
		{
			REPORT
				Tracer et("MatrixInput");
			int n = Storage();
			if (n<=0) Throw(ProgramException("Loading data to zero length matrix"));
			Real* r; r = Store(); *r = (Real)f; n--;
			return MatrixInput(n, r+1);
		}

		MatrixInput GetSubMatrix::operator<<(double f)
		{
			REPORT
				Tracer et("MatrixInput (GetSubMatrix)");
			SetUpLHS();
			if (row_number != 1 || col_skip != 0 || col_number != gm->Ncols())
			{
				Throw(ProgramException("MatrixInput requires complete rows"));
			}
			MatrixRow mr(gm, DirectPart, row_skip);  // to pick up location and length
			int n = mr.Storage();
			if (n<=0)
			{
				Throw(ProgramException("Loading data to zero length row"));
			}
			Real* r; r = mr.Data(); *r = (Real)f; n--;
			if (+(mr.cw*HaveStore))
			{
				Throw(ProgramException("Fails with this matrix type"));
			}
			return MatrixInput(n, r+1);
		}

		MatrixInput MatrixInput::operator<<(float f)
		{
			REPORT
				Tracer et("MatrixInput");
			if (n<=0) Throw(ProgramException("List of values too long"));
			*r = (Real)f; int n1 = n-1; n=0;   // n=0 so we won't trigger exception
			return MatrixInput(n1, r+1);
		}


		MatrixInput GeneralMatrix::operator<<(float f)
		{
			REPORT
				Tracer et("MatrixInput");
			int n = Storage();
			if (n<=0) Throw(ProgramException("Loading data to zero length matrix"));
			Real* r; r = Store(); *r = (Real)f; n--;
			return MatrixInput(n, r+1);
		}

		MatrixInput GetSubMatrix::operator<<(float f)
		{
			REPORT
				Tracer et("MatrixInput (GetSubMatrix)");
			SetUpLHS();
			if (row_number != 1 || col_skip != 0 || col_number != gm->Ncols())
			{
				Throw(ProgramException("MatrixInput requires complete rows"));
			}
			MatrixRow mr(gm, DirectPart, row_skip);  // to pick up location and length
			int n = mr.Storage();
			if (n<=0)
			{
				Throw(ProgramException("Loading data to zero length row"));
			}
			Real* r; r = mr.Data(); *r = (Real)f; n--;
			if (+(mr.cw*HaveStore))
			{
				Throw(ProgramException("Fails with this matrix type"));
			}
			return MatrixInput(n, r+1);
		}
		MatrixInput::~MatrixInput()
		{
			REPORT
				Tracer et("MatrixInput");
			if (n!=0) Throw(ProgramException("A list of values was too short"));
		}

		MatrixInput BandMatrix::operator<<(double)
		{
			Tracer et("MatrixInput");
			bool dummy = true;
			if (dummy)                                   // get rid of warning message
				Throw(ProgramException("Cannot use list read with a BandMatrix"));
			return MatrixInput(0, 0);
		}

		MatrixInput BandMatrix::operator<<(float)
		{
			Tracer et("MatrixInput");
			bool dummy = true;
			if (dummy)                                   // get rid of warning message
				Throw(ProgramException("Cannot use list read with a BandMatrix"));
			return MatrixInput(0, 0);
		}

		void BandMatrix::operator<<(const double*)
		{ Throw(ProgramException("Cannot use array read with a BandMatrix")); }

		void BandMatrix::operator<<(const float*)
		{ Throw(ProgramException("Cannot use array read with a BandMatrix")); }

		void BandMatrix::operator<<(const int*)
		{ Throw(ProgramException("Cannot use array read with a BandMatrix")); }

		void SymmetricBandMatrix::operator<<(const double*)
		{ Throw(ProgramException("Cannot use array read with a BandMatrix")); }

		void SymmetricBandMatrix::operator<<(const float*)
		{ Throw(ProgramException("Cannot use array read with a BandMatrix")); }

		void SymmetricBandMatrix::operator<<(const int*)
		{ Throw(ProgramException("Cannot use array read with a BandMatrix")); }

		// ************************* Reverse order of elements ***********************

		void GeneralMatrix::ReverseElements(GeneralMatrix* gm)
		{
			// reversing into a new matrix
			REPORT
				int n = Storage(); Real* rx = Store() + n; Real* x = gm->Store();
			while (n--) *(--rx) = *(x++);
		}

		void GeneralMatrix::ReverseElements()
		{
			// reversing in place
			REPORT
				int n = Storage(); Real* x = Store(); Real* rx = x + n;
			n /= 2;
			while (n--) { Real t = *(--rx); *rx = *x; *(x++) = t; }
		}


#ifdef use_namespace
	}
#endif

	//$$ newmat4.cpp       Constructors, resize, basic utilities

	// Copyright (C) 1991,2,3,4,8,9: R B Davies

	//#define WANT_STREAM


#ifdef use_namespace
	namespace NEWMAT {
#endif




#ifdef DO_REPORT
#define REPORT { static ExeCounter ExeCount(__LINE__,4); ++ExeCount; }
#else
#define REPORT {}
#endif


#define DO_SEARCH                   // search for LHS of = in RHS

		// **************************** constructors ***************************/

		GeneralMatrix::GeneralMatrix()
		{ store=0; storage=0; nrows_val=0; ncols_val=0; tag_val=-1; }

		GeneralMatrix::GeneralMatrix(ArrayLengthSpecifier s)
		{
			REPORT
				storage=s.Value(); tag_val=-1;
			if (storage)
			{
				store = new Real [storage]; MatrixErrorNoSpace(store);
				MONITOR_REAL_NEW("Make (GenMatrix)",storage,store)
			}
			else store = 0;
		}

		Matrix::Matrix(int m, int n) : GeneralMatrix(m*n)
		{ REPORT nrows_val=m; ncols_val=n; }

		SquareMatrix::SquareMatrix(ArrayLengthSpecifier n)
			: Matrix(n.Value(),n.Value())
		{ REPORT }

		SymmetricMatrix::SymmetricMatrix(ArrayLengthSpecifier n)
			: GeneralMatrix(tristore(n.Value()))
		{ REPORT nrows_val=n.Value(); ncols_val=n.Value(); }

		UpperTriangularMatrix::UpperTriangularMatrix(ArrayLengthSpecifier n)
			: GeneralMatrix(tristore(n.Value()))
		{ REPORT nrows_val=n.Value(); ncols_val=n.Value(); }

		LowerTriangularMatrix::LowerTriangularMatrix(ArrayLengthSpecifier n)
			: GeneralMatrix(tristore(n.Value()))
		{ REPORT nrows_val=n.Value(); ncols_val=n.Value(); }

		DiagonalMatrix::DiagonalMatrix(ArrayLengthSpecifier m) : GeneralMatrix(m)
		{ REPORT nrows_val=m.Value(); ncols_val=m.Value(); }

		Matrix::Matrix(const BaseMatrix& M)
		{
			REPORT // CheckConversion(M);
				// MatrixConversionCheck mcc;
				GeneralMatrix* gmx=((BaseMatrix&)M).Evaluate(MatrixType::Rt);
			GetMatrix(gmx);
		}

		SquareMatrix::SquareMatrix(const BaseMatrix& M) : Matrix(M)
		{
			REPORT
				if (ncols_val != nrows_val)
				{
					Tracer tr("SquareMatrix");
					Throw(NotSquareException(*this));
				}
		}


		SquareMatrix::SquareMatrix(const Matrix& gm)
		{
			REPORT
				if (gm.ncols_val != gm.nrows_val)
				{
					Tracer tr("SquareMatrix(Matrix)");
					Throw(NotSquareException(gm));
				}
				GetMatrix(&gm);
		}


		RowVector::RowVector(const BaseMatrix& M) : Matrix(M)
		{
			REPORT
				if (nrows_val!=1)
				{
					Tracer tr("RowVector");
					Throw(VectorException(*this));
				}
		}

		ColumnVector::ColumnVector(const BaseMatrix& M) : Matrix(M)
		{
			REPORT
				if (ncols_val!=1)
				{
					Tracer tr("ColumnVector");
					Throw(VectorException(*this));
				}
		}

		SymmetricMatrix::SymmetricMatrix(const BaseMatrix& M)
		{
			REPORT  // CheckConversion(M);
				// MatrixConversionCheck mcc;
				GeneralMatrix* gmx=((BaseMatrix&)M).Evaluate(MatrixType::Sm);
			GetMatrix(gmx);
		}

		UpperTriangularMatrix::UpperTriangularMatrix(const BaseMatrix& M)
		{
			REPORT // CheckConversion(M);
				// MatrixConversionCheck mcc;
				GeneralMatrix* gmx=((BaseMatrix&)M).Evaluate(MatrixType::UT);
			GetMatrix(gmx);
		}

		LowerTriangularMatrix::LowerTriangularMatrix(const BaseMatrix& M)
		{
			REPORT // CheckConversion(M);
				// MatrixConversionCheck mcc;
				GeneralMatrix* gmx=((BaseMatrix&)M).Evaluate(MatrixType::LT);
			GetMatrix(gmx);
		}

		DiagonalMatrix::DiagonalMatrix(const BaseMatrix& M)
		{
			REPORT //CheckConversion(M);
				// MatrixConversionCheck mcc;
				GeneralMatrix* gmx=((BaseMatrix&)M).Evaluate(MatrixType::Dg);
			GetMatrix(gmx);
		}

		IdentityMatrix::IdentityMatrix(const BaseMatrix& M)
		{
			REPORT //CheckConversion(M);
				// MatrixConversionCheck mcc;
				GeneralMatrix* gmx=((BaseMatrix&)M).Evaluate(MatrixType::Id);
			GetMatrix(gmx);
		}

		GeneralMatrix::~GeneralMatrix()
		{
			if (store)
			{
				MONITOR_REAL_DELETE("Free (GenMatrix)",storage,store)
					delete [] store;
			}
		}

		CroutMatrix::CroutMatrix(const BaseMatrix& m)
		{
			REPORT
				Tracer tr("CroutMatrix");
			indx = 0;                     // in case of exception at next line
			GeneralMatrix* gm = ((BaseMatrix&)m).Evaluate();
			if (gm->nrows_val!=gm->ncols_val)
			{ gm->tDelete(); Throw(NotSquareException(*gm)); }
			if (gm->type() == MatrixType::Ct)
			{ REPORT  ((CroutMatrix*)gm)->get_aux(*this); GetMatrix(gm); }
			else
			{
				REPORT
					GeneralMatrix* gm1 = gm->Evaluate(MatrixType::Rt);
				GetMatrix(gm1);
				d=true; sing=false;
				indx=new int [nrows_val]; MatrixErrorNoSpace(indx);
				MONITOR_INT_NEW("Index (CroutMat)",nrows_val,indx)
					ludcmp();
			}
		}

		// could we use SetParameters instead of this
		void CroutMatrix::get_aux(CroutMatrix& X)
		{
			X.d = d; X.sing = sing;
			if (tag_val == 0 || tag_val == 1) // reuse the array 
			{ REPORT  X.indx = indx; indx = 0; d = true; sing = true; return; }
			else if (nrows_val == 0)
			{ REPORT indx = 0; d = true; sing = true; return; }
			else                              // copy the array
			{ 
				REPORT
					Tracer tr("CroutMatrix::get_aux");
				int *ix = new int [nrows_val]; MatrixErrorNoSpace(ix);
				MONITOR_INT_NEW("Index (CM::get_aux)", nrows_val, ix)
					int n = nrows_val; int* i = ix; int* j = indx;
				while(n--) *i++ = *j++;
				X.indx = ix;
			}
		}

		CroutMatrix::CroutMatrix(const CroutMatrix& gm)
		{
			REPORT
				Tracer tr("CroutMatrix(const CroutMatrix&)");
			((CroutMatrix&)gm).get_aux(*this);
			GetMatrix(&gm);
		}

		CroutMatrix::~CroutMatrix()
		{
			MONITOR_INT_DELETE("Index (CroutMat)",nrows_val,indx)
				delete [] indx;
		}

		//ReturnMatrix::ReturnMatrix(GeneralMatrix& gmx)
		//{
		//   REPORT
		//   gm = gmx.Image(); gm->ReleaseAndDelete();
		//}


		GeneralMatrix::operator ReturnMatrix() const
		{
			REPORT
				GeneralMatrix* gm = Image(); gm->ReleaseAndDelete();
			return ReturnMatrix(gm);
		}



		ReturnMatrix GeneralMatrix::for_return() const
		{
			REPORT
				GeneralMatrix* gm = Image(); gm->ReleaseAndDelete();
			return ReturnMatrix(gm);
		}

		// ************ Constructors for use with NR in C++ interface ***********

#ifdef SETUP_C_SUBSCRIPTS

		Matrix::Matrix(Real a, int m, int n) : GeneralMatrix(m * n)
		{ REPORT nrows_val=m; ncols_val=n; operator=(a); }

		Matrix::Matrix(const Real* a, int m, int n) : GeneralMatrix(m * n)
		{ REPORT nrows_val=m; ncols_val=n; *this << a; }

#endif



		// ************************** resize matrices ***************************/

		void GeneralMatrix::resize(int nr, int nc, int s)
		{
			REPORT
				if (store)
				{
					MONITOR_REAL_DELETE("Free (ReDimensi)",storage,store)
						delete [] store;
				}
				storage=s; nrows_val=nr; ncols_val=nc; tag_val=-1;
				if (s)
				{
					store = new Real [storage]; MatrixErrorNoSpace(store);
					MONITOR_REAL_NEW("Make (ReDimensi)",storage,store)
				}
				else store = 0;
		}

		void Matrix::resize(int nr, int nc)
		{ REPORT GeneralMatrix::resize(nr,nc,nr*nc); }

		void SquareMatrix::resize(int n)
		{ REPORT GeneralMatrix::resize(n,n,n*n); }

		void SquareMatrix::resize(int nr, int nc)
		{
			REPORT
				Tracer tr("SquareMatrix::resize");
			if (nc != nr) Throw(NotSquareException(*this));
			GeneralMatrix::resize(nr,nc,nr*nc);
		}

		void SymmetricMatrix::resize(int nr)
		{ REPORT GeneralMatrix::resize(nr,nr,tristore(nr)); }

		void UpperTriangularMatrix::resize(int nr)
		{ REPORT GeneralMatrix::resize(nr,nr,tristore(nr)); }

		void LowerTriangularMatrix::resize(int nr)
		{ REPORT GeneralMatrix::resize(nr,nr,tristore(nr)); }

		void DiagonalMatrix::resize(int nr)
		{ REPORT GeneralMatrix::resize(nr,nr,nr); }

		void RowVector::resize(int nc)
		{ REPORT GeneralMatrix::resize(1,nc,nc); }

		void ColumnVector::resize(int nr)
		{ REPORT GeneralMatrix::resize(nr,1,nr); }

		void RowVector::resize(int nr, int nc)
		{
			Tracer tr("RowVector::resize");
			if (nr != 1) Throw(VectorException(*this));
			REPORT GeneralMatrix::resize(1,nc,nc);
		}

		void ColumnVector::resize(int nr, int nc)
		{
			Tracer tr("ColumnVector::resize");
			if (nc != 1) Throw(VectorException(*this));
			REPORT GeneralMatrix::resize(nr,1,nr);
		}

		void IdentityMatrix::resize(int nr)
		{ REPORT GeneralMatrix::resize(nr,nr,1); *store = 1; }


		void Matrix::resize(const GeneralMatrix& A)
		{ REPORT  resize(A.Nrows(), A.Ncols()); }

		void SquareMatrix::resize(const GeneralMatrix& A)
		{
			REPORT
				int n = A.Nrows();
			if (n != A.Ncols())
			{
				Tracer tr("SquareMatrix::resize(GM)");
				Throw(NotSquareException(*this));
			}
			resize(n);
		}

		void nricMatrix::resize(const GeneralMatrix& A)
		{ REPORT  resize(A.Nrows(), A.Ncols()); }

		void ColumnVector::resize(const GeneralMatrix& A)
		{ REPORT  resize(A.Nrows(), A.Ncols()); }

		void RowVector::resize(const GeneralMatrix& A)
		{ REPORT  resize(A.Nrows(), A.Ncols()); }

		void SymmetricMatrix::resize(const GeneralMatrix& A)
		{
			REPORT
				int n = A.Nrows();
			if (n != A.Ncols())
			{
				Tracer tr("SymmetricMatrix::resize(GM)");
				Throw(NotSquareException(*this));
			}
			resize(n);
		}

		void DiagonalMatrix::resize(const GeneralMatrix& A)
		{
			REPORT
				int n = A.Nrows();
			if (n != A.Ncols())
			{
				Tracer tr("DiagonalMatrix::resize(GM)");
				Throw(NotSquareException(*this));
			}
			resize(n);
		}

		void UpperTriangularMatrix::resize(const GeneralMatrix& A)
		{
			REPORT
				int n = A.Nrows();
			if (n != A.Ncols())
			{
				Tracer tr("UpperTriangularMatrix::resize(GM)");
				Throw(NotSquareException(*this));
			}
			resize(n);
		}

		void LowerTriangularMatrix::resize(const GeneralMatrix& A)
		{
			REPORT
				int n = A.Nrows();
			if (n != A.Ncols())
			{
				Tracer tr("LowerTriangularMatrix::resize(GM)");
				Throw(NotSquareException(*this));
			}
			resize(n);
		}

		void IdentityMatrix::resize(const GeneralMatrix& A)
		{
			REPORT
				int n = A.Nrows();
			if (n != A.Ncols())
			{
				Tracer tr("IdentityMatrix::resize(GM)");
				Throw(NotSquareException(*this));
			}
			resize(n);
		}

		void GeneralMatrix::resize(const GeneralMatrix&)
		{
			REPORT
				Tracer tr("GeneralMatrix::resize(GM)");
			Throw(NotDefinedException("resize", "this type of matrix"));
		}

		//*********************** resize_keep *******************************

		void Matrix::resize_keep(int nr, int nc)
		{
			Tracer tr("Matrix::resize_keep");
			if (nr == nrows_val && nc == ncols_val) { REPORT return; }

			if (nr <= nrows_val && nc <= ncols_val)
			{
				REPORT
					Matrix X = submatrix(1,nr,1,nc);
				swap(X);
			}
			else if (nr >= nrows_val && nc >= ncols_val)
			{
				REPORT
					Matrix X(nr, nc); X = 0;
				X.submatrix(1,nrows_val,1,ncols_val) = *this;
				swap(X);
			}
			else
			{
				REPORT
					Matrix X(nr, nc); X = 0;
				if (nr > nrows_val) nr = nrows_val;
				if (nc > ncols_val) nc = ncols_val;
				X.submatrix(1,nr,1,nc) = submatrix(1,nr,1,nc);
				swap(X);
			}
		} 

		void SquareMatrix::resize_keep(int nr)
		{
			Tracer tr("SquareMatrix::resize_keep");
			if (nr < nrows_val)
			{
				REPORT
					SquareMatrix X = sym_submatrix(1,nr);
				swap(X);
			}
			else if (nr > nrows_val)
			{
				REPORT
					SquareMatrix X(nr); X = 0;
				X.sym_submatrix(1,nrows_val) = *this;
				swap(X);
			}
		}

		void SquareMatrix::resize_keep(int nr, int nc)
		{
			Tracer tr("SquareMatrix::resize_keep 2");
			REPORT
				if (nr != nc) Throw(NotSquareException(*this));
			resize_keep(nr);
		}


		void SymmetricMatrix::resize_keep(int nr)
		{
			Tracer tr("SymmetricMatrix::resize_keep");
			if (nr < nrows_val)
			{
				REPORT
					SymmetricMatrix X = sym_submatrix(1,nr);
				swap(X);
			}
			else if (nr > nrows_val)
			{
				REPORT
					SymmetricMatrix X(nr); X = 0;
				X.sym_submatrix(1,nrows_val) = *this;
				swap(X);
			}
		} 

		void UpperTriangularMatrix::resize_keep(int nr)
		{
			Tracer tr("UpperTriangularMatrix::resize_keep");
			if (nr < nrows_val)
			{
				REPORT
					UpperTriangularMatrix X = sym_submatrix(1,nr);
				swap(X);
			}
			else if (nr > nrows_val)
			{
				REPORT
					UpperTriangularMatrix X(nr); X = 0;
				X.sym_submatrix(1,nrows_val) = *this;
				swap(X);
			}
		} 

		void LowerTriangularMatrix::resize_keep(int nr)
		{
			Tracer tr("LowerTriangularMatrix::resize_keep");
			if (nr < nrows_val)
			{
				REPORT
					LowerTriangularMatrix X = sym_submatrix(1,nr);
				swap(X);
			}
			else if (nr > nrows_val)
			{
				REPORT
					LowerTriangularMatrix X(nr); X = 0;
				X.sym_submatrix(1,nrows_val) = *this;
				swap(X);
			}
		} 

		void DiagonalMatrix::resize_keep(int nr)
		{
			Tracer tr("DiagonalMatrix::resize_keep");
			if (nr < nrows_val)
			{
				REPORT
					DiagonalMatrix X = sym_submatrix(1,nr);
				swap(X);
			}
			else if (nr > nrows_val)
			{
				REPORT
					DiagonalMatrix X(nr); X = 0;
				X.sym_submatrix(1,nrows_val) = *this;
				swap(X);
			}
		} 

		void RowVector::resize_keep(int nc)
		{
			Tracer tr("RowVector::resize_keep");
			if (nc < ncols_val)
			{
				REPORT
					RowVector X = columns(1,nc);
				swap(X);
			}
			else if (nc > ncols_val)
			{
				REPORT
					RowVector X(nc); X = 0;
				X.columns(1,ncols_val) = *this;
				swap(X);
			}
		}

		void RowVector::resize_keep(int nr, int nc)
		{
			Tracer tr("RowVector::resize_keep 2");
			REPORT
				if (nr != 1) Throw(VectorException(*this));
			resize_keep(nc);
		}

		void ColumnVector::resize_keep(int nr)
		{
			Tracer tr("ColumnVector::resize_keep");
			if (nr < nrows_val)
			{
				REPORT
					ColumnVector X = rows(1,nr);
				swap(X);
			}
			else if (nr > nrows_val)
			{
				REPORT
					ColumnVector X(nr); X = 0;
				X.rows(1,nrows_val) = *this;
				swap(X);
			}
		} 

		void ColumnVector::resize_keep(int nr, int nc)
		{
			Tracer tr("ColumnVector::resize_keep 2");
			REPORT
				if (nc != 1) Throw(VectorException(*this));
			resize_keep(nr);
		}


		/*
		void GeneralMatrix::resizeForAdd(const GeneralMatrix& A, const GeneralMatrix&)
		{ REPORT resize(A); }

		void GeneralMatrix::resizeForSP(const GeneralMatrix& A, const GeneralMatrix&)
		{ REPORT resize(A); }


		// ************************* SameStorageType ******************************

		// SameStorageType checks A and B have same storage type including bandwidth
		// It does not check same dimensions since we assume this is already done

		bool GeneralMatrix::SameStorageType(const GeneralMatrix& A) const
		{
		REPORT
		return type() == A.type();
		}
		*/

		// ******************* manipulate types, storage **************************/

		int GeneralMatrix::search(const BaseMatrix* s) const
		{ REPORT return (s==this) ? 1 : 0; }

		int GenericMatrix::search(const BaseMatrix* s) const
		{ REPORT return gm->search(s); }

		int MultipliedMatrix::search(const BaseMatrix* s) const
		{ REPORT return bm1->search(s) + bm2->search(s); }

		int ShiftedMatrix::search(const BaseMatrix* s) const
		{ REPORT return bm->search(s); }

		int NegatedMatrix::search(const BaseMatrix* s) const
		{ REPORT return bm->search(s); }

		int ReturnMatrix::search(const BaseMatrix* s) const
		{ REPORT return (s==gm) ? 1 : 0; }

		MatrixType Matrix::type() const { return MatrixType::Rt; }
		MatrixType SquareMatrix::type() const { return MatrixType::Sq; }
		MatrixType SymmetricMatrix::type() const { return MatrixType::Sm; }
		MatrixType UpperTriangularMatrix::type() const { return MatrixType::UT; }
		MatrixType LowerTriangularMatrix::type() const { return MatrixType::LT; }
		MatrixType DiagonalMatrix::type() const { return MatrixType::Dg; }
		MatrixType RowVector::type() const { return MatrixType::RV; }
		MatrixType ColumnVector::type() const { return MatrixType::CV; }
		MatrixType CroutMatrix::type() const { return MatrixType::Ct; }
		MatrixType BandMatrix::type() const { return MatrixType::BM; }
		MatrixType UpperBandMatrix::type() const { return MatrixType::UB; }
		MatrixType LowerBandMatrix::type() const { return MatrixType::LB; }
		MatrixType SymmetricBandMatrix::type() const { return MatrixType::SB; }

		MatrixType IdentityMatrix::type() const { return MatrixType::Id; }



		MatrixBandWidth BaseMatrix::bandwidth() const { REPORT return -1; }
		MatrixBandWidth DiagonalMatrix::bandwidth() const { REPORT return 0; }
		MatrixBandWidth IdentityMatrix::bandwidth() const { REPORT return 0; }

		MatrixBandWidth UpperTriangularMatrix::bandwidth() const
		{ REPORT return MatrixBandWidth(0,-1); }

		MatrixBandWidth LowerTriangularMatrix::bandwidth() const
		{ REPORT return MatrixBandWidth(-1,0); }

		MatrixBandWidth BandMatrix::bandwidth() const
		{ REPORT return MatrixBandWidth(lower_val,upper_val); }

		MatrixBandWidth BandLUMatrix::bandwidth() const
		{ REPORT return MatrixBandWidth(m1,m2); }

		MatrixBandWidth GenericMatrix::bandwidth()const
		{ REPORT return gm->bandwidth(); }

		MatrixBandWidth AddedMatrix::bandwidth() const
		{ REPORT return gm1->bandwidth() + gm2->bandwidth(); }

		MatrixBandWidth SPMatrix::bandwidth() const
		{ REPORT return gm1->bandwidth().minimum(gm2->bandwidth()); }

		MatrixBandWidth KPMatrix::bandwidth() const
		{
			int lower, upper;
			MatrixBandWidth bw1 = gm1->bandwidth(), bw2 = gm2->bandwidth();
			if (bw1.Lower() < 0)
			{
				if (bw2.Lower() < 0) { REPORT lower = -1; }
				else { REPORT lower = bw2.Lower() + (gm1->Nrows() - 1) * gm2->Nrows(); }
			}
			else
			{
				if (bw2.Lower() < 0)
				{ REPORT lower = (1 + bw1.Lower()) * gm2->Nrows() - 1; }
				else { REPORT lower = bw2.Lower() + bw1.Lower() * gm2->Nrows(); }
			}
			if (bw1.Upper() < 0)
			{
				if (bw2.Upper() < 0) { REPORT upper = -1; }
				else { REPORT upper = bw2.Upper() + (gm1->Nrows() - 1) * gm2->Nrows(); }
			}
			else
			{
				if (bw2.Upper() < 0)
				{ REPORT upper = (1 + bw1.Upper()) * gm2->Nrows() - 1; }
				else { REPORT upper = bw2.Upper() + bw1.Upper() * gm2->Nrows(); }
			}
			return MatrixBandWidth(lower, upper);
		}

		MatrixBandWidth MultipliedMatrix::bandwidth() const
		{ REPORT return gm1->bandwidth() * gm2->bandwidth(); }

		MatrixBandWidth ConcatenatedMatrix::bandwidth() const { REPORT return -1; }

		MatrixBandWidth SolvedMatrix::bandwidth() const
		{
			if (+gm1->type() & MatrixType::Diagonal)
			{ REPORT return gm2->bandwidth(); }
			else { REPORT return -1; }
		}

		MatrixBandWidth ScaledMatrix::bandwidth() const
		{ REPORT return gm->bandwidth(); }

		MatrixBandWidth NegatedMatrix::bandwidth() const
		{ REPORT return gm->bandwidth(); }

		MatrixBandWidth TransposedMatrix::bandwidth() const
		{ REPORT return gm->bandwidth().t(); }

		MatrixBandWidth InvertedMatrix::bandwidth() const
		{
			if (+gm->type() & MatrixType::Diagonal)
			{ REPORT return MatrixBandWidth(0,0); }
			else { REPORT return -1; }
		}

		MatrixBandWidth RowedMatrix::bandwidth() const { REPORT return -1; }
		MatrixBandWidth ColedMatrix::bandwidth() const { REPORT return -1; }
		MatrixBandWidth DiagedMatrix::bandwidth() const { REPORT return 0; }
		MatrixBandWidth MatedMatrix::bandwidth() const { REPORT return -1; }
		MatrixBandWidth ReturnMatrix::bandwidth() const
		{ REPORT return gm->bandwidth(); }

		MatrixBandWidth GetSubMatrix::bandwidth() const
		{

			if (row_skip==col_skip && row_number==col_number)
			{ REPORT return gm->bandwidth(); }
			else { REPORT return MatrixBandWidth(-1); }
		}

		// ********************** the memory managment tools **********************/

		//  Rules regarding tDelete, reuse, GetStore, BorrowStore
		//    All matrices processed during expression evaluation must be subject
		//    to exactly one of reuse(), tDelete(), GetStore() or BorrowStore().
		//    If reuse returns true the matrix must be reused.
		//    GetMatrix(gm) always calls gm->GetStore()
		//    gm->Evaluate obeys rules
		//    bm->Evaluate obeys rules for matrices in bm structure

		//  Meaning of tag_val
		//    tag_val = -1          memory cannot be reused (default situation)
		//    tag_val = -2          memory has been borrowed from another matrix
		//                               (don't change values)
		//    tag_val = i > 0       delete or reuse memory after i operations
		//    tag_val = 0           like value 1 but matrix was created by new
		//                               so delete it

		void GeneralMatrix::tDelete()
		{
			if (tag_val<0)
			{
				if (tag_val<-1) { REPORT store = 0; delete this; return; }  // borrowed
				else { REPORT return; }   // not a temporary matrix - leave alone
			}
			if (tag_val==1)
			{
				if (store)
				{
					REPORT  MONITOR_REAL_DELETE("Free   (tDelete)",storage,store)
						delete [] store;
				}
				MiniCleanUp(); return;                           // CleanUp
			}
			if (tag_val==0) { REPORT delete this; return; }

			REPORT tag_val--; return;
		}

		void newmat_block_copy(int n, Real* from, Real* to)
		{
			REPORT
				int i = (n >> 3);
			while (i--)
			{
				*to++ = *from++; *to++ = *from++; *to++ = *from++; *to++ = *from++;
				*to++ = *from++; *to++ = *from++; *to++ = *from++; *to++ = *from++;
			}
			i = n & 7; while (i--) *to++ = *from++;
		}

		bool GeneralMatrix::reuse()
		{
			if (tag_val < -1)                 // borrowed storage
			{
				if (storage)
				{
					REPORT
						Real* s = new Real [storage]; MatrixErrorNoSpace(s);
					MONITOR_REAL_NEW("Make     (reuse)",storage,s)
						newmat_block_copy(storage, store, s); store = s;
				}
				else { REPORT MiniCleanUp(); }                // CleanUp
				tag_val = 0; return true;
			}
			if (tag_val < 0 ) { REPORT return false; }
			if (tag_val <= 1 )  { REPORT return true; }
			REPORT tag_val--; return false;
		}

		Real* GeneralMatrix::GetStore()
		{
			if (tag_val<0 || tag_val>1)
			{
				Real* s;
				if (storage)
				{
					s = new Real [storage]; MatrixErrorNoSpace(s);
					MONITOR_REAL_NEW("Make  (GetStore)",storage,s)
						newmat_block_copy(storage, store, s);
				}
				else s = 0;
				if (tag_val > 1) { REPORT tag_val--; }
				else if (tag_val < -1) { REPORT store = 0; delete this; } // borrowed store
				else { REPORT }
				return s;
			}
			Real* s = store;                             // cleanup - done later
			if (tag_val==0) { REPORT store = 0; delete this; }
			else { REPORT  MiniCleanUp(); }
			return s;
		}

		void GeneralMatrix::GetMatrix(const GeneralMatrix* gmx)
		{
			REPORT  tag_val=-1; nrows_val=gmx->Nrows(); ncols_val=gmx->Ncols();
			storage=gmx->storage; SetParameters(gmx);
			store=((GeneralMatrix*)gmx)->GetStore();
		}

		GeneralMatrix* GeneralMatrix::BorrowStore(GeneralMatrix* gmx, MatrixType mt)
			// Copy storage of *this to storage of *gmx. Then convert to type mt.
			// If mt == 0 just let *gmx point to storage of *this if tag_val==-1.
		{
			if (!mt)
			{
				if (tag_val == -1) { REPORT gmx->tag_val = -2; gmx->store = store; }
				else { REPORT gmx->tag_val = 0; gmx->store = GetStore(); }
			}
			else if (Compare(gmx->type(),mt))
			{ REPORT  gmx->tag_val = 0; gmx->store = GetStore(); }
			else
			{
				REPORT gmx->tag_val = -2; gmx->store = store;
				gmx = gmx->Evaluate(mt); gmx->tag_val = 0; tDelete();
			}

			return gmx;
		}

		void GeneralMatrix::Eq(const BaseMatrix& X, MatrixType mt)
			// Count number of references to this in X.
			// If zero delete storage in this;
			// otherwise tag this to show when storage can be deleted
			// evaluate X and copy to this
		{
#ifdef DO_SEARCH
			int counter=X.search(this);
			if (counter==0)
			{
				REPORT
					if (store)
					{
						MONITOR_REAL_DELETE("Free (operator=)",storage,store)
							REPORT delete [] store; storage = 0; store = 0;
					}
			}
			else { REPORT Release(counter); }
			GeneralMatrix* gmx = ((BaseMatrix&)X).Evaluate(mt);
			if (gmx!=this) { REPORT GetMatrix(gmx); }
			else { REPORT }
			Protect();
#else
			GeneralMatrix* gmx = ((BaseMatrix&)X).Evaluate(mt);
			if (gmx!=this)
			{
				REPORT
					if (store)
					{
						MONITOR_REAL_DELETE("Free (operator=)",storage,store)
							REPORT delete [] store; storage = 0; store = 0;
					}
					GetMatrix(gmx);
			}
			else { REPORT }
			Protect();
#endif
		}

		// version with no conversion
		void GeneralMatrix::Eq(const GeneralMatrix& X)
		{
			GeneralMatrix* gmx = (GeneralMatrix*)&X;
			if (gmx!=this)
			{
				REPORT
					if (store)
					{
						MONITOR_REAL_DELETE("Free (operator=)",storage,store)
							REPORT delete [] store; storage = 0; store = 0;
					}
					GetMatrix(gmx);
			}
			else { REPORT }
			Protect();
		}

		// version to work with operator<<
		void GeneralMatrix::Eq(const BaseMatrix& X, MatrixType mt, bool ldok)
		{
			REPORT
				if (ldok) mt.SetDataLossOK();
			Eq(X, mt);
		}

		void GeneralMatrix::Eq2(const BaseMatrix& X, MatrixType mt)
			// a cut down version of Eq for use with += etc.
			// we know BaseMatrix points to two GeneralMatrix objects,
			// the first being this (may be the same).
			// we know tag_val has been set correctly in each.
		{
			GeneralMatrix* gmx = ((BaseMatrix&)X).Evaluate(mt);
			if (gmx!=this) { REPORT GetMatrix(gmx); }  // simplify GetMatrix ?
			else { REPORT }
			Protect();
		}

		void GeneralMatrix::inject(const GeneralMatrix& X)
			// copy stored values of X; otherwise leave els of *this unchanged
		{
			REPORT
				Tracer tr("inject");
			if (nrows_val != X.nrows_val || ncols_val != X.ncols_val)
				Throw(IncompatibleDimensionsException());
			MatrixRow mr((GeneralMatrix*)&X, LoadOnEntry);
			MatrixRow mrx(this, LoadOnEntry+StoreOnExit+DirectPart);
			int i=nrows_val;
			while (i--) { mrx.Inject(mr); mrx.Next(); mr.Next(); }
		}

		// ************* checking for data loss during conversion *******************/

		bool Compare(const MatrixType& source, MatrixType& destination)
		{
			if (!destination) { destination=source; return true; }
			if (destination==source) return true;
			if (!destination.DataLossOK && !(destination>=source))
				Throw(ProgramException("Illegal Conversion", source, destination));
			return false;
		}

		// ************* Make a copy of a matrix on the heap *********************/

		GeneralMatrix* Matrix::Image() const
		{
			REPORT
				GeneralMatrix* gm = new Matrix(*this); MatrixErrorNoSpace(gm);
			return gm;
		}

		GeneralMatrix* SquareMatrix::Image() const
		{
			REPORT
				GeneralMatrix* gm = new SquareMatrix(*this); MatrixErrorNoSpace(gm);
			return gm;
		}

		GeneralMatrix* SymmetricMatrix::Image() const
		{
			REPORT
				GeneralMatrix* gm = new SymmetricMatrix(*this); MatrixErrorNoSpace(gm);
			return gm;
		}

		GeneralMatrix* UpperTriangularMatrix::Image() const
		{
			REPORT
				GeneralMatrix* gm = new UpperTriangularMatrix(*this);
			MatrixErrorNoSpace(gm); return gm;
		}

		GeneralMatrix* LowerTriangularMatrix::Image() const
		{
			REPORT
				GeneralMatrix* gm = new LowerTriangularMatrix(*this);
			MatrixErrorNoSpace(gm); return gm;
		}

		GeneralMatrix* DiagonalMatrix::Image() const
		{
			REPORT
				GeneralMatrix* gm = new DiagonalMatrix(*this); MatrixErrorNoSpace(gm);
			return gm;
		}

		GeneralMatrix* RowVector::Image() const
		{
			REPORT
				GeneralMatrix* gm = new RowVector(*this); MatrixErrorNoSpace(gm);
			return gm;
		}

		GeneralMatrix* ColumnVector::Image() const
		{
			REPORT
				GeneralMatrix* gm = new ColumnVector(*this); MatrixErrorNoSpace(gm);
			return gm;
		}

		GeneralMatrix* nricMatrix::Image() const
		{
			REPORT
				GeneralMatrix* gm = new nricMatrix(*this); MatrixErrorNoSpace(gm);
			return gm;
		}

		GeneralMatrix* IdentityMatrix::Image() const
		{
			REPORT
				GeneralMatrix* gm = new IdentityMatrix(*this); MatrixErrorNoSpace(gm);
			return gm;
		}

		GeneralMatrix* CroutMatrix::Image() const
		{
			REPORT
				GeneralMatrix* gm = new CroutMatrix(*this); MatrixErrorNoSpace(gm);
			return gm;
		}

		GeneralMatrix* GeneralMatrix::Image() const
		{
			bool dummy = true;
			if (dummy)                                   // get rid of warning message
				Throw(InternalException("Cannot apply Image to this matrix type"));
			return 0;
		}


		// *********************** nricMatrix routines *****************************/

		void nricMatrix::MakeRowPointer()
		{
			REPORT
				if (nrows_val > 0)
				{
					row_pointer = new Real* [nrows_val]; MatrixErrorNoSpace(row_pointer);
					Real* s = Store() - 1; int i = nrows_val; Real** rp = row_pointer;
					if (i) for (;;) { *rp++ = s; if (!(--i)) break; s+=ncols_val; }
				}
				else row_pointer = 0;
		}

		void nricMatrix::DeleteRowPointer()
		{ REPORT if (nrows_val) delete [] row_pointer; }

		void GeneralMatrix::CheckStore() const
		{
			REPORT
				if (!store)
					Throw(ProgramException("NRIC accessing matrix with unset dimensions"));
		}


		// *************************** CleanUp routines *****************************/

		void GeneralMatrix::cleanup()
		{
			// set matrix dimensions to zero, delete storage
			REPORT
				if (store && storage)
				{
					MONITOR_REAL_DELETE("Free (cleanup)    ",storage,store)
						REPORT delete [] store;
				}
				store=0; storage=0; nrows_val=0; ncols_val=0; tag_val = -1;
		}

		void nricMatrix::cleanup()
		{ REPORT DeleteRowPointer(); GeneralMatrix::cleanup(); }

		void nricMatrix::MiniCleanUp()
		{ REPORT DeleteRowPointer(); GeneralMatrix::MiniCleanUp(); }

		void RowVector::cleanup()
		{ REPORT GeneralMatrix::cleanup(); nrows_val=1; }

		void ColumnVector::cleanup()
		{ REPORT GeneralMatrix::cleanup(); ncols_val=1; }

		void CroutMatrix::cleanup()
		{
			REPORT
				if (nrows_val) delete [] indx;
			GeneralMatrix::cleanup();
		}

		void CroutMatrix::MiniCleanUp()
		{
			REPORT
				if (nrows_val) delete [] indx;
			GeneralMatrix::MiniCleanUp();
		}

		void BandLUMatrix::cleanup()
		{
			REPORT
				if (nrows_val) delete [] indx;
			if (storage2) delete [] store2;
			GeneralMatrix::cleanup();
		}

		void BandLUMatrix::MiniCleanUp()
		{
			REPORT
				if (nrows_val) delete [] indx;
			if (storage2) delete [] store2;
			GeneralMatrix::MiniCleanUp();
		}

		// ************************ simple integer array class ***********************

		// construct a new array of length xn. Check that xn is non-negative and
		// that space is available

		SimpleIntArray::SimpleIntArray(int xn) : n(xn)
		{
			if (n < 0) Throw(Logic_error("invalid array length"));
			else if (n == 0) { REPORT  a = 0; }
			else { REPORT  a = new int [n]; if (!a) Throw(Bad_alloc()); }
		}

		// destroy an array - return its space to memory

		SimpleIntArray::~SimpleIntArray() { REPORT  if (a) delete [] a; }

		// access an element of an array; return a "reference" so elements
		// can be modified.
		// check index is within range
		// in this array class the index runs from 0 to n-1

		int& SimpleIntArray::operator[](int i)
		{
			REPORT
				if (i < 0 || i >= n) Throw(Logic_error("array index out of range"));
			return a[i];
		}

		// same thing again but for arrays declared constant so we can't
		// modify its elements

		int SimpleIntArray::operator[](int i) const
		{
			REPORT
				if (i < 0 || i >= n) Throw(Logic_error("array index out of range"));
			return a[i];
		}

		// set all the elements equal to a given value

		void SimpleIntArray::operator=(int ai)
		{ REPORT  for (int i = 0; i < n; i++) a[i] = ai; }

		// set the elements equal to those of another array.
		// now allow length to be changed

		void SimpleIntArray::operator=(const SimpleIntArray& b)
		{
			REPORT
				if (b.n != n) resize(b.n);
			for (int i = 0; i < n; i++) a[i] = b.a[i];
		}

		// construct a new array equal to an existing array
		// check that space is available

		SimpleIntArray::SimpleIntArray(const SimpleIntArray& b) : n(b.n)
		{
			if (n == 0) { REPORT  a = 0; }
			else
			{
				REPORT
					a = new int [n]; if (!a) Throw(Bad_alloc());
				for (int i = 0; i < n; i++) a[i] = b.a[i];
			}
		}

		// change the size of an array; optionally copy data from old array to
		// new array

		void SimpleIntArray::resize(int n1, bool keep)
		{
			if (n1 == n) { REPORT  return; }
			else if (n1 == 0) { REPORT  n = 0; delete [] a; a = 0; }
			else if (n == 0)
			{
				REPORT
					a = new int [n1]; if (!a) Throw(Bad_alloc());
				n = n1;
				if (keep) operator=(0);
			}
			else
			{
				int* a1 = a;
				if (keep)
				{
					REPORT
						int i;
					a = new int [n1]; if (!a) Throw(Bad_alloc());
					if (n > n1) n = n1;
					else for (i = n; i < n1; i++) a[i] = 0;
					for (i = 0; i < n; i++) a[i] = a1[i];
					n = n1; delete [] a1;
				}
				else
				{
					REPORT  n = n1; delete [] a1;
					a = new int [n]; if (!a) Throw(Bad_alloc());
				}
			}
		}

		//************************** swap values ********************************

		// swap values

		void GeneralMatrix::swap(GeneralMatrix& gm)
		{
			REPORT
				int t;
			t = tag_val; tag_val = gm.tag_val; gm.tag_val = t;
			t = nrows_val; nrows_val = gm.nrows_val; gm.nrows_val = t;
			t = ncols_val; ncols_val = gm.ncols_val; gm.ncols_val = t;
			t = storage; storage = gm.storage; gm.storage = t;
			Real* s = store; store = gm.store; gm.store = s;
		}

		void nricMatrix::swap(nricMatrix& gm)
		{
			REPORT
				GeneralMatrix::swap((GeneralMatrix&)gm);
			Real** rp = row_pointer; row_pointer = gm.row_pointer; gm.row_pointer = rp;
		}

		void CroutMatrix::swap(CroutMatrix& gm)
		{
			REPORT
				GeneralMatrix::swap((GeneralMatrix&)gm);
			int* i = indx; indx = gm.indx; gm.indx = i;
			bool b;
			b = d; d = gm.d; gm.d = b;
			b = sing; sing = gm.sing; gm.sing = b;
		}

		void BandMatrix::swap(BandMatrix& gm)
		{
			REPORT
				GeneralMatrix::swap((GeneralMatrix&)gm);
			int i;
			i = lower_val; lower_val = gm.lower_val; gm.lower_val = i;
			i = upper_val; upper_val = gm.upper_val; gm.upper_val = i;
		}

		void SymmetricBandMatrix::swap(SymmetricBandMatrix& gm)
		{
			REPORT
				GeneralMatrix::swap((GeneralMatrix&)gm);
			int i;
			i = lower_val; lower_val = gm.lower_val; gm.lower_val = i;
		}

		void BandLUMatrix::swap(BandLUMatrix& gm)
		{
			REPORT
				GeneralMatrix::swap((GeneralMatrix&)gm);
			int* i = indx; indx = gm.indx; gm.indx = i;
			bool b;
			b = d; d = gm.d; gm.d = b;
			b = sing; sing = gm.sing; gm.sing = b;
			int m;
			m = storage2; storage2 = gm.storage2; gm.storage2 = m;
			m = m1; m1 = gm.m1; gm.m1 = m;
			m = m2; m2 = gm.m2; gm.m2 = m;
			Real* s = store2; store2 = gm.store2; gm.store2 = s;
		}

		void GenericMatrix::swap(GenericMatrix& x)
		{
			REPORT
				GeneralMatrix* tgm = gm; gm = x.gm; x.gm = tgm;
		}

		// ********************** C subscript classes ****************************

		RealStarStar::RealStarStar(Matrix& A)
		{
			REPORT
				Tracer tr("RealStarStar");
			int n = A.ncols();
			int m = A.nrows();
			a = new Real*[m];
			MatrixErrorNoSpace(a);
			Real* d = A.data();
			for (int i = 0; i < m; ++i) a[i] = d + i * n;
		} 

		ConstRealStarStar::ConstRealStarStar(const Matrix& A)
		{
			REPORT
				Tracer tr("ConstRealStarStar");
			int n = A.ncols();
			int m = A.nrows();
			a = new const Real*[m];
			MatrixErrorNoSpace(a);
			const Real* d = A.data();
			for (int i = 0; i < m; ++i) a[i] = d + i * n;
		} 



#ifdef use_namespace
	}
#endif

//$$ newmat3.cpp        Matrix get and restore rows and columns

// Copyright (C) 1991,2,3,4: R B Davies

//#define WANT_STREAM


#ifdef use_namespace
namespace NEWMAT {
#endif


#ifdef DO_REPORT
#define REPORT { static ExeCounter ExeCount(__LINE__,3); ++ExeCount; }
#else
#define REPORT {}
#endif

//#define MONITOR(what,storage,store)
//   { cout << what << " " << storage << " at " << (long)store << "\n"; }

#define MONITOR(what,store,storage) {}


// Control bits codes for GetRow, GetCol, RestoreRow, RestoreCol
//
// LoadOnEntry:
//    Load data into MatrixRow or Col dummy array under GetRow or GetCol
// StoreOnExit:
//    Restore data to original matrix under RestoreRow or RestoreCol
// DirectPart:
//    Load or restore only part directly stored; must be set with StoreOnExit
//    Still have decide how to handle this with symmetric
// StoreHere:
//    used in columns only - store data at supplied storage address;
//    used for GetCol, NextCol & RestoreCol. No need to fill out zeros
// HaveStore:
//    dummy array has been assigned (internal use only).

// For symmetric matrices, treat columns as rows unless StoreHere is set;
// then stick to columns as this will give better performance for doing
// inverses

// How components are used:

// Use rows wherever possible in preference to columns

// Columns without StoreHere are used in in-exact transpose, sum column,
// multiply a column vector, and maybe in future access to column,
// additional multiply functions, add transpose

// Columns with StoreHere are used in exact transpose (not symmetric matrices
// or vectors, load only)

// Columns with MatrixColX (Store to full column) are used in inverse and solve

// Functions required for each matrix class

// GetRow(MatrixRowCol& mrc)
// GetCol(MatrixRowCol& mrc)
// GetCol(MatrixColX& mrc)
// RestoreRow(MatrixRowCol& mrc)
// RestoreCol(MatrixRowCol& mrc)
// RestoreCol(MatrixColX& mrc)
// NextRow(MatrixRowCol& mrc)
// NextCol(MatrixRowCol& mrc)
// NextCol(MatrixColX& mrc)

// The Restore routines assume StoreOnExit has already been checked
// Defaults for the Next routines are given below
// Assume cannot have both !DirectPart && StoreHere for MatrixRowCol routines


// Default NextRow and NextCol:
// will work as a default but need to override NextRow for efficiency

void GeneralMatrix::NextRow(MatrixRowCol& mrc)
{
   REPORT
   if (+(mrc.cw*StoreOnExit)) { REPORT this->RestoreRow(mrc); }
   mrc.rowcol++;
   if (mrc.rowcol<nrows_val) { REPORT this->GetRow(mrc); }
   else { REPORT mrc.cw -= StoreOnExit; }
}

void GeneralMatrix::NextCol(MatrixRowCol& mrc)
{
   REPORT                                                // 423
   if (+(mrc.cw*StoreOnExit)) { REPORT this->RestoreCol(mrc); }
   mrc.rowcol++;
   if (mrc.rowcol<ncols_val) { REPORT this->GetCol(mrc); }
   else { REPORT mrc.cw -= StoreOnExit; }
}

void GeneralMatrix::NextCol(MatrixColX& mrc)
{
   REPORT                                                // 423
   if (+(mrc.cw*StoreOnExit)) { REPORT this->RestoreCol(mrc); }
   mrc.rowcol++;
   if (mrc.rowcol<ncols_val) { REPORT this->GetCol(mrc); }
   else { REPORT mrc.cw -= StoreOnExit; }
}


// routines for matrix

void Matrix::GetRow(MatrixRowCol& mrc)
{
   REPORT
   mrc.skip=0; mrc.storage=mrc.length=ncols_val;
   mrc.data=store+mrc.rowcol*ncols_val;
}


void Matrix::GetCol(MatrixRowCol& mrc)
{
   REPORT
   mrc.skip=0; mrc.storage=mrc.length=nrows_val;
   if ( ncols_val==1 && !(mrc.cw*StoreHere) )      // ColumnVector
      { REPORT mrc.data=store; }
   else
   {
      Real* ColCopy;
      if ( !(mrc.cw*(HaveStore+StoreHere)) )
      {
         REPORT
         ColCopy = new Real [nrows_val]; MatrixErrorNoSpace(ColCopy);
         MONITOR_REAL_NEW("Make (MatGetCol)",nrows_val,ColCopy)
         mrc.data = ColCopy; mrc.cw += HaveStore;
      }
      else { REPORT ColCopy = mrc.data; }
      if (+(mrc.cw*LoadOnEntry))
      {
         REPORT
         Real* Mstore = store+mrc.rowcol; int i=nrows_val;
         //while (i--) { *ColCopy++ = *Mstore; Mstore+=ncols_val; }
         if (i) for (;;)
            { *ColCopy++ = *Mstore; if (!(--i)) break; Mstore+=ncols_val; }
      }
   }
}

void Matrix::GetCol(MatrixColX& mrc)
{
   REPORT
   mrc.skip=0; mrc.storage=nrows_val; mrc.length=nrows_val;
   if (+(mrc.cw*LoadOnEntry))
   {
      REPORT  Real* ColCopy = mrc.data;
      Real* Mstore = store+mrc.rowcol; int i=nrows_val;
      //while (i--) { *ColCopy++ = *Mstore; Mstore+=ncols_val; }
      if (i) for (;;)
          { *ColCopy++ = *Mstore; if (!(--i)) break; Mstore+=ncols_val; }
   }
}

void Matrix::RestoreCol(MatrixRowCol& mrc)
{
   // always check StoreOnExit before calling RestoreCol
   REPORT                                   // 429
   if (+(mrc.cw*HaveStore))
   {
      REPORT                                // 426
      Real* Mstore = store+mrc.rowcol; int i=nrows_val;
      Real* Cstore = mrc.data;
      // while (i--) { *Mstore = *Cstore++; Mstore+=ncols_val; }
      if (i) for (;;)
          { *Mstore = *Cstore++; if (!(--i)) break; Mstore+=ncols_val; }
   }
}

void Matrix::RestoreCol(MatrixColX& mrc)
{
   REPORT
   Real* Mstore = store+mrc.rowcol; int i=nrows_val; Real* Cstore = mrc.data;
   // while (i--) { *Mstore = *Cstore++; Mstore+=ncols_val; }
   if (i) for (;;)
      { *Mstore = *Cstore++; if (!(--i)) break; Mstore+=ncols_val; }
}

void Matrix::NextRow(MatrixRowCol& mrc) { REPORT mrc.IncrMat(); }  // 1808

void Matrix::NextCol(MatrixRowCol& mrc)
{
   REPORT                                        // 632
   if (+(mrc.cw*StoreOnExit)) { REPORT RestoreCol(mrc); }
   mrc.rowcol++;
   if (mrc.rowcol<ncols_val)
   {
      if (+(mrc.cw*LoadOnEntry))
      {
         REPORT
         Real* ColCopy = mrc.data;
         Real* Mstore = store+mrc.rowcol; int i=nrows_val;
         //while (i--) { *ColCopy++ = *Mstore; Mstore+=ncols_val; }
         if (i) for (;;)
            { *ColCopy++ = *Mstore; if (!(--i)) break; Mstore+=ncols_val; }
      }
   }
   else { REPORT mrc.cw -= StoreOnExit; }
}

void Matrix::NextCol(MatrixColX& mrc)
{
   REPORT
   if (+(mrc.cw*StoreOnExit)) { REPORT RestoreCol(mrc); }
   mrc.rowcol++;
   if (mrc.rowcol<ncols_val)
   {
      if (+(mrc.cw*LoadOnEntry))
      {
         REPORT
         Real* ColCopy = mrc.data;
         Real* Mstore = store+mrc.rowcol; int i=nrows_val;
         // while (i--) { *ColCopy++ = *Mstore; Mstore+=ncols_val; }
         if (i) for (;;)
            { *ColCopy++ = *Mstore; if (!(--i)) break; Mstore+=ncols_val; }
      }
   }
   else { REPORT mrc.cw -= StoreOnExit; }
}

// routines for diagonal matrix

void DiagonalMatrix::GetRow(MatrixRowCol& mrc)
{
   REPORT
   mrc.skip=mrc.rowcol; mrc.storage=1;
   mrc.data=store+mrc.skip; mrc.length=ncols_val;
}

void DiagonalMatrix::GetCol(MatrixRowCol& mrc)
{
   REPORT
   mrc.skip=mrc.rowcol; mrc.storage=1; mrc.length=nrows_val;
   if (+(mrc.cw*StoreHere))              // should not happen
      Throw(InternalException("DiagonalMatrix::GetCol(MatrixRowCol&)"));
   else  { REPORT mrc.data=store+mrc.skip; }
                                                      // not accessed
}

void DiagonalMatrix::GetCol(MatrixColX& mrc)
{
   REPORT
   mrc.skip=mrc.rowcol; mrc.storage=1; mrc.length=nrows_val;
   mrc.data = mrc.store+mrc.skip;
   *(mrc.data)=*(store+mrc.skip);
}

void DiagonalMatrix::NextRow(MatrixRowCol& mrc) { REPORT mrc.IncrDiag(); }
                        // 800

void DiagonalMatrix::NextCol(MatrixRowCol& mrc) { REPORT mrc.IncrDiag(); }
                        // not accessed

void DiagonalMatrix::NextCol(MatrixColX& mrc)
{
   REPORT
   if (+(mrc.cw*StoreOnExit))
      { REPORT *(store+mrc.rowcol)=*(mrc.data); }
   mrc.IncrDiag();
   int t1 = +(mrc.cw*LoadOnEntry);
   if (t1 && mrc.rowcol < ncols_val)
      { REPORT *(mrc.data)=*(store+mrc.rowcol); }
}

// routines for upper triangular matrix

void UpperTriangularMatrix::GetRow(MatrixRowCol& mrc)
{
   REPORT
   int row = mrc.rowcol; mrc.skip=row; mrc.length=ncols_val;
   mrc.storage=ncols_val-row; mrc.data=store+(row*(2*ncols_val-row+1))/2;
}


void UpperTriangularMatrix::GetCol(MatrixRowCol& mrc)
{
   REPORT
   mrc.skip=0; int i=mrc.rowcol+1; mrc.storage=i;
   mrc.length=nrows_val; Real* ColCopy;
   if ( !(mrc.cw*(StoreHere+HaveStore)) )
   {
      REPORT                                              // not accessed
      ColCopy = new Real [nrows_val]; MatrixErrorNoSpace(ColCopy);
      MONITOR_REAL_NEW("Make (UT GetCol)",nrows_val,ColCopy)
      mrc.data = ColCopy; mrc.cw += HaveStore;
   }
   else { REPORT ColCopy = mrc.data; }
   if (+(mrc.cw*LoadOnEntry))
   {
      REPORT
      Real* Mstore = store+mrc.rowcol; int j = ncols_val;
      // while (i--) { *ColCopy++ = *Mstore; Mstore += --j; }
      if (i) for (;;)
         { *ColCopy++ = *Mstore; if (!(--i)) break; Mstore += --j; }
   }
}

void UpperTriangularMatrix::GetCol(MatrixColX& mrc)
{
   REPORT
   mrc.skip=0; int i=mrc.rowcol+1; mrc.storage=i;
   mrc.length=nrows_val;
   if (+(mrc.cw*LoadOnEntry))
   {
      REPORT
      Real* ColCopy = mrc.data;
      Real* Mstore = store+mrc.rowcol; int j = ncols_val;
      // while (i--) { *ColCopy++ = *Mstore; Mstore += --j; }
      if (i) for (;;)
         { *ColCopy++ = *Mstore; if (!(--i)) break; Mstore += --j; }
   }
}

void UpperTriangularMatrix::RestoreCol(MatrixRowCol& mrc)
{
  REPORT
  Real* Mstore = store+mrc.rowcol; int i=mrc.rowcol+1; int j = ncols_val;
  Real* Cstore = mrc.data;
  // while (i--) { *Mstore = *Cstore++; Mstore += --j; }
  if (i) for (;;)
     { *Mstore = *Cstore++; if (!(--i)) break; Mstore += --j; }
}

void UpperTriangularMatrix::NextRow(MatrixRowCol& mrc) { REPORT mrc.IncrUT(); }
						      // 722

// routines for lower triangular matrix

void LowerTriangularMatrix::GetRow(MatrixRowCol& mrc)
{
   REPORT
   int row=mrc.rowcol; mrc.skip=0; mrc.storage=row+1; mrc.length=ncols_val;
   mrc.data=store+(row*(row+1))/2;
}

void LowerTriangularMatrix::GetCol(MatrixRowCol& mrc)
{
   REPORT
   int col=mrc.rowcol; mrc.skip=col; mrc.length=nrows_val;
   int i=nrows_val-col; mrc.storage=i; Real* ColCopy;
   if ( +(mrc.cw*(StoreHere+HaveStore)) )
      { REPORT  ColCopy = mrc.data; }
   else
   {
      REPORT                                            // not accessed
      ColCopy = new Real [nrows_val]; MatrixErrorNoSpace(ColCopy);
      MONITOR_REAL_NEW("Make (LT GetCol)",nrows_val,ColCopy)
      mrc.cw += HaveStore; mrc.data = ColCopy;
   }

   if (+(mrc.cw*LoadOnEntry))
   {
      REPORT
      Real* Mstore = store+(col*(col+3))/2;
      // while (i--) { *ColCopy++ = *Mstore; Mstore += ++col; }
      if (i) for (;;)
         { *ColCopy++ = *Mstore; if (!(--i)) break; Mstore += ++col; }
   }
}

void LowerTriangularMatrix::GetCol(MatrixColX& mrc)
{
   REPORT
   int col=mrc.rowcol; mrc.skip=col; mrc.length=nrows_val;
   int i=nrows_val-col; mrc.storage=i; mrc.data = mrc.store + col;

   if (+(mrc.cw*LoadOnEntry))
   {
      REPORT  Real* ColCopy = mrc.data;
      Real* Mstore = store+(col*(col+3))/2;
      // while (i--) { *ColCopy++ = *Mstore; Mstore += ++col; }
      if (i) for (;;)
         { *ColCopy++ = *Mstore; if (!(--i)) break; Mstore += ++col; }
   }
}

void LowerTriangularMatrix::RestoreCol(MatrixRowCol& mrc)
{
   REPORT
   int col=mrc.rowcol; Real* Cstore = mrc.data;
   Real* Mstore = store+(col*(col+3))/2; int i=nrows_val-col;
   //while (i--) { *Mstore = *Cstore++; Mstore += ++col; }
   if (i) for (;;)
      { *Mstore = *Cstore++; if (!(--i)) break; Mstore += ++col; }
}

void LowerTriangularMatrix::NextRow(MatrixRowCol& mrc) { REPORT mrc.IncrLT(); }
					                 //712

// routines for symmetric matrix

void SymmetricMatrix::GetRow(MatrixRowCol& mrc)
{
   REPORT                                                //571
   mrc.skip=0; int row=mrc.rowcol; mrc.length=ncols_val;
   if (+(mrc.cw*DirectPart))
      { REPORT mrc.storage=row+1; mrc.data=store+(row*(row+1))/2; }
   else
   {
      // do not allow StoreOnExit and !DirectPart
      if (+(mrc.cw*StoreOnExit))
         Throw(InternalException("SymmetricMatrix::GetRow(MatrixRowCol&)"));
      mrc.storage=ncols_val; Real* RowCopy;
      if (!(mrc.cw*HaveStore))
      {
         REPORT
         RowCopy = new Real [ncols_val]; MatrixErrorNoSpace(RowCopy);
         MONITOR_REAL_NEW("Make (SymGetRow)",ncols_val,RowCopy)
         mrc.cw += HaveStore; mrc.data = RowCopy;
      }
      else { REPORT RowCopy = mrc.data; }
      if (+(mrc.cw*LoadOnEntry))
      {
         REPORT                                         // 544
         Real* Mstore = store+(row*(row+1))/2; int i = row;
         while (i--) *RowCopy++ = *Mstore++;
         i = ncols_val-row;
         // while (i--) { *RowCopy++ = *Mstore; Mstore += ++row; }
         if (i) for (;;)
            { *RowCopy++ = *Mstore; if (!(--i)) break; Mstore += ++row; }
      }
   }
}

void SymmetricMatrix::GetCol(MatrixRowCol& mrc)
{
   // do not allow StoreHere
   if (+(mrc.cw*StoreHere))
      Throw(InternalException("SymmetricMatrix::GetCol(MatrixRowCol&)"));

   int col=mrc.rowcol; mrc.length=nrows_val;
   REPORT
   mrc.skip=0;
   if (+(mrc.cw*DirectPart))    // actually get row ??
      { REPORT mrc.storage=col+1; mrc.data=store+(col*(col+1))/2; }
   else
   {
      // do not allow StoreOnExit and !DirectPart
      if (+(mrc.cw*StoreOnExit))
         Throw(InternalException("SymmetricMatrix::GetCol(MatrixRowCol&)"));

      mrc.storage=ncols_val; Real* ColCopy;
      if ( +(mrc.cw*HaveStore)) { REPORT ColCopy = mrc.data; }
      else
      {
         REPORT                                      // not accessed
         ColCopy = new Real [ncols_val]; MatrixErrorNoSpace(ColCopy);
         MONITOR_REAL_NEW("Make (SymGetCol)",ncols_val,ColCopy)
         mrc.cw += HaveStore; mrc.data = ColCopy;
      }
      if (+(mrc.cw*LoadOnEntry))
      {
         REPORT
         Real* Mstore = store+(col*(col+1))/2; int i = col;
         while (i--) *ColCopy++ = *Mstore++;
         i = ncols_val-col;
         // while (i--) { *ColCopy++ = *Mstore; Mstore += ++col; }
         if (i) for (;;)
            { *ColCopy++ = *Mstore; if (!(--i)) break; Mstore += ++col; }
      }
   }
}

void SymmetricMatrix::GetCol(MatrixColX& mrc)
{
   int col=mrc.rowcol; mrc.length=nrows_val;
   if (+(mrc.cw*DirectPart))
   {
      REPORT
      mrc.skip=col; int i=nrows_val-col; mrc.storage=i;
      mrc.data = mrc.store+col;
      if (+(mrc.cw*LoadOnEntry))
      {
         REPORT                           // not accessed
         Real* ColCopy = mrc.data;
         Real* Mstore = store+(col*(col+3))/2;
         // while (i--) { *ColCopy++ = *Mstore; Mstore += ++col; }
         if (i) for (;;)
            { *ColCopy++ = *Mstore; if (!(--i)) break; Mstore += ++col; }
      }
   }
   else
   {
      REPORT
      // do not allow StoreOnExit and !DirectPart
      if (+(mrc.cw*StoreOnExit))
         Throw(InternalException("SymmetricMatrix::GetCol(MatrixColX&)"));

      mrc.skip=0; mrc.storage=ncols_val;
      if (+(mrc.cw*LoadOnEntry))
      {
         REPORT
         Real* ColCopy = mrc.data;
         Real* Mstore = store+(col*(col+1))/2; int i = col;
         while (i--) *ColCopy++ = *Mstore++;
         i = ncols_val-col;
         // while (i--) { *ColCopy++ = *Mstore; Mstore += ++col; }
         if (i) for (;;)
            { *ColCopy++ = *Mstore; if (!(--i)) break; Mstore += ++col; }
      }
   }
}

// Do not need RestoreRow because we do not allow !DirectPart && StoreOnExit

void SymmetricMatrix::RestoreCol(MatrixColX& mrc)
{
   REPORT
   // Really do restore column
   int col=mrc.rowcol; Real* Cstore = mrc.data;
   Real* Mstore = store+(col*(col+3))/2; int i = nrows_val-col;
   // while (i--) { *Mstore = *Cstore++; Mstore+= ++col; }
   if (i) for (;;)
      { *Mstore = *Cstore++; if (!(--i)) break; Mstore+= ++col; }

}

// routines for row vector

void RowVector::GetCol(MatrixRowCol& mrc)
{
   REPORT
   // do not allow StoreHere
   if (+(mrc.cw*StoreHere))
      Throw(InternalException("RowVector::GetCol(MatrixRowCol&)"));

   mrc.skip=0; mrc.storage=1; mrc.length=nrows_val;
   mrc.data = store+mrc.rowcol;

}

void RowVector::GetCol(MatrixColX& mrc)
{
   REPORT
   mrc.skip=0; mrc.storage=1; mrc.length=nrows_val;
   if (mrc.cw >= LoadOnEntry)
      { REPORT *(mrc.data) = *(store+mrc.rowcol); }

}

void RowVector::NextCol(MatrixRowCol& mrc)
{ REPORT mrc.rowcol++; mrc.data++; }

void RowVector::NextCol(MatrixColX& mrc)
{
   if (+(mrc.cw*StoreOnExit)) { REPORT *(store+mrc.rowcol)=*(mrc.data); }

   mrc.rowcol++;
   if (mrc.rowcol < ncols_val)
   {
      if (+(mrc.cw*LoadOnEntry)) { REPORT *(mrc.data)=*(store+mrc.rowcol); }
   }
   else { REPORT mrc.cw -= StoreOnExit; }
}

void RowVector::RestoreCol(MatrixColX& mrc)
   { REPORT *(store+mrc.rowcol)=*(mrc.data); }           // not accessed


// routines for band matrices

void BandMatrix::GetRow(MatrixRowCol& mrc)
{
   REPORT
   int r = mrc.rowcol; int w = lower_val+1+upper_val; mrc.length=ncols_val;
   int s = r-lower_val;
   if (s<0) { mrc.data = store+(r*w-s); w += s; s = 0; }
   else mrc.data = store+r*w;
   mrc.skip = s; s += w-ncols_val; if (s>0) w -= s; mrc.storage = w;
}

// should make special versions of this for upper and lower band matrices

void BandMatrix::NextRow(MatrixRowCol& mrc)
{
   REPORT
   int r = ++mrc.rowcol;
   if (r<=lower_val) { mrc.storage++; mrc.data += lower_val+upper_val; }
   else  { mrc.skip++; mrc.data += lower_val+upper_val+1; }
   if (r>=ncols_val-upper_val) mrc.storage--;
}

void BandMatrix::GetCol(MatrixRowCol& mrc)
{
   REPORT
   int c = mrc.rowcol; int n = lower_val+upper_val; int w = n+1;
   mrc.length=nrows_val; Real* ColCopy;
   int b; int s = c-upper_val;
   if (s<=0) { w += s; s = 0; b = c+lower_val; } else b = s*w+n;
   mrc.skip = s; s += w-nrows_val; if (s>0) w -= s; mrc.storage = w;
   if ( +(mrc.cw*(StoreHere+HaveStore)) )
      { REPORT ColCopy = mrc.data; }
   else
   {
      REPORT
      ColCopy = new Real [n+1]; MatrixErrorNoSpace(ColCopy);
      MONITOR_REAL_NEW("Make (BMGetCol)",n+1,ColCopy)
      mrc.cw += HaveStore; mrc.data = ColCopy;
   }

   if (+(mrc.cw*LoadOnEntry))
   {
      REPORT
      Real* Mstore = store+b;
      // while (w--) { *ColCopy++ = *Mstore; Mstore+=n; }
      if (w) for (;;)
         { *ColCopy++ = *Mstore; if (!(--w)) break; Mstore+=n; }
   }
}

void BandMatrix::GetCol(MatrixColX& mrc)
{
   REPORT
   int c = mrc.rowcol; int n = lower_val+upper_val; int w = n+1;
   mrc.length=nrows_val; int b; int s = c-upper_val;
   if (s<=0) { w += s; s = 0; b = c+lower_val; } else b = s*w+n;
   mrc.skip = s; s += w-nrows_val; if (s>0) w -= s; mrc.storage = w;
   mrc.data = mrc.store+mrc.skip;

   if (+(mrc.cw*LoadOnEntry))
   {
      REPORT
      Real* ColCopy = mrc.data; Real* Mstore = store+b;
      // while (w--) { *ColCopy++ = *Mstore; Mstore+=n; }
      if (w) for (;;)
         { *ColCopy++ = *Mstore; if (!(--w)) break; Mstore+=n; }
   }
}

void BandMatrix::RestoreCol(MatrixRowCol& mrc)
{
   REPORT
   int c = mrc.rowcol; int n = lower_val+upper_val; int s = c-upper_val;
   Real* Mstore = store + ((s<=0) ? c+lower_val : s*n+s+n);
   Real* Cstore = mrc.data;
   int w = mrc.storage;
   // while (w--) { *Mstore = *Cstore++; Mstore += n; }
   if (w) for (;;)
      { *Mstore = *Cstore++; if (!(--w)) break; Mstore += n; }
}

// routines for symmetric band matrix

void SymmetricBandMatrix::GetRow(MatrixRowCol& mrc)
{
   REPORT
   int r=mrc.rowcol; int s = r-lower_val; int w1 = lower_val+1; int o = r*w1;
   mrc.length = ncols_val;
   if (s<0) { w1 += s; o -= s; s = 0; }
   mrc.skip = s;

   if (+(mrc.cw*DirectPart))
      { REPORT  mrc.data = store+o; mrc.storage = w1; }
   else
   {
      // do not allow StoreOnExit and !DirectPart
      if (+(mrc.cw*StoreOnExit))
         Throw(InternalException("SymmetricBandMatrix::GetRow(MatrixRowCol&)"));
      int w = w1+lower_val; s += w-ncols_val; Real* RowCopy;
      if (s>0) w -= s; mrc.storage = w; int w2 = w-w1;
      if (!(mrc.cw*HaveStore))
      {
         REPORT
         RowCopy = new Real [2*lower_val+1]; MatrixErrorNoSpace(RowCopy);
         MONITOR_REAL_NEW("Make (SmBGetRow)",2*lower_val+1,RowCopy)
         mrc.cw += HaveStore; mrc.data = RowCopy;
      }
      else { REPORT  RowCopy = mrc.data; }

      if (+(mrc.cw*LoadOnEntry) && ncols_val > 0)
      {
         REPORT
         Real* Mstore = store+o;
         while (w1--) *RowCopy++ = *Mstore++;
         Mstore--;
         while (w2--) { Mstore += lower_val; *RowCopy++ = *Mstore; }
      }
   }
}

void SymmetricBandMatrix::GetCol(MatrixRowCol& mrc)
{
   // do not allow StoreHere
   if (+(mrc.cw*StoreHere))
      Throw(InternalException("SymmetricBandMatrix::GetCol(MatrixRowCol&)"));

   int c=mrc.rowcol; int w1 = lower_val+1; mrc.length=nrows_val;
   REPORT
   int s = c-lower_val; int o = c*w1;
   if (s<0) { w1 += s; o -= s; s = 0; }
   mrc.skip = s;

   if (+(mrc.cw*DirectPart))
   { REPORT  mrc.data = store+o; mrc.storage = w1; }
   else
   {
      // do not allow StoreOnExit and !DirectPart
      if (+(mrc.cw*StoreOnExit))
         Throw(InternalException("SymmetricBandMatrix::GetCol(MatrixRowCol&)"));
      int w = w1+lower_val; s += w-ncols_val; Real* ColCopy;
      if (s>0) w -= s; mrc.storage = w; int w2 = w-w1;

      if ( +(mrc.cw*HaveStore) ) { REPORT ColCopy = mrc.data; }
      else
      {
         REPORT ColCopy = new Real [2*lower_val+1]; MatrixErrorNoSpace(ColCopy);
         MONITOR_REAL_NEW("Make (SmBGetCol)",2*lower_val+1,ColCopy)
         mrc.cw += HaveStore; mrc.data = ColCopy;
      }

      if (+(mrc.cw*LoadOnEntry))
      {
         REPORT
         Real* Mstore = store+o;
         while (w1--) *ColCopy++ = *Mstore++;
         Mstore--;
         while (w2--) { Mstore += lower_val; *ColCopy++ = *Mstore; }
      }
   }
}

void SymmetricBandMatrix::GetCol(MatrixColX& mrc)
{
   int c=mrc.rowcol; int w1 = lower_val+1; mrc.length=nrows_val;
   if (+(mrc.cw*DirectPart))
   {
      REPORT
      int b = c*w1+lower_val;
      mrc.skip = c; c += w1-nrows_val; w1 -= c; mrc.storage = w1;
      Real* ColCopy = mrc.data = mrc.store+mrc.skip;
      if (+(mrc.cw*LoadOnEntry))
      {
         REPORT
         Real* Mstore = store+b;
         // while (w1--) { *ColCopy++ = *Mstore; Mstore += lower; }
         if (w1) for (;;)
            { *ColCopy++ = *Mstore; if (!(--w1)) break; Mstore += lower_val; }
      }
   }
   else
   {
      REPORT
      // do not allow StoreOnExit and !DirectPart
      if (+(mrc.cw*StoreOnExit))
         Throw(InternalException("SymmetricBandMatrix::GetCol(MatrixColX&)"));
      int s = c-lower_val; int o = c*w1;
      if (s<0) { w1 += s; o -= s; s = 0; }
      mrc.skip = s;

      int w = w1+lower_val; s += w-ncols_val;
      if (s>0) w -= s; mrc.storage = w; int w2 = w-w1;

      Real* ColCopy = mrc.data = mrc.store+mrc.skip;

      if (+(mrc.cw*LoadOnEntry))
      {
         REPORT
         Real* Mstore = store+o;
         while (w1--) *ColCopy++ = *Mstore++;
         Mstore--;
         while (w2--) { Mstore += lower_val; *ColCopy++ = *Mstore; }
      }

   }
}

void SymmetricBandMatrix::RestoreCol(MatrixColX& mrc)
{
   REPORT
   int c = mrc.rowcol;
   Real* Mstore = store + c*lower_val+c+lower_val;
   Real* Cstore = mrc.data; int w = mrc.storage;
   // while (w--) { *Mstore = *Cstore++; Mstore += lower_val; }
   if (w) for (;;)
      { *Mstore = *Cstore++; if (!(--w)) break; Mstore += lower_val; }
}

// routines for identity matrix

void IdentityMatrix::GetRow(MatrixRowCol& mrc)
{
   REPORT
   mrc.skip=mrc.rowcol; mrc.storage=1; mrc.data=store; mrc.length=ncols_val;
}

void IdentityMatrix::GetCol(MatrixRowCol& mrc)
{
   REPORT
   mrc.skip=mrc.rowcol; mrc.storage=1; mrc.length=nrows_val;
   if (+(mrc.cw*StoreHere))              // should not happen
      Throw(InternalException("IdentityMatrix::GetCol(MatrixRowCol&)"));
   else  { REPORT mrc.data=store; }
}

void IdentityMatrix::GetCol(MatrixColX& mrc)
{
   REPORT
   mrc.skip=mrc.rowcol; mrc.storage=1; mrc.length=nrows_val;
   mrc.data = mrc.store+mrc.skip; *(mrc.data)=*store;
}

void IdentityMatrix::NextRow(MatrixRowCol& mrc) { REPORT mrc.IncrId(); }

void IdentityMatrix::NextCol(MatrixRowCol& mrc) { REPORT mrc.IncrId(); }

void IdentityMatrix::NextCol(MatrixColX& mrc)
{
   REPORT
   if (+(mrc.cw*StoreOnExit)) { REPORT *store=*(mrc.data); }
   mrc.IncrDiag();            // must increase mrc.data so need IncrDiag
   int t1 = +(mrc.cw*LoadOnEntry);
   if (t1 && mrc.rowcol < ncols_val) { REPORT *(mrc.data)=*store; }
}




// *************************** destructors *******************************

MatrixRowCol::~MatrixRowCol()
{
   if (+(cw*HaveStore))
   {
      MONITOR_REAL_DELETE("Free    (RowCol)",-1,data)  // do not know length
      delete [] data;
   }
}

MatrixRow::~MatrixRow() { if (+(cw*StoreOnExit)) gm->RestoreRow(*this); }

MatrixCol::~MatrixCol() { if (+(cw*StoreOnExit)) gm->RestoreCol(*this); }

MatrixColX::~MatrixColX() { if (+(cw*StoreOnExit)) gm->RestoreCol(*this); }

#ifdef use_namespace
}
#endif



/************************* MatrixType functions *****************************/


// Skew needs more work <<<<<<<<<

// all operations to return MatrixTypes which correspond to valid types
// of matrices.
// Eg: if it has the Diagonal attribute, then it must also have
// Upper, Lower, Band, Square and Symmetric.


MatrixType MatrixType::operator*(const MatrixType& mt) const
{
   REPORT
   int a = attribute & mt.attribute & ~(Symmetric | Skew);
   a |= (a & Diagonal) * 63;                   // recognise diagonal
   return MatrixType(a);
}

MatrixType MatrixType::SP(const MatrixType& mt) const
// elementwise product
// Lower, Upper, Diag, Band if only one is
// Symmetric, Ones, Valid (and Real) if both are
// Need to include Lower & Upper => Diagonal
// Will need to include both Skew => Symmetric
{
   REPORT
   int a = ((attribute | mt.attribute) & ~(Symmetric + Skew + Valid + Ones))
      | (attribute & mt.attribute);
   if ((a & Lower) != 0  &&  (a & Upper) != 0) a |= Diagonal;
   if ((attribute & Skew) != 0)
   {
      if ((mt.attribute & Symmetric) != 0) a |= Skew;  
      if ((mt.attribute & Skew) != 0) { a &= ~Skew; a |= Symmetric; }
   }
   else if ((mt.attribute & Skew) != 0 && (attribute & Symmetric) != 0)
      a |= Skew;  
   a |= (a & Diagonal) * 63;                   // recognise diagonal
   return MatrixType(a);
}

MatrixType MatrixType::KP(const MatrixType& mt) const
// Kronecker product
// Lower, Upper, Diag, Symmetric, Band, Valid if both are
// Band if LHS is band & other is square 
// Ones is complicated so leave this out
{
   REPORT
   int a = (attribute & mt.attribute)  & ~Ones;
   if ((attribute & Band) != 0 && (mt.attribute & Square) != 0)
      a |= Band;
   //int a = ((attribute & mt.attribute) | (attribute & Band)) & ~Ones;

   return MatrixType(a);
}

MatrixType MatrixType::i() const               // type of inverse
{
   REPORT
   int a = attribute & ~(Band+LUDeco);
   a |= (a & Diagonal) * 63;                   // recognise diagonal
   return MatrixType(a);
}

MatrixType MatrixType::t() const
// swap lower and upper attributes
// assume Upper is in bit above Lower
{
   REPORT
   int a = attribute;
   a ^= (((a >> 1) ^ a) & Lower) * 3;
   return MatrixType(a);
}

MatrixType MatrixType::MultRHS() const
{
   REPORT
   // remove symmetric attribute unless diagonal
   return (attribute >= Dg) ? attribute : (attribute & ~Symmetric);
}

// this is used for deciding type of multiplication
bool Rectangular(MatrixType a, MatrixType b, MatrixType c)
{
   REPORT
   return
      ((a.attribute | b.attribute | c.attribute)
      & ~(MatrixType::Valid | MatrixType::Square)) == 0;
}

const char* MatrixType::value() const
{
// make a string with the name of matrix with the given attributes
   switch (attribute)
   {
   case Valid:                              REPORT return "Rect ";
   case Valid+Square:                       REPORT return "Squ  ";
   case Valid+Symmetric+Square:             REPORT return "Sym  ";
   case Valid+Skew+Square:                  REPORT return "Skew ";
   case Valid+Band+Square:                  REPORT return "Band ";
   case Valid+Symmetric+Band+Square:        REPORT return "SmBnd";
   case Valid+Skew+Band+Square:             REPORT return "SkBnd";
   case Valid+Upper+Square:                 REPORT return "UT   ";
   case Valid+Diagonal+Symmetric+Band+Upper+Lower+Square:
                                            REPORT return "Diag ";
   case Valid+Diagonal+Symmetric+Band+Upper+Lower+Ones+Square:
                                            REPORT return "Ident";
   case Valid+Band+Upper+Square:            REPORT return "UpBnd";
   case Valid+Lower+Square:                 REPORT return "LT   ";
   case Valid+Band+Lower+Square:            REPORT return "LwBnd";
   default:
      REPORT
      if (!(attribute & Valid))             return "UnSp ";
      if (attribute & LUDeco)
         return (attribute & Band) ?     "BndLU" : "Crout";
                                            return "?????";
   }
}


GeneralMatrix* MatrixType::New(int nr, int nc, BaseMatrix* bm) const
{
// make a new matrix with the given attributes

   Tracer tr("New"); GeneralMatrix* gm=0;   // initialised to keep gnu happy
   switch (attribute)
   {
   case Valid:
      REPORT
      if (nc==1) { gm = new ColumnVector(nr); break; }
      if (nr==1) { gm = new RowVector(nc); break; }
      gm = new Matrix(nr, nc); break;

   case Valid+Square:
      REPORT
      if (nc!=nr) { Throw(NotSquareException()); }
      gm = new SquareMatrix(nr); break;

   case Valid+Symmetric+Square:
      REPORT gm = new SymmetricMatrix(nr); break;

   case Valid+Band+Square:
      {
         REPORT
         MatrixBandWidth bw = bm->bandwidth();
         gm = new BandMatrix(nr,bw.lower_val,bw.upper_val); break;
      }

   case Valid+Symmetric+Band+Square:
      REPORT gm = new SymmetricBandMatrix(nr,bm->bandwidth().lower_val); break;

   case Valid+Upper+Square:
      REPORT gm = new UpperTriangularMatrix(nr); break;

   case Valid+Diagonal+Symmetric+Band+Upper+Lower+Square:
      REPORT gm = new DiagonalMatrix(nr); break;

   case Valid+Band+Upper+Square:
      REPORT gm = new UpperBandMatrix(nr,bm->bandwidth().upper_val); break;

   case Valid+Lower+Square:
      REPORT gm = new LowerTriangularMatrix(nr); break;

   case Valid+Band+Lower+Square:
      REPORT gm = new LowerBandMatrix(nr,bm->bandwidth().lower_val); break;

   case Valid+Diagonal+Symmetric+Band+Upper+Lower+Ones+Square:
      REPORT gm = new IdentityMatrix(nr); break;

   default:
      Throw(ProgramException("Invalid matrix type"));
   }
   
   MatrixErrorNoSpace(gm); gm->Protect(); return gm;
}


#ifdef use_namespace
}
#endif

//$$ newmat2.cpp      Matrix row and column operations

// Copyright (C) 1991,2,3,4: R B Davies

#define WANT_MATH


#ifdef use_namespace
namespace NEWMAT {
#endif


#ifdef DO_REPORT
#define REPORT { static ExeCounter ExeCount(__LINE__,2); ++ExeCount; }
#else
#define REPORT {}
#endif

	//#define MONITOR(what,storage,store) { cout << what << " " << storage << " at " << (long)store << "\n"; }

#define MONITOR(what,store,storage) {}

	/************************** Matrix Row/Col functions ************************/

	void MatrixRowCol::Add(const MatrixRowCol& mrc)
	{
		// THIS += mrc
		REPORT
			int f = mrc.skip; int l = f + mrc.storage; int lx = skip + storage;
		if (f < skip) f = skip; if (l > lx) l = lx; l -= f;
		if (l<=0) return;
		Real* elx=data+(f-skip); Real* el=mrc.data+(f-mrc.skip);
		while (l--) *elx++ += *el++;
	}

	void MatrixRowCol::AddScaled(const MatrixRowCol& mrc, Real x)
	{
		REPORT
			// THIS += (mrc * x)
			int f = mrc.skip; int l = f + mrc.storage; int lx = skip + storage;
		if (f < skip) f = skip; if (l > lx) l = lx; l -= f;
		if (l<=0) return;
		Real* elx=data+(f-skip); Real* el=mrc.data+(f-mrc.skip);
		while (l--) *elx++ += *el++ * x;
	}

	void MatrixRowCol::Sub(const MatrixRowCol& mrc)
	{
		REPORT
			// THIS -= mrc
			int f = mrc.skip; int l = f + mrc.storage; int lx = skip + storage;
		if (f < skip) f = skip; if (l > lx) l = lx; l -= f;
		if (l<=0) return;
		Real* elx=data+(f-skip); Real* el=mrc.data+(f-mrc.skip);
		while (l--) *elx++ -= *el++;
	}

	void MatrixRowCol::Inject(const MatrixRowCol& mrc)
		// copy stored elements only
	{
		REPORT
			int f = mrc.skip; int l = f + mrc.storage; int lx = skip + storage;
		if (f < skip) f = skip; if (l > lx) l = lx; l -= f;
		if (l<=0) return;
		Real* elx=data+(f-skip); Real* ely=mrc.data+(f-mrc.skip);
		while (l--) *elx++ = *ely++;
	}

	Real DotProd(const MatrixRowCol& mrc1, const MatrixRowCol& mrc2)
	{
		REPORT                                         // not accessed
			int f = mrc1.skip; int f2 = mrc2.skip;
		int l = f + mrc1.storage; int l2 = f2 + mrc2.storage;
		if (f < f2) f = f2; if (l > l2) l = l2; l -= f;
		if (l<=0) return 0.0;

		Real* el1=mrc1.data+(f-mrc1.skip); Real* el2=mrc2.data+(f-mrc2.skip);
		Real sum = 0.0;
		while (l--) sum += *el1++ * *el2++;
		return sum;
	}

	void MatrixRowCol::Add(const MatrixRowCol& mrc1, const MatrixRowCol& mrc2)
	{
		// THIS = mrc1 + mrc2
		int f = skip; int l = skip + storage;
		int f1 = mrc1.skip; int l1 = f1 + mrc1.storage;
		if (f1<f) f1=f; if (l1>l) l1=l;
		int f2 = mrc2.skip; int l2 = f2 + mrc2.storage;
		if (f2<f) f2=f; if (l2>l) l2=l;
		Real* el = data + (f-skip);
		Real* el1 = mrc1.data+(f1-mrc1.skip); Real* el2 = mrc2.data+(f2-mrc2.skip);
		if (f1<f2)
		{
			int i = f1-f; while (i--) *el++ = 0.0;
			if (l1<=f2)                              // disjoint
			{
				REPORT                                // not accessed
					i = l1-f1;     while (i--) *el++ = *el1++;
				i = f2-l1;     while (i--) *el++ = 0.0;
				i = l2-f2;     while (i--) *el++ = *el2++;
				i = l-l2;      while (i--) *el++ = 0.0;
			}
			else
			{
				i = f2-f1;    while (i--) *el++ = *el1++;
				if (l1<=l2)
				{
					REPORT
						i = l1-f2; while (i--) *el++ = *el1++ + *el2++;
					i = l2-l1; while (i--) *el++ = *el2++;
					i = l-l2;  while (i--) *el++ = 0.0;
				}
				else
				{
					REPORT
						i = l2-f2; while (i--) *el++ = *el1++ + *el2++;
					i = l1-l2; while (i--) *el++ = *el1++;
					i = l-l1;  while (i--) *el++ = 0.0;
				}
			}
		}
		else
		{
			int i = f2-f; while (i--) *el++ = 0.0;
			if (l2<=f1)                              // disjoint
			{
				REPORT                                // not accessed
					i = l2-f2;     while (i--) *el++ = *el2++;
				i = f1-l2;     while (i--) *el++ = 0.0;
				i = l1-f1;     while (i--) *el++ = *el1++;
				i = l-l1;      while (i--) *el++ = 0.0;
			}
			else
			{
				i = f1-f2;    while (i--) *el++ = *el2++;
				if (l2<=l1)
				{
					REPORT
						i = l2-f1; while (i--) *el++ = *el1++ + *el2++;
					i = l1-l2; while (i--) *el++ = *el1++;
					i = l-l1;  while (i--) *el++ = 0.0;
				}
				else
				{
					REPORT
						i = l1-f1; while (i--) *el++ = *el1++ + *el2++;
					i = l2-l1; while (i--) *el++ = *el2++;
					i = l-l2;  while (i--) *el++ = 0.0;
				}
			}
		}
	}

	void MatrixRowCol::Sub(const MatrixRowCol& mrc1, const MatrixRowCol& mrc2)
	{
		// THIS = mrc1 - mrc2
		int f = skip; int l = skip + storage;
		int f1 = mrc1.skip; int l1 = f1 + mrc1.storage;
		if (f1<f) f1=f; if (l1>l) l1=l;
		int f2 = mrc2.skip; int l2 = f2 + mrc2.storage;
		if (f2<f) f2=f; if (l2>l) l2=l;
		Real* el = data + (f-skip);
		Real* el1 = mrc1.data+(f1-mrc1.skip); Real* el2 = mrc2.data+(f2-mrc2.skip);
		if (f1<f2)
		{
			int i = f1-f; while (i--) *el++ = 0.0;
			if (l1<=f2)                              // disjoint
			{
				REPORT                                // not accessed
					i = l1-f1;     while (i--) *el++ = *el1++;
				i = f2-l1;     while (i--) *el++ = 0.0;
				i = l2-f2;     while (i--) *el++ = - *el2++;
				i = l-l2;      while (i--) *el++ = 0.0;
			}
			else
			{
				i = f2-f1;    while (i--) *el++ = *el1++;
				if (l1<=l2)
				{
					REPORT
						i = l1-f2; while (i--) *el++ = *el1++ - *el2++;
					i = l2-l1; while (i--) *el++ = - *el2++;
					i = l-l2;  while (i--) *el++ = 0.0;
				}
				else
				{
					REPORT
						i = l2-f2; while (i--) *el++ = *el1++ - *el2++;
					i = l1-l2; while (i--) *el++ = *el1++;
					i = l-l1;  while (i--) *el++ = 0.0;
				}
			}
		}
		else
		{
			int i = f2-f; while (i--) *el++ = 0.0;
			if (l2<=f1)                              // disjoint
			{
				REPORT                                // not accessed
					i = l2-f2;     while (i--) *el++ = - *el2++;
				i = f1-l2;     while (i--) *el++ = 0.0;
				i = l1-f1;     while (i--) *el++ = *el1++;
				i = l-l1;      while (i--) *el++ = 0.0;
			}
			else
			{
				i = f1-f2;    while (i--) *el++ = - *el2++;
				if (l2<=l1)
				{
					REPORT
						i = l2-f1; while (i--) *el++ = *el1++ - *el2++;
					i = l1-l2; while (i--) *el++ = *el1++;
					i = l-l1;  while (i--) *el++ = 0.0;
				}
				else
				{
					REPORT
						i = l1-f1; while (i--) *el++ = *el1++ - *el2++;
					i = l2-l1; while (i--) *el++ = - *el2++;
					i = l-l2;  while (i--) *el++ = 0.0;
				}
			}
		}
	}


	void MatrixRowCol::Add(const MatrixRowCol& mrc1, Real x)
	{
		// THIS = mrc1 + x
		REPORT
			if (!storage) return;
		int f = mrc1.skip; int l = f + mrc1.storage; int lx = skip + storage;
		if (f < skip) { f = skip; if (l < f) l = f; }
		if (l > lx) { l = lx; if (f > lx) f = lx; }

		Real* elx = data; Real* ely = mrc1.data+(f-mrc1.skip);

		int l1 = f-skip;  while (l1--) *elx++ = x;
		l1 = l-f;     while (l1--) *elx++ = *ely++ + x;
		lx -= l;      while (lx--) *elx++ = x;
	}

	void MatrixRowCol::NegAdd(const MatrixRowCol& mrc1, Real x)
	{
		// THIS = x - mrc1
		REPORT
			if (!storage) return;
		int f = mrc1.skip; int l = f + mrc1.storage; int lx = skip + storage;
		if (f < skip) { f = skip; if (l < f) l = f; }
		if (l > lx) { l = lx; if (f > lx) f = lx; }

		Real* elx = data; Real* ely = mrc1.data+(f-mrc1.skip);

		int l1 = f-skip;  while (l1--) *elx++ = x;
		l1 = l-f;     while (l1--) *elx++ = x - *ely++;
		lx -= l;      while (lx--) *elx++ = x;
	}

	void MatrixRowCol::RevSub(const MatrixRowCol& mrc1)
	{
		// THIS = mrc - THIS
		REPORT
			if (!storage) return;
		int f = mrc1.skip; int l = f + mrc1.storage; int lx = skip + storage;
		if (f < skip) { f = skip; if (l < f) l = f; }
		if (l > lx) { l = lx; if (f > lx) f = lx; }

		Real* elx = data; Real* ely = mrc1.data+(f-mrc1.skip);

		int l1 = f-skip;  while (l1--) { *elx = - *elx; elx++; }
		l1 = l-f;     while (l1--) { *elx = *ely++ - *elx; elx++; }
		lx -= l;      while (lx--) { *elx = - *elx; elx++; }
	}

	void MatrixRowCol::ConCat(const MatrixRowCol& mrc1, const MatrixRowCol& mrc2)
	{
		// THIS = mrc1 | mrc2
		REPORT
			int f1 = mrc1.skip; int l1 = f1 + mrc1.storage; int lx = skip + storage;
		if (f1 < skip) { f1 = skip; if (l1 < f1) l1 = f1; }
		if (l1 > lx) { l1 = lx; if (f1 > lx) f1 = lx; }

		Real* elx = data;

		int i = f1-skip;  while (i--) *elx++ =0.0;
		i = l1-f1;
		if (i)                       // in case f1 would take ely out of range
		{ Real* ely = mrc1.data+(f1-mrc1.skip);  while (i--) *elx++ = *ely++; }

		int f2 = mrc2.skip; int l2 = f2 + mrc2.storage; i = mrc1.length;
		int skipx = l1 - i; lx -= i; // addresses rel to second seg, maybe -ve
		if (f2 < skipx) { f2 = skipx; if (l2 < f2) l2 = f2; }
		if (l2 > lx) { l2 = lx; if (f2 > lx) f2 = lx; }

		i = f2-skipx; while (i--) *elx++ = 0.0;
		i = l2-f2;
		if (i)                       // in case f2 would take ely out of range
		{ Real* ely = mrc2.data+(f2-mrc2.skip); while (i--) *elx++ = *ely++; }
		lx -= l2;     while (lx--) *elx++ = 0.0;
	}

	void MatrixRowCol::Multiply(const MatrixRowCol& mrc1)
		// element by element multiply into
	{
		REPORT
			if (!storage) return;
		int f = mrc1.skip; int l = f + mrc1.storage; int lx = skip + storage;
		if (f < skip) { f = skip; if (l < f) l = f; }
		if (l > lx) { l = lx; if (f > lx) f = lx; }

		Real* elx = data; Real* ely = mrc1.data+(f-mrc1.skip);

		int l1 = f-skip;  while (l1--) *elx++ = 0;
		l1 = l-f;     while (l1--) *elx++ *= *ely++;
		lx -= l;      while (lx--) *elx++ = 0;
	}

	void MatrixRowCol::Multiply(const MatrixRowCol& mrc1, const MatrixRowCol& mrc2)
		// element by element multiply
	{
		int f = skip; int l = skip + storage;
		int f1 = mrc1.skip; int l1 = f1 + mrc1.storage;
		if (f1<f) f1=f; if (l1>l) l1=l;
		int f2 = mrc2.skip; int l2 = f2 + mrc2.storage;
		if (f2<f) f2=f; if (l2>l) l2=l;
		Real* el = data + (f-skip); int i;
		if (f1<f2) f1 = f2; if (l1>l2) l1 = l2;
		if (l1<=f1) { REPORT i = l-f; while (i--) *el++ = 0.0; }  // disjoint
		else
		{
			REPORT
				Real* el1 = mrc1.data+(f1-mrc1.skip);
			Real* el2 = mrc2.data+(f1-mrc2.skip);
			i = f1-f ;    while (i--) *el++ = 0.0;
			i = l1-f1;    while (i--) *el++ = *el1++ * *el2++;
			i = l-l1;     while (i--) *el++ = 0.0;
		}
	}

	void MatrixRowCol::KP(const MatrixRowCol& mrc1, const MatrixRowCol& mrc2)
		// row for Kronecker product
	{
		int f = skip; int s = storage; Real* el = data; int i;

		i = mrc1.skip * mrc2.length;
		if (i > f)
		{
			i -= f; f = 0; if (i > s) { i = s; s = 0; }  else s -= i;
			while (i--) *el++ = 0.0;
			if (s == 0) return;
		}
		else f -= i;

		i = mrc1.storage; Real* el1 = mrc1.data;
		int mrc2_skip = mrc2.skip; int mrc2_storage = mrc2.storage;
		int mrc2_length = mrc2.length;
		int mrc2_remain = mrc2_length - mrc2_skip - mrc2_storage;
		while (i--)
		{
			int j; Real* el2 = mrc2.data; Real vel1 = *el1;
			if (f == 0 && mrc2_length <= s)
			{
				j = mrc2_skip; s -= j;    while (j--) *el++ = 0.0;
				j = mrc2_storage; s -= j; while (j--) *el++ = vel1 * *el2++;
				j = mrc2_remain; s -= j;  while (j--) *el++ = 0.0;
			}
			else if (f >= mrc2_length) f -= mrc2_length;
			else
			{
				j = mrc2_skip;
				if (j > f)
				{
					j -= f; f = 0; if (j > s) { j = s; s = 0; } else s -= j;
					while (j--) *el++ = 0.0;
				}
				else f -= j;

				j = mrc2_storage;
				if (j > f)
				{
					j -= f; el2 += f; f = 0; if (j > s) { j = s; s = 0; } else s -= j;
					while (j--) *el++ = vel1 * *el2++;
				}
				else f -= j;

				j = mrc2_remain;
				if (j > f)
				{
					j -= f; f = 0; if (j > s) { j = s; s = 0; } else s -= j;
					while (j--) *el++ = 0.0;
				}
				else f -= j;
			}
			if (s == 0) return;
			++el1;
		}

		i = (mrc1.length - mrc1.skip - mrc1.storage) * mrc2.length;
		if (i > f)
		{
			i -= f; if (i > s) i = s;
			while (i--) *el++ = 0.0;
		}
	}


	void MatrixRowCol::Copy(const MatrixRowCol& mrc1)
	{
		// THIS = mrc1
		REPORT
			if (!storage) return;
		int f = mrc1.skip; int l = f + mrc1.storage; int lx = skip + storage;
		if (f < skip) { f = skip; if (l < f) l = f; }
		if (l > lx) { l = lx; if (f > lx) f = lx; }

		Real* elx = data; Real* ely = 0;

		if (l-f) ely = mrc1.data+(f-mrc1.skip);

		int l1 = f-skip;  while (l1--) *elx++ = 0.0;
		l1 = l-f;     while (l1--) *elx++ = *ely++;
		lx -= l;      while (lx--) *elx++ = 0.0;
	}

	void MatrixRowCol::CopyCheck(const MatrixRowCol& mrc1)
		// Throw an exception if this would lead to a loss of data
	{
		REPORT
			if (!storage) return;
		int f = mrc1.skip; int l = f + mrc1.storage; int lx = skip + storage;
		if (f < skip || l > lx) Throw(ProgramException("Illegal Conversion"));

		Real* elx = data; Real* ely = mrc1.data+(f-mrc1.skip);

		int l1 = f-skip;  while (l1--) *elx++ = 0.0;
		l1 = l-f;     while (l1--) *elx++ = *ely++;
		lx -= l;      while (lx--) *elx++ = 0.0;
	}

	void MatrixRowCol::Check(const MatrixRowCol& mrc1)
		// Throw an exception if +=, -=, copy etc would lead to a loss of data
	{
		REPORT
			int f = mrc1.skip; int l = f + mrc1.storage; int lx = skip + storage;
		if (f < skip || l > lx) Throw(ProgramException("Illegal Conversion"));
	}

	void MatrixRowCol::Check()
		// Throw an exception if +=, -= of constant would lead to a loss of data
		// that is: check full row is present
		// may not be appropriate for symmetric matrices
	{
		REPORT
			if (skip!=0 || storage!=length)
				Throw(ProgramException("Illegal Conversion"));
	}

	void MatrixRowCol::Negate(const MatrixRowCol& mrc1)
	{
		// THIS = -mrc1
		REPORT
			if (!storage) return;
		int f = mrc1.skip; int l = f + mrc1.storage; int lx = skip + storage;
		if (f < skip) { f = skip; if (l < f) l = f; }
		if (l > lx) { l = lx; if (f > lx) f = lx; }

		Real* elx = data; Real* ely = mrc1.data+(f-mrc1.skip);

		int l1 = f-skip;  while (l1--) *elx++ = 0.0;
		l1 = l-f;     while (l1--) *elx++ = - *ely++;
		lx -= l;      while (lx--) *elx++ = 0.0;
	}

	void MatrixRowCol::Multiply(const MatrixRowCol& mrc1, Real s)
	{
		// THIS = mrc1 * s
		REPORT
			if (!storage) return;
		int f = mrc1.skip; int l = f + mrc1.storage; int lx = skip + storage;
		if (f < skip) { f = skip; if (l < f) l = f; }
		if (l > lx) { l = lx; if (f > lx) f = lx; }

		Real* elx = data; Real* ely = mrc1.data+(f-mrc1.skip);

		int l1 = f-skip;  while (l1--) *elx++ = 0.0;
		l1 = l-f;     while (l1--) *elx++ = *ely++ * s;
		lx -= l;      while (lx--) *elx++ = 0.0;
	}

	void DiagonalMatrix::Solver(MatrixColX& mrc, const MatrixColX& mrc1)
	{
		// mrc = mrc / mrc1   (elementwise)
		REPORT
			int f = mrc1.skip; int f0 = mrc.skip;
		int l = f + mrc1.storage; int lx = f0 + mrc.storage;
		if (f < f0) { f = f0; if (l < f) l = f; }
		if (l > lx) { l = lx; if (f > lx) f = lx; }

		Real* elx = mrc.data; Real* eld = store+f;

		int l1 = f-f0;    while (l1--) *elx++ = 0.0;
		l1 = l-f;     while (l1--) *elx++ /= *eld++;
		lx -= l;      while (lx--) *elx++ = 0.0;
		// Solver makes sure input and output point to same memory
	}

	void IdentityMatrix::Solver(MatrixColX& mrc, const MatrixColX& mrc1)
	{
		// mrc = mrc / mrc1   (elementwise)
		REPORT
			int f = mrc1.skip; int f0 = mrc.skip;
		int l = f + mrc1.storage; int lx = f0 + mrc.storage;
		if (f < f0) { f = f0; if (l < f) l = f; }
		if (l > lx) { l = lx; if (f > lx) f = lx; }

		Real* elx = mrc.data; Real eldv = *store;

		int l1 = f-f0;    while (l1--) *elx++ = 0.0;
		l1 = l-f;     while (l1--) *elx++ /= eldv;
		lx -= l;      while (lx--) *elx++ = 0.0;
		// Solver makes sure input and output point to same memory
	}

	void MatrixRowCol::Copy(const double*& r)
	{
		// THIS = *r
		REPORT
			Real* elx = data; const double* ely = r+skip; r += length;
		int l = storage; while (l--) *elx++ = (Real)*ely++;
	}

	void MatrixRowCol::Copy(const float*& r)
	{
		// THIS = *r
		REPORT
			Real* elx = data; const float* ely = r+skip; r += length;
		int l = storage; while (l--) *elx++ = (Real)*ely++;
	}

	void MatrixRowCol::Copy(const int*& r)
	{
		// THIS = *r
		REPORT
			Real* elx = data; const int* ely = r+skip; r += length;
		int l = storage; while (l--) *elx++ = (Real)*ely++;
	}

	void MatrixRowCol::Copy(Real r)
	{
		// THIS = r
		REPORT  Real* elx = data; int l = storage; while (l--) *elx++ = r;
	}

	void MatrixRowCol::Zero()
	{
		// THIS = 0
		REPORT  Real* elx = data; int l = storage; while (l--) *elx++ = 0;
	}

	void MatrixRowCol::Multiply(Real r)
	{
		// THIS *= r
		REPORT  Real* elx = data; int l = storage; while (l--) *elx++ *= r;
	}

	void MatrixRowCol::Add(Real r)
	{
		// THIS += r
		REPORT
			Real* elx = data; int l = storage; while (l--) *elx++ += r;
	}

	Real MatrixRowCol::SumAbsoluteValue()
	{
		REPORT
			Real sum = 0.0; Real* elx = data; int l = storage;
		while (l--) sum += fabs(*elx++);
		return sum;
	}

	// max absolute value of r and elements of row/col
	// we use <= or >= in all of these so we are sure of getting
	// r reset at least once.
	Real MatrixRowCol::MaximumAbsoluteValue1(Real r, int& i)
	{
		REPORT
			Real* elx = data; int l = storage; int li = -1;
		while (l--) { Real f = fabs(*elx++); if (r <= f) { r = f; li = l; } }
		i = (li >= 0) ? storage - li + skip : 0;
		return r;
	}

	// min absolute value of r and elements of row/col
	Real MatrixRowCol::MinimumAbsoluteValue1(Real r, int& i)
	{
		REPORT
			Real* elx = data; int l = storage; int li = -1;
		while (l--) { Real f = fabs(*elx++); if (r >= f) { r = f; li = l; } }
		i = (li >= 0) ? storage - li + skip : 0;
		return r;
	}

	// max value of r and elements of row/col
	Real MatrixRowCol::Maximum1(Real r, int& i)
	{
		REPORT
			Real* elx = data; int l = storage; int li = -1;
		while (l--) { Real f = *elx++; if (r <= f) { r = f; li = l; } }
		i = (li >= 0) ? storage - li + skip : 0;
		return r;
	}

	// min value of r and elements of row/col
	Real MatrixRowCol::Minimum1(Real r, int& i)
	{
		REPORT
			Real* elx = data; int l = storage; int li = -1;
		while (l--) { Real f = *elx++; if (r >= f) { r = f; li = l; } }
		i = (li >= 0) ? storage - li + skip : 0;
		return r;
	}

	Real MatrixRowCol::Sum()
	{
		REPORT
			Real sum = 0.0; Real* elx = data; int l = storage;
		while (l--) sum += *elx++;
		return sum;
	}

	void MatrixRowCol::SubRowCol(MatrixRowCol& mrc, int skip1, int l1) const
	{
		mrc.length = l1;  int d = skip - skip1;
		if (d<0) { mrc.skip = 0; mrc.data = data - d; }
		else  { mrc.skip = d; mrc.data = data; }
		d = skip + storage - skip1;
		d = ((l1 < d) ? l1 : d) - mrc.skip;  mrc.storage = (d < 0) ? 0 : d;
		mrc.cw = 0;
	}

#ifdef use_namespace
}
#endif


//$$ bandmat.cpp                     Band matrix definitions

// Copyright (C) 1991,2,3,4,9: R B Davies

#define WANT_MATH                    // include.h will get math fns

//#define WANT_STREAM


#ifdef use_namespace
namespace NEWMAT {
#endif



#ifdef DO_REPORT
#define REPORT { static ExeCounter ExeCount(__LINE__,10); ++ExeCount; }
#else
#define REPORT {}
#endif

	static inline int my_min(int x, int y) { return x < y ? x : y; }
	static inline int my_max(int x, int y) { return x > y ? x : y; }


	BandMatrix::BandMatrix(const BaseMatrix& M)
	{
		REPORT // CheckConversion(M);
			// MatrixConversionCheck mcc;
			GeneralMatrix* gmx=((BaseMatrix&)M).Evaluate(MatrixType::BM);
		GetMatrix(gmx); CornerClear();
	}

	void BandMatrix::SetParameters(const GeneralMatrix* gmx)
	{
		REPORT
			MatrixBandWidth bw = gmx->bandwidth();
		lower_val = bw.lower_val; upper_val = bw.upper_val;
	}

	void BandMatrix::resize(int n, int lb, int ub)
	{
		REPORT
			Tracer tr("BandMatrix::resize");
		if (lb<0 || ub<0) Throw(ProgramException("Undefined bandwidth"));
		lower_val = (lb<=n) ? lb : n-1; upper_val = (ub<=n) ? ub : n-1;
		GeneralMatrix::resize(n,n,n*(lower_val+1+upper_val)); CornerClear();
	}

	// SimpleAddOK shows when we can add etc two matrices by a simple vector add
	// and when we can add one matrix into another
	// *gm must be the same type as *this
	// return 0 if simple add is OK
	// return 1 if we can add into *gm only
	// return 2 if we can add into *this only
	// return 3 if we can't add either way
	// For SP this will still be valid if we swap 1 and 2

	short BandMatrix::SimpleAddOK(const GeneralMatrix* gm)
	{
		const BandMatrix* bm = (const BandMatrix*)gm;
		if (bm->lower_val == lower_val && bm->upper_val == upper_val)
		{ REPORT return 0; }
		else if (bm->lower_val >= lower_val && bm->upper_val >= upper_val)
		{ REPORT return 1; }
		else if (bm->lower_val <= lower_val && bm->upper_val <= upper_val)
		{ REPORT return 2; }
		else { REPORT return 3; }
	}

	short SymmetricBandMatrix::SimpleAddOK(const GeneralMatrix* gm)
	{
		const SymmetricBandMatrix* bm = (const SymmetricBandMatrix*)gm;
		if (bm->lower_val == lower_val) { REPORT return 0; }
		else if (bm->lower_val > lower_val) { REPORT return 1; }
		else { REPORT return 2; }
	}

	void UpperBandMatrix::resize(int n, int lb, int ub)
	{
		REPORT
			if (lb != 0)
			{
				Tracer tr("UpperBandMatrix::resize");
				Throw(ProgramException("UpperBandMatrix with non-zero lower band" ));
			}
			BandMatrix::resize(n, lb, ub);
	}

	void LowerBandMatrix::resize(int n, int lb, int ub)
	{
		REPORT
			if (ub != 0)
			{
				Tracer tr("LowerBandMatrix::resize");
				Throw(ProgramException("LowerBandMatrix with non-zero upper band" ));
			}
			BandMatrix::resize(n, lb, ub);
	}

	void BandMatrix::resize(const GeneralMatrix& A)
	{
		REPORT
			int n = A.Nrows();
		if (n != A.Ncols())
		{
			Tracer tr("BandMatrix::resize(GM)");
			Throw(NotSquareException(*this));
		}
		MatrixBandWidth mbw = A.bandwidth();
		resize(n, mbw.Lower(), mbw.Upper());
	}
	/*
	bool BandMatrix::SameStorageType(const GeneralMatrix& A) const
	{
	if (type() != A.type()) { REPORT return false; }
	REPORT
	return bandwidth() == A.bandwidth();
	}

	void BandMatrix::resizeForAdd(const GeneralMatrix& A, const GeneralMatrix& B)
	{
	REPORT
	Tracer tr("BandMatrix::resizeForAdd");
	MatrixBandWidth A_BW = A.bandwidth(); MatrixBandWidth B_BW = B.bandwidth();
	if ((A_BW.Lower() < 0) | (A_BW.Upper() < 0) | (B_BW.Lower() < 0)
	| (A_BW.Upper() < 0))
	Throw(ProgramException("Can't resize to BandMatrix" ));
	// already know A and B are square
	resize(A.Nrows(), my_max(A_BW.Lower(), B_BW.Lower()),
	my_max(A_BW.Upper(), B_BW.Upper()));
	}

	void BandMatrix::resizeForSP(const GeneralMatrix& A, const GeneralMatrix& B)
	{
	REPORT
	Tracer tr("BandMatrix::resizeForSP");
	MatrixBandWidth A_BW = A.bandwidth(); MatrixBandWidth B_BW = B.bandwidth();
	if ((A_BW.Lower() < 0) | (A_BW.Upper() < 0) | (B_BW.Lower() < 0)
	| (A_BW.Upper() < 0))
	Throw(ProgramException("Can't resize to BandMatrix" ));
	// already know A and B are square
	resize(A.Nrows(), my_min(A_BW.Lower(), B_BW.Lower()),
	my_min(A_BW.Upper(), B_BW.Upper()));
	}
	*/

	void BandMatrix::operator=(const BaseMatrix& X)
	{
		REPORT // CheckConversion(X);
			// MatrixConversionCheck mcc;
			Eq(X,MatrixType::BM); CornerClear();
	}

	void BandMatrix::CornerClear() const
	{
		// set unused parts of BandMatrix to zero
		REPORT
			int i = lower_val; Real* s = store; int bw = lower_val + 1 + upper_val;
		while (i)
		{ int j = i--; Real* sj = s; s += bw; while (j--) *sj++ = 0.0; }
		i = upper_val; s = store + storage;
		while (i)
		{ int j = i--; Real* sj = s; s -= bw; while (j--) *(--sj) = 0.0; }
	}

	MatrixBandWidth MatrixBandWidth::operator+(const MatrixBandWidth& bw) const
	{
		REPORT
			int l = bw.lower_val; int u = bw.upper_val;
		l = (lower_val < 0 || l < 0) ? -1 : (lower_val > l) ? lower_val : l;
		u = (upper_val < 0 || u < 0) ? -1 : (upper_val > u) ? upper_val : u;
		return MatrixBandWidth(l,u);
	}

	MatrixBandWidth MatrixBandWidth::operator*(const MatrixBandWidth& bw) const
	{
		REPORT
			int l = bw.lower_val; int u = bw.upper_val;
		l = (lower_val < 0 || l < 0) ? -1 : lower_val+l;
		u = (upper_val < 0 || u < 0) ? -1 : upper_val+u;
		return MatrixBandWidth(l,u);
	}

	MatrixBandWidth MatrixBandWidth::minimum(const MatrixBandWidth& bw) const
	{
		REPORT
			int l = bw.lower_val; int u = bw.upper_val;
		if ((lower_val >= 0) && ( (l < 0) || (l > lower_val) )) l = lower_val;
		if ((upper_val >= 0) && ( (u < 0) || (u > upper_val) )) u = upper_val;
		return MatrixBandWidth(l,u);
	}

	UpperBandMatrix::UpperBandMatrix(const BaseMatrix& M)
	{
		REPORT // CheckConversion(M);
			// MatrixConversionCheck mcc;
			GeneralMatrix* gmx=((BaseMatrix&)M).Evaluate(MatrixType::UB);
		GetMatrix(gmx); CornerClear();
	}

	void UpperBandMatrix::operator=(const BaseMatrix& X)
	{
		REPORT // CheckConversion(X);
			// MatrixConversionCheck mcc;
			Eq(X,MatrixType::UB); CornerClear();
	}

	LowerBandMatrix::LowerBandMatrix(const BaseMatrix& M)
	{
		REPORT // CheckConversion(M);
			// MatrixConversionCheck mcc;
			GeneralMatrix* gmx=((BaseMatrix&)M).Evaluate(MatrixType::LB);
		GetMatrix(gmx); CornerClear();
	}

	void LowerBandMatrix::operator=(const BaseMatrix& X)
	{
		REPORT // CheckConversion(X);
			// MatrixConversionCheck mcc;
			Eq(X,MatrixType::LB); CornerClear();
	}

	BandLUMatrix::BandLUMatrix(const BaseMatrix& m)
	{
		REPORT
			Tracer tr("BandLUMatrix");
		storage2 = 0; store2 = 0; indx = 0; // in event of exception during build
		GeneralMatrix* gm = ((BaseMatrix&)m).Evaluate();
		if (gm->nrows() != gm->ncols())
		{ gm->tDelete(); Throw(NotSquareException(*this)); }
		if (gm->type() == MatrixType::BC)
		{ REPORT  ((BandLUMatrix*)gm)->get_aux(*this); GetMatrix(gm); }
		else
		{
			REPORT
				BandMatrix* gm1 = (BandMatrix*)(gm->Evaluate(MatrixType::BM));
			m1 = gm1->lower_val; m2 = gm1->upper_val;
			GetMatrix(gm1);
			d = true; sing = false;
			indx = new int [nrows_val]; MatrixErrorNoSpace(indx);
			MONITOR_INT_NEW("Index (BndLUMat)",nrows_val,indx)
				storage2 = nrows_val * m1;
			store2 = new Real [storage2]; MatrixErrorNoSpace(store2);
			MONITOR_REAL_NEW("Make (BandLUMat)",storage2,store2)
				ludcmp();
		}
	}

	GeneralMatrix* BandLUMatrix::Evaluate(MatrixType mt)
	{
		if (Compare(this->Type(),mt)) { REPORT return this; }
		REPORT
			Tracer et("BandLUMatrix::Evaluate");
		bool dummy = true;
		if (dummy) Throw(ProgramException("Illegal use of BandLUMatrix", *this));
		return this;
	}

	// could we use SetParameters instead of this
	void BandLUMatrix::get_aux(BandLUMatrix& X)
	{
		X.d = d; X.sing = sing; X.storage2 = storage2; X.m1 = m1; X.m2 = m2;   
		if (tag_val == 0 || tag_val == 1) // reuse the array 
		{
			REPORT
				X.indx = indx; indx = 0;
			X.store2 = store2; store2 = 0;
			d = true; sing = true; storage2 = 0; m1 = 0; m2 = 0;
			return;
		}
		else if (nrows_val == 0)
		{
			REPORT
				indx = 0; store2 = 0; storage2 = 0;
			d = true; sing = true; m1 = m2 = 0;
			return;
		}
		else                              // copy the array
		{
			REPORT
				Tracer tr("BandLUMatrix::get_aux");
			int *ix = new int [nrows_val]; MatrixErrorNoSpace(ix);
			MONITOR_INT_NEW("Index (BLUM::get_aux)", nrows_val, ix)
				int n = nrows_val; int* i = ix; int* j = indx;
			while(n--) *i++ = *j++;
			X.indx = ix;
			Real *rx = new Real [storage2]; MatrixErrorNoSpace(indx);
			MONITOR_REAL_NEW("Index (BLUM::get_aux)", storage2, rx)
				newmat_block_copy(storage2, store2, rx);
			X.store2 = rx;
		}
	}

	BandLUMatrix::BandLUMatrix(const BandLUMatrix& gm)
	{
		REPORT
			Tracer tr("BandLUMatrix(const BandLUMatrix&)");
		((BandLUMatrix&)gm).get_aux(*this);
		GetMatrix(&gm);
	}

	void BandLUMatrix::operator=(const BandLUMatrix& gm)
	{
		if (&gm == this) { REPORT tag_val = -1; return; }
		REPORT
			delete [] indx; indx = 0;
		delete [] store2; store2 = 0; storage2 = 0;
		((BandLUMatrix&)gm).get_aux(*this);
		Eq(gm);
	}








	BandLUMatrix::~BandLUMatrix()
	{
		REPORT
			MONITOR_INT_DELETE("Index (BndLUMat)",nrows_val,indx)
			MONITOR_REAL_DELETE("Delete (BndLUMt)",storage2,store2)
			delete [] indx; delete [] store2;
	}

	MatrixType BandLUMatrix::type() const { REPORT return MatrixType::BC; }


	LogAndSign BandLUMatrix::log_determinant() const
	{
		REPORT
			if (sing) return 0.0;
		Real* a = store; int w = m1+1+m2; LogAndSign sum; int i = nrows_val;
		// while (i--) { sum *= *a; a += w; }
		if (i) for (;;) { sum *= *a; if (!(--i)) break; a += w; }
		if (!d) sum.ChangeSign(); return sum;
	}

	GeneralMatrix* BandMatrix::MakeSolver()
	{
		REPORT
			GeneralMatrix* gm = new BandLUMatrix(*this);
		MatrixErrorNoSpace(gm); gm->ReleaseAndDelete(); return gm;
	}


	void BandLUMatrix::ludcmp()
	{
		REPORT
			Real* a = store2; int i = storage2;
		// clear store2 - so unused locations are always zero -
		// required by operator==
		while (i--) *a++ = 0.0;
		a = store;
		i = m1; int j = m2; int k; int n = nrows_val; int w = m1 + 1 + m2;
		while (i)
		{
			Real* ai = a + i;
			k = ++j; while (k--) *a++ = *ai++;
			k = i--; while (k--) *a++ = 0.0;
		}

		a = store; int l = m1;
		for (k=0; k<n; k++)
		{
			Real x = *a; i = k; Real* aj = a;
			if (l < n) l++;
			for (j=k+1; j<l; j++)
			{ aj += w; if (fabs(x) < fabs(*aj)) { x = *aj; i = j; } }
			indx[k] = i;
			if (x==0) { sing = true; return; }
			if (i!=k)
			{
				d = !d; Real* ak = a; Real* ai = store + i * w; j = w;
				while (j--) { x = *ak; *ak++ = *ai; *ai++ = x; }
			}
			aj = a + w; Real* m = store2 + m1 * k;
			for (j=k+1; j<l; j++)
			{
				*m++ = x = *aj / *a; i = w; Real* ak = a;
				while (--i) { Real* aj1 = aj++; *aj1 = *aj - x * *(++ak); }
				*aj++ = 0.0;
			}
			a += w;
		}
	}

	void BandLUMatrix::lubksb(Real* B, int mini)
	{
		REPORT
			Tracer tr("BandLUMatrix::lubksb");
		if (sing) Throw(SingularException(*this));
		int n = nrows_val; int l = m1; int w = m1 + 1 + m2;

		for (int k=0; k<n; k++)
		{
			int i = indx[k];
			if (i!=k) { Real x=B[k]; B[k]=B[i]; B[i]=x; }
			if (l<n) l++;
			Real* m = store2 + k*m1; Real* b = B+k; Real* bi = b;
			for (i=k+1; i<l; i++)  *(++bi) -= *m++ * *b;
		}

		l = -m1;
		for (int i = n-1; i>=mini; i--)
		{
			Real* b = B + i; Real* bk = b; Real x = *bk;
			Real* a = store + w*i; Real y = *a;
			int k = l+m1; while (k--) x -=  *(++a) * *(++bk);
			*b = x / y;
			if (l < m2) l++;
		}
	}

	void BandLUMatrix::Solver(MatrixColX& mcout, const MatrixColX& mcin)
	{
		REPORT
			int i = mcin.skip; Real* el = mcin.data-i; Real* el1=el;
		while (i--) *el++ = 0.0;
		el += mcin.storage; i = nrows_val - mcin.skip - mcin.storage;
		while (i--) *el++ = 0.0;
		lubksb(el1, mcout.skip);
	}

	// Do we need check for entirely zero output?


	void UpperBandMatrix::Solver(MatrixColX& mcout,
		const MatrixColX& mcin)
	{
		REPORT
			int i = mcin.skip-mcout.skip; Real* elx = mcin.data-i;
		while (i-- > 0) *elx++ = 0.0;
		int nr = mcin.skip+mcin.storage;
		elx = mcin.data+mcin.storage; Real* el = elx;
		int j = mcout.skip+mcout.storage-nr; i = nr-mcout.skip;
		while (j-- > 0) *elx++ = 0.0;

		Real* Ael = store + (upper_val+1)*(i-1)+1; j = 0;
		if (i > 0) for(;;)
		{
			elx = el; Real sum = 0.0; int jx = j;
			while (jx--) sum += *(--Ael) * *(--elx);
			elx--; *elx = (*elx - sum) / *(--Ael);
			if (--i <= 0) break;
			if (j<upper_val) Ael -= upper_val - (++j); else el--;
		}
	}

	void LowerBandMatrix::Solver(MatrixColX& mcout,
		const MatrixColX& mcin)
	{
		REPORT
			int i = mcin.skip-mcout.skip; Real* elx = mcin.data-i;
		while (i-- > 0) *elx++ = 0.0;
		int nc = mcin.skip; i = nc+mcin.storage; elx = mcin.data+mcin.storage;
		int nr = mcout.skip+mcout.storage; int j = nr-i; i = nr-nc;
		while (j-- > 0) *elx++ = 0.0;

		Real* el = mcin.data;
		Real* Ael = store + (lower_val+1)*nc + lower_val;
		j = 0;
		if (i > 0) for(;;)
		{
			elx = el; Real sum = 0.0; int jx = j;
			while (jx--) sum += *Ael++ * *elx++;
			*elx = (*elx - sum) / *Ael++;
			if (--i <= 0) break;
			if (j<lower_val) Ael += lower_val - (++j); else el++;
		}
	}


	LogAndSign BandMatrix::log_determinant() const
	{
		REPORT
			BandLUMatrix C(*this); return C.log_determinant();
	}

	LogAndSign LowerBandMatrix::log_determinant() const
	{
		REPORT
			int i = nrows_val; LogAndSign sum;
		Real* s = store + lower_val; int j = lower_val + 1;
		//   while (i--) { sum *= *s; s += j; }
		if (i) for (;;) { sum *= *s; if (!(--i)) break; s += j; }
		((GeneralMatrix&)*this).tDelete(); return sum;
	}

	LogAndSign UpperBandMatrix::log_determinant() const
	{
		REPORT
			int i = nrows_val; LogAndSign sum; Real* s = store; int j = upper_val + 1;
		//   while (i--) { sum *= *s; s += j; }
		if (i) for (;;) { sum *= *s; if (!(--i)) break; s += j; }
		((GeneralMatrix&)*this).tDelete(); return sum;
	}

	GeneralMatrix* SymmetricBandMatrix::MakeSolver()
	{
		REPORT
			GeneralMatrix* gm = new BandLUMatrix(*this);
		MatrixErrorNoSpace(gm); gm->ReleaseAndDelete(); return gm;
	}

	SymmetricBandMatrix::SymmetricBandMatrix(const BaseMatrix& M)
	{
		REPORT  // CheckConversion(M);
			// MatrixConversionCheck mcc;
			GeneralMatrix* gmx=((BaseMatrix&)M).Evaluate(MatrixType::SB);
		GetMatrix(gmx);
	}

	GeneralMatrix* SymmetricBandMatrix::Transpose(TransposedMatrix*, MatrixType mt)
	{ REPORT  return Evaluate(mt); }

	LogAndSign SymmetricBandMatrix::log_determinant() const
	{
		REPORT
			BandLUMatrix C(*this); return C.log_determinant();
	}

	void SymmetricBandMatrix::SetParameters(const GeneralMatrix* gmx)
	{ REPORT lower_val = gmx->bandwidth().lower_val; }

	void SymmetricBandMatrix::resize(int n, int lb)
	{
		REPORT
			Tracer tr("SymmetricBandMatrix::resize");
		if (lb<0) Throw(ProgramException("Undefined bandwidth"));
		lower_val = (lb<=n) ? lb : n-1;
		GeneralMatrix::resize(n,n,n*(lower_val+1));
	}

	void SymmetricBandMatrix::resize(const GeneralMatrix& A)
	{
		REPORT
			int n = A.Nrows();
		if (n != A.Ncols())
		{
			Tracer tr("SymmetricBandMatrix::resize(GM)");
			Throw(NotSquareException(*this));
		}
		MatrixBandWidth mbw = A.bandwidth(); int b = mbw.Lower();
		if (b != mbw.Upper())
		{
			Tracer tr("SymmetricBandMatrix::resize(GM)");
			Throw(ProgramException("Upper and lower band-widths not equal"));
		}
		resize(n, b);
	}
	/*
	bool SymmetricBandMatrix::SameStorageType(const GeneralMatrix& A) const
	{
	if (type() != A.type()) { REPORT return false; }
	REPORT
	return bandwidth() == A.bandwidth();
	}

	void SymmetricBandMatrix::resizeForAdd(const GeneralMatrix& A,
	const GeneralMatrix& B)
	{
	REPORT
	Tracer tr("SymmetricBandMatrix::resizeForAdd");
	MatrixBandWidth A_BW = A.bandwidth(); MatrixBandWidth B_BW = B.bandwidth();
	if ((A_BW.Lower() < 0) | (B_BW.Lower() < 0))
	Throw(ProgramException("Can't resize to SymmetricBandMatrix" ));
	// already know A and B are square
	resize(A.Nrows(), my_max(A_BW.Lower(), B_BW.Lower()));
	}

	void SymmetricBandMatrix::resizeForSP(const GeneralMatrix& A,
	const GeneralMatrix& B)
	{
	REPORT
	Tracer tr("SymmetricBandMatrix::resizeForSP");
	MatrixBandWidth A_BW = A.bandwidth(); MatrixBandWidth B_BW = B.bandwidth();
	if ((A_BW.Lower() < 0) | (B_BW.Lower() < 0))
	Throw(ProgramException("Can't resize to SymmetricBandMatrix" ));
	// already know A and B are square
	resize(A.Nrows(), my_min(A_BW.Lower(), B_BW.Lower()));
	}
	*/

	void SymmetricBandMatrix::operator=(const BaseMatrix& X)
	{
		REPORT // CheckConversion(X);
			// MatrixConversionCheck mcc;
			Eq(X,MatrixType::SB);
	}

	void SymmetricBandMatrix::CornerClear() const
	{
		// set unused parts of BandMatrix to zero
		REPORT
			int i = lower_val; Real* s = store; int bw = lower_val + 1;
		if (i) for(;;)
		{
			int j = i;
			Real* sj = s;
			while (j--) *sj++ = 0.0;
			if (!(--i)) break;
			s += bw;
		}
	}

	MatrixBandWidth SymmetricBandMatrix::bandwidth() const
	{ REPORT return MatrixBandWidth(lower_val,lower_val); }

	GeneralMatrix* BandMatrix::Image() const
	{
		REPORT
			GeneralMatrix* gm = new BandMatrix(*this); MatrixErrorNoSpace(gm);
		return gm;
	}

	GeneralMatrix* UpperBandMatrix::Image() const
	{
		REPORT
			GeneralMatrix* gm = new UpperBandMatrix(*this); MatrixErrorNoSpace(gm);
		return gm;
	}

	GeneralMatrix* LowerBandMatrix::Image() const
	{
		REPORT
			GeneralMatrix* gm = new LowerBandMatrix(*this); MatrixErrorNoSpace(gm);
		return gm;
	}

	GeneralMatrix* SymmetricBandMatrix::Image() const
	{
		REPORT
			GeneralMatrix* gm = new SymmetricBandMatrix(*this); MatrixErrorNoSpace(gm);
		return gm;
	}

	GeneralMatrix* BandLUMatrix::Image() const
	{
		REPORT
			GeneralMatrix* gm = new BandLUMatrix(*this); MatrixErrorNoSpace(gm);
		return gm;
	}


//	inline Real square(Real x) { return x*x; }

	Real SymmetricBandMatrix::sum_square() const
	{
		REPORT
			CornerClear();
		Real sum1=0.0; Real sum2=0.0; Real* s=store; int i=nrows_val;
		int l=lower_val;
		while (i--)
		{ int j = l; while (j--) sum2 += square(*s++); sum1 += square(*s++); }
		((GeneralMatrix&)*this).tDelete(); return sum1 + 2.0 * sum2;
	}

	Real SymmetricBandMatrix::sum_absolute_value() const
	{
		REPORT
			CornerClear();
		Real sum1=0.0; Real sum2=0.0; Real* s=store; int i=nrows_val;
		int l=lower_val;
		while (i--)
		{ int j = l; while (j--) sum2 += fabs(*s++); sum1 += fabs(*s++); }
		((GeneralMatrix&)*this).tDelete(); return sum1 + 2.0 * sum2;
	}

	Real SymmetricBandMatrix::sum() const
	{
		REPORT
			CornerClear();
		Real sum1=0.0; Real sum2=0.0; Real* s=store; int i=nrows_val;
		int l=lower_val;
		while (i--)
		{ int j = l; while (j--) sum2 += *s++; sum1 += *s++; }
		((GeneralMatrix&)*this).tDelete(); return sum1 + 2.0 * sum2;
	}





#ifdef use_namespace
}
#endif

//$$ hholder.cpp                                   QR decomposition

// Copyright (C) 1991,2,3,4: R B Davies

#define WANT_MATH
//#define WANT_STREAM


#ifdef use_namespace
namespace NEWMAT {
#endif

#ifdef DO_REPORT
#define REPORT { static ExeCounter ExeCount(__LINE__,16); ++ExeCount; }
#else
#define REPORT {}
#endif


	/*************************** QR decompositions ***************************/

//	inline Real square(Real x) { return x*x; }

	void QRZT(Matrix& X, LowerTriangularMatrix& L)
	{
		REPORT
			Tracer et("QRZT(1)");
		int n = X.Ncols(); int s = X.Nrows(); L.resize(s);
		if (n == 0 || s == 0) { L = 0.0; return; }
		Real* xi = X.Store(); int k;
		for (int i=0; i<s; i++)
		{
			Real sum = 0.0;
			Real* xi0=xi; k=n; while(k--) { sum += square(*xi++); }
			sum = sqrt(sum);
			if (sum == 0.0)
			{
				REPORT
					k=n; while(k--) { *xi0++ = 0.0; }
				for (int j=i; j<s; j++) L.element(j,i) = 0.0;
			}
			else
			{
				L.element(i,i) = sum;
				Real* xj0=xi0; k=n; while(k--) { *xj0++ /= sum; }
				for (int j=i+1; j<s; j++)
				{
					sum=0.0;
					xi=xi0; Real* xj=xj0; k=n; while(k--) { sum += *xi++ * *xj++; }
					xi=xi0; k=n; while(k--) { *xj0++ -= sum * *xi++; }
					L.element(j,i) = sum;
				}
			}
		}
	}

	void QRZT(const Matrix& X, Matrix& Y, Matrix& M)
	{
		REPORT
			Tracer et("QRZT(2)");
		int n = X.Ncols(); int s = X.Nrows(); int t = Y.Nrows();
		if (Y.Ncols() != n)
		{ Throw(ProgramException("Unequal row lengths",X,Y)); }
		M.resize(t,s);
		Real* xi = X.Store(); int k;
		for (int i=0; i<s; i++)
		{
			Real* xj0 = Y.Store(); Real* xi0 = xi;
			for (int j=0; j<t; j++)
			{
				Real sum=0.0;
				xi=xi0; Real* xj=xj0; k=n; while(k--) { sum += *xi++ * *xj++; }
				xi=xi0; k=n; while(k--) { *xj0++ -= sum * *xi++; }
				M.element(j,i) = sum;
			}
		}
	}

	/*
	void QRZ(Matrix& X, UpperTriangularMatrix& U)
	{
	Tracer et("QRZ(1)");
	int n = X.Nrows(); int s = X.Ncols(); U.resize(s);
	Real* xi0 = X.Store(); int k;
	for (int i=0; i<s; i++)
	{
	Real sum = 0.0;
	Real* xi = xi0; k=n; while(k--) { sum += square(*xi); xi+=s; }
	sum = sqrt(sum);
	U.element(i,i) = sum;
	if (sum==0.0) Throw(SingularException(U));
	Real* xj0=xi0; k=n; while(k--) { *xj0 /= sum; xj0+=s; }
	xj0 = xi0;
	for (int j=i+1; j<s; j++)
	{
	sum=0.0;
	xi=xi0; k=n; xj0++; Real* xj=xj0;
	while(k--) { sum += *xi * *xj; xi+=s; xj+=s; }
	xi=xi0; k=n; xj=xj0;
	while(k--) { *xj -= sum * *xi; xj+=s; xi+=s; }
	U.element(i,j) = sum;
	}
	xi0++;
	}
	}
	*/

	void QRZ(Matrix& X, UpperTriangularMatrix& U)
	{
		REPORT
			Tracer et("QRZ(1)");
		int n = X.Nrows(); int s = X.Ncols(); U.resize(s); U = 0.0;
		if (n == 0 || s == 0) return;
		Real* xi0 = X.Store(); Real* u0 = U.Store(); Real* u;
		int j, k; int J = s; int i = s;
		while (i--)
		{
			Real* xj0 = xi0; Real* xi = xi0; k = n;
			if (k) for (;;)
			{
				u = u0; Real Xi = *xi; Real* xj = xj0;
				j = J; while(j--) *u++ += Xi * *xj++;
				if (!(--k)) break;
				xi += s; xj0 += s;
			}

			Real sum = sqrt(*u0); *u0 = sum; u = u0+1;
			if (sum == 0.0)
			{
				REPORT
					j = J - 1; while(j--) *u++ = 0.0;

				xj0 = xi0++; k = n;
				if (k) for (;;)
				{
					*xj0 = 0.0;
					if (!(--k)) break;
					xj0 += s;
				}
				u0 += J--;
			}
			else
			{
				int J1 = J-1; j = J1; while(j--) *u++ /= sum;

				xj0 = xi0; xi = xi0++; k = n;
				if (k) for (;;)
				{
					u = u0+1; Real Xi = *xi; Real* xj = xj0;
					Xi /= sum; *xj++ = Xi;
					j = J1; while(j--) *xj++ -= *u++ * Xi;
					if (!(--k)) break;
					xi += s; xj0 += s;
				}
				u0 += J--;
			}
		}
	}

	void QRZ(const Matrix& X, Matrix& Y, Matrix& M)
	{
		REPORT
			Tracer et("QRZ(2)");
		int n = X.Nrows(); int s = X.Ncols(); int t = Y.Ncols();
		if (Y.Nrows() != n)
		{ Throw(ProgramException("Unequal column lengths",X,Y)); }
		M.resize(s,t); M = 0;Real* m0 = M.Store(); Real* m;
		Real* xi0 = X.Store();
		int j, k; int i = s;
		while (i--)
		{
			Real* xj0 = Y.Store(); Real* xi = xi0; k = n;
			if (k) for (;;)
			{
				m = m0; Real Xi = *xi; Real* xj = xj0;
				j = t; while(j--) *m++ += Xi * *xj++;
				if (!(--k)) break;
				xi += s; xj0 += t;
			}

			xj0 = Y.Store(); xi = xi0++; k = n;
			if (k) for (;;)
			{
				m = m0; Real Xi = *xi; Real* xj = xj0;
				j = t; while(j--) *xj++ -= *m++ * Xi;
				if (!(--k)) break;
				xi += s; xj0 += t;
			}
			m0 += t;
		}
	}

	/*

	void QRZ(const Matrix& X, Matrix& Y, Matrix& M)
	{
	Tracer et("QRZ(2)");
	int n = X.Nrows(); int s = X.Ncols(); int t = Y.Ncols();
	if (Y.Nrows() != n)
	{ Throw(ProgramException("Unequal column lengths",X,Y)); }
	M.resize(s,t);
	Real* xi0 = X.Store(); int k;
	for (int i=0; i<s; i++)
	{
	Real* xj0 = Y.Store();
	for (int j=0; j<t; j++)
	{
	Real sum=0.0;
	Real* xi=xi0; Real* xj=xj0; k=n;
	while(k--) { sum += *xi * *xj; xi+=s; xj+=t; }
	xi=xi0; k=n; xj=xj0++;
	while(k--) { *xj -= sum * *xi; xj+=t; xi+=s; }
	M.element(i,j) = sum;
	}
	xi0++;
	}
	}
	*/

	void updateQRZT(Matrix& X, LowerTriangularMatrix& L)
	{
		REPORT
			Tracer et("updateQRZT");
		int n = X.Ncols(); int s = X.Nrows();
		if (s != L.Nrows())
			Throw(ProgramException("Incompatible dimensions",X,L)); 
		if (n == 0 || s == 0) return;
		Real* xi = X.Store(); int k;
		for (int i=0; i<s; i++)
		{
			Real r = L.element(i,i); 
			Real sum = 0.0;
			Real* xi0=xi; k=n; while(k--) { sum += square(*xi++); }
			sum = sqrt(sum + square(r));
			if (sum == 0.0)
			{
				REPORT
					k=n; while(k--) { *xi0++ = 0.0; }
				for (int j=i; j<s; j++) L.element(j,i) = 0.0;
			}
			else
			{
				Real frs = fabs(r) + sum;
				Real a0 = sqrt(frs / sum); Real alpha = a0 / frs;
				if (r <= 0) { REPORT L.element(i,i) = sum; alpha = -alpha; }
				else { REPORT L.element(i,i) = -sum; }
				Real* xj0=xi0; k=n; while(k--) { *xj0++ *= alpha; }
				for (int j=i+1; j<s; j++)
				{
					sum = 0.0;
					xi=xi0; Real* xj=xj0; k=n; while(k--) { sum += *xi++ * *xj++; }
					sum += a0 * L.element(j,i);
					xi=xi0; k=n; while(k--) { *xj0++ -= sum * *xi++; }
					L.element(j,i) -= sum * a0;
				}
			}
		}
	}

	void updateQRZ(Matrix& X, UpperTriangularMatrix& U)
	{
		REPORT
			Tracer et("updateQRZ");
		int n = X.Nrows(); int s = X.Ncols();
		if (s != U.Ncols())
			Throw(ProgramException("Incompatible dimensions",X,U));
		if (n == 0 || s == 0) return; 
		Real* xi0 = X.Store(); Real* u0 = U.Store(); Real* u;
		RowVector V(s); Real* v0 = V.Store(); Real* v; V = 0.0;
		int j, k; int J = s; int i = s;
		while (i--)
		{
			Real* xj0 = xi0; Real* xi = xi0; k = n;
			if (k) for (;;)
			{
				v = v0; Real Xi = *xi; Real* xj = xj0;
				j = J; while(j--) *v++ += Xi * *xj++;
				if (!(--k)) break;
				xi += s; xj0 += s;
			}

			Real r = *u0;
			Real sum = sqrt(*v0 + square(r));

			if (sum == 0.0)
			{
				REPORT
					u = u0; v = v0;
				j = J; while(j--) { *u++ = 0.0; *v++ = 0.0; }
				xj0 = xi0++; k = n;
				if (k) for (;;)
				{
					*xj0 = 0.0;
					if (!(--k)) break;
					xj0 += s;
				}
				u0 += J--;
			}
			else
			{
				Real frs = fabs(r) + sum;
				Real a0 = sqrt(frs / sum); Real alpha = a0 / frs;
				if (r <= 0) { REPORT alpha = -alpha; *u0 = sum; }
				else { REPORT *u0 = -sum; }

				j = J - 1; v = v0 + 1; u = u0 + 1;     
				while (j--)
				{ *v = a0 * *u + alpha * *v; *u -= a0 * *v; ++v; ++u; }

				xj0 = xi0; xi = xi0++; k = n;
				if (k) for (;;)
				{
					v = v0 + 1; Real Xi = *xi; Real* xj = xj0;
					Xi *= alpha; *xj++ = Xi;
					j = J - 1; while(j--) *xj++ -= *v++ * Xi;
					if (!(--k)) break;
					xi += s; xj0 += s;
				}

				j = J; v = v0;
				while (j--) *v++ = 0.0;

				u0 += J--;
			}
		}
	}

	// Matrix A's first n columns are orthonormal
	// so A.Columns(1,n).t() * A.Columns(1,n) is the identity matrix.
	// Fill out the remaining columns of A to make them orthonormal
	// so A.t() * A is the identity matrix 
	void extend_orthonormal(Matrix& A, int n)
	{
		REPORT
			Tracer et("extend_orthonormal");
		int nr = A.nrows(); int nc = A.ncols();
		if (nc > nr) Throw(IncompatibleDimensionsException(A));
		if (n > nc) Throw(IncompatibleDimensionsException(A));
		ColumnVector SSR;
		{ Matrix A1 = A.Columns(1,n); SSR = A1.sum_square_rows(); }
		for (int i = n; i < nc; ++i)
		{
			// pick row with smallest SSQ
			int k; SSR.minimum1(k);
			// orthogonalise column with 1 at element k, 0 elsewhere
			// next line is rather inefficient
			ColumnVector X = - A.Columns(1, i) * A.SubMatrix(k, k, 1, i).t();
			X(k) += 1.0;
			// normalise
			X /= sqrt(X.SumSquare());
			// update row sums of squares
			for (k = 1; k <= nr; ++k) SSR(k) += square(X(k));
			// load new column into matrix
			A.Column(i+1) = X;
		}
	}





#ifdef use_namespace
}
#endif

//$$myexcept.cpp                        Exception handler

// Copyright (C) 1993,4,6: R B Davies


#define WANT_STREAM                    // include.h will get stream fns
#define WANT_STRING



#ifdef use_namespace
namespace RBD_COMMON {
#endif


	//#define REG_DEREG                    // for print out uses of new/delete
	//#define CLEAN_LIST                   // to print entries being added to
	// or deleted from cleanup list

#ifdef SimulateExceptions

	void Throw()
	{
		for (Janitor* jan = JumpBase::jl->janitor; jan; jan = jan->NextJanitor)
			jan->CleanUp();
		JumpItem* jx = JumpBase::jl->ji;    // previous jumpbase;
		if ( !jx ) { MatrixTerminate(); }         // jl was initial JumpItem
		JumpBase::jl = jx;                  // drop down a level; cannot be in front
		// of previous line
		Tracer::last = JumpBase::jl->trace;
		longjmp(JumpBase::jl->env, 1);
	}

#endif                                 // end of simulate exceptions


	unsigned long BaseException::Select;
	char* BaseException::what_error;
	int BaseException::SoFar;
	int BaseException::LastOne;

	BaseException::BaseException(const char* a_what)
	{
		Select++; SoFar = 0;
		if (!what_error)                   // make space for exception message
		{
			LastOne = 511;
			what_error = new char[512];
			if (!what_error)                // fail to make space
			{
				LastOne = 0;
				what_error = (char *)"No heap space for exception message\n";
			}
		}
		AddMessage("\n\nAn exception has been thrown\n");
		AddMessage(a_what);
		if (a_what) Tracer::AddTrace();
	}

	void BaseException::AddMessage(const char* a_what)
	{
		if (a_what)
		{
			int l = strlen(a_what); int r = LastOne - SoFar;
			if (l < r) { strcpy(what_error+SoFar, a_what); SoFar += l; }
			else if (r > 0)
			{
				strncpy(what_error+SoFar, a_what, r);
				what_error[LastOne] = 0;
				SoFar = LastOne;
			}
		}
	}

	void BaseException::AddInt(int value)
	{
		bool negative;
		if (value == 0) { AddMessage("0"); return; }
		else if (value < 0) { value = -value; negative = true; }
		else negative = false;
		int n = 0; int v = value;        // how many digits will we need?
		while (v > 0) { v /= 10; n++; }
		if (negative) n++;
		if (LastOne-SoFar < n) { AddMessage("***"); return; }

		SoFar += n; n = SoFar; what_error[n] = 0;
		while (value > 0)
		{
			int nv = value / 10; int rm = value - nv * 10;  value = nv;
			what_error[--n] = (char)(rm + '0');
		}
		if (negative) what_error[--n] = '-';
		return;
	}

	void Tracer::PrintTrace()
	{
		std::cout << "\n";
		for (Tracer* et = last; et; et=et->previous)
			std::cout << "  * " << et->entry << "\n";
	}

	void Tracer::AddTrace()
	{
		if (last)
		{
			BaseException::AddMessage("Trace: ");
			BaseException::AddMessage(last->entry);
			for (Tracer* et = last->previous; et; et=et->previous)
			{
				BaseException::AddMessage("; ");
				BaseException::AddMessage(et->entry);
			}
			BaseException::AddMessage(".\n");
		}
	}

#ifdef SimulateExceptions


	Janitor::Janitor()
	{
		if (do_not_link)
		{
			do_not_link = false; NextJanitor = 0; OnStack = false;
#ifdef CLEAN_LIST
			cout << "Not added to clean-list " << (unsigned long)this << "\n";
#endif
		}
		else
		{
			OnStack = true;
#ifdef CLEAN_LIST
			cout << "Add to       clean-list " << (unsigned long)this << "\n";
#endif
			NextJanitor = JumpBase::jl->janitor; JumpBase::jl->janitor=this;
		}
	}

	Janitor::~Janitor()
	{
		// expect the item to be deleted to be first on list
		// but must be prepared to search list
		if (OnStack)
		{
#ifdef CLEAN_LIST
			cout << "Delete from  clean-list " << (unsigned long)this << "\n";
#endif
			Janitor* lastjan = JumpBase::jl->janitor;
			if (this == lastjan) JumpBase::jl->janitor = NextJanitor;
			else
			{
				for (Janitor* jan = lastjan->NextJanitor; jan;
					jan = lastjan->NextJanitor)
				{
					if (jan==this)
					{ lastjan->NextJanitor = jan->NextJanitor; return; }
					lastjan=jan;
				}

				Throw(BaseException(
					"Cannot resolve memory linked list\nSee notes in myexcept.cpp for details\n"
					));


				// This message occurs when a call to ~Janitor() occurs, apparently
				// without a corresponding call to Janitor(). This could happen if my
				// way of deciding whether a constructor is being called by new
				// fails.

				// It may happen if you are using my simulated exceptions and also have
				// your compiler s exceptions turned on.

				// It can also happen if you have a class derived from Janitor
				// which does not include a copy constructor [ eg X(const &X) ].
				// Possibly also if delete is applied an object on the stack (ie not
				// called by new). Otherwise, it is a bug in myexcept or your compiler.
				// If you do not #define TEMPS_DESTROYED_QUICKLY you will get this
				// error with Microsoft C 7.0. There are probably situations where
				// you will get this when you do define TEMPS_DESTROYED_QUICKLY. This
				// is a bug in MSC. Beware of "operator" statements for defining
				// conversions; particularly for converting from a Base class to a
				// Derived class.

				// You may get away with simply deleting this error message and Throw
				// statement if you can not find a better way of overcoming the
				// problem. In any case please tell me if you get this error message,
				// particularly for compilers apart from Microsoft C 7.0.


			}
		}
	}

	JumpItem* JumpBase::jl;              // will be set to zero
	jmp_buf JumpBase::env;
	bool Janitor::do_not_link;           // will be set to false


	int JanitorInitializer::ref_count;

	JanitorInitializer::JanitorInitializer()
	{
		if (ref_count++ == 0) new JumpItem;
		// need JumpItem at head of list
	}

#endif                              // end of SimulateExceptions

	Tracer* Tracer::last;               // will be set to zero


	void MatrixTerminate()
	{
		cout << "\n\nThere has been an exception with no handler - exiting";
		const char* what = BaseException::what();
		if (what) cout << what << "\n";
		//exit(1);
	}



#ifdef DO_FREE_CHECK
	// Routines for tracing whether new and delete calls are balanced

	FreeCheckLink::FreeCheckLink() : next(FreeCheck::next)
	{ FreeCheck::next = this; }

	FCLClass::FCLClass(void* t, char* name) : ClassName(name) { ClassStore=t; }

	FCLRealArray::FCLRealArray(void* t, char* o, int s)
		: Operation(o), size(s) { ClassStore=t; }

	FCLIntArray::FCLIntArray(void* t, char* o, int s)
		: Operation(o), size(s) { ClassStore=t; }

	FreeCheckLink* FreeCheck::next;
	int FreeCheck::BadDelete;

	void FCLClass::Report()
	{ cout << "   " << ClassName << "   " << (unsigned long)ClassStore << "\n"; }

	void FCLRealArray::Report()
	{
		cout << "   " << Operation << "   " << (unsigned long)ClassStore <<
			"   " << size << "\n";
	}

	void FCLIntArray::Report()
	{
		cout << "   " << Operation << "   " << (unsigned long)ClassStore <<
			"   " << size << "\n";
	}

	void FreeCheck::Register(void* t, char* name)
	{
		FCLClass* f = new FCLClass(t,name);
		if (!f) { cout << "Out of memory in FreeCheck\n"; exit(1); }
#ifdef REG_DEREG
		cout << "Registering   " << name << "   " << (unsigned long)t << "\n";
#endif
	}

	void FreeCheck::RegisterR(void* t, char* o, int s)
	{
		FCLRealArray* f = new FCLRealArray(t,o,s);
		if (!f) { cout << "Out of memory in FreeCheck\n"; exit(1); }
#ifdef REG_DEREG
		cout << o << "   " << s << "   " << (unsigned long)t << "\n";
#endif
	}

	void FreeCheck::RegisterI(void* t, char* o, int s)
	{
		FCLIntArray* f = new FCLIntArray(t,o,s);
		if (!f) { cout << "Out of memory in FreeCheck\n"; exit(1); }
#ifdef REG_DEREG
		cout << o << "   " << s << "   " << (unsigned long)t << "\n";
#endif
	}

	void FreeCheck::DeRegister(void* t, char* name)
	{
		FreeCheckLink* last = 0;
#ifdef REG_DEREG
		cout << "Deregistering " << name << "   " << (unsigned long)t << "\n";
#endif
		for (FreeCheckLink* fcl = next; fcl; fcl = fcl->next)
		{
			if (fcl->ClassStore==t)
			{
				if (last) last->next = fcl->next; else next = fcl->next;
				delete fcl; return;
			}
			last = fcl;
		}
		cout << "\nRequest to delete non-existent object of class and location:\n";
		cout << "   " << name << "   " << (unsigned long)t << "\n";
		BadDelete++;
		Tracer::PrintTrace();
		cout << "\n";
	}

	void FreeCheck::DeRegisterR(void* t, char* o, int s)
	{
		FreeCheckLink* last = 0;
#ifdef REG_DEREG
		cout << o << "   " << s << "   " << (unsigned long)t << "\n";
#endif
		for (FreeCheckLink* fcl = next; fcl; fcl = fcl->next)
		{
			if (fcl->ClassStore==t)
			{
				if (last) last->next = fcl->next; else next = fcl->next;
				if (s >= 0 && ((FCLRealArray*)fcl)->size != s)
				{
					cout << "\nArray sizes do not agree:\n";
					cout << "   " << o << "   " << (unsigned long)t
						<< "   " << ((FCLRealArray*)fcl)->size << "   " << s << "\n";
					Tracer::PrintTrace();
					cout << "\n";
				}
				delete fcl; return;
			}
			last = fcl;
		}
		cout << "\nRequest to delete non-existent real array:\n";
		cout << "   " << o << "   " << (unsigned long)t << "   " << s << "\n";
		BadDelete++;
		Tracer::PrintTrace();
		cout << "\n";
	}

	void FreeCheck::DeRegisterI(void* t, char* o, int s)
	{
		FreeCheckLink* last = 0;
#ifdef REG_DEREG
		cout << o << "   " << s << "   " << (unsigned long)t << "\n";
#endif
		for (FreeCheckLink* fcl = next; fcl; fcl = fcl->next)
		{
			if (fcl->ClassStore==t)
			{
				if (last) last->next = fcl->next; else next = fcl->next;
				if (s >= 0 && ((FCLIntArray*)fcl)->size != s)
				{
					cout << "\nArray sizes do not agree:\n";
					cout << "   " << o << "   " << (unsigned long)t
						<< "   " << ((FCLIntArray*)fcl)->size << "   " << s << "\n";
					Tracer::PrintTrace();
					cout << "\n";
				}
				delete fcl; return;
			}
			last = fcl;
		}
		cout << "\nRequest to delete non-existent int array:\n";
		cout << "   " << o << "   " << (unsigned long)t << "   " << s << "\n";
		BadDelete++;
		Tracer::PrintTrace();
		cout << "\n";
	}

	void FreeCheck::Status()
	{
		if (next)
		{
			cout << "\nObjects of the following classes remain undeleted:\n";
			for (FreeCheckLink* fcl = next; fcl; fcl = fcl->next) fcl->Report();
			cout << "\n";
		}
		else cout << "\nNo objects remain undeleted\n\n";
		if (BadDelete)
		{
			cout << "\nThere were " << BadDelete << 
				" requests to delete non-existent items\n\n";
		}
	}

#endif                            // end of DO_FREE_CHECK

	// derived exception bodies

	Logic_error::Logic_error(const char* a_what) : BaseException()
	{
		Select = BaseException::Select;
		AddMessage("Logic error:- "); AddMessage(a_what);
		if (a_what) Tracer::AddTrace();
	}

	Runtime_error::Runtime_error(const char* a_what)
		: BaseException()
	{
		Select = BaseException::Select;
		AddMessage("Runtime error:- "); AddMessage(a_what);
		if (a_what) Tracer::AddTrace();
	}

	Domain_error::Domain_error(const char* a_what) : Logic_error()
	{
		Select = BaseException::Select;
		AddMessage("domain error\n"); AddMessage(a_what);
		if (a_what) Tracer::AddTrace();
	}

	Invalid_argument::Invalid_argument(const char* a_what) : Logic_error()
	{
		Select = BaseException::Select;
		AddMessage("invalid argument\n"); AddMessage(a_what);
		if (a_what) Tracer::AddTrace();
	}

	Length_error::Length_error(const char* a_what) : Logic_error()
	{
		Select = BaseException::Select;
		AddMessage("length error\n"); AddMessage(a_what);
		if (a_what) Tracer::AddTrace();
	}

	Out_of_range::Out_of_range(const char* a_what) : Logic_error()
	{
		Select = BaseException::Select;
		AddMessage("out of range\n"); AddMessage(a_what);
		if (a_what) Tracer::AddTrace();
	}

	//Bad_cast::Bad_cast(const char* a_what) : Logic_error()
	//{
	//   Select = BaseException::Select;
	//   AddMessage("bad cast\n"); AddMessage(a_what);
	//   if (a_what) Tracer::AddTrace();
	//}

	//Bad_typeid::Bad_typeid(const char* a_what) : Logic_error()
	//{
	//   Select = BaseException::Select;
	//   AddMessage("bad type id.\n"); AddMessage(a_what);
	//   if (a_what) Tracer::AddTrace();
	//}

	Range_error::Range_error(const char* a_what) : Runtime_error()
	{
		Select = BaseException::Select;
		AddMessage("range error\n"); AddMessage(a_what);
		if (a_what) Tracer::AddTrace();
	}

	Overflow_error::Overflow_error(const char* a_what) : Runtime_error()
	{
		Select = BaseException::Select;
		AddMessage("overflow error\n"); AddMessage(a_what);
		if (a_what) Tracer::AddTrace();
	}

	Bad_alloc::Bad_alloc(const char* a_what) : BaseException()
	{
		Select = BaseException::Select;
		AddMessage("bad allocation\n"); AddMessage(a_what);
		if (a_what) Tracer::AddTrace();
	}




	unsigned long Logic_error::Select;
	unsigned long Runtime_error::Select;
	unsigned long Domain_error::Select;
	unsigned long Invalid_argument::Select;
	unsigned long Length_error::Select;
	unsigned long Out_of_range::Select;
	//unsigned long Bad_cast::Select;
	//unsigned long Bad_typeid::Select;
	unsigned long Range_error::Select;
	unsigned long Overflow_error::Select;
	unsigned long Bad_alloc::Select;

#ifdef use_namespace
}
#endif

//$$ newmatex.cpp                    Exception handler

// Copyright (C) 1992,3,4,7: R B Davies

#define WANT_STREAM                  // include.h will get stream fns


#ifdef use_namespace
namespace NEWMAT {
#endif

	unsigned long OverflowException::Select;
	unsigned long SingularException::Select;
	unsigned long NPDException::Select;
	unsigned long ConvergenceException::Select;
	unsigned long ProgramException::Select;
	unsigned long IndexException::Select;
	unsigned long VectorException::Select;
	unsigned long NotSquareException::Select;
	unsigned long SubMatrixDimensionException::Select;
	unsigned long IncompatibleDimensionsException::Select;
	unsigned long NotDefinedException::Select;
	unsigned long CannotBuildException::Select;
	unsigned long InternalException::Select;



	static void MatrixDetails(const GeneralMatrix& A)
		// write matrix details to Exception buffer
	{
		MatrixBandWidth bw = A.bandwidth();
		int ubw = bw.upper_val; int lbw = bw.lower_val;
		BaseException::AddMessage("MatrixType = ");
		BaseException::AddMessage(A.Type().Value());
		BaseException::AddMessage("  # Rows = "); BaseException::AddInt(A.Nrows());
		BaseException::AddMessage("; # Cols = "); BaseException::AddInt(A.Ncols());
		if (lbw >=0)
		{
			BaseException::AddMessage("; lower BW = ");
			BaseException::AddInt(lbw);
		}
		if (ubw >=0)
		{
			BaseException::AddMessage("; upper BW = ");
			BaseException::AddInt(ubw);
		}
		BaseException::AddMessage("\n");
	}

	NPDException::NPDException(const GeneralMatrix& A)
		: Runtime_error()
	{
		Select = BaseException::Select;
		AddMessage("detected by Newmat: matrix not positive definite\n\n");
		MatrixDetails(A);
		Tracer::AddTrace();
	}

	SingularException::SingularException(const GeneralMatrix& A)
		: Runtime_error()
	{
		Select = BaseException::Select;
		AddMessage("detected by Newmat: matrix is singular\n\n");
		MatrixDetails(A);
		Tracer::AddTrace();
	}

	ConvergenceException::ConvergenceException(const GeneralMatrix& A)
		: Runtime_error()
	{
		Select = BaseException::Select;
		AddMessage("detected by Newmat: process fails to converge\n\n");
		MatrixDetails(A);
		Tracer::AddTrace();
	}

	ConvergenceException::ConvergenceException(const char* c) : Runtime_error()
	{
		Select = BaseException::Select;
		AddMessage("detected by Newmat: ");
		AddMessage(c); AddMessage("\n\n");
		if (c) Tracer::AddTrace();
	}

	OverflowException::OverflowException(const char* c) : Runtime_error()
	{
		Select = BaseException::Select;
		AddMessage("detected by Newmat: ");
		AddMessage(c); AddMessage("\n\n");
		if (c) Tracer::AddTrace();
	}

	ProgramException::ProgramException(const char* c) : Logic_error()
	{
		Select = BaseException::Select;
		AddMessage("detected by Newmat: ");
		AddMessage(c); AddMessage("\n\n");
		if (c) Tracer::AddTrace();
	}

	ProgramException::ProgramException(const char* c, const GeneralMatrix& A)
		: Logic_error()
	{
		Select = BaseException::Select;
		AddMessage("detected by Newmat: ");
		AddMessage(c); AddMessage("\n\n");
		MatrixDetails(A);
		if (c) Tracer::AddTrace();
	}

	ProgramException::ProgramException(const char* c, const GeneralMatrix& A,
		const GeneralMatrix& B) : Logic_error()
	{
		Select = BaseException::Select;
		AddMessage("detected by Newmat: ");
		AddMessage(c); AddMessage("\n\n");
		MatrixDetails(A); MatrixDetails(B);
		if (c) Tracer::AddTrace();
	}

	ProgramException::ProgramException(const char* c, MatrixType a, MatrixType b)
		: Logic_error()
	{
		Select = BaseException::Select;
		AddMessage("detected by Newmat: ");
		AddMessage(c); AddMessage("\nMatrixTypes = ");
		AddMessage(a.Value()); AddMessage("; ");
		AddMessage(b.Value()); AddMessage("\n\n");
		if (c) Tracer::AddTrace();
	}

	VectorException::VectorException() : Logic_error()
	{
		Select = BaseException::Select;
		AddMessage("detected by Newmat: cannot convert matrix to vector\n\n");
		Tracer::AddTrace();
	}

	VectorException::VectorException(const GeneralMatrix& A)
		: Logic_error()
	{
		Select = BaseException::Select;
		AddMessage("detected by Newmat: cannot convert matrix to vector\n\n");
		MatrixDetails(A);
		Tracer::AddTrace();
	}

	NotSquareException::NotSquareException(const GeneralMatrix& A)
		: Logic_error()
	{
		Select = BaseException::Select;
		AddMessage("detected by Newmat: matrix is not square\n\n");
		MatrixDetails(A);
		Tracer::AddTrace();
	}

	NotSquareException::NotSquareException()
		: Logic_error()
	{
		Select = BaseException::Select;
		AddMessage("detected by Newmat: matrix is not square\n\n");
		Tracer::AddTrace();
	}

	SubMatrixDimensionException::SubMatrixDimensionException()
		: Logic_error()
	{
		Select = BaseException::Select;
		AddMessage("detected by Newmat: incompatible submatrix dimension\n\n");
		Tracer::AddTrace();
	}

	IncompatibleDimensionsException::IncompatibleDimensionsException()
		: Logic_error()
	{
		Select = BaseException::Select;
		AddMessage("detected by Newmat: incompatible dimensions\n\n");
		Tracer::AddTrace();
	}

	IncompatibleDimensionsException::IncompatibleDimensionsException
		(const GeneralMatrix& A, const GeneralMatrix& B)
		: Logic_error()
	{
		Select = BaseException::Select;
		AddMessage("detected by Newmat: incompatible dimensions\n\n");
		MatrixDetails(A); MatrixDetails(B);
		Tracer::AddTrace();
	}

	IncompatibleDimensionsException::IncompatibleDimensionsException
		(const GeneralMatrix& A)
		: Logic_error()
	{
		Select = BaseException::Select;
		AddMessage("detected by Newmat: incompatible dimensions\n\n");
		MatrixDetails(A);
		Tracer::AddTrace();
	}

	NotDefinedException::NotDefinedException(const char* op, const char* matrix)
		: Logic_error()
	{
		Select = BaseException::Select;
		AddMessage("detected by Newmat: ");
		AddMessage(op);
		AddMessage(" not defined for ");
		AddMessage(matrix);
		AddMessage("\n\n");
		Tracer::AddTrace();
	}

	CannotBuildException::CannotBuildException(const char* matrix)
		: Logic_error()
	{
		Select = BaseException::Select;
		AddMessage("detected by Newmat: cannot build matrix type ");
		AddMessage(matrix); AddMessage("\n\n");
		Tracer::AddTrace();
	}

	IndexException::IndexException(int i, const GeneralMatrix& A)
		: Logic_error()
	{
		Select = BaseException::Select;
		AddMessage("detected by Newmat: index error: requested index = ");
		AddInt(i); AddMessage("\n\n");
		MatrixDetails(A);
		Tracer::AddTrace();
	}

	IndexException::IndexException(int i, int j, const GeneralMatrix& A)
		: Logic_error()
	{
		Select = BaseException::Select;
		AddMessage("detected by Newmat: index error: requested indices = ");
		AddInt(i); AddMessage(", "); AddInt(j);
		AddMessage("\n\n");
		MatrixDetails(A);
		Tracer::AddTrace();
	}


	IndexException::IndexException(int i, const GeneralMatrix& A, bool)
		: Logic_error()
	{
		Select = BaseException::Select;
		AddMessage("detected by Newmat: element error: requested index (wrt 0) = ");
		AddInt(i);
		AddMessage("\n\n");
		MatrixDetails(A);
		Tracer::AddTrace();
	}

	IndexException::IndexException(int i, int j, const GeneralMatrix& A, bool)
		: Logic_error()
	{
		Select = BaseException::Select;
		AddMessage(
			"detected by Newmat: element error: requested indices (wrt 0) = ");
		AddInt(i); AddMessage(", "); AddInt(j);
		AddMessage("\n\n");
		MatrixDetails(A);
		Tracer::AddTrace();
	}

	InternalException::InternalException(const char* c) : Logic_error()
	{
		Select = BaseException::Select;
		AddMessage("internal error detected by Newmat: please inform author\n");
		AddMessage(c); AddMessage("\n\n");
		Tracer::AddTrace();
	}




	/************************* ExeCounter functions *****************************/

#ifdef DO_REPORT

	int ExeCounter::nreports;                      // will be set to zero

	ExeCounter::ExeCounter(int xl, int xf) : line(xl), fileid(xf), nexe(0) {}

	ExeCounter::~ExeCounter()
	{
		nreports++;
		cout << "REPORT  " << setw(6) << nreports << "  "
			<< setw(6) << fileid << "  " << setw(6) << line
			<< "  " << setw(6) << nexe << "\n";
	}

#endif

	/**************************** error handler *******************************/

	void MatrixErrorNoSpace(const void* v) { if (!v) Throw(Bad_alloc()); }
	// throw exception if v is null




	/************************* miscellanous errors ***************************/


	void CroutMatrix::GetRow(MatrixRowCol&)
	{ Throw(NotDefinedException("GetRow","Crout")); }
	void CroutMatrix::GetCol(MatrixRowCol&)
	{ Throw(NotDefinedException("GetCol","Crout")); }
	void BandLUMatrix::GetRow(MatrixRowCol&)
	{ Throw(NotDefinedException("GetRow","BandLUMatrix")); }
	void BandLUMatrix::GetCol(MatrixRowCol&)
	{ Throw(NotDefinedException("GetCol","BandLUMatrix")); }
	void BaseMatrix::IEQND() const
	{ Throw(NotDefinedException("inequalities", "matrices")); }


#ifdef use_namespace
}
#endif

//$$ submat.cpp                         submatrices

// Copyright (C) 1991,2,3,4: R B Davies


#ifdef use_namespace
namespace NEWMAT {
#endif

#ifdef DO_REPORT
#define REPORT { static ExeCounter ExeCount(__LINE__,11); ++ExeCount; }
#else
#define REPORT {}
#endif


	/****************************** submatrices *********************************/

	GetSubMatrix BaseMatrix::submatrix(int first_row, int last_row, int first_col,
		int last_col) const
	{
		REPORT
			Tracer tr("submatrix");
		int a = first_row - 1; int b = last_row - first_row + 1;
		int c = first_col - 1; int d = last_col - first_col + 1;
		if (a<0 || b<0 || c<0 || d<0) Throw(SubMatrixDimensionException());
		// allow zero rows or columns
		return GetSubMatrix(this, a, b, c, d, false);
	}

	GetSubMatrix BaseMatrix::sym_submatrix(int first_row, int last_row) const
	{
		REPORT
			Tracer tr("sym_submatrix");
		int a = first_row - 1; int b = last_row - first_row + 1;
		if (a<0 || b<0) Throw(SubMatrixDimensionException());
		// allow zero rows or columns
		return GetSubMatrix( this, a, b, a, b, true);
	}

	GetSubMatrix BaseMatrix::row(int first_row) const
	{
		REPORT
			Tracer tr("SubMatrix(row)");
		int a = first_row - 1;
		if (a<0) Throw(SubMatrixDimensionException());
		return GetSubMatrix(this, a, 1, 0, -1, false);
	}

	GetSubMatrix BaseMatrix::rows(int first_row, int last_row) const
	{
		REPORT
			Tracer tr("SubMatrix(rows)");
		int a = first_row - 1; int b = last_row - first_row + 1;
		if (a<0 || b<0) Throw(SubMatrixDimensionException());
		// allow zero rows or columns
		return GetSubMatrix(this, a, b, 0, -1, false);
	}

	GetSubMatrix BaseMatrix::column(int first_col) const
	{
		REPORT
			Tracer tr("SubMatrix(column)");
		int c = first_col - 1;
		if (c<0) Throw(SubMatrixDimensionException());
		return GetSubMatrix(this, 0, -1, c, 1, false);
	}

	GetSubMatrix BaseMatrix::columns(int first_col, int last_col) const
	{
		REPORT
			Tracer tr("SubMatrix(columns)");
		int c = first_col - 1; int d = last_col - first_col + 1;
		if (c<0 || d<0) Throw(SubMatrixDimensionException());
		// allow zero rows or columns
		return GetSubMatrix(this, 0, -1, c, d, false);
	}

	void GetSubMatrix::SetUpLHS()
	{
		REPORT
			Tracer tr("SubMatrix(LHS)");
		const BaseMatrix* bm1 = bm;
		GeneralMatrix* gm1 = ((BaseMatrix*&)bm)->Evaluate();
		if ((BaseMatrix*)gm1!=bm1)
			Throw(ProgramException("Invalid LHS"));
		if (row_number < 0) row_number = gm1->Nrows();
		if (col_number < 0) col_number = gm1->Ncols();
		if (row_skip+row_number > gm1->Nrows()
			|| col_skip+col_number > gm1->Ncols())
			Throw(SubMatrixDimensionException());
	}

	void GetSubMatrix::operator<<(const BaseMatrix& bmx)
	{
		REPORT
			Tracer tr("SubMatrix(<<)"); GeneralMatrix* gmx = 0;
		Try
		{
			SetUpLHS(); gmx = ((BaseMatrix&)bmx).Evaluate();
			if (row_number != gmx->Nrows() || col_number != gmx->Ncols())
				Throw(IncompatibleDimensionsException());
			MatrixRow mrx(gmx, LoadOnEntry);
			MatrixRow mr(gm, LoadOnEntry+StoreOnExit+DirectPart, row_skip);
			// do need LoadOnEntry
			MatrixRowCol sub; int i = row_number;
			while (i--)
			{
				mr.SubRowCol(sub, col_skip, col_number);   // put values in sub
				sub.Copy(mrx); mr.Next(); mrx.Next();
			}
			gmx->tDelete();
		}

		CatchAll
		{
			if (gmx) gmx->tDelete();
			ReThrow;
		}
	}

	void GetSubMatrix::operator=(const BaseMatrix& bmx)
	{
		REPORT
			Tracer tr("SubMatrix(=)"); GeneralMatrix* gmx = 0;
		// MatrixConversionCheck mcc;         // Check for loss of info
		Try
		{
			SetUpLHS(); gmx = ((BaseMatrix&)bmx).Evaluate();
			if (row_number != gmx->Nrows() || col_number != gmx->Ncols())
				Throw(IncompatibleDimensionsException());
			LoadAndStoreFlag lasf =
				(  row_skip == col_skip
				&& gm->type().is_symmetric()
				&& gmx->type().is_symmetric() )
				? LoadOnEntry+DirectPart
				: LoadOnEntry;
			MatrixRow mrx(gmx, lasf);
			MatrixRow mr(gm, LoadOnEntry+StoreOnExit+DirectPart, row_skip);
			// do need LoadOnEntry
			MatrixRowCol sub; int i = row_number;
			while (i--)
			{
				mr.SubRowCol(sub, col_skip, col_number);   // put values in sub
				sub.CopyCheck(mrx); mr.Next(); mrx.Next();
			}
			gmx->tDelete();
		}

		CatchAll
		{
			if (gmx) gmx->tDelete();
			ReThrow;
		}
	}

	void GetSubMatrix::operator<<(const double* r)
	{
		REPORT
			Tracer tr("SubMatrix(<<double*)");
		SetUpLHS();
		if (row_skip+row_number > gm->Nrows() || col_skip+col_number > gm->Ncols())
			Throw(SubMatrixDimensionException());
		MatrixRow mr(gm, LoadOnEntry+StoreOnExit+DirectPart, row_skip);
		// do need LoadOnEntry
		MatrixRowCol sub; int i = row_number;
		while (i--)
		{
			mr.SubRowCol(sub, col_skip, col_number);   // put values in sub
			sub.Copy(r); mr.Next();
		}
	}

	void GetSubMatrix::operator<<(const float* r)
	{
		REPORT
			Tracer tr("SubMatrix(<<float*)");
		SetUpLHS();
		if (row_skip+row_number > gm->Nrows() || col_skip+col_number > gm->Ncols())
			Throw(SubMatrixDimensionException());
		MatrixRow mr(gm, LoadOnEntry+StoreOnExit+DirectPart, row_skip);
		// do need LoadOnEntry
		MatrixRowCol sub; int i = row_number;
		while (i--)
		{
			mr.SubRowCol(sub, col_skip, col_number);   // put values in sub
			sub.Copy(r); mr.Next();
		}
	}

	void GetSubMatrix::operator<<(const int* r)
	{
		REPORT
			Tracer tr("SubMatrix(<<int*)");
		SetUpLHS();
		if (row_skip+row_number > gm->Nrows() || col_skip+col_number > gm->Ncols())
			Throw(SubMatrixDimensionException());
		MatrixRow mr(gm, LoadOnEntry+StoreOnExit+DirectPart, row_skip);
		// do need LoadOnEntry
		MatrixRowCol sub; int i = row_number;
		while (i--)
		{
			mr.SubRowCol(sub, col_skip, col_number);   // put values in sub
			sub.Copy(r); mr.Next();
		}
	}

	void GetSubMatrix::operator=(Real r)
	{
		REPORT
			Tracer tr("SubMatrix(=Real)");
		SetUpLHS();
		MatrixRow mr(gm, LoadOnEntry+StoreOnExit+DirectPart, row_skip);
		// do need LoadOnEntry
		MatrixRowCol sub; int i = row_number;
		while (i--)
		{
			mr.SubRowCol(sub, col_skip, col_number);   // put values in sub
			sub.Copy(r); mr.Next();
		}
	}

	void GetSubMatrix::inject(const GeneralMatrix& gmx)
	{
		REPORT
			Tracer tr("SubMatrix(inject)");
		SetUpLHS();
		if (row_number != gmx.Nrows() || col_number != gmx.Ncols())
			Throw(IncompatibleDimensionsException());
		MatrixRow mrx((GeneralMatrix*)(&gmx), LoadOnEntry);
		MatrixRow mr(gm, LoadOnEntry+StoreOnExit+DirectPart, row_skip);
		// do need LoadOnEntry
		MatrixRowCol sub; int i = row_number;
		while (i--)
		{
			mr.SubRowCol(sub, col_skip, col_number);   // put values in sub
			sub.Inject(mrx); mr.Next(); mrx.Next();
		}
	}

	void GetSubMatrix::operator+=(const BaseMatrix& bmx)
	{
		REPORT
			Tracer tr("SubMatrix(+=)"); GeneralMatrix* gmx = 0;
		// MatrixConversionCheck mcc;         // Check for loss of info
		Try
		{
			SetUpLHS(); gmx = ((BaseMatrix&)bmx).Evaluate();
			if (row_number != gmx->Nrows() || col_number != gmx->Ncols())
				Throw(IncompatibleDimensionsException());
			MatrixRow mrx(gmx, LoadOnEntry);
			MatrixRow mr(gm, LoadOnEntry+StoreOnExit+DirectPart, row_skip);
			// do need LoadOnEntry
			MatrixRowCol sub; int i = row_number;
			while (i--)
			{
				mr.SubRowCol(sub, col_skip, col_number);   // put values in sub
				sub.Check(mrx);                            // check for loss of info
				sub.Add(mrx); mr.Next(); mrx.Next();
			}
			gmx->tDelete();
		}

		CatchAll
		{
			if (gmx) gmx->tDelete();
			ReThrow;
		}
	}

	void GetSubMatrix::operator-=(const BaseMatrix& bmx)
	{
		REPORT
			Tracer tr("SubMatrix(-=)"); GeneralMatrix* gmx = 0;
		// MatrixConversionCheck mcc;         // Check for loss of info
		Try
		{
			SetUpLHS(); gmx = ((BaseMatrix&)bmx).Evaluate();
			if (row_number != gmx->Nrows() || col_number != gmx->Ncols())
				Throw(IncompatibleDimensionsException());
			MatrixRow mrx(gmx, LoadOnEntry);
			MatrixRow mr(gm, LoadOnEntry+StoreOnExit+DirectPart, row_skip);
			// do need LoadOnEntry
			MatrixRowCol sub; int i = row_number;
			while (i--)
			{
				mr.SubRowCol(sub, col_skip, col_number);   // put values in sub
				sub.Check(mrx);                            // check for loss of info
				sub.Sub(mrx); mr.Next(); mrx.Next();
			}
			gmx->tDelete();
		}

		CatchAll
		{
			if (gmx) gmx->tDelete();
			ReThrow;
		}
	}

	void GetSubMatrix::operator+=(Real r)
	{
		REPORT
			Tracer tr("SubMatrix(+= or -= Real)");
		// MatrixConversionCheck mcc;         // Check for loss of info
		Try
		{
			SetUpLHS();
			MatrixRow mr(gm, LoadOnEntry+StoreOnExit+DirectPart, row_skip);
			// do need LoadOnEntry
			MatrixRowCol sub; int i = row_number;
			while (i--)
			{
				mr.SubRowCol(sub, col_skip, col_number);   // put values in sub
				sub.Check();                               // check for loss of info
				sub.Add(r); mr.Next();
			}
		}

		CatchAll
		{
			ReThrow;
		}
	}

	void GetSubMatrix::operator*=(Real r)
	{
		REPORT
			Tracer tr("SubMatrix(*= or /= Real)");
		// MatrixConversionCheck mcc;         // Check for loss of info
		Try
		{
			SetUpLHS();
			MatrixRow mr(gm, LoadOnEntry+StoreOnExit+DirectPart, row_skip);
			// do need LoadOnEntry
			MatrixRowCol sub; int i = row_number;
			while (i--)
			{
				mr.SubRowCol(sub, col_skip, col_number);   // put values in sub
				sub.Multiply(r); mr.Next();
			}
		}

		CatchAll
		{
			ReThrow;
		}
	}

#ifdef use_namespace
}
#endif

