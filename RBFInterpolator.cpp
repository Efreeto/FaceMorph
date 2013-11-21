//////////////////////////////////////////////////////////////////////////
//
// RBFInterpolator : interpolation by radial basis functions
// 
// 2009 Karsten Noe
//
// Read the blog at cg.alexandra.dk for more information
//
//////////////////////////////////////////////////////////////////////////

#include "RBFInterpolator.h"
#include <math.h>

RBFInterpolator::RBFInterpolator()
{
	successfullyInitialized = false;
}

RBFInterpolator::RBFInterpolator(vector<real> x, vector<real> y, vector<real> z, vector<real> f)
{
	successfullyInitialized = false; // default value for if we end init prematurely.

	M = f.size();

	// all four input vectors must have the same length.
	if ( x.size() != M || y.size() != M || z.size() != M )
		return;	

	ColumnVector F = ColumnVector(M + 4);
	P = Matrix(M, 3);

	Matrix G(M + 4,M + 4);

	// copy function values
	for (unsigned int i = 1; i <= M; i++)
		F(i) = f[i-1];
	
	F(M+1) = 0;  F(M+2) = 0;  F(M+3) = 0;  F(M+4) = 0;

	// fill xyz coordinates into P 
	for (unsigned int i = 1; i <= M; i++)
	{
		P(i,1) = x[i-1];
		P(i,2) = y[i-1];
		P(i,3) = z[i-1];
	}

	// the matrix below is symmetric, so I could save some calculations Hmmm. must be a todo
	for (unsigned int i = 1; i <= M; i++)
	for (unsigned int j = 1; j <= M; j++)
	{
		real dx = x[i-1] - x[j-1];
		real dy = y[i-1] - y[j-1];
		real dz = z[i-1] - z[j-1];

		real distance_squared = dx*dx + dy*dy + dz*dz;

		G(i,j) = g(distance_squared);
	}

	//Set last 4 columns of G
	for (unsigned int i = 1; i <= M; i++)
	{
		G( i, M+1 ) = 1;
		G( i, M+2 ) = x[i-1];
		G( i, M+3 ) = y[i-1];
		G( i, M+4 ) = z[i-1];
	}

	for (unsigned int i = M+1; i <= M+4; i++)
	for (unsigned int j = M+1; j <= M+4; j++)
		G( i, j ) = 0;

	//Set last 4 rows of G
	for (unsigned int j = 1; j <= M; j++)
	{
		G( M+1, j ) = 1;
		G( M+2, j ) = x[j-1];
		G( M+3, j ) = y[j-1];
		G( M+4, j ) = z[j-1];
	}

	Try 
	{ 
		Ginv = G.i(); 

		A = Ginv*F;
		successfullyInitialized = true;
	}
    CatchAll { cout << BaseException::what() << endl; }

}

RBFInterpolator::~RBFInterpolator()
{

}

real RBFInterpolator::interpolate(real x, real y, real z)
{
	if (!successfullyInitialized)
		return 0.0f;

	real sum = 0.0f;

	// RBF part
	for (unsigned int i = 1; i <= M; i++)
	{
		real dx = x - P(i,1);
		real dy = y - P(i,2);
		real dz = z - P(i,3);

		real distance_squared = dx*dx + dy*dy + dz*dz;

		sum += A(i) * g(distance_squared);
	}
	
	//affine part
	sum += A(M+1) + A(M+2)*x + A(M+3)*y + A(M+4)*z;

	return sum;
}

//note: assuming the input is t squared
real RBFInterpolator::g(real t_squared)
{	
	return sqrt(log10(t_squared + 1.0f));
}

void RBFInterpolator::UpdateFunctionValues(vector<real> f)
{
	successfullyInitialized = false;

	ColumnVector F(M+4);

	// copy function values
	for (unsigned int i = 1; i <= M; i++)
		F(i) = f[i-1];
	
	F(M+1) = 0;  F(M+2) = 0;  F(M+3) = 0;  F(M+4) = 0;

	Try 
	{ 
		A = Ginv*F;
		successfullyInitialized = true;
	}
    CatchAll { cout << BaseException::what() << endl; }

}
