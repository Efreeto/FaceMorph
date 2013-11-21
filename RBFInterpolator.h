//////////////////////////////////////////////////////////////////////////
//
// RBFInterpolator : interpolation by radial basis functions
// 
// 2009 Karsten Noe
//
// Read the blog at cg.alexandra.dk for more information
//
//////////////////////////////////////////////////////////////////////////


#ifndef RBFInterpolator_H
#define RBFInterpolator_H

#define WANT_STREAM                  
#define WANT_MATH                    

#include "newmat.h"                  // use matrix operations from the newmat library by R B Davies
									 // http://www.robertnz.net/nm_intro.htm

#include <vector>

typedef float real;

using namespace std;

class RBFInterpolator
{
public:
	RBFInterpolator();
	~RBFInterpolator();

	//create an interpolation function f that obeys F_i = f(x_i, y_i, z_i)
	RBFInterpolator(vector<real> x, vector<real> y, vector<real> z, vector<real> F);

	//specify new function values F_i while keeping the same 
	void UpdateFunctionValues(vector<real> F);

	//evaluate the interpolation function f at the 3D position (x,y,z)
	real interpolate(real x, real y, real z);

private:

	// read the tutorial on cg.alexandra.dk to make sense of these variable names ;-)
	real g(real t_squared);
	bool successfullyInitialized;
	ColumnVector A;
	Matrix P;		//contains positions
	Matrix Ginv;
	unsigned int M; 
};

#endif