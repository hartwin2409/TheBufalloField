#pragma once


// c-runtime
#define _USE_MATH_DEFINES
#include <math.h>

namespace gaussian{


	// gaussian function of f(x) = ( 1 / A) * exp( - 0.5 * ( ( x - mean) / radius)**2)
	// with
	//      A =  1 / ( radius * sqrt( 2 * M_PI))
	inline double Distribution1D( double x, double sigma = 1.0, double mean = 0.0){
		double A = 1 / ( sigma * sqrt( 2 * M_PI));
		double f = ( x - mean) /  sigma;
		return A * exp( -0.5 * f *  f);
	}


	// gaussian gradient function of f(x) = ( 1 / A) * ( (m-x) radius**2) * exp( - 0.5 * ( ( x - mean) / radius)**2)
	// with
	//      A =  1 / ( radius * sqrt( 2 * M_PI))
	inline double DistributionGradient1D( double x, double sigma = 1.0, double mean = 0.0){
		double A = (mean - x) / ( sigma * sigma * sigma * sqrt( 2 * M_PI));
		double f = ( x - mean) /  sigma;
		return A * exp( -0.5 * f *  f);
	}


	// gaussian function of f(x) = ( 1 / A) * exp( - 0.5 * ( ( x*x + y*y + z*z) / ( radius**2))
	// with
	//      A =  1 / ( radius **3 * pow( 2 * M_PI, 3/2))
	inline double Distribution3D( double x, double y, double z, double sigma = 1.0){
		//double meanX = meanY = meanZ = 0.0;
		double A = 1 / ( sigma * sigma * sigma * pow( 2 * M_PI, 3.0 / 2.0));
		double f = ( x * x + y * y + z * z) / ( sigma * sigma);
		return A * exp( -0.5 * f);
	}


	// gaussian function of f(x) = ( 1 / A) * exp( - 0.5 * ( ( x*x + y*y + z*z) / ( radius**2))
	// with
	//      A =  1 / ( radius **3 * pow( 2 * M_PI, 3/2))
	inline double DistributionGradientX3D( double x, double y, double z, double sigma = 1.0){
		//double meanX = meanY = meanZ = 0.0;
		double A = x / ( pow( sigma, 5) * pow( 2 * M_PI, 3.0 / 2.0));
		double f = ( x * x + y * y + z * z) / ( sigma * sigma);
		return A * exp( -0.5 * f);
	}


	// gaussian function of f(x) = ( 1 / A) * exp( - 0.5 * ( ( radius) / radius)**2)
	// with
	//      A =  1 / ( radius **3 * pow( 2 * M_PI, 3/2))
	inline double Distribution3D( double radius, double sigma = 1.0){
		//double meanX = meanY = meanZ = 0.0;
		double A = 1 / ( sigma * sigma * sigma * pow( 2 * M_PI, 3.0 / 2.0));
		double f = ( radius * radius) / ( sigma * sigma);
		return A * exp( - 0.5 * f);
	}


	// inverse function of f(x) = ( 1 / A) * exp( - 0.5 * ( ( x - mean) / radius)**2)
	// with
	//      A =  1 / ( radius * sqrt( 2 * M_PI))
	// gives back only the positive branch
	double Inverse1D( double y, double sigma = 1.0, double mean = 0.0){
		if( y <= 0)
			return Inverse1D( 1.e-27);
		double A = sigma * sqrt( 2 * M_PI);
		double f = 2 * log( 1 / ( A * y));
		return sigma * sqrt( f) + mean;     // positive branch
		// return -radius * sqrt( _f) + mean  // is the negative branch
	}


	// inverse function of f(x) = ( 1 / A) * exp( - 0.5 * ( ( x - mean) / radius)**2)
	// with
	//      A =  1 / ( radius * sqrt( 2 * M_PI))
	// gives back only the positive branch
	double Inverse3D( double y, double sigma = 1.0, double mean = 0.0){
		if( y <= 0)
			return Inverse1D( 1.e-27);
		double A = sigma * sigma * sigma * pow( 2 * M_PI, 3.0 / 2.0);
		double f = 2 * log( 1 / ( A * y));
		return sigma * sqrt( f) + mean;     // positive branch
		// return -radius * sqrt( _f) + mean  // is the negative branch
	}


}
