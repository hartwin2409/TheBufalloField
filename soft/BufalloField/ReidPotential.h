#pragma once

// std
#include <stdexcept>

/*
	Literature:

		Reid, Roderick V. (1968). "Local phenomenological nucleon–nucleon potentials". Annals of Physics. 50 (3): 411–448.
		Bibcode:1968AnPhy..50..411R. doi:10.1016/0003-4916(68)90126-7.

		Nuclear magnetic shielding in the hydrogen molecule
		Reid, Roderick V., Jr. 
		1974AnPhy..84..535R 1974/05 

*/

namespace physics{
namespace ReidPotential{

	constexpr double m =  0.7 / physics::femtoMeter;

	/*
		model of strong nuclear force 1968 see https://en.wikipedia.org/wiki/Nuclear_force

		Vreid( r) = [ -10.463 * exp( -m * r) - 1650.6 * exp( -4 * m * r) + 6484.2 * exp(-7 * m * r) ] / ( m * r)

		r : distance between the two nucleii
		m : 0.7 /  fm

		The nuclear force is powerfully attractive between nucleons at distances of about
		1 femtometre (fm, or 1.0 × 10−15 metres), but it rapidly decreases to insignificance
		at distances beyond about 2.5 fm.
		At distances less than 0.7 fm, the nuclear force becomes repulsive.
		This repulsive component is responsible for the physical size of nuclei,
		since the nucleons can come no closer than the force allows. 
	*/


	/**
		Helper function Hf(g, a, r)

		f(r) = g * exp( a * m * r) / ( m * r)
	*/
	static double HelperFct( double g, double a, double _radius){
		return g *  exp( a * m * _radius) / ( m * _radius);
	}


	/**
	f(r)' = ( g * exp( a * m * r) / ( m * r))'
	= g * [ ( a * m) * exp( a * m * r) / ( m * r) + exp( a * m * r)( -1 / ( m * r * r)]
	= g * [ ( a * m)                   / ( m * r) +                ( -1 / ( m * r * r)] * exp( a * m * r)
	= g * [ ( a * m * r) / ( m * r * r) - 1 / ( m * r * r)] * exp( a * m * r)
	= g * [ ( a * m * r) / m            - 1 /   m         ] * exp( a * m * r) / ( r * r)
	= g * [   a * r     - 1 /  m         ] * exp( a * m * r) / ( r * r)
	f'(r) = g *( r* a   - 1/m) exp( a * m * r) / ( r * r)
	*/
	static double HelperDeviationFct( double g, double a,
		_In_ double _radius     ///< distance between two nucleii
	){
		return g * ( _radius * a  - 1 / m) * exp( a * m * _radius) / ( _radius * _radius);
	}


	/**
		Potential of strong nuclear force in [MeV]
	*/
	static double PotentialMeV(
		_In_ double _radius    ///< distance between two nucleii
	){
		double result2 =
			+ HelperFct( -10.463, -1, _radius)
			+ HelperFct( -1650.6, -4, _radius)
			+ HelperFct(  6484.2, -7, _radius);
		return result2;
	}


	/**
	Potential of strong nuclear force in [J]
	*/
	static double Potential(
		_In_ double _radius    ///< distance between two nucleii
	){
		return PotentialMeV( _radius) / physics::MeV;
	}
	/**
		Strong nuclear force between two nucleons like neutron or proton in [N]
			Neutron1                 Neutron2
			. x x x                    x x x   
			x       x                x       x 
			x   +   x                x   +   x 
			x       x                x       x 
			. x x x                    x x x   
			.   |                        |
			.   +- - - distance - - - - -+
	*/
	static double StrongNuclearForce(
		double _radius               // distance between two nucleii
	){
		if( _radius < -0.01 * physics::femtoMeter)
			return NAN;
		double mev =
			- HelperDeviationFct( -10.463,  -1, _radius)
			- HelperDeviationFct( -1650.6,  -4, _radius)
			- HelperDeviationFct(  6484.2,  -7, _radius);
		return mev * MeV;
	}


	static double StrongNuclearForceAttractivePart(
		double _radius               // distance between two nucleii
	){
		if( _radius < -0.01 * physics::femtoMeter)
			return NAN;
		double mev =
			- HelperDeviationFct( -10.463,  -1, _radius)
			- HelperDeviationFct( -1650.6,  -4, _radius);
		return mev * MeV;
	}


	static double StrongNuclearForceRepellingPart(
		double _radius               // distance between two nucleii
	){
		if( _radius < -0.01 * physics::femtoMeter)
			return NAN;
		double mev =
			- HelperDeviationFct(  6484.2,  -7, _radius);
		return mev * MeV;
	}


	static double GetMaximumForce(){
		double start = 0.5 * physics::femtoMeter;
		double delta = 0.1 * physics::femtoMeter;

		double f0 = 0;
		double f1 = 0;


		do{
			double d0    = start;
			double d1    = d0 + delta;

			f0 = StrongNuclearForce( d0);
			f1 = StrongNuclearForce( d1);

			do{
				double d2 = d1 + delta;
				double f2 = StrongNuclearForce( d2);

				if( f0 <= f1 && f2 <= f1){
					break;
				}
				d0 = d1; f0 = f1;
				d1 = d2; f1 = f2;

				if( d2 > 2){
					throw  std::runtime_error( " reid potential maximum force failed");
				}
			} while( 1);

			delta /= 10;
			start = d0;
		} while( fabs( delta) >= 0.01 * physics::femtoMeter);
		return f1;
	}


}}