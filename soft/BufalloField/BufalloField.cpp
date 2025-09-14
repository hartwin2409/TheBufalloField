/*
   Calculates the strong forces between two neutrons.
*/

// std
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
namespace fs = std::filesystem;
#include <string>
#include <sstream>
#include <vector>

#include <assert.h>
#define _USE_MATH_DEFINES
#include <math.h>
#define HAVE_OMP
#ifdef HAVE_OMP
#include <omp.h>
#endif

namespace physics{
	constexpr double femtoMeter             = 1e-15;                // [m]
	constexpr double gravitationConstant    = 6.6743e-11;           // ± 0.00085e-11 in [N * m^2 / kg^2] bzw.[m^3/( kg * s^2)], force between two masses m1 and m2 in distance  r see: https://en.wikipedia.org/wiki/Gravitational_constant
	constexpr double lightVelocity          = 299792458;            // [m / secs] <-- exact: 
	constexpr double eV                     = 1.602176634e-19;      ///< [J] exact, see https://de.wikipedia.org/wiki/Elektronenvolt
	constexpr double MeV                    = eV * 1e6;             ///< [J] mega electron volt


	namespace gravitation{
		// gravitational acceleration of 1 kg in the domain of an
		// other mass in [m /( sec * sec)]
		static double Acceleration( double distance, double mass){ 
			return gravitationConstant  *  mass  / ( distance * distance);
		}

		// return Curvature of space and time  of a central mass in [1/m]
		static double SpaceTimeCurvature( double distance, double mass){
			double _acc = Acceleration( distance, mass);
			return _acc / (lightVelocity * lightVelocity); // in [ m/(s^2) /( m^2 / s^2)] = 1 / m
		}
	}

	namespace neutron{
		constexpr double mass           = 1.67492749804e-27;    // [kg]
		constexpr double radius          = 1.6e-15 / 2;          // [m] see: https://en.wikipedia.org/wiki/Neutron
	}


	namespace ReidPotential{

		constexpr double m =  0.7 / physics::femtoMeter;

		/*
		model of strong nuclear force 1968 see https://en.wikipedia.org/wiki/Nuclear_force

		Vreid( r) = [ -10.463 * exp( -m * r) - 1650.6 * exp( -4 * m * r) + 6484.2 * exp(-7 * m * r) ] / ( m * r)

		r : distance between the two nucleii from center to center
		m : 0.7 /  fm
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

		// Potential of strong nuclear force in [MeV]
		static double PotentialMeV(
			_In_ double _radius    ///< distance between two nucleii
		){
			double result2 =
				+ HelperFct( -10.463, -1, _radius)
				+ HelperFct( -1650.6, -4, _radius)
				+ HelperFct(  6484.2, -7, _radius);
			return result2;
		}

		//Potential of strong nuclear force in [J]
		static double Potential(
			_In_ double _radius    ///< distance between two nucleii
		){
			return PotentialMeV( _radius) / physics::MeV;
		}

		//Strong nuclear force between two nucleons like neutron or proton in [N]
		static double StrongNuclearForce(
			double _radius               // distance between two nucleii from center to center
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
	}
}

#include "Distribution.h"
#include "geom.h"


class Schwarzschildmetrik{
public:
	double centralMass;
	double distanceNN;
	double radiusSchwarzschild; // Schwarzschildradius in [m]


	Schwarzschildmetrik( double _distanceNN, double _centralMass)
	: distanceNN( _distanceNN)
	, centralMass( _centralMass)
	{
		radiusSchwarzschild = 2.0 * physics::gravitationConstant * centralMass / ( physics::lightVelocity * physics::lightVelocity); // Schwarzschildradius in [m]
	}

	double ProperTimeDifference(  // [1]
		double deltaRadius                                  // [m]   symmetrical to distanceNN
	){
		// timely curvature
		double dr_2 = deltaRadius / 2.0;
		double diff = 1.0 / ( distanceNN - dr_2) - 1.0 / ( distanceNN + dr_2);  // in [1/m]
		return (radiusSchwarzschild / 2.0) * diff;                              // in [ m * 1/m] = [1]
	}


	inline double Radius(){
		return radiusSchwarzschild;
	}


	// returns the difference of time1 and time2 separated spatially by delta radius in [1]
	//    t2/t1 = ProperTimeRatio(dr) = 1 + f(dr)     with t2 = t1 + eps
	//    1 + eps/t1 = 1 + f(dr)
	//    eps/t1 = f(dr)
	double ProperTimeRatioFraction( // [1]
		double deltaRadius          // [m]  symmetrical to distanceNN
	){
		double dr_2 = deltaRadius / 2.0;
		double rad1 =  distanceNN - dr_2;
		double rad2 =  distanceNN + dr_2;
		return radiusSchwarzschild * 0.5 * ( 1 / rad1 - 1 / rad2);
	}

	double ProperTimeRatioFraction2(  // [1]
		double deltaRadius          // [m]  symmetrical to distanceNN
	){
		double dr_2 = deltaRadius / 2.0;
		double rad1 = distanceNN - dr_2;
		double rad2 = distanceNN + dr_2;
		double timeRatioExact = sqrt( ( 1 - radiusSchwarzschild/ rad2)/( 1 - radiusSchwarzschild/ rad1));
		double eps = radiusSchwarzschild * ( 1 - rad1 / rad2) / ( rad1 - radiusSchwarzschild);
		double exact = sqrt( 1 + eps);
		double p1  = eps / 2;
		double p2  = - eps * eps / 8.0;
		double p3  = - eps * eps *eps / 16.0;
		return p1 + p2 + p3;
	}
};


class SphereAspectRatio
{
public:
	double radiusOfSphere    = physics::neutron::radius;   // unit
	double distanceNN        = physics::neutron::radius * 100;
	double centralMass       = physics::neutron::mass;
	double radiusOfCurvature = 1 / physics::gravitation::SpaceTimeCurvature( distanceNN, centralMass); // distance to curvature centre in [m]


	SphereAspectRatio( double _radiusOfSphere, double _distanceNN,double _centralMass)
		: radiusOfSphere( _radiusOfSphere)
		, distanceNN( _distanceNN)
		, centralMass( _centralMass)
	{
		radiusOfCurvature = 1 / physics::gravitation::SpaceTimeCurvature( distanceNN, centralMass);
	}


	double Curvature(){
		return 1 / radiusOfCurvature;
	}


	double RingArea( double innerRadius, double outerRadius){
		return geom::circle::Area( outerRadius) - geom::circle::Area( innerRadius);
	}


	void SetCurvature(
		double _curvature // in [ 1 / m]
	){
		radiusOfCurvature = 1 / _curvature;
	}


	void SetRadiusOfSphere(
		double _radiusOfSphere // in [ m]
	){
		radiusOfSphere = _radiusOfSphere;
	}


	/**
	Calculate cross volume with ring and areas of bottom and top face of cone
	*/
	void TubeParams(
		_Out_ double& volume,              // cone volume                    in [m^3]
		_Out_ double& faceArea,            // cross area of tube             in [m*m]
		_Out_ double& timeAreaDelta,       // different proper times (caused by curvature) converted into area delta  in [ m^2]
		_In_  double  height,              // height of the cone             in [m]
		_In_  double  rInner,              // inner radius of cone at face 1 in [m]
		_In_  double  rOuter,              // outer radius of cone at face 1 in [m]
		_In_  double  deltaZCurvature,     // ratio of curvature per delta z in [m] ??? should be [m*s]
		_In_  double  deltaR               // thickness of ring              in [m]
	){
		faceArea = geom::circle::Area( rOuter) - geom::circle::Area( rInner);
		timeAreaDelta = 2 * (rInner + deltaR / 2) * M_PI * deltaZCurvature; // approximation of very thin ring in [ m^2 * s]
		double _volOuter = geom::cone::Volume( height, rOuter, rOuter);
		double _volInner = geom::cone::Volume( height, rInner, rInner);
		volume = _volOuter - _volInner;
	}


	/**
	Calculate the difference of faces of a sphere.
	Face 1 points to the center of space-time-curvature
	Face 2 points into the opposite direction
	*/
	double FaceAreaDelta(){         // in [ m^2]
		double timeAreaDelta = 0;                                  // delta area caused by curvature in [ m^2]
		double deltaR        = radiusOfSphere / 1000;              // radius step                    in [m]
		int    numSteps      = ( int) ( radiusOfSphere / deltaR);  // number of infinitesimal tubes

		// approximate sphere by infinitesimal tubes in x-direction
		for( int i = 0; i < numSteps; i++){
			double r = ( i + 0.5) * deltaR;                             // middle radius of tube in [m]
			double h = sqrt( radiusOfSphere * radiusOfSphere - r * r);  // half height of tube   in [m]

			double rInner = r - deltaR / 2;                             // inner radius of tube  in [m]
			double rOuter = r + deltaR / 2;                             // outer radius of tube  in [m]

			double _deltaZCurvature = deltaR * 2 * h / radiusOfCurvature; // in [m]

			double _vol     = 0;                     // volume of tube
			double _fA      = 0;                     // cross area of tube
			double _fADelta = 0;                     // difference of front areas of tube caused by curvature
			double _height  = h * 2;                 // height of the tube
			TubeParams( _vol, _fA, _fADelta, _height, rInner, rOuter, _deltaZCurvature, deltaR);

			timeAreaDelta  += _fADelta;
		}
		return  timeAreaDelta;
	}


	/**
	Calculate the difference of times on the spheres surface.
	Face 1 points to the point of space-time-curvature
	Face 2 points into the opposite direction
	*/
	double TimelyFaceAreaDelta( // return in [ m^2]
	){
		double deltaArea        = 0;                                  // in [ m* m]
		Schwarzschildmetrik ssm = { distanceNN, centralMass};
		double deltaR           = radiusOfSphere / 1000;              // radius step
		int    numSteps         = ( int) ( radiusOfSphere / deltaR);  // number of rings
		double radSphere2       = radiusOfSphere * radiusOfSphere;

		// approximate sphere by tubes in x-direction
		for( int i = 0; i < numSteps; i++){
			double r1 = ( i      ) * deltaR;           // inner  radius of ring
			double r  = ( i + 0.5) * deltaR;           // middle radius of ring
			double r2 = ( i + 1.0) * deltaR;           // outer  radius of ring

			double H                 = sqrt( radSphere2 - r * r);            // half height of cone in x-direction
			double timeRatioFraction = ssm.ProperTimeRatioFraction( H * 2);  // in [1]
			double _fRA              = RingArea( r1, r2);                    // area of the ring
			deltaArea               += _fRA * timeRatioFraction;             // time dilation recalculated as area difference in [m^2 * s]
		}
		return deltaArea;
	}


	// new aspect ratio calculation, sphere is approximated by rods
	double AnalyticalFaceAreaDelta(){ // return in [ m^2]
		// Bowl radius	:   r
		// Bowl volume	:   Vb = 4 * π * r^3/ 3
		// Bowl cross area	:   Ab  =  π * r^2

		// Approximate bowl volume  by rods pointing into the direction of space curvature centre.

		// Rod:
		// Length     :   Li
		// Rod volume :   Vr= Li * Δs
		// Count      :   n
		// Face-Area  :   Δs = Ab / n

		// The face of the rod pointing into the direction of space curvature centre is smaller than the other because of the curvature:

		// Curvature factor	: f
		// Face-Area-Delta 	: εi = f * Li * Δs

		double _f  = 1.0 / radiusOfCurvature; // in [1 / m] ??? must be [ s / m]

		// Sum of all Face-Area-Deltas = Space curvature of bowl

		//	Sum( εi)  = Sum( f * Li * Δs)
		//			  = f * Sum( Vr)
		//			  = f * Vb
		//			  = f * 4 * π * r^3 / 3	
		double _result = _f * geom::sphere::Volume( radiusOfSphere);
		return _result;
	}

};



class RadialPoint{
public:
	double r;    // radius
	double v;    // value

	RadialPoint( double _r = 0, double _v = 0)
		: r( _r)
		, v( _v)
	{}
};


class BufalloNucleonBase{
public:
	distribution::Base& distrib;
	std::string         name             = "Base";                    // name of this class instance
	double              centralMass      = physics::neutron::mass;    // mass of mass in the center of gravitational field
	double              distanceNN       = physics::neutron::radius * 100; // distance between nucleon an the gravitational center
	double              nucleonMass      = physics::neutron::mass;    // mass of nucleon
	double              nucleonRadius    = physics::neutron::radius;  // radius of nucleon
	double              normVolume       = 0.0;                       // volume of nucleon within nucleon radius
	double              massDensity      = 0.0;                       // mean mass density within nucleon radius
	double              massEnergy       = 0.0;                       // mass energy in [J]=[Nm]
	double              massPressure     = 0.0;                       // in [Pa] = [J/m^3], Energy meanDensity
	double              gravForce        = 0.0;                       // in [N] 


	inline double ProbabilityDensity(         // return probability meanDensity in [1/m^3]
		double _radius                        // radius in [m]
	){
		return  distrib.Distribution( _radius);
	}


	inline double ProbabilityGradient(         // return probability meanDensity in [1/m^3]
		double radius, double delta            // radius in [m]
	){
		return (ProbabilityDensity( radius + delta) - ProbabilityDensity( radius - delta)) / delta;
	}


	//inline double CouplingFactorB_RadRad( // return coupling factor of nucleon B
	//	double radiusA,                   // radius of the shielding nucleon A
	//	double radiusB                    // radius of nucleon B
	//){
	//	double pdA = ProbabilityDensity( radiusA);
	//	double pdB = ProbabilityDensity( radiusB);
	//	if( pdB == 0)
	//		return 0; // avoid denominator of 0
	//	return pdB / ( pdA + pdB);
	//}


	inline double CouplingFactorB_RadProb( // return coupling factor of nucleon B
		double radiusA,                    // radius of the shielding nucleon A
		double probB                       // prbability density of nucleon B
	){
		double pdA = ProbabilityDensity( radiusA);
		if( probB == 0)
			return 0; // avoid denominator of 0
		return probB / ( pdA + probB);
	}


	// properties
	BufalloNucleonBase( std::string _name, double _nucleonMass, double _nucleonRadius, distribution::Base& _distrib)
		: name( _name)
		, nucleonMass( _nucleonMass)
		, nucleonRadius( _nucleonRadius)
		, distrib( _distrib)
	{
		normVolume   = geom::sphere::Volume( nucleonRadius);
		massDensity  = nucleonMass / normVolume;
		massEnergy   = nucleonMass * physics::lightVelocity * physics::lightVelocity;
		massPressure = massEnergy / normVolume;
		gravForce    =  physics::gravitation::Acceleration( distanceNN, nucleonMass) * nucleonMass;
	}


	// accelerational force of nucleon caused by central gravitational mass
	double GravitationalForce(){ // in [N]
		return  physics::gravitation::Acceleration( distanceNN, centralMass) * nucleonMass;
	}


	virtual double StrongForce(
		double _distance        ///< distance nucleon centers in [m]
	) = 0;

	virtual double PressureEnergy() = 0;
};


/*
The hard sphere nucleon has  a unique mean density.
Its radius is that of a neutron

The nucleon is approximated by a single sphere.
*/
class BufalloNucleonHardSphere : public BufalloNucleonBase{
public:

	distribution::HardSphere hardSphere           = { physics::neutron::radius}; // meanDensity distribution of nucleon
	double                   bufalloFieldPressure = 0;                           // pressure of the bufallo field in [N/m^2] = [Pa]
	double                   bufalloConstant      = 0;                           // in [1]
	double                   areaTimeDelta        = 0;                           // in [m^2]
	double                   pressureEnergy       = 0;


	BufalloNucleonHardSphere(
		std::string    _name            = "unknown",
		double         _nucleonMass     = physics::neutron::mass,
		double         _nucleonRadius   = physics::neutron::radius
	)
	:hardSphere( _nucleonRadius)
	, BufalloNucleonBase( _name, _nucleonMass, _nucleonRadius, hardSphere)
	{
		/*
		The bufallo field puts some pressure on the neutron.
		Time runs faster on the surfaces of the neutrons facing outwards.
		Therefore the pressure caused by momentum exchange over time is higher and presses the neutrons together.

		How great must the pressure be so that it exerts the same force on the neutron
		as gravity?
		*/
		// We consider two neutrons located at a distance equal to 100 times their radii
		// Therefore we have a flat space time metric and can use the Newtonian approximation
		SphereAspectRatio sar = { nucleonRadius, distanceNN, centralMass};

		areaTimeDelta     = sar.AnalyticalFaceAreaDelta();               // in [ m^2]
		bufalloConstant   = gravForce / ( areaTimeDelta * massPressure); // in N/(m^2 * J / m^3) = N / ( N m / m) = 1
		double tfad       = sar.TimelyFaceAreaDelta();                   // in [ m^2]
		double bufalloFrequencyTfad = gravForce / ( tfad * massPressure); // in N / ( m^2 * J / m^3) = N / ( m * N / m) = 1

		bufalloFieldPressure =  gravForce / areaTimeDelta;               // in [N/m^2]

		// validate
		pressureEnergy = PressureEnergy();
		assert( fabs( pressureEnergy - massEnergy) < 0.0001);
	}


	double PressureEnergy() override{
		return massPressure * normVolume;
	}


	// Magnitude of force in [N] pressing two touching nucleons of same size together
	double StrongForce(
		double _distance        ///< distance nucleon centers in [m]
	) override{
		double _strongForce = 0;
		double _crossPos = _distance / 2.0;
		// position of cross section
		if( nucleonRadius >= _crossPos){
			double _crossRad2 = nucleonRadius * nucleonRadius - _crossPos * _crossPos;
			assert( _crossRad2 >= 0);
			_strongForce = _crossRad2 * M_PI * bufalloFieldPressure; // in [ m^2 * Pa] = [ m^2 * N / m^2 ] = [N]
		}
		double couplingFactor = 0.5;
		return couplingFactor * _strongForce;
	}
};


#define SIGN(v)  ( v >= 0 ? 1 : -1)


class BufalloNucleonCartesian : public BufalloNucleonBase{
public:

	double maxXYZ               = 3 * nucleonRadius;    // grid range in x-, y-, and z-direction
	double distanceNN           = nucleonRadius * 100;  // distance between two nucleons
	int    numSteps             = 250;                  // number of steps in grid range
	double step                 = maxXYZ / numSteps;    // cartesian grid step
	double bufalloPressure      = 0.0;                  // in [N/m^2]
	double timeAreaDelta        = 0.0;                  // in [m^2]
	double bufalloConstant      = 0.0;                  // in [1]


	BufalloNucleonCartesian( std::string _name, distribution::Base& _distrib)
		: BufalloNucleonBase( _name, physics::neutron::mass, physics::neutron::radius, _distrib)
	{
		/*
		The bufallo field puts some pressure on the neutron.
		Time runs faster on the surfaces of the neutrons facing outwards.
		Therefore the pressure caused by momentum exchange over time is higher and presses the neutrons together.

		How great must the pressure be so that it exerts the same force on the neutron
		as gravity?

		[kg*(m/s)/(m*m*s)] = [ kg/m*s^2] = [ (kg*m/s^2)/ m^2] = [N/m^2] =[Pa]
		*/
		double pressureEnergy = PressureEnergy();

		timeAreaDelta   = TimeAreaDelta();                                // in [m^2]
		bufalloConstant = gravForce / ( timeAreaDelta * pressureEnergy);  // in N / ( m^2  * J/m^3) = N / ( m^2  * N * m / m^3) = 1
		bufalloPressure = gravForce / timeAreaDelta;                      // in N / m^2

		// debug
		assert( fabs( pressureEnergy / massEnergy - 1) < 0.01);
	}


	// Returns the time area delta in [ m^2]
	double TimeAreaDelta(){
		// sample quarter: octant of +x and -x
		// Assume space curvature center in  negative x - direction
		// only the gradient in x direction is considered, because of symmetry reasons the forces in z and y direction cancel out
		// calculate the weighted volume difference between the left (-x) and right(+x) side of the nucleon
		double sumAreaDelta        = 0;
		double maxRadius2          = maxXYZ * maxXYZ;
		double infArea             = step * step;              // area of a single side of the infinitesimal volume
		double infVol              = infArea * step;           // in [m^3]
		Schwarzschildmetrik ssm    = { distanceNN, centralMass};

		// omp parallel for collapse(3) - combines the three loops (zi, yi, xi) into a single iteration domain
		// reduction(+:sumAreaDelta)    - makes sure that every thread has its own local copy sumAreaDelta, which are combined at the end
		// schedule(dynamic)            -  can help to balance load imbalances caused by inner conditions  (if(radius2 <= maxRadius2)), else schedule(static).
		//#pragma omp parallel for collapse(3) reduction(+:sumForce) schedule(dynamic)
#ifdef HAVE_OMP
#pragma omp parallel for reduction(+:sumAreaDelta)
#endif
		for( int zi = 0; zi <= numSteps; zi++){
			double z = (zi + 0.5) * step;
			double radZ2 = z * z;
			for( int yi = 0; yi <= numSteps; yi++){
				double y = ( yi + 0.5) * step;
				double _radZY2 = radZ2 + y * y;
				for( int xi = 0; xi <= numSteps; xi++){
					double x  = (xi + 0.5) * step;
					double x2 = x + step;
					double radius2 = _radZY2 + x * x;
					if( radius2 <= maxRadius2){
						double radiusP1 = sqrt( radius2);                 // radius of position 1 to calculate the gradient
						double radiusP2 = sqrt( _radZY2 + x2 * x2);       // radius of position 2 to calculate the gradient
						double pdA      = distrib.Distribution( radiusP1);
						double pdB      = distrib.Distribution( radiusP2);
						double pdDelta  = ( - pdB + pdA) * SIGN(x);                           // negative Gradient in [ 1 / m^3]
						double deltaX   = x + step/2;
						double tad      = ssm.ProperTimeRatioFraction( deltaX * 2) * infArea; // time area delta in [m^2]

						sumAreaDelta += pdDelta * tad; // in [ 1 *  m^2] = [ m^2]
					}
				}
			}
		}
		return sumAreaDelta;
	}


	inline int GetSteps(){ return numSteps;}


	double StrongForce(
		double distance        ///< distance nucleon centers in [m]
	) override{
		// sample quarter of nucleon B: ( -x...+x, 0..+y, 0...+z)
		// center of nucleon B is( 0, 0, 0)
		// center of nucleon A is( -distance, 0, 0)
		// Assume space curvature center in  negative x - direction
		// calculate the weighted volume difference between the left (-x) and right(+x) side of the nucleon
		double sumForce     = 0;
		double maxRadius2   = maxXYZ * maxXYZ;
		double infArea      = step * step;         // infinitesimal area of a single side of the infinitesimal volume
		double infVol       = infArea * step;      // infinitesimal volume in [m^3]

		// omp parallel for collapse(3) - combines the three loops (zi, yi, xi) into a single iteration domain
		// reduction(+:sumAreaDelta)    - makes sure that every thread has its own local copy sumAreaDelta, which are combined at the end
		// schedule(dynamic)            -  can help to balance load imbalances caused by inner conditions  (if(radius2 <= maxRadius2)), else schedule(static).
		//#pragma omp parallel for collapse(3) reduction(+:sumForce) schedule(dynamic)
#ifdef HAVE_OMP
#pragma omp parallel for reduction(+:sumForce)
#endif
		for( int zi = 0; zi <= numSteps; zi++){
			double z = (zi + 0.5) * step;          // z coordinate of the center of the infinitesimal volume
			double radZ2 = z * z;
			for( int yi = 0; yi <= numSteps; yi++){
				double y = ( yi + 0.5) * step;               // y coordinate of the center of the infinitesimal volume
				double _radZY2 = radZ2 + y * y;
				for( int xi = -numSteps; xi <= numSteps; xi++){
					double x1 = xi * step;
					double x =  (xi + 0.5) * step;                       // x coordinate of the center of the infinitesimal volume
					double x2 = ( xi + 1.0) * step;
					double radius2 = _radZY2 + x * x;                    // squared radius from center of nucleon B to center of infinitesimal volume
					if( radius2 <= maxRadius2){
						double radiusP1 = sqrt( _radZY2 + x1 * x1);         // radius of position x1 to calculate the difference
						double radiusP2 = sqrt( _radZY2 + x2 * x2);         // radius of position x2 to calculate the difference
						double pd1      = distrib.Distribution( radiusP1);
						double pd2      = distrib.Distribution( radiusP2);
						double pdDelta  = - pd2 + pd1;                      // negative difference at infinitesimal volume in [ 1 / m^3]
						double radiusB  = sqrt( radius2);                   // radius to center of nucleon B
						double diffX    = distance + x;
						double radiusA  = sqrt( _radZY2 + diffX * diffX);   // radius to center of nucleon A

						double cfB = CouplingFactorB_RadProb( radiusA, (pd2 + pd1) / 2); // couplingfactor in [1]
						double deltaRho = pdDelta * massEnergy;             // Energy density difference at infinitesimal volume in [J/m^3]
						sumForce += cfB * infArea * deltaRho;               // [ m^2 * (N * m / m^3)] = [N]
					}
				}
			}
		}
		sumForce *= 4;   // *4 becaus only a quarter ( -x...+x, 0..+y, 0...+z) has been sampled

		return sumForce;
	}


	inline double LocalDensity( double _radius){
		return ProbabilityDensity( _radius) * nucleonMass; // [kg/m m*m]
	}


	inline double LocalPressure( double _radius){
		return LocalDensity( _radius) * physics::lightVelocity * physics::lightVelocity; // [N/m*m]
	}


	double PressureEnergy() override{
		double _sumOctant = 0;
		double _massOctant = 0;
		double _infVol = step * step * step;
		for( int zi = 0; zi <= numSteps; zi++){
			double z = ( zi + 0.5) * step;
			for( int yi = 0; yi <= numSteps; yi++){
				double y = ( yi  + 0.5) * step;
				for( int xi = 0; xi <= numSteps; xi++){
					double x = ( xi  + 0.5) * step;
					double _radius = sqrt( z * z + y * y + x * x);
					_sumOctant += LocalPressure( _radius) * _infVol;
					_massOctant += LocalDensity( _radius) * _infVol;
				}
			}
		}
		double mass = _massOctant * 8;
		double _pressEnerg = _sumOctant * 8;
		return _pressEnerg;
	}
};


void ProgressOutputToTerminal( double distance){
	int _idx = (int)( distance / 0.025 * physics::femtoMeter);
	if( ( _idx % 100) == 0){
		std::cout << "dist=" << std::setprecision( 3) << std::showpoint << std::setw( 9) << distance << " --> " << 4.0 * physics::femtoMeter << "     \r";
	}
}


// Calcluates the strong forces from bufallo field pressure
void StrongForces(){
	const double fm              = physics::femtoMeter;

	std::cout << "Strong Forces"     << std::endl;

	distribution::HardSphere distrHard    = { physics::neutron::radius};
	double                   bindingFactor = 0.37;
	distribution::Fine       distrFine     = { physics::neutron::radius, bindingFactor};
	BufalloNucleonHardSphere hardSphere( "Hard");
	BufalloNucleonCartesian cartFine( "Fine", distrFine);

	std::string filename = "StrongForces.csv";
	std::ofstream out( filename);
	out << "Distance[fm];Reid;Hard;Fine;ReidRepel*0.55;ReidAttrac;Fbufallo" << std::endl;
	for( double distance = 0.025 * fm; distance < 4.0 * fm; distance += 0.025 * fm){

		ProgressOutputToTerminal( distance);

		double reid        = physics::ReidPotential::StrongNuclearForce( distance);
		double hardForce   = hardSphere.StrongForce( distance);
		double fineForce   = cartFine.StrongForce( distance);
		double reidRepel   = 0.55 * physics::ReidPotential::StrongNuclearForceRepellingPart( distance);
		double reidAttrac  = physics::ReidPotential::StrongNuclearForceAttractivePart( distance);
		out
			<< distance / physics::femtoMeter << ";"
			<<  reid                          << ";"
			<< -hardForce                     << ";"
			<< -fineForce                     << ";"
			<<  reidRepel                     << ";"
			<<  reidAttrac                    << ";"
			<<  fineForce + reidRepel         << std::endl;
	}
	out.close();
	std::cout << "\n";
	std::cout << "results saved to \'" << filename << "\'\n";
}


int main(){
	StrongForces();
}