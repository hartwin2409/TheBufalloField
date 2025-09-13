#pragma once

#include "gaussian.h"


namespace distribution{

	class Base{
	public:
		enum class Type{
			uniform,
			gauss,
			fine
		};
		virtual double Distribution( double radius) = 0;
		virtual double Distribution( double x, double y, double z) = 0;
		virtual double Inverse( double y) = 0;
		double DeltaX( double x1, double x2, double y, double z){
			return ( Distribution( x2, y, z) - Distribution( x1, y, z)) / ( x2 - x1);
		};
		virtual double Sigma() = 0;
		Type type = Type::uniform;

		Base( Type _type)
		: type( _type)
		{}
	};


	// Distribution of a spherical gaussian distribution
	// 
	class Gauss3D : public Base{
	public:
		double sigma = 1.0;        // ´standard deviation of gaussian distribution

		Gauss3D( double _sigma)
		: Base( Type::gauss)
		, sigma( _sigma)
		{}


		inline double Distribution( double radius) override{
			return gaussian::Distribution3D( radius, sigma);
		};


		double Distribution( double x, double y, double z) override{
			return gaussian::Distribution3D( x, y, z, sigma);
		};


		double Inverse( double y) override{
			return gaussian::Inverse3D( y, sigma);
		}

		inline double Sigma() override{
			return sigma;
		}
	};


	// Distribution of a sphere with radius.
	// Inside the radius there is a homogeneous probability density
	// The probability volume of the sphere is 1.
	class HardSphere : public Base{
	public:
		double radius  = 1.0;
		double density = 0;        // probability density in [ 1 / m^3]


		HardSphere( double _sigma)
			: Base( Type::uniform)
			, radius( _sigma)
		{
			density = 1 / ( radius * radius * radius * 4 * M_PI / 3.0);
		}


		inline double Distribution( double _radius) override{
			return _radius <= radius ? 1.0 * density: 0.0;
		};


		double Distribution( double x, double y, double z) override{
			double radius = sqrt( x * x + y * y + z * z);
			return Distribution( radius);
		};


		double Inverse( double y) override{
			if( y <= 0)
				return 1.e100;
			return y < density ? radius : 0;
		}

		inline double Sigma() override{
			return radius;
		}
	};


	// Distribution of a sphere with fine structure of binding and rest energy distributions.
	// The probability volume of the sphere is 1.
	//    _radiusBindingMass = physics::neutron::radius * 0.70; // lower limit of nucleon radius lean factor = 0.57;
	//    _restMass          = physics::neutron::mass *  1.0 / ( 80.0 + 1.0) ^= 0.012 ^=  1.2%
	//    _bindingMass       = physics::neutron::mass * 80.0 / ( 80.0 + 1.0) ^= 0.988 ^= 98.8%
	// see: NeutronFineStructure.odt
	class Fine : public Base{
	public:
		double densityBindings = 0;        // probability density in [ 1 / m^3]
		double densityRest     = 0;        // probability density in [ 1 / m^3]
		double sigmaBindings   = 0.7;
		double sigmaRest       = 1.0;
		double restShare       = 0.012;    // share of rest energy
		double bindingsShare   = 0.988;    // share of binding energy
		double bindingFactor   = 1.0;      // pre factor of sigma to calculate the gauss distribution of the binding energy
		double restFactor      = 1.0;      // pre factor of sigma to calculate the gauss distribution of the rest energy


		Fine( double _sigma, double _bindingFactor)
		: Base( Type::fine)
		, bindingFactor( _bindingFactor)
		{
			sigmaBindings   = bindingFactor * _sigma;
			sigmaRest       = restFactor    * _sigma;
			densityBindings = 0.988 / ( sigmaBindings * sigmaBindings * sigmaBindings * 4 * M_PI / 3.0);
			densityRest     = 0.012 / ( sigmaRest     * sigmaRest     * sigmaRest     * 4 * M_PI / 3.0);
		}


		inline double Distribution( double _radius) override{
			return bindingsShare * gaussian::Distribution3D( _radius, sigmaBindings)
			+      restShare     * gaussian::Distribution3D( _radius, sigmaRest);
		};


		double Distribution( double x, double y, double z) override{
			double radius = sqrt( x * x + y * y + z * z);
			return Distribution( radius);
		};


		double Inverse( double y) override{
			return bindingsShare * gaussian::Inverse3D( y, sigmaBindings);
			+      restShare     * gaussian::Inverse3D( y, sigmaRest);
		}


		inline double Sigma() override{
			return sigmaRest > sigmaBindings ? sigmaRest : sigmaBindings;
		}

	};

}