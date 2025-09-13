#pragma once

namespace geom{
	namespace circle{


		static inline  double Circumference( double _radius){
			return _radius * 2 * M_PI;
		}


		static inline double Diameter( double _radius){
			return _radius * 2;
		}


		static inline double Area( double _radius){
			return  M_PI * _radius * _radius;
		}

	}


	namespace cone{ // truncated cone

		//see: "geom/_docs/Truncated Cone.odt"

		static inline double Volume(
			_In_ double _height,
			_In_ double _radiusTop,
			_In_ double _radiusBottom
		){
			return _height * M_PI * ( _radiusBottom * _radiusBottom + _radiusBottom * _radiusTop  + _radiusTop * _radiusTop) / 3;
		}


		static inline  double LateralArea( // Mantelfläche
			_In_ double _height,
			_In_ double _radiusTop,
			_In_ double _radiusBottom
		){
			double deltaRadius = _radiusBottom - _radiusTop;
			double m           = sqrt( deltaRadius * deltaRadius + _height * _height);
			return ( _radiusTop + _radiusBottom) * m * M_PI;
		}

	}

	namespace sphere{
		static inline double Volume(  double radius){
			return 4 * M_PI * radius * radius * radius / 3.0;
		}

		static inline double SurfaceArea( double radius){
			return 4 * M_PI * radius * radius;
		}
	}

}
