#ifndef _GEOAC_PARAMETERS_CPP_
#define _GEOAC_PARAMETERS_CPP_

using namespace std;

#include "GeoAc.Parameters.h"

//Initial angles at the source
	double	GeoAc_theta;	// inclination angle
	double 	GeoAc_phi;      // azimuth angle

// Parameters for configuring solver
    int     GeoAc_EqCnt;		// Number of equations to solve for
    int     GeoAc_dim;          // Number of dimensions

    bool    GeoAc_CalcAmp;      // Is amplitude to be calculated?
    bool    GeoAc_AtmoStrat;	// Is the medium stratified?

    double GeoAc_ds_min =  0.001;      // Smallest possible ds for solver
    double GeoAc_ds_max =  0.5;        // Largest possible ds for solver

// Parameters to define propagation region
    double	GeoAc_ray_limit = 	5000.0;	// s for limiting ray length
    double 	GeoAc_vert_limit = 	200.0;      // Altitude at which to stop ray tracing
    double 	GeoAc_range_limit = 2000.0;     // Maximum propgaation range from the source

// Mathematical constants
    double Pi =  3.141592653589793238462643;
    double gam = 1.4;
    double R =   287.05;
    double g =   9.8;

#endif /* _GEOAC_PARAMETERS_H_ */
