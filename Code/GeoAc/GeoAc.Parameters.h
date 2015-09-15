#ifndef _GEOAC_PARAMETERS_H_
#define _GEOAC_PARAMETERS_H_

using namespace std;

//Initial angles at the source
	extern double   GeoAc_theta;	// inclination angle of ray path
	extern double   GeoAc_phi;      // azimuth angle of ray path

// Parameters for configuring solver
    extern int      GeoAc_EqCnt;		// Number of equations to solve
    extern int      GeoAc_dim;          // Number of dimensions (2 or 3)

    extern bool     GeoAc_CalcAmp;      // Is amplitude to be calculated?
    extern bool     GeoAc_AtmoStrat;	// Is the medium stratified?

    extern double   GeoAc_ds_min;       // Smallest possible ds for solver
    extern double   GeoAc_ds_max;       // Largest possible ds for solver

// Parameters to define propagation region
    extern double   GeoAc_ray_limit;    // Ray length at which to stop ray tracing
    extern double   GeoAc_vert_limit;   // Altitude at which to stop ray tracing
    extern double 	GeoAc_range_limit;  // Range at which to stop ray tracing

    extern double 	GeoAc_x_min_limit;  // E-W range at which to stop ray tracing in 3D.RngDep
    extern double 	GeoAc_x_max_limit;  // E-W range at which to stop ray tracing in 3D.RngDep
    extern double 	GeoAc_y_min_limit;  // N-S range at which to stop ray tracing in 3D.RngDep
    extern double 	GeoAc_y_max_limit;  // N-S range at which to stop ray tracing in 3D.RngDep

    extern double 	GeoAc_lat_min_limit; // Minimum latitude of propagation region in Global.RngDep
    extern double 	GeoAc_lat_max_limit; // Maximum latitude of propagation region in Global.RngDep
    extern double 	GeoAc_lon_min_limit; // Minimum longitude of propagation region in Global.RngDep
    extern double 	GeoAc_lon_max_limit; // Maximum longitude of propagation region in Global.RngDep

// Mathematical constants
    extern double   Pi;
    extern double   gam;
    extern double   R;
    extern double   g;

#endif /* _GEOAC_PARAMETERS_H_ */
