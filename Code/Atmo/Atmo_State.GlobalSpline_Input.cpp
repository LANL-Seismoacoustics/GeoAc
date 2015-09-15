# ifndef SPLINED_ATMO_CPP_
# define SPLINED_ATMO_CPP_

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "GeoAc.Parameters.h"
#include "Atmo_State.h"

using namespace std;

//-----------------------------------------------//
//---------Topographical Ground Function---------//
//-----------------------------------------------//
double r_earth = 6370.0;

double GroundTopography(double lat, double lon){
    return r_earth + z_grnd;
}


//----------------------------------//
//--------Define GSL Splines--------//
//----------------------------------//

int atmo_spline_length;
double atmo_spline_rmin, atmo_spline_rmax;

double wind_FD_cutoff = 0.5;
double wind_FD_width = 0.2;

gsl_spline *splined_c;		gsl_interp_accel *acc_c;
gsl_spline *splined_u;		gsl_interp_accel *acc_u;
gsl_spline *splined_v;		gsl_interp_accel *acc_v;
gsl_spline *splined_w;		gsl_interp_accel *acc_w;
gsl_spline *splined_rho;	gsl_interp_accel *acc_rho;

//----------------------------------------------------------//
//--------Set Limits for Breaking Out of Ray Tracing--------//
//----------------------------------------------------------//
void GeoAc_SetAtmo(){
    // In this case, the cublic spline set up defines the limits, this repeats it
    GeoAc_vert_limit = 		atmo_spline_rmax;
	GeoAc_range_limit = 	5000.0;
    
	GeoAc_grazing_vert_limit = 		0.001;
	GeoAc_grazing_nu_tolerance = 	0.00001;
}

//----------------------------------------//
//---Functions to deal with file input----//
//----------------------------------------//

int  file_length(string file_name){
	ifstream file_in;	file_in.open(file_name.c_str() );

	if(!file_in.is_open()){
		cout << "Error opening file, check file name";
		return 0;
	} else {
		int count = 0;
		string line;
	
		while(!file_in.eof()) {
			getline (file_in, line);
			count++;
		}
		file_in.close();
		return count - 1;
	}
}

int  file_width(string file_name){
	double temp;
	int known_length = file_length(file_name);
	int count = 0;

	ifstream file_in;
	file_in.open(file_name.c_str() );

	if(!file_in.is_open()){
		cout << "Error opening file, check file name";
		return 0;
	} else {
		file_in.close();
		file_in.open(file_name.c_str() );
		

		while(!file_in.eof()) {
			file_in >> temp;
			count++;
		}
	
		file_in.close();
		return (count - 1)/known_length;
	}
}

double** file_to_array_build(string file_name){
	double** array;
	int input_rows = file_length(file_name);
	int input_cols = file_width(file_name);
	
	ifstream file_in;
	file_in.open(file_name.c_str() );

	if(!file_in.is_open()){
		cout << "Error opening file, check file name";
	} else {
		array = new double* [input_rows];	
		for (int i = 0; i < input_rows; i++){     //null vector then read data from file into array
			array[i] = new double [input_cols];
			for(int j = 0; j < input_cols; j++){
				array[i][j] = 0;
				file_in >> array[i][j];
			}
		}
		file_in.close();
	}
	return array;
}

//------------------------------------------------------------------------------//
//---The cubic spline routine from GSL is used to fit input data from a file----//
//------------------------------------------------------------------------------//

void load_g2s(string file_name){
	// Load a g2s file with columns describing (r-r_g, T(r) u(r), v(r), rho(r), p(r))
	atmo_spline_length = file_length(file_name);

	splined_u = 	gsl_spline_alloc (gsl_interp_cspline, atmo_spline_length);	
	splined_v = 	gsl_spline_alloc (gsl_interp_cspline, atmo_spline_length);	
	splined_w = 	gsl_spline_alloc (gsl_interp_cspline, atmo_spline_length);	
	splined_c = 	gsl_spline_alloc (gsl_interp_cspline, atmo_spline_length);	
	splined_rho = 	gsl_spline_alloc (gsl_interp_cspline, atmo_spline_length);	

	acc_u = 	gsl_interp_accel_alloc ();
	acc_v = 	gsl_interp_accel_alloc ();
	acc_w = 	gsl_interp_accel_alloc ();
	acc_c = 	gsl_interp_accel_alloc ();
	acc_rho =	gsl_interp_accel_alloc ();

	double r_array[atmo_spline_length];
	double u_array[atmo_spline_length],	v_array[atmo_spline_length], w_array[atmo_spline_length];	
	double c_array[atmo_spline_length], rho_array[atmo_spline_length];
	
	double** IO_array;
	IO_array = file_to_array_build(file_name);

	r_array[0]	= R_Ground0;
	u_array[0] 	= 0.0;
	v_array[0] 	= 0.0;
	w_array[0] 	= 0.0;
	c_array[0]	= sqrt(R*gam*IO_array[0][1])/1000.0;
	rho_array[0]= IO_array[0][5];

	for(int i = 1; i < atmo_spline_length; i++){
		r_array[i] = 	r_earth + IO_array[i][0];
		u_array[i] = 	IO_array[i][2]/(1.0 + exp((wind_FD_cutoff - IO_array[i][0])/wind_FD_width))/1000.0;
		v_array[i] = 	IO_array[i][3]/(1.0 + exp((wind_FD_cutoff - IO_array[i][0])/wind_FD_width))/1000.0;
		w_array[i] = 	0.0;
		c_array[i] =	sqrt(R*gam*IO_array[i][1])/1000.0;
		rho_array[i]= 	IO_array[i][5];
	}

	
	atmo_spline_rmin = r_array[0];
	atmo_spline_rmax = r_array[atmo_spline_length - 1];
	GeoAc_vert_limit = r_array[atmo_spline_length - 1];
	
	gsl_spline_init (splined_u, r_array, u_array, atmo_spline_length);
	gsl_spline_init (splined_v, r_array, v_array, atmo_spline_length);
	gsl_spline_init (splined_w, r_array, w_array, atmo_spline_length);
	gsl_spline_init (splined_c, r_array, c_array, atmo_spline_length);
	gsl_spline_init (splined_rho, r_array, rho_array, atmo_spline_length);

	for(int i = 0; i < atmo_spline_length; i++){
		delete IO_array[i];
	}
	delete IO_array;
}


void load_ecmwf(string file_name){
	// Load a g2s file with columns describing (r-r_g, T(r) u(r), v(r), rho(r), p(r))
	atmo_spline_length = file_length(file_name);
    
	splined_u = 	gsl_spline_alloc (gsl_interp_cspline, atmo_spline_length);
	splined_v = 	gsl_spline_alloc (gsl_interp_cspline, atmo_spline_length);
	splined_w = 	gsl_spline_alloc (gsl_interp_cspline, atmo_spline_length);
	splined_c = 	gsl_spline_alloc (gsl_interp_cspline, atmo_spline_length);
	splined_rho = 	gsl_spline_alloc (gsl_interp_cspline, atmo_spline_length);
    
	acc_u = 	gsl_interp_accel_alloc ();
	acc_v = 	gsl_interp_accel_alloc ();
	acc_w = 	gsl_interp_accel_alloc ();
	acc_c = 	gsl_interp_accel_alloc ();
	acc_rho =	gsl_interp_accel_alloc ();
    
	double r_array[atmo_spline_length];
	double u_array[atmo_spline_length],	v_array[atmo_spline_length], w_array[atmo_spline_length];
	double c_array[atmo_spline_length], rho_array[atmo_spline_length];
	
	double** IO_array;
	IO_array = file_to_array_build(file_name);
    
	r_array[0]	= R_Ground0;
	u_array[0] 	= 0.0;
	v_array[0] 	= 0.0;
	w_array[0] 	= 0.0;
	c_array[0]	= sqrt(R*gam*IO_array[0][4])/1000.0;
	rho_array[0]= IO_array[0][5];
    
	for(int i = 1; i < atmo_spline_length; i++){
		r_array[i] = 	R_Ground0 + IO_array[i][0];
		u_array[i] = 	IO_array[i][1]/(1.0 + exp((wind_FD_cutoff - IO_array[i][0])/wind_FD_width))/1000.0;
		v_array[i] = 	IO_array[i][2]/(1.0 + exp((wind_FD_cutoff - IO_array[i][0])/wind_FD_width))/1000.0;
		w_array[i] = 	0.0;
		c_array[i] =	sqrt(R*gam*IO_array[i][4])/1000.0;
		rho_array[i]= 	IO_array[i][5];
	}
    
	
	atmo_spline_rmin = r_array[0];
	atmo_spline_rmax = r_array[atmo_spline_length - 1];
	GeoAc_vert_limit = r_array[atmo_spline_length - 1];
	
	gsl_spline_init (splined_u, r_array, u_array, atmo_spline_length);
	gsl_spline_init (splined_v, r_array, v_array, atmo_spline_length);
	gsl_spline_init (splined_w, r_array, w_array, atmo_spline_length);
	gsl_spline_init (splined_c, r_array, c_array, atmo_spline_length);
	gsl_spline_init (splined_rho, r_array, rho_array, atmo_spline_length);
    
	for(int i = 0; i < atmo_spline_length; i++){
		delete IO_array[i];
	}
	delete IO_array;
}


void clear_splines(){	
 	gsl_spline_free (splined_u);	gsl_interp_accel_free (acc_u);
	gsl_spline_free (splined_v);	gsl_interp_accel_free (acc_v);
	gsl_spline_free (splined_w);	gsl_interp_accel_free (acc_w);
	gsl_spline_free (splined_c);	gsl_interp_accel_free (acc_c);
	gsl_spline_free (splined_rho);	gsl_interp_accel_free (acc_rho);
}

//----------------------------------------------------------------//
//---Define the propagation medium parameters from the splines----//
//----------------------------------------------------------------//

double c(double r, double theta, double phi) {
	if(r > atmo_spline_rmax){
		return gsl_spline_eval(splined_c,atmo_spline_rmax, acc_c)
			+ gsl_spline_eval_deriv(splined_c,atmo_spline_rmax, acc_c)*(r - atmo_spline_rmax)
			+ 0.5*gsl_spline_eval_deriv2(splined_c, atmo_spline_rmax, acc_c)*pow(r - atmo_spline_rmax,2);
	}else if( r < atmo_spline_rmin){
		return gsl_spline_eval(splined_c,atmo_spline_rmin, acc_c)
			+ gsl_spline_eval_deriv(splined_c,atmo_spline_rmin, acc_c)*(r - atmo_spline_rmin)
			+ 0.5*gsl_spline_eval_deriv2(splined_c, atmo_spline_rmin, acc_c)*pow(r - atmo_spline_rmin,2);
	}else{
		return gsl_spline_eval(splined_c, r, acc_c);
	}
}

double u(double r, double theta, double phi) {
	if(r > atmo_spline_rmax){
		return gsl_spline_eval(splined_u,atmo_spline_rmax, acc_u)
			+ gsl_spline_eval_deriv(splined_u,atmo_spline_rmax, acc_u)*(r - atmo_spline_rmax)
			+ 0.5*gsl_spline_eval_deriv2(splined_u, atmo_spline_rmax, acc_u)*pow(r - atmo_spline_rmax,2);
	}else if( r < atmo_spline_rmin){
		return gsl_spline_eval(splined_u,atmo_spline_rmin, acc_u)
			+ gsl_spline_eval_deriv(splined_u,atmo_spline_rmin, acc_u)*(r - atmo_spline_rmin)
			+ 0.5*gsl_spline_eval_deriv2(splined_u, atmo_spline_rmin, acc_u)*pow(r - atmo_spline_rmin,2);
	}else{
		return gsl_spline_eval(splined_u, r, acc_u);
	}
}

double v(double r, double theta, double phi) {
	if(r > atmo_spline_rmax){
		return gsl_spline_eval(splined_v,atmo_spline_rmax, acc_v)
			+ gsl_spline_eval_deriv(splined_v,atmo_spline_rmax, acc_v)*(r - atmo_spline_rmax)
			+ 0.5*gsl_spline_eval_deriv2(splined_v, atmo_spline_rmax, acc_v)*pow(r - atmo_spline_rmax,2);
	}else if( r < atmo_spline_rmin){
		return gsl_spline_eval(splined_v,atmo_spline_rmin, acc_v)
			+ gsl_spline_eval_deriv(splined_v,atmo_spline_rmin, acc_v)*(r - atmo_spline_rmin)
			+ 0.5*gsl_spline_eval_deriv2(splined_v, atmo_spline_rmin, acc_v)*pow(r - atmo_spline_rmin,2);
	}else{
		return gsl_spline_eval(splined_v, r, acc_v);
	}
}

double w(double r, double theta, double phi) {
	if(r > atmo_spline_rmax){
		return gsl_spline_eval(splined_w,atmo_spline_rmax, acc_w)
			+ gsl_spline_eval_deriv(splined_w,atmo_spline_rmax, acc_w)*(r - atmo_spline_rmax)
			+ 0.5*gsl_spline_eval_deriv2(splined_w, atmo_spline_rmax, acc_w)*pow(r - atmo_spline_rmax,2);
	}else if( r < atmo_spline_rmin){
		return gsl_spline_eval(splined_w,atmo_spline_rmin, acc_w)
			+ gsl_spline_eval_deriv(splined_w,atmo_spline_rmin, acc_w)*(r - atmo_spline_rmin)
			+ 0.5*gsl_spline_eval_deriv2(splined_w, atmo_spline_rmin, acc_w)*pow(r - atmo_spline_rmin,2);
	}else{
		return gsl_spline_eval(splined_w, r, acc_w);
	}
}

double rho(double r, double theta, double phi) {
	if(r > atmo_spline_rmax){
		return gsl_spline_eval(splined_rho,atmo_spline_rmax, acc_rho)
			+ gsl_spline_eval_deriv(splined_rho,atmo_spline_rmax, acc_rho)*(r - atmo_spline_rmax)
			+ 0.5*gsl_spline_eval_deriv2(splined_rho, atmo_spline_rmax, acc_rho)*pow(r - atmo_spline_rmax,2);
	}else if( r < atmo_spline_rmin){
		return gsl_spline_eval(splined_rho,atmo_spline_rmin, acc_rho)
			+ gsl_spline_eval_deriv(splined_rho,atmo_spline_rmin, acc_rho)*(r - atmo_spline_rmin)
			+ 0.5*gsl_spline_eval_deriv2(splined_rho, atmo_spline_rmin, acc_rho)*pow(r - atmo_spline_rmin,2);
	}else{
		return gsl_spline_eval(splined_rho, r, acc_rho);
	}
}



double c_diff(double r, double theta, double phi, int dim) {
	double result = 0.0;
	if(dim==0){
		if(r > atmo_spline_rmax){
			result = gsl_spline_eval_deriv(splined_c,atmo_spline_rmax, acc_c)
				+ gsl_spline_eval_deriv2(splined_c, atmo_spline_rmax, acc_c)*(r - atmo_spline_rmax);
		} else if ( r < atmo_spline_rmin){
			result = gsl_spline_eval_deriv(splined_c, atmo_spline_rmin, acc_c)
				+ gsl_spline_eval_deriv2(splined_c, atmo_spline_rmin, acc_c)*(r - atmo_spline_rmin);
		} else{
			result = gsl_spline_eval_deriv(splined_c, r, acc_c);
		}
	}
	return result;
}

double u_diff(double r, double theta, double phi, int dim) {
	double result = 0.0;
	if(dim==0){
		if(r > atmo_spline_rmax){
			result = gsl_spline_eval_deriv(splined_u,atmo_spline_rmax, acc_u)
						+ gsl_spline_eval_deriv2(splined_u, atmo_spline_rmax, acc_u)*(r - atmo_spline_rmax);
		} else if ( r < atmo_spline_rmin){
			result = gsl_spline_eval_deriv(splined_u,atmo_spline_rmin, acc_u)
						+ gsl_spline_eval_deriv2(splined_u, atmo_spline_rmin, acc_u)*(r - atmo_spline_rmin);
		} else{
			result = gsl_spline_eval_deriv(splined_u, r, acc_u);
		}
	}
	return result;
}

double v_diff(double r, double theta, double phi, int dim) {
	double result = 0.0;
	if(dim==0){
		if(r > atmo_spline_rmax){
			result = gsl_spline_eval_deriv(splined_v,atmo_spline_rmax, acc_v)
						+ gsl_spline_eval_deriv2(splined_v, atmo_spline_rmax, acc_v)*(r - atmo_spline_rmax);
		} else if ( r < atmo_spline_rmin){
			result = gsl_spline_eval_deriv(splined_v,atmo_spline_rmin, acc_v) 
						+ gsl_spline_eval_deriv2(splined_v, atmo_spline_rmin, acc_v)*(r - atmo_spline_rmin);
		} else{
			result = gsl_spline_eval_deriv(splined_v, r, acc_v);
		}
	}
	return result;
}

double w_diff(double r, double theta, double phi, int dim) {
	double result = 0.0;
	if(dim==0){
		if(r > atmo_spline_rmax){
			result = gsl_spline_eval_deriv(splined_w,atmo_spline_rmax, acc_w)
						+ gsl_spline_eval_deriv2(splined_w, atmo_spline_rmax, acc_w)*(r - atmo_spline_rmax);
		} else if ( r < atmo_spline_rmin){
			result = gsl_spline_eval_deriv(splined_w,atmo_spline_rmin, acc_w)
						+ gsl_spline_eval_deriv2(splined_w, atmo_spline_rmin, acc_w)*(r - atmo_spline_rmin);
		} else{
			result = gsl_spline_eval_deriv(splined_w, r, acc_w);
		}
	}
	return result;
}


double c_ddiff(double r, double theta, double phi, int dim1, int dim2){
	double result;
	if(dim1==0 && dim2==0){
		if(r > atmo_spline_rmax){
			result = gsl_spline_eval_deriv2(splined_c, atmo_spline_rmax, acc_c);
		} else if ( r < atmo_spline_rmin){
			result = gsl_spline_eval_deriv2(splined_c, atmo_spline_rmin, acc_c);
		} else{
			return gsl_spline_eval_deriv2(splined_c, r, acc_c);
		}
	} else {
		result = 0.0;
	}
	return result;
}

double u_ddiff(double r, double theta, double phi, int dim1, int dim2) {
	double result;
	if(dim1==0 && dim2==0){
		if(r > atmo_spline_rmax){
			result = gsl_spline_eval_deriv2(splined_u, atmo_spline_rmax, acc_u);
		} else if ( r < atmo_spline_rmin){
			result = gsl_spline_eval_deriv2(splined_u, atmo_spline_rmin, acc_u);
		} else{
			return gsl_spline_eval_deriv2(splined_u, r, acc_u);
		}
	} else {
		result = 0.0;
	}
	return result;
}
double v_ddiff(double r, double theta, double phi, int dim1, int dim2) {
	double result;
	if(dim1==0 && dim2==0){
		if(r > atmo_spline_rmax){
			result = gsl_spline_eval_deriv2(splined_v, atmo_spline_rmax, acc_v);
		} else if ( r < atmo_spline_rmin){
			result = gsl_spline_eval_deriv2(splined_v, atmo_spline_rmin, acc_v);
		} else{
			return gsl_spline_eval_deriv2(splined_v, r, acc_v);
		}
	} else {
		result = 0.0;
	}
	return result;
}
double w_ddiff(double r, double theta, double phi, int dim1, int dim2) {
	double result;
	if(dim1==0 && dim2==0){
		if(r > atmo_spline_rmax){
			result = gsl_spline_eval_deriv2(splined_w, atmo_spline_rmax, acc_w);
		} else if ( r < atmo_spline_rmin){
			result = gsl_spline_eval_deriv2(splined_w, atmo_spline_rmin, acc_w);
		} else{
			return gsl_spline_eval_deriv2(splined_w, r, acc_w);
		}
	} else {
		result = 0.0;
	}
	return result;
}

#endif  /* SPLINED_ATMO_CPP_*/
