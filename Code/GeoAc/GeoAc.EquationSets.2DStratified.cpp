#ifndef GEOAC_EQSETS_2DSTRAT_CPP_
#define GEOAC_EQSETS_2DSTRAT_CPP_

#include <math.h>
#include <iostream>
#include <sstream>

#include "GeoAc.Parameters.h"
#include "Atmo_State.h"

using namespace std;

//----------------------------------------------//
//-------Propagate in 2D, assume stratified-----//
//----------------------------------------------//
void GeoAc_SetSystem(){
    GeoAc_dim = 2;				// Dimensions
    GeoAc_AtmoStrat = true;     // Is the medium stratified?
}

//-----------------------------------------------------------//
//-------Structure containing source functions which are-----//
//-------called multiple times in the solver-----------------//
//-----------------------------------------------------------//
struct GeoAc_Sources_Struct{
	// A structure containing the source terms for geometric acoustics in 2D
	double c_eff;			// Effective sound speed
	double c_eff_0;			// Effective sound speed at source
	double c_eff_diff;		// dc/dz = c'
	double c_eff_ddiff;		// d^2 c/d z^2 = c''
};

	struct GeoAc_Sources_Struct GeoAc_Sources = {0.0, 0.0, 0.0, 0.0};

//----------------------------------------------------------------------//
//-------Fill in solution[0][n] with the appropriate initial values-----//
//----------------------------------------------------------------------//
void GeoAc_SetInitialConditions(double ** & solution, double r0, double z0){
	GeoAc_Sources.c_eff_0 = c(0.0,0.0,z0) + u(0.0,0.0,z0)*cos(GeoAc_phi) + v(0.0,0.0,z0)*sin(GeoAc_phi);

	for(int index = 0; index < GeoAc_EqCnt; index++){
		switch(index){
			case(0):			// r(0) = 0
                solution[0][index] = r0;
                break;
                
			case(1):			// z(0) = 0
                solution[0][index] = z0;
                break;
                
			case(3):			// dr dtheta(0) = 0
			case(4):			// dz dtheta(0) = 0
				solution[0][index] =  0.0;
				break;

			case(2):			// nu_z(0) = sin(theta)
				solution[0][index] =  sin(GeoAc_theta);
				break;
			
            case(5):			// dnu_z dtheta(0) = cos(theta)
				solution[0][index] =  cos(GeoAc_theta);
				break;
			
            default:
				cout << "Unexpected index in Initial_Cond.  Model includes 6 variables." << '\n';
		}
	}
}


//--------------------------------------------------------------------------//
//-------Taylor series fit to more accurately deterine intercept values-----//
//--------------------------------------------------------------------------//
void GeoAc_ApproximateIntercept(double ** solution, int k, double* & prev){
	double dz_k = solution[k][1] - solution[k-1][1];    // set dz for step = z_k - z_{k-1}
	double dz_grnd = solution[k-1][1] - z_grnd;         // set dz from z_{k-1} to ground

	for(int index = 0; index < GeoAc_EqCnt; index++){
		prev[index] = solution[k-1][index] + (solution[k-1][index] - solution[k][index])/dz_k*dz_grnd
						+ 1.0/2.0*(solution[k][index] + solution[k-2][index] - 2.0*solution[k-1][index])/pow(dz_k,2.0)*pow(dz_grnd,2.0);
		
	}
}

//-------------------------------------------------------------------------//
//-------Fill in solution[0][n] with the appropriate reflection values-----//
//-------------------------------------------------------------------------//
void GeoAc_SetReflectionConditions(double** & solution, int k_end){
	double* prev = new double [GeoAc_EqCnt];
	GeoAc_ApproximateIntercept(solution, k_end, prev);

    double c_eff_diff = c_diff(0.0,0.0,z_grnd,2) + u_diff(0.0,0.0,z_grnd,2)*cos(GeoAc_phi) + v_diff(0.0,0.0,z_grnd,2)*sin(GeoAc_phi);
    double dnuz_ds = - GeoAc_Sources.c_eff_0 / pow(c(0.0,0.0,z_grnd),2) * c_eff_diff;
    
	for(int index = 0; index < GeoAc_EqCnt; index++){
		switch(index){
			case(0):			// r(0) = r(s0)
			case(3):			// R(0) = R(s0)
				solution[0][index] = prev[index]; 
				break;

			case(1):			// z(0) = 0
				solution[0][index] = z_grnd;
				break;
			case(2):			// nu_z(0)=	-nu_z(s0)
			case(4):			// Z(0) = 	-Z(s0)
				solution[0][index] = -prev[index];
				break;
			case(5):			// mu_zt(0) =  	-mu_z(s0) + 2 d nu_z/ ds * ds0/dtheta
				solution[0][index] = -prev[index] + 2.0 * dnuz_ds * prev[4]/(c(0.0,0.0,z_grnd)/GeoAc_Sources.c_eff_0 * prev[2]);
				break;
			default:
				cout << "Unexpected index in Initial_Cond.  Model includes 6 variables." << '\n';
		}
	}
	delete [] prev;
}


//-----------------------------------------------------------------------------------//
//-------Vary the solver step size, currently uses smaller steps near the ground-----//
//-----------------------------------------------------------------------------------//
double GeoAc_Set_ds(double* current_values){
	double z = current_values[1];
	double result = 0.05 - 0.049 * exp(-(z - z_grnd)/0.75);

	result = min(result, GeoAc_ds_max);
	result = max(result, GeoAc_ds_min);
	return result;
}

//---------------------------------------//
//-------Update the source functions-----//
//---------------------------------------//
void GeoAc_UpdateSources(double ray_length, double* current_values){
	double z = current_values[1];			// Extract current z value
	
	GeoAc_Sources.c_eff = c(0,0,z) + u(0,0,z)*cos(GeoAc_phi) + v(0,0,z)*sin(GeoAc_phi);	// Update c_{eff}
	GeoAc_Sources.c_eff_diff = c_diff(0,0,z,2) + u_diff(0,0,z,2)*cos(GeoAc_phi)			// Update d/dz c_{eff}
									+ v_diff(0,0,z,2)*sin(GeoAc_phi);							

	if(GeoAc_CalcAmp){
		GeoAc_Sources.c_eff_ddiff = c_ddiff(0,0,z,2,2) + u_ddiff(0,0,z,2,2)*cos(GeoAc_phi) // Update 2nd derivative
										+ v_ddiff(0,0,z,2,2)*sin(GeoAc_phi);
	}

}

//-----------------------------------------------------------//
//-------Evaluate the Source Equation For Specific Index-----//
//-----------------------------------------------------------//
double GeoAc_EvalSrcEq(double ray_length, double* current_values, int Eq_Number){
	double result;
	
	double nu_z =	current_values [2];		double dzt = current_values [4];
	double mu_z =  current_values [5];
	double c = GeoAc_Sources.c_eff;			double c0 = GeoAc_Sources.c_eff_0;
	double dc = GeoAc_Sources.c_eff_diff;	double ddc = GeoAc_Sources.c_eff_ddiff;
	
	switch(Eq_Number){
		case(0): // dr/ds
			result = c/c0*cos(GeoAc_theta);
			break;
		case(1): // dz/ds
			result = c/c0*nu_z;
			break;
		case(2): // d(zeta)/ds
			result = -c0/pow(c,2)*dc;
			break;
		case(3): // dR/ds
			result = dc*dzt/c0*cos(GeoAc_theta) - c/c0*sin(GeoAc_theta);
			break;
		case(4): // dZ/ds
			result = dc*dzt/c0*nu_z + c/c0*mu_z;
			break;
		case(5): // d(eta)/ds
			result = (2*pow(dc/c,2) - ddc/c)*c0/c*dzt;
			break;
	}
	return result;
}

//-------------------------------------------------------------------//
//-------Calculate the Hamiltonian (Eikonal) To Check For Errors-----//
//-------------------------------------------------------------------//
double GeoAc_EvalHamiltonian(double ray_length, double* current_values){
	double nu = sqrt( pow(cos(GeoAc_theta),2) + pow(current_values[2], 2));
	return nu - GeoAc_Sources.c_eff_0/GeoAc_Sources.c_eff;
}

//--------------------------------------------------------------------------//
//-------Check if ray has left propagation region or returned to ground-----//
//--------------------------------------------------------------------------//
bool GeoAc_BreakCheck(double ** solution, int index){
	double r 	= solution[index][0];
	double z 	= solution[index][1];
	
	bool check = false;	
	if(z > GeoAc_vert_limit)    check = true;
	if(r > GeoAc_range_limit)   check = true;
	
	return check;
}
bool GeoAc_GroundCheck(double ** solution, int index){
	double z 	= solution[index][1];
	double nu_z	= solution[index][2];

	bool check = false;	
	if(z < z_grnd) check = true;
	
	return check;
}

//----------------------------------------------------------------------------------//
//-------Calculate the travel time from source to location or between locations-----//
//----------------------------------------------------------------------------------//
double GeoAc_TravelTime(double ** solution, int index){
	double dr, dz, ds, z_avg, c_eff;
	double traveltime = 0.0;

	for (int n = 0; n < index; n++){
		dr = solution[n+1][0] - solution[n][0];
		dz = solution[n+1][1] - solution[n][1];
		z_avg = solution[n][1] + dz/2.0;

		c_eff = c(0,0,z_avg) + u(0,0,z_avg)*cos(GeoAc_phi) + v(0,0,z_avg)*sin(GeoAc_phi);

		ds = sqrt(pow(dr,2) + pow(dz,2));

		traveltime += ds/c_eff;
	}
	return traveltime;	
}

void GeoAc_TravelTimeSegment(double & time, double ** solution, int start, int end){
	double dr, dz, ds, z_avg, c_eff;
    
	for (int n = start; n < end; n++){
		dr = solution[n+1][0] - solution[n][0];
		dz = solution[n+1][1] - solution[n][1];
		z_avg = solution[n][1] + dz/2.0;
        
		c_eff = c(0,0,z_avg) + u(0,0,z_avg)*cos(GeoAc_phi) + v(0,0,z_avg)*sin(GeoAc_phi);
        
		ds = sqrt(pow(dr,2) + pow(dz,2));
        
		time += ds/c_eff;
	}
}

//--------------------------------------------------------------------------//
//------Integrate the Sutherland Bass Attenuation Through the Ray Path------//
//--------------------------------------------------------------------------//
double GeoAc_SB_Atten(double ** solution, int end, double freq){
	double dr, dz, ds, x, y, z;
    double atten = 0.0;
	for (int n = 0; n < end; n++){
		
		dr = solution[n+1][0] - solution[n][0];		// Calculate ds = sqrt(dx^2 + dy^2 + dz^2)
		dz = solution[n+1][1] - solution[n][1];
		ds = sqrt(dr*dr + dz*dz);
        
		x = (solution[n][0] + dr/2.0)*cos(GeoAc_phi);				// Calculate (x,y,z) at the midpoint between [n] and [n+1]
		y = (solution[n][0] + dr/2.0)*sin(GeoAc_phi);
		z = solution[n][1] + dz/2.0;
        
		atten += SuthBass_Alpha(x, y, z, freq)*ds;	// Add contribution to the attenuation
	}
    return atten;
}

void GeoAc_SB_AttenSegment(double & atten, double ** solution, int start, int end, double freq){
	double dr, dz, ds, x, y, z;
	for (int n = start; n < end; n++){
		
		dr = solution[n+1][0] - solution[n][0];		// Calculate ds = sqrt(dx^2 + dy^2 + dz^2)
		dz = solution[n+1][1] - solution[n][1];
		ds = sqrt(dr*dr + dz*dz);
        
		x = (solution[n][0] + dr/2.0)*cos(GeoAc_phi);				// Calculate (x,y,z) at the midpoint between [n] and [n+1]
		y = (solution[n][0] + dr/2.0)*sin(GeoAc_phi);
		z = solution[n][1] + dz/2.0;
        
		atten += SuthBass_Alpha(x, y, z, freq)*ds;	// Add contribution to the attenuation
	}
}

//-----------------------------------------------------------------------------------//
//-------Calculate the Jacobian determinant and from it the ampltude coefficient-----//
//-----------------------------------------------------------------------------------//
double GeoAc_Jacobian(double ** solution, int index){
	double r = solution[index][0];						
	double z = solution[index][1];
	double drds = c(0.0,0.0,z)/GeoAc_Sources.c_eff_0*cos(GeoAc_theta);
	double dzds = c(0.0,0.0,z)/GeoAc_Sources.c_eff_0*solution[index][2];
	double drdtheta = solution[index][3];				
	double dzdtheta = solution[index][4];
		
	return r * (drds * dzdtheta - dzds * drdtheta);
}


double GeoAc_Amplitude(double ** solution, int index){
	double z = solution[index][1];
	double D = GeoAc_Jacobian(solution, index);
	
	double Amp_Num = rho(0.0,0.0,z)*c(0.0,0.0,z)*cos(GeoAc_theta);
	double Amp_Den = rho(0.0,0.0,z_grnd)*GeoAc_Sources.c_eff_0*D;
	
	
	return 1.0/(4.0*Pi)*sqrt(fabs(Amp_Num/Amp_Den));	
}

//--------------------------------------------------------------//
//---------Count the caustics encountered by monitoring---------//
//----how many times the Jacobian determinant changes sign------//
//--------------------------------------------------------------//
int GeoAc_CausticCnt(double ** solution, int start, int end){
	int count = 0;
	double current, prev = GeoAc_Jacobian(solution, start);
	for(int n = start+1; n < end; n++){
		current = GeoAc_Jacobian(solution,n);
		if(current*prev < 0.0){
			count++;
		}
		prev = current;
	}
	return count;	
}

#endif /* GEOAC_EQSETS_2DSTRAT_CPP_ */
