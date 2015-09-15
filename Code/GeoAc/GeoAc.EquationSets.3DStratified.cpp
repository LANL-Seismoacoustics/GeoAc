#ifndef GEOAC_EQSETS_3DSTRAT_H_
#define GEOAC_EQSETS_3DSTRAT_H_

#include <math.h>
#include <iostream>
#include "GeoAc.Parameters.h"
#include "Atmo_State.h"

using namespace std;

//----------------------------------------------//
//-------Propagate in 3D, assume stratified-----//
//----------------------------------------------//
void GeoAc_SetSystem(){
    GeoAc_dim = 3;				// Dimensions
    GeoAc_AtmoStrat = true;     // Is the medium stratified?
}

//-----------------------------------------------------------//
//-------Structure containing source functions which are-----//
//-------called multiple times in the solver-----------------//
//-----------------------------------------------------------//
struct GeoAc_Sources_Struct{
    double src_loc[3];  // Source location
	double c0;          // Thermodynamic sound speed at source
    
    double nu0_xy[2];       // Horizontal components of the eikonal vector remain constant for stratified atmosphere
    double mu0_xy[2][2];    // Launch angle derivatives of the eikonal vector (also constant)

    double c;		// Thermodynamic sound speed
	double dc;		// dc/dz = c'
	double ddc;		// d^2 c/d z^2 = c''

	double u;		// Easterly wind speed
	double du;		// du/dz = u'
	double ddu;		// d^2 u/d z^2 = u''

	double v;		// Northerly wind speed
	double dv;		// dv/dz = v'
	double ddv;		// d^2 v/d z^2 = v''

	double w;		// Vertical wind speed
	double dw;		// dw/dz = w'
	double ddw;		// d^2 w/d z^2 = w''
	
	double nu_mag;      // Eikonal vector magnitude
    double dnu_mag[2];  // Angular derivatives of nu_mag
	
	double c_prop[3];	// Propagation velocity cp = c {nu_x,nu_y,nu_z} /|nu| + {u,v,w}
	double c_prop_mag;	// Propagation velocity magnitude
    
    double dc_prop[3][2];   // Angular derivatives of c_prop (components)
    double dc_prop_mag[2];	// Angular derivatives of c_prop (magnitudes)
};

    struct GeoAc_Sources_Struct GeoAc_Sources = {{0.0, 0.0, 0.0}, 0.0,
                                                 {0.0, 0.0}, {{0.0, 0.0}, {0.0, 0.0}},
                                                 0.0, 0.0, 0.0,
												 0.0, 0.0, 0.0,
												 0.0, 0.0, 0.0,
												 0.0, 0.0, 0.0,
                                                 0.0, {0.0, 0.0},
                                                 {0.0, 0.0, 0.0}, 0.0,
                                                 {{0.0, 0.0}, {0.0, 0.0},{0.0, 0.0}}, {0.0, 0.0}};

//----------------------------------------------------------------------//
//-------Fill in solution[0][n] with the appropriate initial values-----//
//----------------------------------------------------------------------//
void GeoAc_SetInitialConditions(double ** & solution, double x0, double y0, double z0){
    GeoAc_Sources.src_loc[0] = x0;
    GeoAc_Sources.src_loc[1] = y0;
    GeoAc_Sources.src_loc[2] = z0;
	GeoAc_Sources.c0 = c(x0, y0, z0);
    
    double M_Comps[3] = { u(x0, y0, z0) / GeoAc_Sources.c0,   v(x0, y0, z0) / GeoAc_Sources.c0,     w(x0, y0, z0) / GeoAc_Sources.c0};
    double nu0[3] =     { cos(GeoAc_theta) * cos(GeoAc_phi),  cos(GeoAc_theta) * sin(GeoAc_phi),    sin(GeoAc_theta)};
    double mu0_th[3] =  {-sin(GeoAc_theta) * cos(GeoAc_phi), -sin(GeoAc_theta) * sin(GeoAc_phi),    cos(GeoAc_theta)};
    double mu0_ph[3] =  {-cos(GeoAc_theta) * sin(GeoAc_phi),  cos(GeoAc_theta) * cos(GeoAc_phi),    0.0};
    
    double M = 1.0 + (nu0[0] * M_Comps[0] + nu0[1] * M_Comps[1] + nu0[2] * M_Comps[2]);
    double dM_th = mu0_th[0] * M_Comps[0] + mu0_th[1] * M_Comps[1] + mu0_th[2] * M_Comps[2];
    double dM_ph = mu0_ph[0] * M_Comps[0] + mu0_ph[1] * M_Comps[1] + mu0_ph[2] * M_Comps[2];
    
    GeoAc_Sources.nu0_xy[0] = nu0[0] / M;
    GeoAc_Sources.nu0_xy[1] = nu0[1] / M;
    GeoAc_Sources.mu0_xy[0][0] = mu0_th[0] / M - nu0[0] / pow(M,2.0) * dM_th;
    GeoAc_Sources.mu0_xy[1][0] = mu0_th[1] / M - nu0[1] / pow(M,2.0) * dM_th;
    GeoAc_Sources.mu0_xy[0][1] = mu0_ph[0] / M - nu0[0] / pow(M,2.0) * dM_ph;
    GeoAc_Sources.mu0_xy[1][1] = mu0_ph[1] / M - nu0[1] / pow(M,2.0) * dM_ph;

	for(int index = 0; index < GeoAc_EqCnt; index++){
		switch(index){
			case(0):			// x(0) = x0
                solution[0][index] = x0;
                break;
			
            case(1):			// y(0) = y0
                solution[0][index] = y0;
                break;
			
            case(2):			// z(0) = z0
                solution[0][index] = z0;
                break;
                
			case(4):			// Xt(0) = 0
			case(5):			// Yt(0) = 0
			case(6):			// Zt(0) = 0
                
			case(8):			// Xp(0) = 0
			case(9):			// Yp(0) = 0
			case(10):			// Zp(0) = 0
				solution[0][index] =  0.0;
				break;
                
			case(3):			// nu_z(0) = nu0/M
				solution[0][index] =  nu0[2] / M;
				break;

			case(7):			// mu_zt(0) = dnu0/M + sin(theta)/M^2 * dM
				solution[0][index] = mu0_th[2] / M - nu0[2] / pow(M,2.0) * dM_th;
				break;
            
            case(11):			// mu_zp(0) = dnu0/M + sin(theta)/M^2 * dM
                solution[0][index] = mu0_ph[2] / M - nu0[2] / pow(M,2.0) * dM_ph;
                break;
                
			default:
				cout << "Unexpected index in Initial_Cond.  Model includes 12 variables." << '\n';
		}
	}
}

//--------------------------------------------------------------------------//
//-------Taylor series fit to more accurately deterine intercept values-----//
//--------------------------------------------------------------------------//
void GeoAc_ApproximateIntercept(double ** solution, int k, double* & prev){
	double result;
	double dz_k = solution[k][2] - solution[k-1][2];	// set dz for step = z_k - z_{k-1}
	double dz_grnd = solution[k-1][2] - z_grnd;         // set dz from z_{k-1} to ground (assuming z = 0 is ground)

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

	double x = prev[0], y = prev[1];
	double dnuz_ds = - 1.0/c(x,y,z_grnd) * (GeoAc_Sources.c0 / c(x,y,z_grnd) * c_diff(x,y,z_grnd,2)
                                            + GeoAc_Sources.nu0_xy[0] * u_diff(x,y,z_grnd,2)
											 + GeoAc_Sources.nu0_xy[1] * v_diff(x,y,z_grnd,2)
											  + prev[3] * w_diff(x,y,z_grnd,2));

	for(int index = 0; index < GeoAc_EqCnt; index++){
		switch(index){
			case(0):			// x(0) = x(s0)
			case(1):			// y(0) = y(s0)
			case(2):			// z(0) = 0 = +/-z(s0)
			case(4):			// Xt(0) = Xt(s0)
			case(5):			// Yt(0) = Yt(s0)
			case(8):			// Xp(0) = Xp(s0)
			case(9):			// Yp(0) = Yp(s0)
				solution[0][index] = prev[index]; 
				break;
			case(3):			// nu_z(0) = -nu_z(s0) = sin(theta)
			case(6):			// Zt(0) = -Zt(s0)
			case(10):			// Zp(0) = -Zp(s0)
				solution[0][index] = -prev[index];
				break;
			case(7):			// mu_zt(0) = -mu_zt(s0) + 2 dnu_z/ds * ds_0/dtheta
			case(11):			// mu_zp(0) = -mu_zp(s0) + 2 dnu_z/ds * ds_0/dphi
				solution[0][index] = -prev[index] + 2.0 * dnuz_ds * prev[index - 1]/(c(x,y,z_grnd) / GeoAc_Sources.c0 * prev[3]);
				break;
			default:
				cout << "Unexpected index in Initial_Cond.  Model includes 12 variables." << '\n';
		}
	}
	delete [] prev;
}

//-----------------------------------------------------------------------------------//
//-------Vary the solver step size, currently uses smaller steps near the ground-----//
//-----------------------------------------------------------------------------------//
double GeoAc_Set_ds(double* current_values){
	double z = current_values[2];
	double result = 0.05 - 0.049 * exp(-(z - z_grnd)/0.75);

	result = min(result, GeoAc_ds_max);
	result = max(result, GeoAc_ds_min);
	return result;
}

//---------------------------------------//
//-------Update the source functions-----//
//---------------------------------------//
void GeoAc_UpdateSources(double ray_length, double* current_values){
	double x = current_values[0];
	double y = current_values[1];
	double z = current_values[2];
	double nu[3] = {GeoAc_Sources.nu0_xy[0], GeoAc_Sources.nu0_xy[1], current_values[3]};

    GeoAc_Sources.c = c(x,y,z);		GeoAc_Sources.dc = c_diff(x,y,z,2);
    GeoAc_Sources.u = u(x,y,z);		GeoAc_Sources.du = u_diff(x,y,z,2);
    GeoAc_Sources.v = v(x,y,z);		GeoAc_Sources.dv = v_diff(x,y,z,2);
    GeoAc_Sources.w = w(x,y,z);		GeoAc_Sources.dw = w_diff(x,y,z,2);

    GeoAc_Sources.nu_mag = 	GeoAc_Sources.c0 / GeoAc_Sources.c * (1.0 - (nu[0] * GeoAc_Sources.u + nu[1] * GeoAc_Sources.v + nu[2] * GeoAc_Sources.w) / GeoAc_Sources.c0);

    GeoAc_Sources.c_prop[0] = GeoAc_Sources.c * nu[0] / GeoAc_Sources.nu_mag + GeoAc_Sources.u;
    GeoAc_Sources.c_prop[1] = GeoAc_Sources.c * nu[1] / GeoAc_Sources.nu_mag + GeoAc_Sources.v;
    GeoAc_Sources.c_prop[2] = GeoAc_Sources.c * nu[2] / GeoAc_Sources.nu_mag + GeoAc_Sources.w;
    GeoAc_Sources.c_prop_mag = sqrt(pow(GeoAc_Sources.c_prop[0],2) + pow(GeoAc_Sources.c_prop[1],2) + pow(GeoAc_Sources.c_prop[2],2));

	double dwinds[3], mu_th[3], mu_ph[3], Zth, Zph, dnu_th, dnu_ph;
	if(GeoAc_CalcAmp){
		GeoAc_Sources.ddc = c_ddiff(x,y,z,2,2);     GeoAc_Sources.ddu = u_ddiff(x,y,z,2,2);
		GeoAc_Sources.ddv = v_ddiff(x,y,z,2,2);     GeoAc_Sources.ddw = w_ddiff(x,y,z,2,2);

		mu_th[0] = GeoAc_Sources.mu0_xy[0][0];  mu_th[1] = GeoAc_Sources.mu0_xy[1][0]; 	mu_th[2] = current_values[7];   Zth = current_values[6];
		mu_ph[0] = GeoAc_Sources.mu0_xy[0][1];	mu_ph[1] = GeoAc_Sources.mu0_xy[1][1];  mu_ph[2] = current_values[11];  Zph = current_values[10];

		GeoAc_Sources.dnu_mag[0] = (nu[0] * mu_th[0] + nu[1] * mu_th[1] + nu[2] * mu_th[2]) / GeoAc_Sources.nu_mag;
		GeoAc_Sources.dnu_mag[1] = (nu[0] * mu_ph[0] + nu[1] * mu_ph[1] + nu[2] * mu_ph[2]) / GeoAc_Sources.nu_mag;   
		
		dwinds[0] = u_diff(x,y,z,2);
        dwinds[1] = v_diff(x,y,z,2);
        dwinds[2] = w_diff(x,y,z,2);
        for(int n = 0; n < 3; n++){
            GeoAc_Sources.dc_prop[n][0] = nu[n] / GeoAc_Sources.nu_mag * GeoAc_Sources.dc * Zth + GeoAc_Sources.c * mu_th[n] / GeoAc_Sources.nu_mag
                        - GeoAc_Sources.c * nu[n] / pow(GeoAc_Sources.nu_mag,2) * GeoAc_Sources.dnu_mag[0] + dwinds[n] * Zth;
            GeoAc_Sources.dc_prop[n][1] = nu[n] / GeoAc_Sources.nu_mag * GeoAc_Sources.dc * Zph + GeoAc_Sources.c * mu_ph[n] / GeoAc_Sources.nu_mag
                        - GeoAc_Sources.c * nu[n] / pow(GeoAc_Sources.nu_mag,2) * GeoAc_Sources.dnu_mag[1] + dwinds[n] * Zph;
        }

        GeoAc_Sources.dc_prop_mag[0] = (GeoAc_Sources.c_prop[0] * GeoAc_Sources.dc_prop[0][0] + GeoAc_Sources.c_prop[1] * GeoAc_Sources.dc_prop[1][0] + GeoAc_Sources.c_prop[2] * GeoAc_Sources.dc_prop[2][0]) / GeoAc_Sources.c_prop_mag;
        GeoAc_Sources.dc_prop_mag[1] = (GeoAc_Sources.c_prop[0] * GeoAc_Sources.dc_prop[0][1] + GeoAc_Sources.c_prop[1] * GeoAc_Sources.dc_prop[1][1] + GeoAc_Sources.c_prop[2] * GeoAc_Sources.dc_prop[2][1]) / GeoAc_Sources.c_prop_mag;

	}
}

//-----------------------------------------------------------//
//-------Evaluate the Source Equation For Specific Index-----//
//-----------------------------------------------------------//
double GeoAc_EvalSrcEq(double ray_length, double* current_values, int Eq_Number){
	double result;

    double cp_mag = GeoAc_Sources.c_prop_mag;
	double nu[3], mu[3];
    
	switch(Eq_Number){
		case(0): // dx/ds
		case(1): // dy/ds
		case(2): // dz/ds
			result = GeoAc_Sources.c_prop[Eq_Number] / cp_mag;
			break;
            
		case(3): // d nu_z /ds
            nu[0] = GeoAc_Sources.nu0_xy[0];
            nu[1] = GeoAc_Sources.nu0_xy[1];
            nu[2] = current_values[3];
            
            result = -1.0 / cp_mag * (GeoAc_Sources.nu_mag * GeoAc_Sources.dc
                                            + (nu[0] * GeoAc_Sources.du + nu[1] * GeoAc_Sources.dv + nu[2] * GeoAc_Sources.dw));
			break;

		case(4): // dXt/ds
		case(5): // dYt/ds
		case(6): // dZt/ds
			result = GeoAc_Sources.dc_prop[Eq_Number - 4][0] / cp_mag
                        - GeoAc_Sources.c_prop[Eq_Number - 4] / pow(cp_mag,2) * GeoAc_Sources.dc_prop_mag[0];
			break;


		case(7): // d mu_zt /ds
            nu[0] = GeoAc_Sources.nu0_xy[0];    mu[0] = GeoAc_Sources.mu0_xy[0][0];
            nu[1] = GeoAc_Sources.nu0_xy[1];    mu[1] = GeoAc_Sources.mu0_xy[1][0];
            nu[2] = current_values[3];          mu[2] = current_values[7];
    
			result = 1.0 / pow(cp_mag, 2) * (GeoAc_Sources.nu_mag * GeoAc_Sources.dc + (nu[0] * GeoAc_Sources.du + nu[1] * GeoAc_Sources.dv + nu[2] *GeoAc_Sources.dw)) * GeoAc_Sources.dc_prop_mag[0]
						-1.0 / cp_mag * (GeoAc_Sources.dnu_mag[0] * GeoAc_Sources.dc + (mu[0] * GeoAc_Sources.du + mu[1] * GeoAc_Sources.dv + mu[2] * GeoAc_Sources.dw
                                                + (GeoAc_Sources.nu_mag * GeoAc_Sources.ddc + nu[0] * GeoAc_Sources.ddu + nu[1] * GeoAc_Sources.ddv + nu[2] * GeoAc_Sources.ddw) * current_values[6]));
			break;


		case(8):  // dXp/ds
		case(9):  // dYp/ds
		case(10): // dZp/ds
			result =  GeoAc_Sources.dc_prop[Eq_Number - 8][1] / cp_mag
                        - GeoAc_Sources.c_prop[Eq_Number - 8] / pow(cp_mag,2) * GeoAc_Sources.dc_prop_mag[1];
			break;

		case(11): // d mu_zp/ds
            nu[0] = GeoAc_Sources.nu0_xy[0];    mu[0] = GeoAc_Sources.mu0_xy[0][1];
            nu[1] = GeoAc_Sources.nu0_xy[1];    mu[1] = GeoAc_Sources.mu0_xy[1][1];
            nu[2] = current_values[3];          mu[2] = current_values[11];
            
			result = 1.0 / pow(cp_mag, 2) * (GeoAc_Sources.nu_mag * GeoAc_Sources.dc + (nu[0] * GeoAc_Sources.du + nu[1] * GeoAc_Sources.dv + nu[2] *GeoAc_Sources.dw)) * GeoAc_Sources.dc_prop_mag[1]
            -1.0 / cp_mag * (GeoAc_Sources.dnu_mag[1] * GeoAc_Sources.dc + (mu[0] * GeoAc_Sources.du + mu[1] * GeoAc_Sources.dv + mu[2] * GeoAc_Sources.dw
                                                                            + (GeoAc_Sources.nu_mag * GeoAc_Sources.ddc + nu[0] * GeoAc_Sources.ddu + nu[1] * GeoAc_Sources.ddv + nu[2] * GeoAc_Sources.ddw) * current_values[10]));
			break;
	}
	return result;
}

//-------------------------------------------------------------------//
//-------Calculate the Hamiltonian (Eikonal) To Check For Errors-----//
//-------------------------------------------------------------------//
double GeoAc_EvalHamiltonian(double ** solution, int k){
	double x = solution[k][0], y = solution[k][1], z = solution[k][2];
    double nu[3] = {GeoAc_Sources.nu0_xy[0], GeoAc_Sources.nu0_xy[1], solution[k][3]};
    	
    return fabs(sqrt(nu[0] * nu[0] + nu[1] * nu[1] + nu[2] * nu[2])
                - GeoAc_Sources.c0 / c(x,y,z)*(1.0 - (nu[0] * u(x,y,z) + nu[1] * v(x,y,z) + nu[2] * w(x,y,z)) / GeoAc_Sources.c0));
}


//--------------------------------------------------------------------------//
//-------Check if ray has left propagation region or returned to ground-----//
//--------------------------------------------------------------------------//
bool GeoAc_BreakCheck(double ** solution, int index){
	double r 	= sqrt(pow(solution[index][0],2) + pow(solution[index][1],2));
	double z 	= solution[index][2];
	
	bool check = false;	
	if(z > GeoAc_vert_limit)    check = true;
	if(r > GeoAc_range_limit)   check = true;
	
	return check;
}

bool GeoAc_GroundCheck(double ** solution, int index){
	bool check = false;
	if(solution[index][2] < z_grnd) check = true;
	
	return check;
}

//----------------------------------------------------------------------------------//
//-------Calculate the travel time from source to location or between locations-----//
//----------------------------------------------------------------------------------//
double GeoAc_TravelTime(double ** solution, int index){
	double dx, dy, dz, ds, x, y, z, nu[3], nu_mag, c_prop[3], c_prop_mag;
	double traveltime = 0;
	
    nu[0] = GeoAc_Sources.nu0_xy[0];
    nu[1] = GeoAc_Sources.nu0_xy[1];
    
	for (int n = 0; n < index; n++){
		
		dx = solution[n+1][0] - solution[n][0];		// Calculate ds = sqrt(dx^2 + dy^2 + dz^2)
		dy = solution[n+1][1] - solution[n][1];
		dz = solution[n+1][2] - solution[n][2];
		ds = sqrt(dx*dx + dy*dy + dz*dz);

		x = solution[n][0] + dx/2.0;				// Calculate (x,y,z) at the midpoint between [n] and [n+1]
		y = solution[n][1] + dy/2.0;
		z = solution[n][2] + dz/2.0;

		nu[2] = solution[n][3] + (solution[n+1][3] - solution[n][3])/2.0;		// Calculate |c_prop|
		nu_mag = (c(0,0,0) - nu[0]*u(x,y,z) - nu[1]*v(x,y,z))/c(x,y,z);
		c_prop[0] = c(x,y,z)*nu[0]/nu_mag + u(x,y,z); 
		c_prop[1] = c(x,y,z)*nu[1]/nu_mag + v(x,y,z);
		c_prop[2] = c(x,y,z)*nu[2]/nu_mag;
		c_prop_mag = sqrt(pow(c_prop[0],2) + pow(c_prop[1],2) + pow(c_prop[2],2));
	
		traveltime += ds/c_prop_mag;	// Add contribution to the travel time
	}
	return traveltime;	
}

void GeoAc_TravelTimeSegment(double & time, double ** solution, int start, int end){
	double dx, dy, dz, ds, x, y, z, nu[3], nu_mag, c_prop[3], c_prop_mag;
	double traveltime = 0;
	
    nu[0] = GeoAc_Sources.nu0_xy[0];
    nu[1] = GeoAc_Sources.nu0_xy[1];
    
	for (int n = start; n < end; n++){
		
		dx = solution[n+1][0] - solution[n][0];		// Calculate ds = sqrt(dx^2 + dy^2 + dz^2)
		dy = solution[n+1][1] - solution[n][1];
		dz = solution[n+1][2] - solution[n][2];
		ds = sqrt(dx*dx + dy*dy + dz*dz);
        
		x = solution[n][0] + dx/2.0;				// Calculate (x,y,z) at the midpoint between [n] and [n+1]
		y = solution[n][1] + dy/2.0;
		z = solution[n][2] + dz/2.0;
        
		nu[2] = solution[n][3] + (solution[n+1][3] - solution[n][3])/2.0;		// Calculate |c_prop|
		nu_mag = (c(0,0,0) - nu[0]*u(x,y,z) - nu[1]*v(x,y,z))/c(x,y,z);
		c_prop[0] = c(x,y,z)*nu[0]/nu_mag + u(x,y,z);
		c_prop[1] = c(x,y,z)*nu[1]/nu_mag + v(x,y,z);
		c_prop[2] = c(x,y,z)*nu[2]/nu_mag;
		c_prop_mag = sqrt(pow(c_prop[0],2) + pow(c_prop[1],2) + pow(c_prop[2],2));
        
		time += ds/c_prop_mag;	// Add contribution to the travel time
	}
}

//-----------------------------------------------------------------------------------//
//-------Calculate the Jacobian determinant and from it the ampltude coefficient-----//
//-----------------------------------------------------------------------------------//
double GeoAc_Jacobian(double ** solution, int index){
	double x = solution[index][0], y = solution[index][1], z = solution[index][2];
    double x0 = GeoAc_Sources.src_loc[0], y0 = GeoAc_Sources.src_loc[1], z0 = GeoAc_Sources.src_loc[2];
    
	double nu[3] = 		{GeoAc_Sources.nu0_xy[0], GeoAc_Sources.nu0_xy[1], solution[index][3]};
	double nu_mag = 	(c(x0,y0,z0) - nu[0]*u(x,y,z) - nu[1]*v(x,y,z))/c(x,y,z);
	double c_prop[3] = 	{c(x,y,z) * nu[0] / nu_mag + u(x,y,z),
                         c(x,y,z) * nu[1] / nu_mag + v(x,y,z),
                         c(x,y,z) * nu[2] / nu_mag + w(x,y,z)};
	double c_prop_mag = sqrt(pow(c_prop[0],2) + pow(c_prop[1],2) + pow(c_prop[2],2));

	double dxds = c_prop[0]/c_prop_mag,		dyds = c_prop[1]/c_prop_mag,	dzds = c_prop[2]/c_prop_mag;	
	double dxdtheta = solution[index][4],	dydtheta = solution[index][5], 	dzdtheta = solution[index][6];
	double dxdphi = solution[index][8], 	dydphi = solution[index][9], 	dzdphi = solution[index][10];				
		
	return dxds*(dydtheta*dzdphi - dydphi*dzdtheta)
			 - dxdtheta*(dyds*dzdphi - dzds*dydphi)
				 + dxdphi*(dyds*dzdtheta - dzds*dydtheta);
}


double GeoAc_Amplitude(double ** solution, int index){
	double x = solution[index][0], y = solution[index][1], z = solution[index][2];
    double x0 = GeoAc_Sources.src_loc[0], y0 = GeoAc_Sources.src_loc[1], z0 = GeoAc_Sources.src_loc[2];
    
	double nu[3] = 		{GeoAc_Sources.nu0_xy[0], GeoAc_Sources.nu0_xy[1], solution[index][3]};
    double nu_mag = 	(c(x0,y0,z0) - nu[0]*u(x,y,z) - nu[1]*v(x,y,z))/c(x,y,z);
    double nu_mag0 =    1.0 - (nu[0] * u(x0,y0,z0) - nu[1] * v(x0,y0,z0))/c(x0,y0,z0);
	
	double c_prop[3] = 	{c(x,y,z) * nu[0] / nu_mag + u(x,y,z),        c(x,y,z) * nu[1] / nu_mag + v(x,y,z),        c(x,y,z) * nu[2] / nu_mag};
	double c_prop0[3] = {c(x0,y0,z0) * nu[0] / nu_mag0 + u(x0,y0,z0), c(x0,y0,z0) * nu[1] / nu_mag0 + v(x0,y0,z0), c(x0,y0,z0) * sqrt(1.0 - pow(nu[0] / nu_mag0,2) - pow(nu[1] / nu_mag0,2))};

	double  c_prop_mag =  sqrt(pow(c_prop[0],2) + pow(c_prop[1],2) + pow(c_prop[2],2)),
            c_prop_mag0 = sqrt(pow(c_prop0[0],2) + pow(c_prop0[1],2) + pow(c_prop0[2],2));
	double  D = GeoAc_Jacobian(solution, index);
	
	double  Amp_Num = rho(x,y,z) * nu_mag * pow(c(x,y,z),3) * c_prop_mag0 * cos(GeoAc_theta);
	double  Amp_Den = rho(x0,y0,z0) * nu_mag0* pow(c(x0,y0,z0),3) * c_prop_mag * D;
	
	
	return 1.0 / (4.0 * Pi) * sqrt(fabs(Amp_Num / Amp_Den));
}

//--------------------------------------------------------------------------//
//------Integrate the Sutherland Bass Attenuation Through the Ray Path------//
//--------------------------------------------------------------------------//
double GeoAc_SB_Atten(double ** solution, int end, double freq){
	double dx, dy, dz, ds, x, y, z;
    double atten = 0.0;
	for (int n = 0; n < end; n++){
		
		dx = solution[n+1][0] - solution[n][0];		// Calculate ds = sqrt(dx^2 + dy^2 + dz^2)
		dy = solution[n+1][1] - solution[n][1];
		dz = solution[n+1][2] - solution[n][2];
		ds = sqrt(dx*dx + dy*dy + dz*dz);
        
		x = solution[n][0] + dx/2.0;				// Calculate (x,y,z) at the midpoint between [n] and [n+1]
		y = solution[n][1] + dy/2.0;
		z = solution[n][2] + dz/2.0;
        
		atten += SuthBass_Alpha(x, y, z, freq)*ds;	// Add contribution to the attenuation
	}
    return atten;
}

void GeoAc_SB_AttenSegment(double & atten, double ** solution, int start, int end, double freq){
	double dx, dy, dz, ds, x, y, z;
	for (int n = start; n < end; n++){
		
		dx = solution[n+1][0] - solution[n][0];		// Calculate ds = sqrt(dx^2 + dy^2 + dz^2)
		dy = solution[n+1][1] - solution[n][1];
		dz = solution[n+1][2] - solution[n][2];
		ds = sqrt(dx*dx + dy*dy + dz*dz);
        
		x = solution[n][0] + dx/2.0;				// Calculate (x,y,z) at the midpoint between [n] and [n+1]
		y = solution[n][1] + dy/2.0;
		z = solution[n][2] + dz/2.0;
                
		atten += SuthBass_Alpha(x, y, z, freq)*ds;	// Add contribution to the attenuation
	}
}


//--------------------------------------------------------------//
//---------Count the caustics encountered by monitoring---------//
//----how many times the Jacobian determinant changes sign------//
//--------------------------------------------------------------//
int GeoAc_CausticCnt(double ** solution, int start, int index){
	int count = 0;
	double current, prev = GeoAc_Jacobian(solution, 1);
	for(int n = 2; n < index; n++){
		current = GeoAc_Jacobian(solution,n);
		if(current*prev < 0.0){
			count++;
		}
		prev = current;
	}
	return count;	
}

#endif /* GEOAC_EQSETS_3DSTRAT_H_ */
