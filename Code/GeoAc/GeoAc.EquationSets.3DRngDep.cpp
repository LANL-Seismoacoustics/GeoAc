#ifndef GEOAC_EQSETS_3DRNGDEP_H_
#define GEOAC_EQSETS_3DRNGDEP_H_

#include <math.h>
#include <iostream>

#include "GeoAc.Parameters.h"
#include "Atmo_State.h"
#include "G2S_MultiDimSpline3D.h"

using namespace std;

//--------------------------------------------------//
//-------Propagate in 3D, assume non-stratified-----//
//--------------------------------------------------//
void GeoAc_SetSystem(){
    GeoAc_dim = 3;				// Dimensions
    GeoAc_AtmoStrat = false;     // Is the medium stratified?
}

//-----------------------------------------------------------//
//-------Structure containing source functions which are-----//
//-------called multiple times in the solver-----------------//
//-----------------------------------------------------------//
struct GeoAc_Sources_Struct{
    
    double src_loc[3];  // Location of the source (x0, y0, z0)
	double c0;			// Thermodynamic sound speed at the source
    
	double c;			// Thermodynamic sound speed
	double dc[5];		// Derivatives with respect to x, y, z, launch inclination, and launch azimuth
	double ddc[3][2];	// Derivatives with respect to x, y, and z; then with respect to either inclination or azimuth
    
	double u;			// E-W wind component
	double du[5];		// Derivatives with respect to x, y, z, launch inclination, and launch azimuth
	double ddu[3][2];	// Derivatives with respect to x, y, and z; then with respect to either inclination or azimuth
    
	double v;			// N-S wind component
	double dv[5];		// Derivatives with respect to x, y, z, launch inclination, and launch azimuth
	double ddv[3][2];	// Derivatives with respect to x, y, and z; then with respect to either inclination or azimuth
    
	double w;			// Vertical wind component
	double dw[5];		// Derivatives with respect to x, y, z, launch inclination, and launch azimuth
	double ddw[3][2];	// Derivatives with respect to x, y, and z; then with respect to either inclination or azimuth
    
    double nu0;         // Eikonal vector magnitude at the source
    double nu_mag;		// Eikonal vector magnitude
	double dnu_mag[2];	// Derivative of eikonal vector magnitude with respect to inclination and azimuth angles
    
	double c_gr[3];         // Group velocity cg = c {nu_x,nu_y,nu_z}/|nu} + {u,v,w}
	double c_gr_mag;        // Group velocity magnitude |cp| = c sqrt(1 + 2 nu/|nu| dot v/c + v^2/c^2)
    double dc_gr[3][2];     // Derivatives of c_gr components with respect to lt and lp
	double dc_gr_mag[2];	// Derivatives of c_gr magnitude with respect to lt and lp
    
};

struct GeoAc_Sources_Struct GeoAc_Sources = {
    {0.0, 0.0, 0.0}, 0.0,
    0.0, {0.0, 0.0, 0.0, 0.0, 0.0}, {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}},
    0.0, {0.0, 0.0, 0.0, 0.0, 0.0}, {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}},
    0.0, {0.0, 0.0, 0.0, 0.0, 0.0}, {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}},
    0.0, {0.0, 0.0, 0.0, 0.0, 0.0}, {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}},
    0.0, 0.0, {0.0, 0.0},
    {0.0, 0.0, 0.0}, 0.0, {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}, {0.0, 0.0}
};

//----------------------------------------------------------------------//
//-------Fill in solution[0][n] with the appropriate initial values-----//
//----------------------------------------------------------------------//
void GeoAc_SetInitialConditions(double ** & solution, double x0, double y0, double z0){
    GeoAc_Sources.src_loc[0] = x0;
    GeoAc_Sources.src_loc[1] = y0;
    GeoAc_Sources.src_loc[2] = z0;
	GeoAc_Sources.c0 = c(x0, y0, z0);
    
    double MachComps[3] =   { u(x0, y0, z0)/GeoAc_Sources.c0,   v(x0, y0, z0)/GeoAc_Sources.c0,     w(x0, y0, z0)/GeoAc_Sources.c0};
    double nu0[3] =         { cos(GeoAc_theta)*cos(GeoAc_phi),  cos(GeoAc_theta)*sin(GeoAc_phi),    sin(GeoAc_theta)};
    double mu0_th[3] =      {-sin(GeoAc_theta)*cos(GeoAc_phi), -sin(GeoAc_theta)*sin(GeoAc_phi),    cos(GeoAc_theta)};
    double mu0_ph[3] =      {-cos(GeoAc_theta)*sin(GeoAc_phi),  cos(GeoAc_theta)*cos(GeoAc_phi),    0.0};
    
    double MachScalar = 1.0 + (nu0[0]*MachComps[0] + nu0[1]*MachComps[1] + nu0[2]*MachComps[2]);
    GeoAc_Sources.nu0 = 1.0/MachScalar;
    
	for(int index = 0; index < GeoAc_EqCnt; index++){
		switch(index){
			case(0):			// x(0) = x0
                solution[0][index] =  x0;
				break;
                
            case(1):			// y(0) = y0
				solution[0][index] =  y0;
				break;
                
            case(2):			// z(0) = z0
                solution[0][index] =  z0;
				break;
                
			case(3):			// nu_x(0)
			case(4):			// nu_y(0)
			case(5):            // nu_z(0)
				solution[0][index] = nu0[index - 3]/MachScalar;
				break;
                
            case(6):			// Xt(0) = 0
			case(7):			// Yt(0) = 0
			case(8):			// Zt(0) = 0
			case(12):			// Xp(0) = 0
			case(13):			// Yp(0) = 0
			case(14):			// Zp(0) = 0
				solution[0][index] =  0.0;
				break;
                
			case(9):			// mu_x_th(0) = d/d lt (nu_x)
			case(10):			// mu_y_th(0) = d/d lt (nu_y)
			case(11):			// mu_z_th(0) = d/d lt (nu_z)
                solution[0][index] =  mu0_th[index - 9]/MachScalar - nu0[index - 9]/pow(MachScalar,2.0) * (mu0_th[0]*MachComps[0] + mu0_th[1]*MachComps[1] + mu0_th[2]*MachComps[2]);
				break;
                
			case(15):			// mu_x_ph(0) = d/d lp (nu_x)
            case(16):			// mu_y_ph(0) = d/d lp (nu_y)
            case(17):			// mu_z_ph(0) = d/d lp (nu_z)
                solution[0][index] =  mu0_ph[index - 15]/MachScalar - nu0[index - 15]/pow(MachScalar,2.0) * (mu0_ph[0]*MachComps[0] + mu0_ph[1]*MachComps[1] + mu0_ph[2]*MachComps[2]);
                break;
                
			default:
                cout << "Unexpected index in Initial_Cond.  Model includes 18 variables." << '\n';
		}
	}

    for(int n = 0; n < 3; n++){
        Temp_Spline.accel[n] =  0;
        Windu_Spline.accel[n] = 0;
        Windv_Spline.accel[n] = 0;
    }
    
}


//--------------------------------------------------------------------------//
//-------Taylor series fit to more accurately deterine intercept values-----//
//--------------------------------------------------------------------------//
void GeoAc_ApproximateIntercept(double ** solution, int k, double* & prev){
	double result;
	double dz_k = solution[k][2] - solution[k-1][2];    // set dz for step = z_k - z_{k-1}
	double dz_grnd = solution[k-1][2] - z_grnd;         // set dz from z_{k-1} to ground (assuming z = 0 is ground)
    
	for(int index = 0; index < GeoAc_EqCnt; index++){
		prev[index] = solution[k-1][index] 	+ (solution[k-1][index] - solution[k][index])/dz_k*dz_grnd
                            + 1.0/2.0*(solution[k][index] + solution[k-2][index] - 2.0*solution[k-1][index])/pow(dz_k,2.0)*pow(dz_grnd,2.0);
		
	}
}

//-------------------------------------------------------------------------//
//-------Fill in solution[0][n] with the appropriate reflection values-----//
//-------------------------------------------------------------------------//
void GeoAc_SetReflectionConditions(double** & solution, int k_end){
	double* prev = new double [GeoAc_EqCnt];
	GeoAc_ApproximateIntercept(solution, k_end, prev);

    double c_grnd = c(prev[0], prev[1], z_grnd);
	double dnuz_ds = - 1.0/c_grnd * (GeoAc_Sources.c0/c_grnd * c_diff(prev[0], prev[1],z_grnd, 2)
                                                   + prev[3] * u_diff(prev[0], prev[1],z_grnd, 2)
                                                        + prev[4] * v_diff(prev[0],prev[1],z_grnd, 2)
                                                            + prev[5] * w_diff(prev[0],prev[1],z_grnd, 2));
	
	for(int index = 0; index < GeoAc_EqCnt; index++){
		switch(index){
			case(0):			// x(0) = 		x(s0)
			case(1):			// y(0) = 		y(s0)
			case(3):			// nux(0) = 	nux(s0)
			case(4):			// nuy(0) = 	nuy(s0)
			case(6):			// Xt(0) = 		Xt(s0)
			case(7):			// Yt(0) = 		Yt(s0)
			case(9):			// mu_xt(0) = 	mu_xt(s0)
            case(10):			// mu_yt(0) = 	mu_yt(s0)
			case(12):			// Xp(0) = 		Xp(s0)
            case(13):			// Yp(0) = 		Yp(s0)
			case(15):			// mu_xp(0) = 	mu_xp(s0)
            case(16):			// mu_yp(0) = 	mu_yp(s0)
				solution[0][index] = prev[index];
				break;
                
			case(2):			// z(0) = 		0.0 = +/-z(s0)
				solution[0][index] = z_grnd;
				break;
                
			case(5):			// nu_z = 		-nu_z(s0)
			case(8):			// Zt(0) = 		-Zt(s0)
			case(14):			// Zp(0) = 		-Zp(s0)
				solution[0][index] = -prev[index];
				break;
                
			case(11):			// mu_zt(0) =  	-mu_zt(s0) + 2 d nu_z/ ds * ds0/dtheta
 			case(17):			// mu_zp(0) = 	-mu_zp(s0) + 2 d nu_z/ ds * ds0/dphi
				solution[0][index] =  -prev[index] + 2.0 * dnuz_ds*prev[index - 3]/(c_grnd/GeoAc_Sources.c0 * prev[5]);
				break;
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
    // Extract ray location and eikonal vector components
    double x = current_values[0],		y = current_values[1], 	z = current_values[2];
	double nu[3] = {current_values[3], 	current_values[4], 		current_values[5]};
    
    double temp;
    double dtemp[3];
    
    if(!GeoAc_CalcAmp){
        Eval_Spline_AllOrder1(x, y, z, Temp_Spline, temp, dtemp[0], dtemp[1], dtemp[2]);
        for(int n = 0; n < 3; n++){
            Windu_Spline.accel[n] = Temp_Spline.accel[n];
            Windv_Spline.accel[n] = Temp_Spline.accel[n];
        }
        
        Eval_Spline_AllOrder1(x, y, z, Windu_Spline, GeoAc_Sources.u, GeoAc_Sources.du[0], GeoAc_Sources.du[1], GeoAc_Sources.du[2]);
        Eval_Spline_AllOrder1(x, y, z, Windv_Spline, GeoAc_Sources.v, GeoAc_Sources.dv[0], GeoAc_Sources.dv[1], GeoAc_Sources.dv[2]);
        GeoAc_Sources.w = w(x,y,z);
 
        GeoAc_Sources.c = sqrt(gamR * temp);
        for(int n = 0; n < 3; n++){
            GeoAc_Sources.dc[n] = gamR / (2.0 * GeoAc_Sources.c) * dtemp[n];
            GeoAc_Sources.dw[n] = w_diff(x,y,z,n);
        }
        
        GeoAc_Sources.nu_mag = 	 sqrt(nu[0]*nu[0] + nu[1]*nu[1] + nu[2]*nu[2]);
        
        GeoAc_Sources.c_gr[0] =  GeoAc_Sources.c*nu[0]/GeoAc_Sources.nu_mag + GeoAc_Sources.u;
        GeoAc_Sources.c_gr[1] =  GeoAc_Sources.c*nu[1]/GeoAc_Sources.nu_mag + GeoAc_Sources.v;
        GeoAc_Sources.c_gr[2] =  GeoAc_Sources.c*nu[2]/GeoAc_Sources.nu_mag + GeoAc_Sources.w;
        
        GeoAc_Sources.c_gr_mag = sqrt(pow(GeoAc_Sources.c_gr[0],2) + pow(GeoAc_Sources.c_gr[1],2) + pow(GeoAc_Sources.c_gr[2],2));
    }

    double X_th[3], X_ph[3], mu_th[3], mu_ph[3];
    double ddtemp[3][3];
    double ddWindu[3][3];
    double ddWindv[3][3];
    
    if(GeoAc_CalcAmp){
        X_th[0]  = current_values[6];       X_th[1]  = current_values[7];       X_th[2]  = current_values[8];
		mu_th[0] = current_values[9];		mu_th[1] = current_values[10];		mu_th[2] = current_values[11];
        X_ph[0]  = current_values[12];      X_ph[1]  = current_values[13];      X_ph[2]  = current_values[14];
		mu_ph[0] = current_values[15];		mu_ph[1] = current_values[16];		mu_ph[2] = current_values[17];
        
        Eval_Spline_AllOrder2(x, y, z, Temp_Spline, temp, dtemp[0], dtemp[1], dtemp[2], ddtemp[0][0], ddtemp[1][1], ddtemp[2][2], ddtemp[0][1], ddtemp[0][2], ddtemp[1][2]);
        for(int n = 0; n < 3; n++){
            Windu_Spline.accel[n] = Temp_Spline.accel[n];
            Windv_Spline.accel[n] = Temp_Spline.accel[n];
        }

        Eval_Spline_AllOrder2(x, y, z, Windu_Spline, GeoAc_Sources.u, GeoAc_Sources.du[0], GeoAc_Sources.du[1], GeoAc_Sources.du[2], ddWindu[0][0], ddWindu[1][1], ddWindu[2][2], ddWindu[0][1], ddWindu[0][2], ddWindu[1][2]);
        Eval_Spline_AllOrder2(x, y, z, Windv_Spline, GeoAc_Sources.v, GeoAc_Sources.dv[0], GeoAc_Sources.dv[1], GeoAc_Sources.dv[2], ddWindv[0][0], ddWindv[1][1], ddWindv[2][2], ddWindv[0][1], ddWindv[0][2], ddWindv[1][2]);
        GeoAc_Sources.w = w(x,y,z);

        ddtemp[1][0] = ddtemp[0][1];    ddtemp[2][0] = ddtemp[0][2];    ddtemp[2][1] = ddtemp[1][2];
        ddWindu[1][0] = ddWindu[0][1];  ddWindu[2][0] = ddWindu[0][2];  ddWindu[2][1] = ddWindu[1][2];
        ddWindv[1][0] = ddWindv[0][1];  ddWindv[2][0] = ddWindv[0][2];  ddWindv[2][1] = ddWindv[1][2];
        
        GeoAc_Sources.c = sqrt(gamR * temp);
        for(int n = 0; n < 3; n++){
            GeoAc_Sources.dc[n] = gamR / (2.0 * GeoAc_Sources.c) * dtemp[n];
            GeoAc_Sources.dw[n] = w_diff(x,y,z,n);
            
            GeoAc_Sources.ddc[n][0] = 0.0;            GeoAc_Sources.ddc[n][1] = 0.0;
            GeoAc_Sources.ddu[n][0] = 0.0;            GeoAc_Sources.ddu[n][1] = 0.0;
            GeoAc_Sources.ddv[n][0] = 0.0;            GeoAc_Sources.ddv[n][1] = 0.0;
            GeoAc_Sources.ddw[n][0] = 0.0;            GeoAc_Sources.ddw[n][1] = 0.0;
            for(int m = 0; m < 3; m++){
                GeoAc_Sources.ddc[n][0] += X_th[m]*(gamR/(2.0*GeoAc_Sources.c) * ddtemp[n][m] - pow(gamR,2)/(4.0 * pow(GeoAc_Sources.c,3)) * dtemp[n]*dtemp[m]);
                GeoAc_Sources.ddc[n][1] += X_ph[m]*(gamR/(2.0*GeoAc_Sources.c) * ddtemp[n][m] - pow(gamR,2)/(4.0 * pow(GeoAc_Sources.c,3)) * dtemp[n]*dtemp[m]);
                
                GeoAc_Sources.ddu[n][0] += X_th[m]*ddWindu[n][m];       GeoAc_Sources.ddu[n][1] += X_ph[m]*ddWindu[n][m];
                GeoAc_Sources.ddv[n][0] += X_th[m]*ddWindv[n][m];       GeoAc_Sources.ddv[n][1] += X_ph[m]*ddWindv[n][m];
                GeoAc_Sources.ddw[n][0] += X_th[m]*w_ddiff(x,y,z,n,m);  GeoAc_Sources.ddw[n][1] += X_ph[m]*w_ddiff(x,y,z,n,m);
            }
        }

        GeoAc_Sources.dc[3] = 0.0;      GeoAc_Sources.dc[4] = 0.0;
        GeoAc_Sources.du[3] = 0.0;      GeoAc_Sources.du[4] = 0.0;
        GeoAc_Sources.dv[3] = 0.0;      GeoAc_Sources.dv[4] = 0.0;
        GeoAc_Sources.dw[3] = 0.0;      GeoAc_Sources.dw[4] = 0.0;
        for(int n = 0; n < 3; n++){
            GeoAc_Sources.dc[3] += X_th[n]*GeoAc_Sources.dc[n];     GeoAc_Sources.dc[4] += X_ph[n]*GeoAc_Sources.dc[n];
            GeoAc_Sources.du[3] += X_th[n]*GeoAc_Sources.du[n];     GeoAc_Sources.du[4] += X_ph[n]*GeoAc_Sources.du[n];
            GeoAc_Sources.dv[3] += X_th[n]*GeoAc_Sources.dv[n];     GeoAc_Sources.dv[4] += X_ph[n]*GeoAc_Sources.dv[n];
            GeoAc_Sources.dw[3] += X_th[n]*GeoAc_Sources.dw[n];     GeoAc_Sources.dw[4] += X_ph[n]*GeoAc_Sources.dw[n];
        }
        
        GeoAc_Sources.nu_mag = 	 sqrt(nu[0]*nu[0] + nu[1]*nu[1] + nu[2]*nu[2]);
		GeoAc_Sources.dnu_mag[0] = (nu[0]*mu_th[0] + nu[1]*mu_th[1] + nu[2]*mu_th[2])/GeoAc_Sources.nu_mag;
		GeoAc_Sources.dnu_mag[1] = (nu[0]*mu_ph[0] + nu[1]*mu_ph[1] + nu[2]*mu_ph[2])/GeoAc_Sources.nu_mag;
        
        GeoAc_Sources.c_gr[0] =  GeoAc_Sources.c*nu[0]/GeoAc_Sources.nu_mag + GeoAc_Sources.u;
        GeoAc_Sources.c_gr[1] =  GeoAc_Sources.c*nu[1]/GeoAc_Sources.nu_mag + GeoAc_Sources.v;
        GeoAc_Sources.c_gr[2] =  GeoAc_Sources.c*nu[2]/GeoAc_Sources.nu_mag + GeoAc_Sources.w;
        GeoAc_Sources.c_gr_mag = sqrt(pow(GeoAc_Sources.c_gr[0],2) + pow(GeoAc_Sources.c_gr[1],2) + pow(GeoAc_Sources.c_gr[2],2));
        
        GeoAc_Sources.dc_gr[0][0] = nu[0]/GeoAc_Sources.nu_mag*GeoAc_Sources.dc[3] + GeoAc_Sources.c*mu_th[0]/GeoAc_Sources.nu_mag - GeoAc_Sources.c*nu[0]/pow(GeoAc_Sources.nu_mag,2) * GeoAc_Sources.dnu_mag[0] + GeoAc_Sources.du[3];
        GeoAc_Sources.dc_gr[1][0] = nu[1]/GeoAc_Sources.nu_mag*GeoAc_Sources.dc[3] + GeoAc_Sources.c*mu_th[1]/GeoAc_Sources.nu_mag - GeoAc_Sources.c*nu[1]/pow(GeoAc_Sources.nu_mag,2) * GeoAc_Sources.dnu_mag[0] + GeoAc_Sources.dv[3];
        GeoAc_Sources.dc_gr[2][0] = nu[2]/GeoAc_Sources.nu_mag*GeoAc_Sources.dc[3] + GeoAc_Sources.c*mu_th[2]/GeoAc_Sources.nu_mag - GeoAc_Sources.c*nu[2]/pow(GeoAc_Sources.nu_mag,2) * GeoAc_Sources.dnu_mag[0] + GeoAc_Sources.dw[3];
        GeoAc_Sources.dc_gr_mag[0] = (GeoAc_Sources.c_gr[0]*GeoAc_Sources.dc_gr[0][0] + GeoAc_Sources.c_gr[1]*GeoAc_Sources.dc_gr[1][0] + GeoAc_Sources.c_gr[2]*GeoAc_Sources.dc_gr[2][0])/GeoAc_Sources.c_gr_mag;
        
        GeoAc_Sources.dc_gr[0][1] = nu[0]/GeoAc_Sources.nu_mag*GeoAc_Sources.dc[4] + GeoAc_Sources.c*mu_ph[0]/GeoAc_Sources.nu_mag - GeoAc_Sources.c*nu[0]/pow(GeoAc_Sources.nu_mag,2) * GeoAc_Sources.dnu_mag[1] + GeoAc_Sources.du[4];
        GeoAc_Sources.dc_gr[1][1] = nu[1]/GeoAc_Sources.nu_mag*GeoAc_Sources.dc[4] + GeoAc_Sources.c*mu_ph[1]/GeoAc_Sources.nu_mag - GeoAc_Sources.c*nu[1]/pow(GeoAc_Sources.nu_mag,2) * GeoAc_Sources.dnu_mag[1] + GeoAc_Sources.dv[4];
        GeoAc_Sources.dc_gr[2][1] = nu[2]/GeoAc_Sources.nu_mag*GeoAc_Sources.dc[4] + GeoAc_Sources.c*mu_ph[2]/GeoAc_Sources.nu_mag - GeoAc_Sources.c*nu[2]/pow(GeoAc_Sources.nu_mag,2) * GeoAc_Sources.dnu_mag[1] + GeoAc_Sources.dw[4];        
        GeoAc_Sources.dc_gr_mag[1] = (GeoAc_Sources.c_gr[0]*GeoAc_Sources.dc_gr[0][1] + GeoAc_Sources.c_gr[1]*GeoAc_Sources.dc_gr[1][1] + GeoAc_Sources.c_gr[2]*GeoAc_Sources.dc_gr[2][1])/GeoAc_Sources.c_gr_mag;
	}
}

//-----------------------------------------------------------//
//-------Evaluate the Source Equation For Specific Index-----//
//-----------------------------------------------------------//
double GeoAc_EvalSrcEq(double ray_length, double* current_values, int Eq_Number){
	double result;
    
	// Set variables used in all equations
    double nu[3] =      {current_values[3], 	current_values[4], 		current_values[5]};
    double mu_th[3] =   {current_values[9], 	current_values[10],		current_values[11]};
    double mu_ph[3] =   {current_values[15], 	current_values[16],		current_values[17]};
    
	switch(Eq_Number){
		case(0):	// d x / d s
        case(1):    // d y / d s
        case(2):    // d z / d s
			result = GeoAc_Sources.c_gr[Eq_Number]/GeoAc_Sources.c_gr_mag;
			break;
            
		case(3): 	// d nu_x /ds
		case(4): 	// d nu_y / ds
		case(5): 	// d nu_z / ds
			result = -1.0/GeoAc_Sources.c_gr_mag*(GeoAc_Sources.nu_mag*GeoAc_Sources.dc[Eq_Number-3]
                                                        + nu[0]*GeoAc_Sources.du[Eq_Number-3] + nu[1]*GeoAc_Sources.dv[Eq_Number-3] + nu[2]*GeoAc_Sources.dw[Eq_Number-3]);
			break;
            
            
		case(6):	// d X_th/ds
        case(7):	// d Y_th/ds
        case(8): 	// d Z_th/ds
			result = GeoAc_Sources.dc_gr[Eq_Number-6][0]/GeoAc_Sources.c_gr_mag
                        - GeoAc_Sources.c_gr[Eq_Number-6]/pow(GeoAc_Sources.c_gr_mag,2) * GeoAc_Sources.dc_gr_mag[0];
			break;
            
            
		case(9): 	// d mu_x_th /ds
		case(10): 	// d mu_y_th /ds
		case(11): 	// d mu_z_th /ds
			result = 1.0/pow(GeoAc_Sources.c_gr_mag,2) * GeoAc_Sources.dc_gr_mag[0]*(GeoAc_Sources.nu_mag*GeoAc_Sources.dc[Eq_Number-9]
                                                                                        + nu[0]*GeoAc_Sources.du[Eq_Number-9] + nu[1]*GeoAc_Sources.dv[Eq_Number-9] + nu[2]*GeoAc_Sources.dw[Eq_Number-9])
                            - 1.0/GeoAc_Sources.c_gr_mag*(GeoAc_Sources.dnu_mag[0]*GeoAc_Sources.dc[Eq_Number-9] + GeoAc_Sources.nu_mag*GeoAc_Sources.ddc[Eq_Number-9][0]
                                                            + mu_th[0]*GeoAc_Sources.du[Eq_Number-9] + mu_th[1]*GeoAc_Sources.dv[Eq_Number-9] + mu_th[2]*GeoAc_Sources.dw[Eq_Number-9]
                                                                    + nu[0]*GeoAc_Sources.ddu[Eq_Number-9][0] + nu[1]*GeoAc_Sources.ddv[Eq_Number-9][0] + nu[2]*GeoAc_Sources.ddw[Eq_Number-9][0]);
            
			break;
            
            
		case(12):  	// dX_ph/ds
		case(13):  	// dY_ph/ds
        case(14):	// dZ_ph/ds
			result = GeoAc_Sources.dc_gr[Eq_Number-12][1]/GeoAc_Sources.c_gr_mag
                        - GeoAc_Sources.c_gr[Eq_Number-12]/pow(GeoAc_Sources.c_gr_mag,2) * GeoAc_Sources.dc_gr_mag[1];
			break;
            
		case(15): 	// d mu_x_ph/ds
		case(16): 	// d mu_y_ph/ds
		case(17): 	// d mu_z_ph/ds
			result = 1.0/pow(GeoAc_Sources.c_gr_mag,2) * GeoAc_Sources.dc_gr_mag[1]*(GeoAc_Sources.nu_mag*GeoAc_Sources.dc[Eq_Number-15]
                                                                                            + nu[0]*GeoAc_Sources.du[Eq_Number-15] + nu[1]*GeoAc_Sources.dv[Eq_Number-15] + nu[2]*GeoAc_Sources.dw[Eq_Number-15])
                            - 1.0/GeoAc_Sources.c_gr_mag*(GeoAc_Sources.dnu_mag[1]*GeoAc_Sources.dc[Eq_Number-15] + GeoAc_Sources.nu_mag*GeoAc_Sources.ddc[Eq_Number-15][1]
                                                            + mu_ph[0]*GeoAc_Sources.du[Eq_Number-15] + mu_ph[1]*GeoAc_Sources.dv[Eq_Number-15] + mu_ph[2]*GeoAc_Sources.dw[Eq_Number-15]
                                                                + nu[0]*GeoAc_Sources.ddu[Eq_Number-15][1] + nu[1]*GeoAc_Sources.ddv[Eq_Number-15][1] + nu[2]*GeoAc_Sources.ddw[Eq_Number-15][1]);
			break;
            
	}
	return result;
}

//-------------------------------------------------------------------//
//-------Calculate the Hamiltonian (Eikonal) To Check For Errors-----//
//-------------------------------------------------------------------//
double GeoAc_EvalHamiltonian(double ray_length, double* current_values, double c0){
	double x = current_values[0],  y = current_values[1], z = current_values[2];
	double nu[3] = {current_values[3], current_values[4], current_values[5]};
	
	return sqrt(nu[0]*nu[0] + nu[1]*nu[1] + nu[2]*nu[2]) - c0/c(x,y,z)
                    + (u(x,y,z)*nu[0] + v(x,y,z)*nu[1] + w(x,y,z)*nu[2])/c(x,y,z);
}

double GeoAc_EvalHamiltonian(double** solution, int index){
	double x = solution[index][0],          y = solution[index][1],         z = solution[index][2];
	double x0 = GeoAc_Sources.src_loc[0], 	y0 = GeoAc_Sources.src_loc[1], 	z0 = GeoAc_Sources.src_loc[2];
	double nu[3] = {solution[index][3],     solution[index][4],             solution[index][5]};
    double c0 = GeoAc_Sources.c0;
    
	return sqrt(nu[0]*nu[0] + nu[1]*nu[1] + nu[2]*nu[2]) - c0/c(x,y,z)*(1.0 - (u(x,y,z)*nu[0] + v(x,y,z)*nu[1] + w(x,y,z)*nu[2])/c0);
}

double GeoAc_EvalHamiltonian_Deriv(double** solution, int index){
	double  x =           solution[index][0],y = solution[index][1], z = solution[index][2],
            X_th[3] =    {solution[index][6],    solution[index][7],     solution[index][8]},
            X_ph[3] =    {solution[index][12],   solution[index][13],    solution[index][14]},
            nu[3] =      {solution[index][3],    solution[index][4],     solution[index][5]},
            mu_th[3] =   {solution[index][9],    solution[index][10],    solution[index][11]},
            mu_ph[3] =   {solution[index][15],   solution[index][16],    solution[index][17]};
    
    double mag_nu = sqrt(nu[0]*nu[0]+nu[1]*nu[1]+nu[2]*nu[2]);
    
    double dc_dlt = X_th[0]*c_diff(x,y,z,0) + X_th[1]*c_diff(x,y,z,1) + X_th[2]*c_diff(x,y,z,2);
    double dc_dlp = X_ph[0]*c_diff(x,y,z,0) + X_ph[1]*c_diff(x,y,z,1) + X_ph[2]*c_diff(x,y,z,2);
    
    double du_dlt = X_th[0]*u_diff(x,y,z,0) + X_th[1]*u_diff(x,y,z,1) + X_th[2]*u_diff(x,y,z,2);
    double du_dlp = X_ph[0]*u_diff(x,y,z,0) + X_ph[1]*u_diff(x,y,z,1) + X_ph[2]*u_diff(x,y,z,2);
    
    double dv_dlt = X_th[0]*v_diff(x,y,z,0) + X_th[1]*v_diff(x,y,z,1) + X_th[2]*v_diff(x,y,z,2);
    double dv_dlp = X_ph[0]*v_diff(x,y,z,0) + X_ph[1]*v_diff(x,y,z,1) + X_ph[2]*v_diff(x,y,z,2);
    
    double dw_dlt = X_th[0]*w_diff(x,y,z,0) + X_th[1]*w_diff(x,y,z,1) + X_th[2]*w_diff(x,y,z,2);
    double dw_dlp = X_ph[0]*w_diff(x,y,z,0) + X_ph[1]*w_diff(x,y,z,1) + X_ph[2]*w_diff(x,y,z,2);
    
    
    double Resid_th = (nu[0]*mu_th[0] + nu[1]*mu_th[1] + nu[2]*mu_th[2])/mag_nu + mag_nu/c(x,y,z) * dc_dlt
                            + 1.0/c(x,y,z) * (mu_th[0]*u(x,y,z) + mu_th[1]*v(x,y,z) + mu_th[2]*w(x,y,z) + nu[0]*du_dlt + nu[1]*dv_dlt + nu[2]*dw_dlt);
    double Resid_ph = (nu[0]*mu_ph[0] + nu[1]*mu_ph[1] + nu[2]*mu_ph[2])/mag_nu + mag_nu/c(x,y,z) * dc_dlp
                            + 1.0/c(x,y,z) * (mu_ph[0]*u(x,y,z) + mu_ph[1]*v(x,y,z) + mu_ph[2]*w(x,y,z) + nu[0]*du_dlp + nu[1]*dv_dlp + nu[2]*dw_dlp);
    
    return sqrt(Resid_th*Resid_th + Resid_ph*Resid_ph);    
}


//---------------------------------------------------------//
//--------Check if the ray has left the propagation--------//
//-------------region or returned to the ground------------//
//---------------------------------------------------------//
bool GeoAc_BreakCheck(double ** solution, int index){
	double x 	= solution[index][0];
    double y    = solution[index][1];
    double z 	= solution[index][2];
    
	bool check = false;
    
	if(x > GeoAc_x_max_limit) check = true;
    if(x < GeoAc_x_min_limit) check = true;
    if(y > GeoAc_y_max_limit) check = true;
    if(y < GeoAc_y_min_limit) check = true;
    if(z > GeoAc_vert_limit) check = true;
	
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
    double SndSpd;
	double traveltime = 0;
	
	for (int n = 0; n < index; n++){
		
		dx = solution[n+1][0] - solution[n][0];		// Calculate ds = sqrt(dx^2 + dy^2 + dz^2)
		dy = solution[n+1][1] - solution[n][1];
		dz = solution[n+1][2] - solution[n][2];
		ds = sqrt(dx*dx + dy*dy + dz*dz);
        
		x = solution[n][0] + dx/2.0;				// Calculate (x,y,z) at the midpoint between [n] and [n+1]
		y = solution[n][1] + dy/2.0;
		z = solution[n][2] + dz/2.0;
        
        // Calculate |c_prop|
		nu[0] = solution[n][3] + (solution[n+1][3] - solution[n][3])/2.0;
		nu[1] = solution[n][4] + (solution[n+1][4] - solution[n][4])/2.0;
		nu[2] = solution[n][5] + (solution[n+1][5] - solution[n][5])/2.0;
		nu_mag = sqrt(nu[0]*nu[0] + nu[1]*nu[1] + nu[2]*nu[2]);
        
        SndSpd = c(x,y,z);
		c_prop[0] = SndSpd*nu[0]/nu_mag + u(x,y,z);
		c_prop[1] = SndSpd*nu[1]/nu_mag + v(x,y,z);
		c_prop[2] = SndSpd*nu[2]/nu_mag + w(x,y,z);
		c_prop_mag = sqrt(pow(c_prop[0],2) + pow(c_prop[1],2) + pow(c_prop[2],2));
        
		traveltime += ds/c_prop_mag;	// Add contribution to the travel time
	}
	return traveltime;
}


void GeoAc_TravelTimeSegment(double & time, double ** solution, int start, int end){
	double dx, dy, dz, ds, x, y, z, nu[3], nu_mag, c_prop[3], c_prop_mag;
    double SndSpd;
    
	for (int n = start; n < end; n++){
        
		dx = solution[n+1][0] - solution[n][0];		// Calculate ds = sqrt(dx^2 + dy^2 + dz^2)
		dy = solution[n+1][1] - solution[n][1];
		dz = solution[n+1][2] - solution[n][2];
		ds = sqrt(dx*dx + dy*dy + dz*dz);
        
		x = solution[n][0] + dx/2.0;				// Calculate (x,y,z) at the midpoint between [n] and [n+1]
		y = solution[n][1] + dy/2.0;
		z = solution[n][2] + dz/2.0;
        
        
		// Calculate |c_prop|
		nu[0] = solution[n][3] + (solution[n+1][3] - solution[n][3])/2.0;
		nu[1] = solution[n][4] + (solution[n+1][4] - solution[n][4])/2.0;
		nu[2] = solution[n][5] + (solution[n+1][5] - solution[n][5])/2.0;
		nu_mag = sqrt(nu[0]*nu[0] + nu[1]*nu[1] + nu[2]*nu[2]);
        
        SndSpd = c(x,y,z);
		c_prop[0] = SndSpd*nu[0]/nu_mag + u(x,y,z);
		c_prop[1] = SndSpd*nu[1]/nu_mag + v(x,y,z);
		c_prop[2] = SndSpd*nu[2]/nu_mag + w(x,y,z);
		c_prop_mag = sqrt(pow(c_prop[0],2) + pow(c_prop[1],2) + pow(c_prop[2],2));
        
		time += ds/c_prop_mag;	// Add contribution to the travel time
	}
}

//-----------------------------------------------------------------------------------//
//-------Calculate the Jacobian determinant and from it the ampltude coefficient-----//
//-----------------------------------------------------------------------------------//
double GeoAc_Jacobian(double ** solution, int index){
	double x = solution[index][0], y = solution[index][1], z = solution[index][2];
	double nu[3] = 	{solution[index][3], solution[index][4], solution[index][5]};
	double nu_mag = sqrt(nu[0]*nu[0] + nu[1]*nu[1] + nu[2]*nu[2]);
    double SndSpd = c(x,y,z);

	double c_prop[3] = 	{SndSpd*nu[0]/nu_mag + u(x,y,z), 	SndSpd*nu[1]/nu_mag + v(x,y,z), 	SndSpd*nu[2]/nu_mag + w(x,y,z)};
	double c_prop_mag = sqrt(pow(c_prop[0],2) + pow(c_prop[1],2) + pow(c_prop[2],2));
    
	double dxds = c_prop[0]/c_prop_mag,		dyds = c_prop[1]/c_prop_mag,	dzds = c_prop[2]/c_prop_mag;
	double dxdtheta = solution[index][6],	dydtheta = solution[index][7], 	dzdtheta = solution[index][8];
	double dxdphi = solution[index][12], 	dydphi = solution[index][13], 	dzdphi = solution[index][14];
    
	return	dxds*(dydtheta*dzdphi - dydphi*dzdtheta)
                - dxdtheta*(dyds*dzdphi - dzds*dydphi)
                    + dxdphi*(dyds*dzdtheta - dzds*dydtheta);
}


double GeoAc_Amplitude(double ** solution, int index){
	double x =      solution[index][0],         y = solution[index][1],         z = solution[index][2];
    double x0 =     GeoAc_Sources.src_loc[0],   y0 = GeoAc_Sources.src_loc[1],  z0 = GeoAc_Sources.src_loc[2];
	double nu[3] = {solution[index][3],         solution[index][4],             solution[index][5]};
    
	double 	c0 = GeoAc_Sources.c0,
            SndSpd = c(x,y,z),
            Windu = u(x,y,z),
            Windv = v(x,y,z),
            Windw = w(x,y,z),
            Windu0 = u(x0,y0,z0),
            Windv0 = v(x0,y0,z0),
            Windw0 = w(x0,y0,z0),
            nu_mag = (c0 - nu[0]*Windu - nu[1]*Windv - nu[2]*Windw)/SndSpd,
            nu_mag0 = 1.0 - nu[0]*Windu0/c0 - nu[1]*Windv0/c0 - nu[2]*Windw0/c0,
            c_prop[3] =  {SndSpd*nu[0]/nu_mag + Windu,                  SndSpd*nu[1]/nu_mag + Windv,                    SndSpd*nu[2]/nu_mag + Windw},
            c_prop0[3] = {c0*cos(GeoAc_theta)*cos(GeoAc_phi) + Windu0, 	c0*cos(GeoAc_theta)*sin(GeoAc_phi) + Windv0, 	c0*sin(GeoAc_theta) + Windw0};
    
	double  c_prop_mag = sqrt(pow(c_prop[0],2) + pow(c_prop[1],2) + pow(c_prop[2],2)),
    c_prop_mag0 = sqrt(pow(c_prop0[0],2) + pow(c_prop0[1],2) + pow(c_prop0[2],2));
	
	double D = GeoAc_Jacobian(solution, index);
	double Amp_Num = rho(x,y,z) * nu_mag * pow(SndSpd,3) * c_prop_mag0 * cos(GeoAc_theta);
	double Amp_Den = rho(x0,y0,z0) * nu_mag0* pow(c0,3) * c_prop_mag * D;
    
	return 1.0/(4.0*Pi)*sqrt(fabs(Amp_Num/Amp_Den));
}


//--------------------------------------------------------------------------//
//------Integrate the Sutherland Bass Attenuation Through the Ray Path------//
//--------------------------------------------------------------------------//
double GeoAc_SB_Atten(double ** solution, int end, double freq){
	double dx, dy, dz, ds, x, y, z;
    double atten = 0.0;
	for (int n = 0; n < end; n++){
		
		dx = solution[n+1][0] - solution[n][0];                     // Calculate ds = sqrt(dx^2 + dy^2 + dz^2)
		dy = solution[n+1][1] - solution[n][1];
        dz = solution[n+1][2] - solution[n][2];
		
        ds = sqrt(dx*dx + dy*dy + dz*dz);
        
		x = solution[n][0] + dx/2.0;
        y = solution[n][1] + dy/2.0;
		z = solution[n][2] + dz/2.0;
        
		atten += SuthBass_Alpha(x, y, z, freq)*ds;	// Add contribution to the travel time
	}
    return atten;
}

void GeoAc_SB_AttenSegment(double & atten, double ** solution, int start, int end, double freq){
	double dx, dy, dz, ds, x, y, z;
	for (int n = start; n < end; n++){
		
		dx = solution[n+1][0] - solution[n][0];                     // Calculate ds = sqrt(dx^2 + dy^2 + dz^2)
		dy = solution[n+1][1] - solution[n][1];
        dz = solution[n+1][2] - solution[n][2];
		
        ds = sqrt(dx*dx + dy*dy + dz*dz);
        
		x = solution[n][0] + dx/2.0;
        y = solution[n][1] + dy/2.0;
		z = solution[n][2] + dz/2.0;
        
		atten += SuthBass_Alpha(x, y, z, freq)*ds;	// Add contribution to the travel time
	}
}

//--------------------------------------------------------------//
//---------Count the caustics encountered by monitoring---------//
//----how many times the Jacobian determinant changes sign------//
//--------------------------------------------------------------//
int GeoAc_CausticCnt(double ** solution, int start, int end){
	int count = 0;
	double current, prev = GeoAc_Jacobian(solution, 1);
	
    for(int n = start; n < end; n++){
		current = GeoAc_Jacobian(solution,n);
		if(current*prev < 0.0) count++;
		prev = current;
	}

	return count;
}


#endif /* GEOAC_EQSETS_3DRNGDEP_H_ */
