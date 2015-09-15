#ifndef GEOAC_EQSETS_GLOBAL_H_
#define GEOAC_EQSETS_GLOBAL_H_

#include <math.h>
#include <iostream>

#include "GeoAc.Parameters.h"
#include "Atmo_State.h"
#include "G2S_GlobalMultiDimSpline3D.h"

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

    double src_loc[3];  // Location of the source (r, theta, phi)
	double c0;			// Thermodynamic sound speed at source 

	double c;			// Thermodynamic sound speed
	double dc[5];		// Dc_r, Dc_t, Dc_p, dc_lt, dc_lp (derivatives with respect to radius, theta, phi, launch theta, launch phi)
	double ddc[3][2];	// ddc_rlt, ddc_rlp; ddc_tlt, ddc_tlp; ddc_plt, ddc_plp (derivatives with respect to r and lt, r and lp; etc.)
    
	double w;			// Vertical wind speed (r component, positive values upward)
	double dw[5];		// dw_r, dw_t, dw_p, dw_lt, dw_lp
    double ddw[3][2];	// ddw_rlt, ddw_rlp; ddw_tlt, ddw_tlp; ddw_plt, ddw_plp
    
	double v;			// N-S wind speed (theta component, positive values toward south pole)
	double dv[5];		// dv_r, dv_t, dv_p, dv_lt, dv_lp
    double ddv[3][2];	// ddv_rlt, ddv_rlp; ddv_tlt, ddv_tlp; ddv_plt, ddv_plp
    
	double u;			// E-W wind speed (phi component, positive values to the west)
	double du[5];		// du_r, du_t, du_p, du_lt, du_lp
    double ddu[3][2];	// ddu_rlt, ddu_rlp; ddu_tlt, ddu_tlp; ddu_plt, ddu_plp
    
    double nu0;         // Eikonal vector magnitude at the source
    double nu_mag;		// Eikonal vector magnitude
	double dnu_mag[2];	// Derivative of |nu| with respect to lt, lp

	double c_gr[3];         // Group velocity cg = c {nu_x,nu_y,nu_z}/|nu} + {u,v,w}
	double c_gr_mag;        // Group velocity magnitude |cp| = c sqrt(1 + 2 nu/|nu| dot v/c + v^2/c^2)
    double dc_gr[3][2];     // Derivatives of c_gr components with respect to lt and lp
	double dc_gr_mag[2];	// Derivatives of c_gr magnitude with respect to lt and lp
    
    double GeoCoeff[3];         // Scalar coefficients for coordiante transformation (spherical -> 1.0, 1.0/r, 1.0/(r*sin(theta)))
    double d_GeoCoeff[3][2];    // Angular derivatives of the coefficients for a coordiante transformation
    
    double GeoTerms[3];         // Extra terms which keep the eikonal vector direction constant as the unit vectors vary with position
    double d_GeoTerms[3][2];    // Angular derivatives of the extra terms keeping the eikonal vector direction constant
};

struct GeoAc_Sources_Struct GeoAc_Sources = {
    {0.0, 0.0, 0.0}, 0.0,
    0.0, {0.0, 0.0, 0.0, 0.0, 0.0}, {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}},
    0.0, {0.0, 0.0, 0.0, 0.0, 0.0}, {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}},
    0.0, {0.0, 0.0, 0.0, 0.0, 0.0}, {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}},
    0.0, {0.0, 0.0, 0.0, 0.0, 0.0}, {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}},
    0.0, 0.0, {0.0, 0.0},
    {0.0, 0.0, 0.0}, 0.0, {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}, {0.0, 0.0},
    {0.0, 0.0, 0.0}, {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}},
    {0.0, 0.0, 0.0}, {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}
};

//----------------------------------------------------------------------//
//-------Fill in solution[0][n] with the appropriate initial values-----//
//----------------------------------------------------------------------//
void GeoAc_SetInitialConditions(double ** & solution, double r0, double theta0, double phi0){
    GeoAc_Sources.src_loc[0] = r0 + r_earth;
    GeoAc_Sources.src_loc[1] = theta0;
    GeoAc_Sources.src_loc[2] = phi0;    
	GeoAc_Sources.c0 = c(r0 + r_earth, theta0, phi0);
    
    double MachComps[3] = { w(r0 + r_earth, theta0, phi0)/GeoAc_Sources.c0,
                            v(r0 + r_earth, theta0, phi0)/GeoAc_Sources.c0,
                            u(r0 + r_earth, theta0, phi0)/GeoAc_Sources.c0};
    double nu0[3] =    {sin(GeoAc_theta),  cos(GeoAc_theta)*sin(GeoAc_phi),    cos(GeoAc_theta)*cos(GeoAc_phi)};
    double mu0_lt[3] = {cos(GeoAc_theta), -sin(GeoAc_theta)*sin(GeoAc_phi),   -sin(GeoAc_theta)*cos(GeoAc_phi)};
    double mu0_lp[3] = {0.0,               cos(GeoAc_theta)*cos(GeoAc_phi),   -cos(GeoAc_theta)*sin(GeoAc_phi)};
    double MachScalar = 1.0 + (nu0[0]*MachComps[0] + nu0[1]*MachComps[1] + nu0[2]*MachComps[2]);
    GeoAc_Sources.nu0 = 1.0/MachScalar;
    
	for(int index = 0; index < GeoAc_EqCnt; index++){
		switch(index){
			case(0):			// r(0) = r0 relative to earth radius
				solution[0][index] =  r0 + r_earth;
				break;
                
            case(1):			// t(0) = theta0
				solution[0][index] =  theta0;
				break;
                
            case(2):			// p(0) = phi0
                solution[0][index] =  phi0;
				break;
                
			case(3):			// nu_r(0)
			case(4):			// nu_t(0)
			case(5):            // nu_p(0)
				solution[0][index] = nu0[index - 3]/MachScalar;
				break;
                
            case(6):			// Rt(0) = 0
			case(7):			// Tt(0) = 0
			case(8):			// Pt(0) = 0
			case(12):			// Rp(0) = 0
			case(13):			// Tp(0) = 0
			case(14):			// Pp(0) = 0
				solution[0][index] =  0.0;
				break;
                
			case(9):			// mu_r_lt(0) = d/d lt (nu_r)
			case(10):			// mu_t_lt(0) = d/d lt (nu_t)
			case(11):			// mu_p_lt(0) = d/d lt (nu_p)
                solution[0][index] =  mu0_lt[index - 9]/MachScalar - nu0[index - 9]/pow(MachScalar,2.0) * (mu0_lt[0]*MachComps[0] + mu0_lt[1]*MachComps[1] + mu0_lt[2]*MachComps[2]);
				break;
                
			case(15):			// mu_r_lp(0) = d/d lp (nu_r)
            case(16):			// mu_t_lp(0) = d/d lp (nu_t)
            case(17):			// mu_p_lp(0) = d/d lp (nu_p)
                solution[0][index] =  mu0_lp[index - 15]/MachScalar - nu0[index - 15]/pow(MachScalar,2.0) * (mu0_lp[0]*MachComps[0] + mu0_lp[1]*MachComps[1] + mu0_lp[2]*MachComps[2]);
                break;

			default:
                cout << "Unexpected index in Initial_Cond.  Model includes 18 variables." << '\n';
		}
	}
}
//--------------------------------------------------------------------------//
//-------Taylor series fit to more accurately deterine intercept values-----//
//--------------------------------------------------------------------------//
void GeoAc_ApproximateIntercept(double ** solution, int k, double* & prev){
	double result;
	double dr_k = solution[k][0] - solution[k-1][0];        // set dr for step = r_k - r_{k-1}
	double dr_grnd = solution[k-1][0] - (r_earth + z_grnd);	// set dr from r_{k-1} to ground
    
	for(int index = 0; index < GeoAc_EqCnt; index++){
		prev[index] = solution[k-1][index] 	+ (solution[k-1][index] - solution[k][index])/dr_k*dr_grnd;
                        + 1.0/2.0*(solution[k][index] + solution[k-2][index] - 2.0*solution[k-1][index])/pow(dr_k,2.0)*pow(dr_grnd,2.0);
	}
}

//-------------------------------------------------------------------------//
//-------Fill in solution[0][n] with the appropriate reflection values-----//
//-------------------------------------------------------------------------//
void GeoAc_SetReflectionConditions(double** & solution, int k_end){
	double* prev = new double [GeoAc_EqCnt];
	GeoAc_ApproximateIntercept(solution, k_end, prev);
    
    double c_ref = c(prev[0], prev[1], prev[2]);
	double dnu_r_ds = - 1.0/c_ref * (GeoAc_Sources.c0/c_ref * c_diff(prev[0],prev[1],prev[2],0)
                                        + prev[3] * w_diff(prev[0],prev[1],prev[2],0)
                                            + prev[4] * v_diff(prev[0],prev[1],prev[2],0)
                                                + prev[5] * u_diff(prev[0],prev[1],prev[2],0)
                                                    + c_ref/prev[0] * (pow(prev[4],2) + pow(prev[5],2)));
    
	for(int index = 0; index < GeoAc_EqCnt; index++){
		switch(index){
			case(0):			// r(s0+) = 		r_ground (theta(s0-),phi(s0-))
                solution[0][index] = r_earth + z_grnd;
                break;
                
			case(1):			// theta(s0+) = 	theta(s0-)
			case(2):			// phi(s0+) =       phi(s0-)

			case(4):			// nu_theta(s0+) = 	nu_theta(s0-)
			case(5):			// nu_phi(s0+) = 	nu_phi(s0-)
			
			case(7):			// Theta_lt(s0+) = 	Theta_lt(s0-)
			case(8):			// Phi_lt(s0+) =     Phi_lt(s0-)
			
			case(10):			// mu_t_lt(s0+) = 	mu_t_lt(s0-)
			case(11):			// mu_p_lt(s0+) =  	mu_p_lt(s0-)
                
			case(13):			// Theta_lp(s0+) = 	Theta_lp(s0-)
			case(14):			// Phi_lp(s0+) = 	Phi_lp(s0-)
			
			case(16):			// mu_t_lp(s0+) = 	mu_t_lp(s0-)
 			case(17):			// mu_p_lp(s0+) = 	mu_p_lp(s0-)
                solution[0][index] =  prev[index];
                break;
		
			case(3):			// nu_r(s0+) =      -nu_r(s0-)
            case(6):			// R_lt(s0+) = 		-R_lt(s0-)
            case(12):			// R_lp(s0+) = 		-R_lp(s0-)
                solution[0][index] =  -prev[index];
                break;
                
            case(9):			// mu_r_lt(s0+) = 	-mu_r_lt(s0-) + 2 d nu_r/ds * ds0/dlt
            case(15):			// mu_r_lp(s0+) = 	mu_r_lp(s0-)  2 d nu_r/ds * sd0/dlp
                solution[0][index] =  -prev[index] + 2.0*dnu_r_ds * prev[index - 3]/( c_ref / GeoAc_Sources.c0 * prev[3]);
                break;


        
        
        }
	}
	delete [] prev;
}

//-----------------------------------------------------------------------------------//
//-------Vary the solver step size, currently uses smaller steps near the ground-----//
//-----------------------------------------------------------------------------------//
double GeoAc_Set_ds(double* current_values){
	double r = current_values[0];
	double result = 0.05 - 0.049 * exp(-(r - (r_earth + z_grnd))/0.75);
    
	result = min(result, GeoAc_ds_max);
	result = max(result, GeoAc_ds_min);
	return result;
}

//---------------------------------------//
//-------Update the source functions-----//
//---------------------------------------//
void GeoAc_UpdateSources(double ray_length, double* current_values){
    // Extract ray location and eikonal vector components
    double r = current_values[0],		theta = current_values[1], 	phi = current_values[2];
	double nu[3] = {current_values[3], 	current_values[4], 		current_values[5]};
    
    double temp;
    double dtemp[3];
    
    if(!GeoAc_CalcAmp){
        Eval_Spline_AllOrder1(r, theta, phi, Temp_Spline, temp, dtemp[0], dtemp[1], dtemp[2]);
        for(int n = 0; n < 3; n++){
            Windu_Spline.accel[n] = Temp_Spline.accel[n];
            Windv_Spline.accel[n] = Temp_Spline.accel[n];
        }
        
        Eval_Spline_AllOrder1(r, theta, phi, Windu_Spline, GeoAc_Sources.u, GeoAc_Sources.du[0], GeoAc_Sources.du[1], GeoAc_Sources.du[2]);
        Eval_Spline_AllOrder1(r, theta, phi, Windv_Spline, GeoAc_Sources.v, GeoAc_Sources.dv[0], GeoAc_Sources.dv[1], GeoAc_Sources.dv[2]);
        GeoAc_Sources.w = w(r, theta, phi);
        
        GeoAc_Sources.c = sqrt(gamR * temp);
        for(int n = 0; n < 3; n++){
            GeoAc_Sources.dc[n] = gamR / (2.0 * GeoAc_Sources.c) * dtemp[n];
            GeoAc_Sources.dw[n] = w_diff(r,theta,phi,n);
        }
        
        GeoAc_Sources.nu_mag = 	 sqrt(nu[0]*nu[0] + nu[1]*nu[1] + nu[2]*nu[2]);
        
        GeoAc_Sources.c_gr[0] =  GeoAc_Sources.c*nu[0]/GeoAc_Sources.nu_mag + GeoAc_Sources.w;
        GeoAc_Sources.c_gr[1] =  GeoAc_Sources.c*nu[1]/GeoAc_Sources.nu_mag + GeoAc_Sources.v;
        GeoAc_Sources.c_gr[2] =  GeoAc_Sources.c*nu[2]/GeoAc_Sources.nu_mag + GeoAc_Sources.u;
        
        GeoAc_Sources.c_gr_mag = sqrt(pow(GeoAc_Sources.c_gr[0],2) + pow(GeoAc_Sources.c_gr[1],2) + pow(GeoAc_Sources.c_gr[2],2));
        
        GeoAc_Sources.GeoCoeff[0] = 1.0;
        GeoAc_Sources.GeoCoeff[1] = 1.0/r;
        GeoAc_Sources.GeoCoeff[2] = 1.0/(r*cos(theta));
        
        GeoAc_Sources.GeoTerms[0] = 0.0;
        GeoAc_Sources.GeoTerms[1] = (nu[0]*GeoAc_Sources.v - nu[1]*GeoAc_Sources.w);
        GeoAc_Sources.GeoTerms[2] = (nu[0]*GeoAc_Sources.u - nu[2]*GeoAc_Sources.w)*cos(theta) + (nu[1]*GeoAc_Sources.u - nu[2]*GeoAc_Sources.v)*sin(theta);
        
        GeoAc_Sources.GeoTerms[0] += 1.0/r * (nu[1]*GeoAc_Sources.c_gr[1] + nu[2]*GeoAc_Sources.c_gr[2]);
        GeoAc_Sources.GeoTerms[1] += -nu[0]*GeoAc_Sources.c_gr[1] + nu[2]*GeoAc_Sources.c_gr[2]*tan(theta);
        GeoAc_Sources.GeoTerms[2] += -GeoAc_Sources.c_gr[2]*(nu[0]*cos(theta) + nu[1]*sin(theta));
    }
    
    double R_th[3], R_ph[3], mu_th[3], mu_ph[3];
    double ddtemp[3][3];
    double ddWindu[3][3];
    double ddWindv[3][3];
    
    if(GeoAc_CalcAmp){
        R_th[0]  = current_values[6];       R_th[1]  = current_values[7];       R_th[2]  = current_values[8];
		mu_th[0] = current_values[9];		mu_th[1] = current_values[10];		mu_th[2] = current_values[11];
        R_ph[0]  = current_values[12];      R_ph[1]  = current_values[13];      R_ph[2]  = current_values[14];
		mu_ph[0] = current_values[15];		mu_ph[1] = current_values[16];		mu_ph[2] = current_values[17];
        
        Eval_Spline_AllOrder2(r, theta, phi, Temp_Spline, temp, dtemp[0], dtemp[1], dtemp[2], ddtemp[0][0], ddtemp[1][1], ddtemp[2][2], ddtemp[0][1], ddtemp[0][2], ddtemp[1][2]);
        for(int n = 0; n < 3; n++){
            Windu_Spline.accel[n] = Temp_Spline.accel[n];
            Windv_Spline.accel[n] = Temp_Spline.accel[n];
        }
        
        Eval_Spline_AllOrder2(r, theta, phi, Windu_Spline, GeoAc_Sources.u, GeoAc_Sources.du[0], GeoAc_Sources.du[1], GeoAc_Sources.du[2], ddWindu[0][0], ddWindu[1][1], ddWindu[2][2], ddWindu[0][1], ddWindu[0][2], ddWindu[1][2]);
        Eval_Spline_AllOrder2(r, theta, phi, Windv_Spline, GeoAc_Sources.v, GeoAc_Sources.dv[0], GeoAc_Sources.dv[1], GeoAc_Sources.dv[2], ddWindv[0][0], ddWindv[1][1], ddWindv[2][2], ddWindv[0][1], ddWindv[0][2], ddWindv[1][2]);
        GeoAc_Sources.w = w(r, theta, phi);
        
        ddtemp[1][0] = ddtemp[0][1];    ddtemp[2][0] = ddtemp[0][2];    ddtemp[2][1] = ddtemp[1][2];
        ddWindu[1][0] = ddWindu[0][1];  ddWindu[2][0] = ddWindu[0][2];  ddWindu[2][1] = ddWindu[1][2];
        ddWindv[1][0] = ddWindv[0][1];  ddWindv[2][0] = ddWindv[0][2];  ddWindv[2][1] = ddWindv[1][2];
        
        GeoAc_Sources.c = sqrt(gamR * temp);
        for(int n = 0; n < 3; n++){
            GeoAc_Sources.dc[n] = gamR / (2.0 * GeoAc_Sources.c) * dtemp[n];
            GeoAc_Sources.dw[n] = w_diff(r,theta,phi,n);
            
            GeoAc_Sources.ddc[n][0] = 0.0;            GeoAc_Sources.ddc[n][1] = 0.0;
            GeoAc_Sources.ddu[n][0] = 0.0;            GeoAc_Sources.ddu[n][1] = 0.0;
            GeoAc_Sources.ddv[n][0] = 0.0;            GeoAc_Sources.ddv[n][1] = 0.0;
            GeoAc_Sources.ddw[n][0] = 0.0;            GeoAc_Sources.ddw[n][1] = 0.0;
            for(int m = 0; m < 3; m++){
                GeoAc_Sources.ddc[n][0] += R_th[m]*(gamR/(2.0*GeoAc_Sources.c) * ddtemp[n][m] - pow(gamR,2)/(4.0 * pow(GeoAc_Sources.c,3)) * dtemp[n]*dtemp[m]);
                GeoAc_Sources.ddc[n][1] += R_ph[m]*(gamR/(2.0*GeoAc_Sources.c) * ddtemp[n][m] - pow(gamR,2)/(4.0 * pow(GeoAc_Sources.c,3)) * dtemp[n]*dtemp[m]);
                
                GeoAc_Sources.ddu[n][0] += R_th[m]*ddWindu[n][m];               GeoAc_Sources.ddu[n][1] += R_ph[m]*ddWindu[n][m];
                GeoAc_Sources.ddv[n][0] += R_th[m]*ddWindv[n][m];               GeoAc_Sources.ddv[n][1] += R_ph[m]*ddWindv[n][m];
                GeoAc_Sources.ddw[n][0] += R_th[m]*w_ddiff(r,theta,phi,n,m);    GeoAc_Sources.ddw[n][1] += R_ph[m]*w_ddiff(r,theta,phi,n,m);
            }
        }
        
        GeoAc_Sources.dc[3] = 0.0;      GeoAc_Sources.dc[4] = 0.0;
        GeoAc_Sources.du[3] = 0.0;      GeoAc_Sources.du[4] = 0.0;
        GeoAc_Sources.dv[3] = 0.0;      GeoAc_Sources.dv[4] = 0.0;
        GeoAc_Sources.dw[3] = 0.0;      GeoAc_Sources.dw[4] = 0.0;
        for(int n = 0; n < 3; n++){
            GeoAc_Sources.dc[3] += R_th[n]*GeoAc_Sources.dc[n];     GeoAc_Sources.dc[4] += R_ph[n]*GeoAc_Sources.dc[n];
            GeoAc_Sources.du[3] += R_th[n]*GeoAc_Sources.du[n];     GeoAc_Sources.du[4] += R_ph[n]*GeoAc_Sources.du[n];
            GeoAc_Sources.dv[3] += R_th[n]*GeoAc_Sources.dv[n];     GeoAc_Sources.dv[4] += R_ph[n]*GeoAc_Sources.dv[n];
            GeoAc_Sources.dw[3] += R_th[n]*GeoAc_Sources.dw[n];     GeoAc_Sources.dw[4] += R_ph[n]*GeoAc_Sources.dw[n];
        }
        
        GeoAc_Sources.nu_mag = 	 sqrt(nu[0]*nu[0] + nu[1]*nu[1] + nu[2]*nu[2]);
		GeoAc_Sources.dnu_mag[0] = (nu[0]*mu_th[0] + nu[1]*mu_th[1] + nu[2]*mu_th[2])/GeoAc_Sources.nu_mag;
		GeoAc_Sources.dnu_mag[1] = (nu[0]*mu_ph[0] + nu[1]*mu_ph[1] + nu[2]*mu_ph[2])/GeoAc_Sources.nu_mag;
        
        GeoAc_Sources.c_gr[0] =  GeoAc_Sources.c*nu[0]/GeoAc_Sources.nu_mag + GeoAc_Sources.w;
        GeoAc_Sources.c_gr[1] =  GeoAc_Sources.c*nu[1]/GeoAc_Sources.nu_mag + GeoAc_Sources.v;
        GeoAc_Sources.c_gr[2] =  GeoAc_Sources.c*nu[2]/GeoAc_Sources.nu_mag + GeoAc_Sources.u;
        GeoAc_Sources.c_gr_mag = sqrt(pow(GeoAc_Sources.c_gr[0],2) + pow(GeoAc_Sources.c_gr[1],2) + pow(GeoAc_Sources.c_gr[2],2));
        
        GeoAc_Sources.dc_gr[0][0] = nu[0]/GeoAc_Sources.nu_mag*GeoAc_Sources.dc[3] + GeoAc_Sources.c*mu_th[0]/GeoAc_Sources.nu_mag - GeoAc_Sources.c*nu[0]/pow(GeoAc_Sources.nu_mag,2) * GeoAc_Sources.dnu_mag[0] + GeoAc_Sources.dw[3];
        GeoAc_Sources.dc_gr[1][0] = nu[1]/GeoAc_Sources.nu_mag*GeoAc_Sources.dc[3] + GeoAc_Sources.c*mu_th[1]/GeoAc_Sources.nu_mag - GeoAc_Sources.c*nu[1]/pow(GeoAc_Sources.nu_mag,2) * GeoAc_Sources.dnu_mag[0] + GeoAc_Sources.dv[3];
        GeoAc_Sources.dc_gr[2][0] = nu[2]/GeoAc_Sources.nu_mag*GeoAc_Sources.dc[3] + GeoAc_Sources.c*mu_th[2]/GeoAc_Sources.nu_mag - GeoAc_Sources.c*nu[2]/pow(GeoAc_Sources.nu_mag,2) * GeoAc_Sources.dnu_mag[0] + GeoAc_Sources.du[3];
        GeoAc_Sources.dc_gr_mag[0] = (GeoAc_Sources.c_gr[0]*GeoAc_Sources.dc_gr[0][0] + GeoAc_Sources.c_gr[1]*GeoAc_Sources.dc_gr[1][0] + GeoAc_Sources.c_gr[2]*GeoAc_Sources.dc_gr[2][0])/GeoAc_Sources.c_gr_mag;
        
        GeoAc_Sources.dc_gr[0][1] = nu[0]/GeoAc_Sources.nu_mag*GeoAc_Sources.dc[4] + GeoAc_Sources.c*mu_ph[0]/GeoAc_Sources.nu_mag - GeoAc_Sources.c*nu[0]/pow(GeoAc_Sources.nu_mag,2) * GeoAc_Sources.dnu_mag[1] + GeoAc_Sources.dw[4];
        GeoAc_Sources.dc_gr[1][1] = nu[1]/GeoAc_Sources.nu_mag*GeoAc_Sources.dc[4] + GeoAc_Sources.c*mu_ph[1]/GeoAc_Sources.nu_mag - GeoAc_Sources.c*nu[1]/pow(GeoAc_Sources.nu_mag,2) * GeoAc_Sources.dnu_mag[1] + GeoAc_Sources.dv[4];
        GeoAc_Sources.dc_gr[2][1] = nu[2]/GeoAc_Sources.nu_mag*GeoAc_Sources.dc[4] + GeoAc_Sources.c*mu_ph[2]/GeoAc_Sources.nu_mag - GeoAc_Sources.c*nu[2]/pow(GeoAc_Sources.nu_mag,2) * GeoAc_Sources.dnu_mag[1] + GeoAc_Sources.du[4];
        GeoAc_Sources.dc_gr_mag[1] = (GeoAc_Sources.c_gr[0]*GeoAc_Sources.dc_gr[0][1] + GeoAc_Sources.c_gr[1]*GeoAc_Sources.dc_gr[1][1] + GeoAc_Sources.c_gr[2]*GeoAc_Sources.dc_gr[2][1])/GeoAc_Sources.c_gr_mag;
        
        GeoAc_Sources.GeoCoeff[0] = 1.0;
        GeoAc_Sources.GeoCoeff[1] = 1.0/r;
        GeoAc_Sources.GeoCoeff[2] = 1.0/(r*cos(theta));
        
        GeoAc_Sources.d_GeoCoeff[0][0] = 0.0;
        GeoAc_Sources.d_GeoCoeff[1][0] = -R_th[0]/(pow(r,2));
        GeoAc_Sources.d_GeoCoeff[2][0] = -R_th[0]/(pow(r,2)*cos(theta)) + sin(theta)*R_th[1]/(r*pow(cos(theta),2));
        
        GeoAc_Sources.d_GeoCoeff[0][1] = 0.0;
        GeoAc_Sources.d_GeoCoeff[1][1] = -R_ph[0]/(pow(r,2));
        GeoAc_Sources.d_GeoCoeff[2][1] = -R_ph[0]/(pow(r,2)*cos(theta)) + sin(theta)*R_ph[1]/(r*pow(cos(theta),2));
        
        GeoAc_Sources.GeoTerms[0] = 0.0;
        GeoAc_Sources.GeoTerms[1] = (nu[0]*GeoAc_Sources.v - nu[1]*GeoAc_Sources.w);
        GeoAc_Sources.GeoTerms[2] = (nu[0]*GeoAc_Sources.u - nu[2]*GeoAc_Sources.w)*cos(theta) + (nu[1]*GeoAc_Sources.u - nu[2]*GeoAc_Sources.v)*sin(theta);
        
        GeoAc_Sources.GeoTerms[0] += 1.0/r * (nu[1]*GeoAc_Sources.c_gr[1] + nu[2]*GeoAc_Sources.c_gr[2]);
        GeoAc_Sources.GeoTerms[1] += -nu[0]*GeoAc_Sources.c_gr[1] + nu[2]*GeoAc_Sources.c_gr[2]*tan(theta);
        GeoAc_Sources.GeoTerms[2] += -GeoAc_Sources.c_gr[2]*(nu[0]*cos(theta) + nu[1]*sin(theta));
        
        GeoAc_Sources.d_GeoTerms[0][0] = 0.0;
        GeoAc_Sources.d_GeoTerms[1][0] = (mu_th[0]*GeoAc_Sources.v + nu[0]*GeoAc_Sources.dv[3] - mu_th[1]*GeoAc_Sources.w - nu[1] * GeoAc_Sources.dw[3]);
        GeoAc_Sources.d_GeoTerms[2][0] = (mu_th[0]*GeoAc_Sources.u + nu[0]*GeoAc_Sources.du[3] - mu_th[2]*GeoAc_Sources.w - nu[2] * GeoAc_Sources.dw[3])*cos(theta) - (nu[0]*GeoAc_Sources.u - nu[2]*GeoAc_Sources.w)*R_th[1]*sin(theta)
                                            + (mu_th[1]*GeoAc_Sources.u + nu[1]*GeoAc_Sources.du[3] - mu_th[2]*GeoAc_Sources.v - nu[2] * GeoAc_Sources.dv[3])*sin(theta) + (nu[1]*GeoAc_Sources.u - nu[2]*GeoAc_Sources.v)*R_th[1]*cos(theta);
        
        GeoAc_Sources.d_GeoTerms[0][0] += -R_th[0]/pow(r,2)*(nu[1]*GeoAc_Sources.c_gr[1] + nu[2]*GeoAc_Sources.c_gr[2])
                                            + 1.0/r*(mu_th[1]*GeoAc_Sources.c_gr[1] + nu[1]*GeoAc_Sources.dc_gr[1][0] + mu_th[2]*GeoAc_Sources.c_gr[2] + nu[2]*GeoAc_Sources.dc_gr[2][0]);
        GeoAc_Sources.d_GeoTerms[1][0] += -mu_th[0]*GeoAc_Sources.c_gr[1] - nu[0]*GeoAc_Sources.dc_gr[1][0] + mu_th[2]*GeoAc_Sources.c_gr[2]*tan(theta) + nu[2]*GeoAc_Sources.dc_gr[2][0]*tan(theta) + nu[2]*GeoAc_Sources.c_gr[2]*R_th[1]/pow(cos(theta),2);
        GeoAc_Sources.d_GeoTerms[2][0] += -GeoAc_Sources.dc_gr[2][0]*(nu[0]*cos(theta) + nu[1]*sin(theta)) - GeoAc_Sources.c_gr[2]*(mu_th[0]*cos(theta) - nu[0]*R_th[1]*sin(theta) + mu_th[1]*sin(theta) + nu[1]*R_th[1]*cos(theta));
        
        GeoAc_Sources.d_GeoTerms[0][1] = 0.0;
        GeoAc_Sources.d_GeoTerms[1][1] = (mu_ph[0]*GeoAc_Sources.v + nu[0]*GeoAc_Sources.dv[4] - mu_ph[1]*GeoAc_Sources.w - nu[1] * GeoAc_Sources.dw[4]);
        GeoAc_Sources.d_GeoTerms[2][1] = (mu_ph[0]*GeoAc_Sources.u + nu[0]*GeoAc_Sources.du[4] - mu_ph[2]*GeoAc_Sources.w - nu[2] * GeoAc_Sources.dw[4])*cos(theta) - (nu[0]*GeoAc_Sources.u - nu[2]*GeoAc_Sources.w)*R_ph[1]*sin(theta)
                                            + (mu_ph[1]*GeoAc_Sources.u + nu[1]*GeoAc_Sources.du[4] - mu_ph[2]*GeoAc_Sources.v - nu[2] * GeoAc_Sources.dv[4])*sin(theta) + (nu[1]*GeoAc_Sources.u - nu[2]*GeoAc_Sources.v)*R_ph[1]*cos(theta);
        
        GeoAc_Sources.d_GeoTerms[0][1] += -R_ph[0]/pow(r,2)*(nu[1]*GeoAc_Sources.c_gr[1] + nu[2]*GeoAc_Sources.c_gr[2])
                                            + 1.0/r*(mu_ph[1]*GeoAc_Sources.c_gr[1] + nu[1]*GeoAc_Sources.dc_gr[1][1] + mu_ph[2]*GeoAc_Sources.c_gr[2] + nu[2]*GeoAc_Sources.dc_gr[2][1]);
        GeoAc_Sources.d_GeoTerms[1][1] += -mu_ph[0]*GeoAc_Sources.c_gr[1] - nu[0]*GeoAc_Sources.dc_gr[1][1] + mu_ph[2]*GeoAc_Sources.c_gr[2]*tan(theta) + nu[2]*GeoAc_Sources.dc_gr[2][1]*tan(theta) + nu[2]*GeoAc_Sources.c_gr[2]*R_ph[1]/pow(cos(theta),2);
        GeoAc_Sources.d_GeoTerms[2][1] += -GeoAc_Sources.dc_gr[2][1]*(nu[0]*cos(theta) + nu[1]*sin(theta)) - GeoAc_Sources.c_gr[2]*(mu_ph[0]*cos(theta) - nu[0]*R_ph[1]*sin(theta) + mu_ph[1]*sin(theta) + nu[1]*R_ph[1]*cos(theta));
	}
}
//-----------------------------------------------------------//
//-------Evaluate the Source Equation For Specific Index-----//
//-----------------------------------------------------------//
double GeoAc_EvalSrcEq(double ray_length, double* current_values, int Eq_Number){
	double result;
    
	// Set variables used in all equations
    	double nu[3] =      {current_values[3], 	current_values[4], 		current_values[5]};
    	double mu_lt[3] =   {current_values[9], 	current_values[10],		current_values[11]};
    	double mu_lp[3] =   {current_values[15], 	current_values[16],		current_values[17]};
    
	switch(Eq_Number){
		case(0):	// d r / d s
        case(1):    // d theta / d s
        case(2):    // d phi / d s
			result = GeoAc_Sources.GeoCoeff[Eq_Number]*GeoAc_Sources.c_gr[Eq_Number]/GeoAc_Sources.c_gr_mag;
			break;
            
		case(3): 	// d nu_r /ds
		case(4): 	// d nu_theta / ds
		case(5): 	// d nu_phi / ds
			result = -GeoAc_Sources.GeoCoeff[Eq_Number-3]/GeoAc_Sources.c_gr_mag*(GeoAc_Sources.nu_mag*GeoAc_Sources.dc[Eq_Number-3]
                        + nu[0]*GeoAc_Sources.dw[Eq_Number-3] + nu[1]*GeoAc_Sources.dv[Eq_Number-3] + nu[2]*GeoAc_Sources.du[Eq_Number-3] + GeoAc_Sources.GeoTerms[Eq_Number-3]);
			break;
            
            
		case(6):	// d R_lt/ds
        case(7):	// d Theta_lt/ds
        case(8): 	// d Phi_lt/ds
			result = GeoAc_Sources.d_GeoCoeff[Eq_Number-6][0]*GeoAc_Sources.c_gr[Eq_Number-6]/GeoAc_Sources.c_gr_mag
                         + GeoAc_Sources.GeoCoeff[Eq_Number-6]*GeoAc_Sources.dc_gr[Eq_Number-6][0]/GeoAc_Sources.c_gr_mag
                            - GeoAc_Sources.GeoCoeff[Eq_Number-6]*GeoAc_Sources.c_gr[Eq_Number-6]/pow(GeoAc_Sources.c_gr_mag,2) * GeoAc_Sources.dc_gr_mag[0];
			break;
            
            
		case(9): 	// d mu_r_lt /ds
		case(10): 	// d mu_theta_lt /ds
		case(11): 	// d mu_phi_lt /ds
			result = -GeoAc_Sources.d_GeoCoeff[Eq_Number-9][0]/GeoAc_Sources.c_gr_mag*(GeoAc_Sources.nu_mag*GeoAc_Sources.dc[Eq_Number-9]
                        + nu[0]*GeoAc_Sources.dw[Eq_Number-9] + nu[1]*GeoAc_Sources.dv[Eq_Number-9] + nu[2]*GeoAc_Sources.du[Eq_Number-9] + GeoAc_Sources.GeoTerms[Eq_Number-9])
                    + GeoAc_Sources.GeoCoeff[Eq_Number-9]/pow(GeoAc_Sources.c_gr_mag,2) * GeoAc_Sources.dc_gr_mag[0]*(GeoAc_Sources.nu_mag*GeoAc_Sources.dc[Eq_Number-9]
                        + nu[0]*GeoAc_Sources.dw[Eq_Number-9] + nu[1]*GeoAc_Sources.dv[Eq_Number-9] + nu[2]*GeoAc_Sources.du[Eq_Number-9])
                    - GeoAc_Sources.GeoCoeff[Eq_Number-9]/GeoAc_Sources.c_gr_mag*(GeoAc_Sources.dnu_mag[0]*GeoAc_Sources.dc[Eq_Number-9] + GeoAc_Sources.nu_mag*GeoAc_Sources.ddc[Eq_Number-9][0]
                        + mu_lt[0]*GeoAc_Sources.dw[Eq_Number-9] + mu_lt[1]*GeoAc_Sources.dv[Eq_Number-9] + mu_lt[2]*GeoAc_Sources.du[Eq_Number-9]
                            + nu[0]*GeoAc_Sources.ddw[Eq_Number-9][0] + nu[1]*GeoAc_Sources.ddv[Eq_Number-9][0] + nu[2]*GeoAc_Sources.ddu[Eq_Number-9][0] + GeoAc_Sources.d_GeoTerms[Eq_Number-9][0]);
			break;
            
            
		case(12):  	// dR_lp/ds
		case(13):  	// dTheta_lp/ds
        case(14):	// dPhi_lp/ds
			result = GeoAc_Sources.d_GeoCoeff[Eq_Number-12][1]*GeoAc_Sources.c_gr[Eq_Number-12]/GeoAc_Sources.c_gr_mag
                        + GeoAc_Sources.GeoCoeff[Eq_Number-12]*GeoAc_Sources.dc_gr[Eq_Number-12][1]/GeoAc_Sources.c_gr_mag
                            - GeoAc_Sources.GeoCoeff[Eq_Number-12]*GeoAc_Sources.c_gr[Eq_Number-12]/pow(GeoAc_Sources.c_gr_mag,2) * GeoAc_Sources.dc_gr_mag[1];
			break;
            
		case(15): 	// d mu_r_lp/ds
		case(16): 	// d mu_t_lp/ds
		case(17): 	// d mu_p_lp/ds
			result = -GeoAc_Sources.d_GeoCoeff[Eq_Number-15][1]/GeoAc_Sources.c_gr_mag*(GeoAc_Sources.nu_mag*GeoAc_Sources.dc[Eq_Number-15]
                        + nu[0]*GeoAc_Sources.dw[Eq_Number-15] + nu[1]*GeoAc_Sources.dv[Eq_Number-15] + nu[2]*GeoAc_Sources.du[Eq_Number-15] + GeoAc_Sources.GeoTerms[Eq_Number-15])
                    + GeoAc_Sources.GeoCoeff[Eq_Number-15]/pow(GeoAc_Sources.c_gr_mag,2) * GeoAc_Sources.dc_gr_mag[1]*(GeoAc_Sources.nu_mag*GeoAc_Sources.dc[Eq_Number-15]
                        + nu[0]*GeoAc_Sources.dw[Eq_Number-15] + nu[1]*GeoAc_Sources.dv[Eq_Number-15] + nu[2]*GeoAc_Sources.du[Eq_Number-15])
                    - GeoAc_Sources.GeoCoeff[Eq_Number-15]/GeoAc_Sources.c_gr_mag*(GeoAc_Sources.dnu_mag[1]*GeoAc_Sources.dc[Eq_Number-15] + GeoAc_Sources.nu_mag*GeoAc_Sources.ddc[Eq_Number-15][1]
                        + mu_lp[0]*GeoAc_Sources.dw[Eq_Number-15] + mu_lp[1]*GeoAc_Sources.dv[Eq_Number-15] + mu_lp[2]*GeoAc_Sources.du[Eq_Number-15]
                            + nu[0]*GeoAc_Sources.ddw[Eq_Number-15][1] + nu[1]*GeoAc_Sources.ddv[Eq_Number-15][1] + nu[2]*GeoAc_Sources.ddu[Eq_Number-15][1] + GeoAc_Sources.d_GeoTerms[Eq_Number-15][1]);
            break;
	}
	return result;
}






//-------------------------------------------------------------------//
//-------Calculate the Hamiltonian (Eikonal) To Check For Errors-----//
//-------------------------------------------------------------------//
double GeoAc_EvalHamiltonian(double ray_length, double* current_values, double c0){
	double r = current_values[0],  theta = current_values[1], phi = current_values[2];
	double nu[3] = {current_values[3], current_values[4], current_values[5]};
	
	return sqrt(nu[0]*nu[0] + nu[1]*nu[1] + nu[2]*nu[2]) - c0/c(r,theta,phi)
    + (w(r,theta,phi)*nu[0] + v(r,theta,phi)*nu[1] + u(r,theta,phi)*nu[2])/c(r, theta, phi);
}

double GeoAc_EvalHamiltonian_solution(double** solution, int index){
	double nu[3] = {solution[index][3], solution[index][4], solution[index][5]};
	double r = solution[index][0],	theta = solution[index][1],	phi = solution[index][2];
	double r0 = GeoAc_Sources.src_loc[0], 	theta0 = GeoAc_Sources.src_loc[1], 	phi0 = GeoAc_Sources.src_loc[2];
    
	return sqrt(nu[0]*nu[0] + nu[1]*nu[1] + nu[2]*nu[2]) - c(r0,theta0,phi0)/c(r,theta,phi)*(1.0 - (w(r,theta,phi)*nu[0] + v(r,theta,phi)*nu[1] + u(r,theta,phi)*nu[2])/c(r0, theta0, phi0));
}

double GeoAc_EvalHamiltonian_Deriv(double** solution, int index){
	double  r = solution[index][0],
            theta = solution[index][1],
            phi = solution[index][2],
            R_lt[3] =    {solution[index][6],    solution[index][7],     solution[index][8]},
            R_lp[3] =    {solution[index][12],   solution[index][13],    solution[index][14]},
            nu[3] =      {solution[index][3],    solution[index][4],     solution[index][5]},
            mu_lt[3] =   {solution[index][9],    solution[index][10],    solution[index][11]},
            mu_lp[3] =   {solution[index][15],   solution[index][16],    solution[index][17]};

    double mag_nu = sqrt(nu[0]*nu[0]+nu[1]*nu[1]+nu[2]*nu[2]);
    
    double dc_dlt = R_lt[0]*c_diff(r,theta,phi,0) + R_lt[1]*c_diff(r,theta,phi,1) + R_lt[2]*c_diff(r,theta,phi,2);
    double dc_dlp = R_lp[0]*c_diff(r,theta,phi,0) + R_lp[1]*c_diff(r,theta,phi,1) + R_lp[2]*c_diff(r,theta,phi,2);

    double dw_dlt = R_lt[0]*w_diff(r,theta,phi,0) + R_lt[1]*w_diff(r,theta,phi,1) + R_lt[2]*w_diff(r,theta,phi,2);
    double dw_dlp = R_lp[0]*w_diff(r,theta,phi,0) + R_lp[1]*w_diff(r,theta,phi,1) + R_lp[2]*w_diff(r,theta,phi,2);

    double dv_dlt = R_lt[0]*v_diff(r,theta,phi,0) + R_lt[1]*v_diff(r,theta,phi,1) + R_lt[2]*v_diff(r,theta,phi,2);
    double dv_dlp = R_lp[0]*v_diff(r,theta,phi,0) + R_lp[1]*v_diff(r,theta,phi,1) + R_lp[2]*v_diff(r,theta,phi,2);

    double du_dlt = R_lt[0]*u_diff(r,theta,phi,0) + R_lt[1]*u_diff(r,theta,phi,1) + R_lt[2]*u_diff(r,theta,phi,2);
    double du_dlp = R_lp[0]*u_diff(r,theta,phi,0) + R_lp[1]*u_diff(r,theta,phi,1) + R_lp[2]*u_diff(r,theta,phi,2);

    
    double Resid_lt = (nu[0]*mu_lt[0] + nu[1]*mu_lt[1] + nu[2]*mu_lt[2])/mag_nu + mag_nu/c(r,theta,phi) * dc_dlt
                        + 1.0/c(r,theta,phi) * (mu_lt[0]*w(r,theta,phi) + mu_lt[1]*v(r,theta,phi) + mu_lt[2]*u(r,theta,phi) + nu[0]*dw_dlt + nu[1]*dv_dlt + nu[2]*du_dlt);
    double Resid_lp = (nu[0]*mu_lp[0] + nu[1]*mu_lp[1] + nu[2]*mu_lp[2])/mag_nu + mag_nu/c(r,theta,phi) * dc_dlp
                        + 1.0/c(r,theta,phi) * (mu_lp[0]*w(r,theta,phi) + mu_lp[1]*v(r,theta,phi) + mu_lp[2]*u(r,theta,phi) + nu[0]*dw_dlp + nu[1]*dv_dlp + nu[2]*du_dlp);
    
    
    return sqrt(Resid_lt*Resid_lt + Resid_lp*Resid_lp);

}




//--------------------------------------------------------------------------//
//-------Check if ray has left propagation region or returned to ground-----//
//--------------------------------------------------------------------------//
bool GeoAc_BreakCheck(double ** solution, int index){
    bool check = false;
    
	double alt = solution[index][0];
    double lat = solution[index][1];
    double lon = solution[index][2];

    if(alt > GeoAc_vert_limit) check = true;
    if(lat < GeoAc_lat_min_limit || lat > GeoAc_lat_max_limit) check = true;
    if(lon < GeoAc_lon_min_limit || lon > GeoAc_lon_max_limit) check = true;
	
	return check;
}

bool GeoAc_GroundCheck(double ** solution, int index){
	double alt 	= solution[index][0];
	
	bool check = false;
	if(alt < (r_earth + z_grnd)) check = true;
	
	return check;
}

//----------------------------------------------------------------------------------//
//-------Calculate the travel time from source to location or between locations-----//
//----------------------------------------------------------------------------------//
double GeoAc_TravelTime(double ** solution, int index){
    double dr, dt, dp, ds, r, t, p, nu[3], nu_mag, c_prop[3], c_prop_mag;
	double traveltime = 0.0;
	
	for (int n = 0; n < index; n++){
		
		dr = solution[n+1][0] - solution[n][0];		// Calculate ds = sqrt(dr^2 + r^2 dt^2 + r^2 sin^2(t) dp^2)
		dt = solution[n+1][1] - solution[n][1];
		dp = solution[n+1][2] - solution[n][2];
        
		r = solution[n][0] + dr/2.0;				// Calculate (r,t,p) at the midpoint between [n] and [n+1]
		t = solution[n][1] + dt/2.0;
		p = solution[n][2] + dp/2.0;

        ds = sqrt(pow(dr,2) + pow(r*dt,2) + pow(r*cos(t)*dp,2));

		// Calculate |c_prop|
		nu[0] = solution[n][3] + (solution[n+1][3] - solution[n][3])/2.0;
		nu[1] = solution[n][4] + (solution[n+1][4] - solution[n][4])/2.0;
		nu[2] = solution[n][5] + (solution[n+1][5] - solution[n][5])/2.0;
		nu_mag = sqrt(nu[0]*nu[0] + nu[1]*nu[1] + nu[2]*nu[2]);
        
		c_prop[0] = c(r,t,p)*nu[0]/nu_mag + w(r,t,p);
		c_prop[1] = c(r,t,p)*nu[1]/nu_mag + v(r,t,p);
		c_prop[2] = c(r,t,p)*nu[2]/nu_mag + u(r,t,p);
		c_prop_mag = sqrt(pow(c_prop[0],2) + pow(c_prop[1],2) + pow(c_prop[2],2));
        
		traveltime += ds/c_prop_mag;	// Add contribution to the travel time
	}
    
	return traveltime;
}


void GeoAc_TravelTimeSegment(double & time, double ** solution, int start, int end){
    double dr, dt, dp, ds, r, t, p, nu[3], nu_mag, c_prop[3], c_prop_mag;
    
	for (int n = start; n < end; n++){
        
		dr = solution[n+1][0] - solution[n][0];		// Calculate ds = sqrt(dx^2 + dy^2 + dz^2)
		dt = solution[n+1][1] - solution[n][1];
		dp = solution[n+1][2] - solution[n][2];
        
		r = solution[n][0] + dr/2.0;				// Calculate (x,y,z) at the midpoint between [n] and [n+1]
		t = solution[n][1] + dt/2.0;
		p = solution[n][2] + dp/2.0;
        
        ds = sqrt(pow(dr,2) + pow(r*dt,2) + pow(r*cos(t)*dp,2));
        
		// Calculate |c_prop|
		nu[0] = solution[n][3] + (solution[n+1][3] - solution[n][3])/2.0;
		nu[1] = solution[n][4] + (solution[n+1][4] - solution[n][4])/2.0;
		nu[2] = solution[n][5] + (solution[n+1][5] - solution[n][5])/2.0;
		nu_mag = sqrt(nu[0]*nu[0] + nu[1]*nu[1] + nu[2]*nu[2]);
        
		c_prop[0] = c(r,t,p)*nu[0]/nu_mag + w(r,t,p);
		c_prop[1] = c(r,t,p)*nu[1]/nu_mag + v(r,t,p);
		c_prop[2] = c(r,t,p)*nu[2]/nu_mag + u(r,t,p);
		c_prop_mag = sqrt(pow(c_prop[0],2) + pow(c_prop[1],2) + pow(c_prop[2],2));
		time += ds/c_prop_mag;	// Add contribution to the travel time
	}
}

//-----------------------------------------------------------------------------------//
//-------Calculate the Jacobian determinant and from it the ampltude coefficient-----//
//-----------------------------------------------------------------------------------//
double GeoAc_Jacobian(double ** solution, int index){
	double r = solution[index][0], theta = solution[index][1], phi = solution[index][2];
    
	double nu[3] = 		{solution[index][3], solution[index][4], solution[index][5]};
	double nu_mag = 	sqrt(nu[0]*nu[0] + nu[1]*nu[1] + nu[2]*nu[2]);
	double c_prop[3] = 	{c(r,theta,phi)*nu[0]/nu_mag + w(r,theta,phi), 	c(r,theta,phi)*nu[1]/nu_mag + v(r,theta,phi), 	c(r,theta,phi)*nu[2]/nu_mag + u(r,theta,phi)};
	double c_prop_mag = sqrt(pow(c_prop[0],2) + pow(c_prop[1],2) + pow(c_prop[2],2));
    
	double dr_ds = c_prop[0]/c_prop_mag,	dt_ds = 1.0/r*c_prop[1]/c_prop_mag,	dp_ds = 1.0/(r*sin(theta))*c_prop[2]/c_prop_mag;
	double dr_dlt = solution[index][6],     dt_dlt = solution[index][7],        dp_dlt = solution[index][8];
	double dr_dlp = solution[index][12], 	dt_dlp = solution[index][13],       dp_dlp = solution[index][14];
    
	return	pow(r,2)*cos(theta)*(dr_ds*(dt_dlt*dp_dlp - dt_dlp*dp_dlt) - dr_dlt*(dt_ds*dp_dlp - dp_ds*dt_dlp) + dr_dlp*(dt_ds*dp_dlt - dp_ds*dt_dlt));
}


double GeoAc_Amplitude(double ** solution, int index){
    double r0 = GeoAc_Sources.src_loc[0], theta0 = GeoAc_Sources.src_loc[1], phi0 = GeoAc_Sources.src_loc[2];
	double r = solution[index][0], theta = solution[index][1], phi = solution[index][2];
	double nu[3] = {solution[index][3], solution[index][4], solution[index][5]};
    double nu0[3] = {sin(GeoAc_theta),  cos(GeoAc_theta)*sin(GeoAc_phi),    cos(GeoAc_theta)*cos(GeoAc_phi)};
    
    double  nu_mag = (GeoAc_Sources.c0 - nu[0]*w(r,theta,phi) - nu[1]*v(r,theta,phi) - nu[2]*u(r,theta,phi))/c(r,theta,phi),
            nu_mag0 = GeoAc_Sources.nu0,
            c_prop[3] =  {c(r,theta,phi)*nu[0]/nu_mag + w(r,theta,phi),			c(r,theta,phi)*nu[1]/nu_mag + v(r,theta,phi), 			c(r,theta,phi)*nu[2]/nu_mag + u(r, theta, phi)},
            c_prop0[3] = {GeoAc_Sources.c0*nu0[0]/nu_mag0 + w(r0,theta0,phi0),	GeoAc_Sources.c0*nu0[1]/nu_mag + v(r0,theta0,phi0),		GeoAc_Sources.c0*nu0[2]/nu_mag + u(r0, theta0, phi0)};
    
	double  c_prop_mag =  sqrt(pow(c_prop[0],2) +  pow(c_prop[1],2) +  pow(c_prop[2],2)),
            c_prop_mag0 = sqrt(pow(c_prop0[0],2) + pow(c_prop0[1],2) + pow(c_prop0[2],2));
	
	double D = GeoAc_Jacobian(solution, index);	
	double Amp_Num = rho(r,theta,phi) *  nu_mag * pow(c(r,theta,phi),3)  *  c_prop_mag0 * cos(GeoAc_theta);
	double Amp_Den = rho(r0,theta0,phi0)*nu_mag0* pow(c(r0,theta0,phi0),3)* c_prop_mag  * D;
    
	return 1.0/(4.0*Pi)*sqrt(fabs(Amp_Num/Amp_Den));	
}

//--------------------------------------------------------------------------//
//------Integrate the Sutherland Bass Attenuation Through the Ray Path------//
//--------------------------------------------------------------------------//
double GeoAc_SB_Atten(double ** solution, int end, double freq){
	double dr, dt, dp, ds, r, t, p;
    double atten = 0.0;
	for (int n = 0; n < end; n++){
		
		dr = solution[n+1][0] - solution[n][0];		// Calculate ds = sqrt(dx^2 + dy^2 + dz^2)
		dt = solution[n+1][1] - solution[n][1];
		dp = solution[n+1][2] - solution[n][2];
        
		r = solution[n][0] + dr/2.0;				// Calculate (x,y,z) at the midpoint between [n] and [n+1]
		t = solution[n][1] + dt/2.0;
		p = solution[n][2] + dp/2.0;
        
        ds = sqrt(pow(dr,2) + pow(r*dt,2) + pow(r*sin(t)*dp,2));
        
		atten += SuthBass_Alpha(r, t, p, freq)*ds;	// Add contribution to the travel time
	}
    return atten;
}

void GeoAc_SB_AttenSegment(double & atten, double ** solution, int start, int end, double freq){
	double dr, dt, dp, ds, r, t, p;
	for (int n = start; n < end; n++){
		
		dr = solution[n+1][0] - solution[n][0];		// Calculate ds = sqrt(dx^2 + dy^2 + dz^2)
		dt = solution[n+1][1] - solution[n][1];
		dp = solution[n+1][2] - solution[n][2];
        
		r = solution[n][0] + dr/2.0;				// Calculate (x,y,z) at the midpoint between [n] and [n+1]
		t = solution[n][1] + dt/2.0;
		p = solution[n][2] + dp/2.0;
        
        ds = sqrt(pow(dr,2) + pow(r*dt,2) + pow(r*sin(t)*dp,2));

		atten += SuthBass_Alpha(r, t, p, freq)*ds;	// Add contribution to the travel time
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
		if(current*prev < 0.0){
			count++;
		}
		prev = current;
	}
	return count;	
}

#endif /* GEOAC_EQSETS_GLOBAL_H_ */
