# ifndef G2S_GLIOBAL_RNGDEP_ATMOSPHERE_CPP_
# define G2S_GLIOBAL_RNGDEP_ATMOSPHERE_CPP_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>

#include "G2S_GlobalMultiDimSpline3D.h"
#include "Atmo_State.h"
#include "GeoAc.Parameters.h"

using namespace std;

//-----------------------------------------------//
//---------Define the Propagation Region---------//
//-----------------------------------------------//
double r_min, r_max;
double t_min, t_max;
double p_min, p_max;

void GeoAc_SetPropRegion(){
    r_min = Windu_Spline.r_vals[0], r_max = Windu_Spline.r_vals[Windu_Spline.length_r-1];
    t_min = Windu_Spline.t_vals[0], t_max = Windu_Spline.t_vals[Windu_Spline.length_t-1];
    p_min = Windu_Spline.p_vals[0], p_max = Windu_Spline.p_vals[Windu_Spline.length_p-1];
    
    GeoAc_vert_limit  =	r_max;
    GeoAc_lat_min_limit = t_min;      GeoAc_lat_max_limit = t_max;
    GeoAc_lon_min_limit = p_min;      GeoAc_lon_max_limit = p_max;
}


//-----------------------------------------------//
//---------Topographical Ground Function---------//
//-----------------------------------------------//
double r_earth = 6370.0;
double z_grnd = 0.0;


//----------------------------------------//
//------Parameters for Interpolation------//
//----------------------------------------//
int r_cnt;  // Number of r (altitude) data points
int t_cnt;  // Number of latitude (theta -> t) data points
int p_cnt;  // Number of longitude (phi -> p) data points

double* r_vals;     // r_i elements (r_cnt length)
double* t_vals;     // t_j elements (t_cnt length)
double* p_vals;     // p_k elements (p_cnt length)

double*** T_vals;           // Temperature at (r_i, t_j, p_k) (r_cnt x t_cnt x p_cnt)
double*** T_slopes;         // Slopes for vertical splines of temperature at each t[], p[] node
double*** T_slopes_dt;      // Slopes for d/dp of vertical splines of temperature at each t[], p[] node
double*** T_slopes_dp;      // Slopes for d/dt of vertical splines of temperature at each t[], p[] node

double*** u_vals;           // E-W winds at (r_i, t_j, p_k) (r_cnt x t_cnt x p_cnt)
double*** u_slopes;         // Slopes for vertical splines of temperature at each t[], p[] node
double*** u_slopes_dt;      // Slopes for d/dt of vertical splines of temperature at each t[], p[] node
double*** u_slopes_dp;      // Slopes for d/dp of vertical splines of temperature at each t[], p[] node

double*** v_vals;           // N-S winds at (r_i, t_j, p_k) (r_cnt x t_cnt x p_cnt)
double*** v_slopes;         // Slopes for vertical splines of temperature at each t[], p[] node
double*** v_slopes_dt;      // Slopes for d/dt of vertical splines of temperature at each t[], p[] node
double*** v_slopes_dp;      // Slopes for d/dp of vertical splines of temperature at each t[], p[] node

double*** rho_vals;         // Density at (r_i, t_j, p_k) (r_cnt x t_cnt x p_cnt)
double*** rho_slopes;         // Slopes for vertical splines of temperature at each t[], p[] node
double*** rho_slopes_dt;      // Slopes for d/dt of vertical splines of temperature at each t[], p[] node
double*** rho_slopes_dp;      // Slopes for d/dp of vertical splines of temperature at each t[], p[] node


//----------------------------------------//
//----------File IO Manipulation----------//
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

int file_width(string file_name){
	double temp;
	int known_length = file_length(file_name);
	int count = 0;
    ifstream file_in;   file_in.open(file_name.c_str() );
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

bool Check_G2S_Format(string file_name){
    int width = file_width(file_name);
    if(width==6){ return true;}
    else{         return false;}
}


//----------------------------------------//
//---------G2S Array Manipulation---------//
//----------------------------------------//
void SetUp_G2S_Arrays(char* file_prefix, char* loclat_file, char* loclon_file){
    char output_buffer [50];
    
    t_cnt = file_length(loclat_file);
    p_cnt = file_length(loclon_file);
    
    // Open the first file and determine its length
    sprintf(output_buffer, "%s%i.met", file_prefix, 0);
    if(Check_G2S_Format(output_buffer)) r_cnt = file_length(output_buffer);
    
    // Set the sizes of the various arrays
    r_vals = new double [r_cnt];
    t_vals = new double [t_cnt];
    p_vals = new double [p_cnt];
    
    T_vals = new double** [r_cnt];   rho_vals = new double** [r_cnt];
    u_vals = new double** [r_cnt];   v_vals = new double** [r_cnt];
    
    for(int nr = 0; nr < r_cnt; nr++){
        T_vals[nr] = new double* [t_cnt];    rho_vals[nr] = new double* [t_cnt];
        u_vals[nr] = new double* [t_cnt];    v_vals[nr] = new double* [t_cnt];
        for(int nt = 0; nt < t_cnt; nt++){
            T_vals[nr][nt] = new double [p_cnt];    rho_vals[nr][nt] = new double [p_cnt];
            u_vals[nr][nt] = new double [p_cnt];    v_vals[nr][nt] = new double [p_cnt];
        }
    }
}

void Load_G2S_Multi(char* file_prefix, char* lat_file, char* lon_file, char* option){
    ifstream file_in; double temp;
    char output_buffer [50];
    
    file_in.open(lat_file);
    for (int i = 0; i < t_cnt; i++){
        file_in >> t_vals[i];
        t_vals[i]*=Pi/180.0;
    }
    file_in.close();
    
    file_in.open(lon_file);
    for (int i = 0; i < p_cnt; i++){
        file_in >> p_vals[i];
        p_vals[i]*=Pi/180.0;
    }
    file_in.close();
    
    // Copy the G2S data into the various arrays
    for(int np = 0; np < p_cnt; np++){
    for(int nt = 0; nt < t_cnt; nt++){
        sprintf(output_buffer, "%s%i.met", file_prefix, nt + np*t_cnt);
        ifstream file_in; file_in.open(output_buffer);
        
        if (strncmp(option, "zTuvdp",6) == 0){
            for (int nr = 0; nr < r_cnt; nr++){
                file_in >> r_vals[nr];           // Extract r_i value
                file_in >> T_vals[nr][nt][np];   // Extract T(r_i)
                file_in >> u_vals[nr][nt][np];   // Extract u(r_i)
                file_in >> v_vals[nr][nt][np];   // Extract v(r_i)
                file_in >> rho_vals[nr][nt][np]; // Extract density(r_i)
                file_in >> temp; // Extract rho(r_i)
                
                // Convert to global radius and convert winds m/s -> km/s and scale near the ground to guarantee u(r_g), v(r_g) = 0
                r_vals[nr]+= r_earth;
                u_vals[nr][nt][np]*=(2.0 / (1.0 + exp(-(r_vals[nr] - r_earth - z_grnd)/0.2)) - 1.0) / 1000.0;
                v_vals[nr][nt][np]*=(2.0 / (1.0 + exp(-(r_vals[nr] - r_earth - z_grnd)/0.2)) - 1.0) / 1000.0;
            }
        } else if (strncmp(option, "zuvwTdp",7) == 0){
            for(int nr = 0; nr < r_cnt; nr++){
                file_in >> r_vals[nr];           // Extract r_i value
                file_in >> u_vals[nr][nt][np];   // Extract u(r_i)
                file_in >> v_vals[nr][nt][np];   // Extract v(r_i)
                file_in >> temp;                 // Extract w(r_i) but don't store
                file_in >> T_vals[nr][nt][np];   // Extract T(r_i)
                file_in >> rho_vals[nr][nt][np]; // Extract density(r_i)
                file_in >> temp;                 // Extract p(r_i) but don't store
                
                // Convert to global radius and convert winds m/s -> km/s and scale near the ground to guarantee u(r_g), v(r_g) = 0
                r_vals[nr]+= r_earth;
                u_vals[nr][nt][np]*=(2.0 / (1.0 + exp(-(r_vals[nr] - r_earth - z_grnd)/0.2)) - 1.0) / 1000.0;
                v_vals[nr][nt][np]*=(2.0 / (1.0 + exp(-(r_vals[nr] - r_earth - z_grnd)/0.2)) - 1.0) / 1000.0;
            }
        } else {
            cout << "Unrecognized profile option: " << option << ".  Valid options are: zTuvdp and zuvwTdp" << '\n';
        }
        file_in.close();
    }}
}

void Clear_G2S_Arrays(){
    // Clear the arrays
    delete r_vals;
    delete t_vals;
    delete p_vals;
    
    for(int nr = 0; nr < r_cnt; nr++){
        for(int nt = 0; nt < t_cnt; nt++){
            delete T_vals[nr][nt];  delete rho_vals[nr][nt];
            delete u_vals[nr][nt];  delete v_vals[nr][nt];
        }
        delete T_vals[nr];  delete rho_vals[nr];
        delete u_vals[nr];  delete v_vals[nr];
    }
    delete T_vals;    delete rho_vals;
    delete u_vals;    delete v_vals;
}

//-------------------------------------//
//----------Three-Dimensional----------//
//-------Interpolation Structure-------//
//-------------------------------------//
double BiCubic_ConversionMatrix [16][16] = {
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{-3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{ 2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
	{ 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0},
	{ 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0},
	{-3, 0, 3, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0},
	{ 0, 0, 0, 0,-3, 0, 3, 0, 0, 0, 0, 0,-2, 0,-1, 0},
	{ 9,-9,-9, 9, 6, 3,-6,-3, 6,-6, 3,-3, 4, 2, 2, 1},
	{-6, 6, 6,-6,-3,-3, 3, 3,-4, 4,-2, 2,-2,-2,-1,-1},
	{ 2, 0,-2, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0},
	{ 0, 0, 0, 0, 2, 0,-2, 0, 0, 0, 0, 0, 1, 0, 1, 0},
	{-6, 6, 6,-6,-4,-2, 4, 2,-3, 3,-3, 3,-2,-1,-2,-1},
	{ 4,-4,-4, 4, 2, 2,-2,-2, 2,-2, 2,-2, 1, 1, 1, 1}
};

//----------------------------------------------//
//------------Functions to Build and------------//
//-------Clear the Input and Slope Arrays-------//
//----------------------------------------------//


void BuildInputArrays(struct MultiDimSpline_3DGlobal & Spline){
    Spline.r_vals = new double [Spline.length_r];
    Spline.t_vals = new double [Spline.length_t];
    Spline.p_vals = new double [Spline.length_p];
    
    Spline.f_vals = new double** [Spline.length_r];
    
    for(int i = 0; i < Spline.length_r; i++){
        Spline.f_vals[i] = new double* [Spline.length_t];
        for(int j = 0; j < Spline.length_t; j++){
            Spline.f_vals[i][j] = new double [Spline.length_p];
        }
    }
}
void BuildSlopesArrays(struct MultiDimSpline_3DGlobal & Spline){
    Spline.f_slopes = new double** [Spline.length_r];
    Spline.dfdt_slopes = new double** [Spline.length_r];
    Spline.dfdp_slopes = new double** [Spline.length_r];
    
    for(int i = 0; i < Spline.length_r; i++){
        Spline.f_slopes[i] = new double* [Spline.length_t];
        Spline.dfdt_slopes[i] = new double* [Spline.length_t];
        Spline.dfdp_slopes[i] = new double* [Spline.length_t];
        for(int j = 0; j < Spline.length_t; j++){
            Spline.f_slopes[i][j] = new double [Spline.length_p];
            Spline.dfdt_slopes[i][j] = new double [Spline.length_p];
            Spline.dfdp_slopes[i][j] = new double [Spline.length_p];
        }
    }
}

void ClearInputArrays(struct MultiDimSpline_3DGlobal & Spline){
    delete Spline.r_vals;
    delete Spline.t_vals;
    delete Spline.p_vals;
    
    for(int i = 0; i < Spline.length_r; i++){
        for(int j = 0; j < Spline.length_t; j++){
            delete Spline.f_vals[i][j];
        }
        delete Spline.f_vals[i];
    }
    delete Spline.f_vals;
    
}

void ClearSlopesArrays(struct MultiDimSpline_3DGlobal & Spline){
    for(int i = 0; i < Spline.length_r; i++){
        for(int j = 0; j < Spline.length_t; j++){
            delete Spline.f_slopes[i][j];
            delete Spline.dfdt_slopes[i][j];
            delete Spline.dfdp_slopes[i][j];
        }
        delete Spline.f_slopes[i];
        delete Spline.dfdt_slopes[i];
        delete Spline.dfdp_slopes[i];
    }
    delete Spline.f_slopes;
    delete Spline.dfdt_slopes;
    delete Spline.dfdp_slopes;
    
}


//--------------------------------------//
//-----------Functions to Set-----------//
//-------the Interpolation Slopes-------//
//--------------------------------------//
void Set_Slopes_Multi(struct MultiDimSpline_3DGlobal & Spline){
    double ai, bi, ci, di;
    
    double new_c[Spline.length_r-1];
    double new_d[Spline.length_r];
    
    // Set f_slopes
    for(int mt = 0; mt < Spline.length_t; mt++){
    for(int mp = 0; mp < Spline.length_p; mp++){
            
        bi = 2.0 / (Spline.r_vals[1] - Spline.r_vals[0]);
        ci = 1.0 / (Spline.r_vals[1] - Spline.r_vals[0]);
        di = 3.0 * (Spline.f_vals[1][mt][mp] - Spline.f_vals[0][mt][mp]) / pow(Spline.r_vals[1] - Spline.r_vals[0], 2);
            
        new_c[0] = ci/bi;
        new_d[0] = di/bi;
            
        for(int i = 1; i < Spline.length_r - 1; i++) {
            ai = 1.0/(Spline.r_vals[i] - Spline.r_vals[i-1]);
            bi = 2.0 * (1.0/(Spline.r_vals[i] - Spline.r_vals[i-1]) + 1.0/(Spline.r_vals[i+1] - Spline.r_vals[i]));
            ci = 1.0/(Spline.r_vals[i+1] - Spline.r_vals[i]);
            di = 3.0 * ((Spline.f_vals[i][mt][mp] - Spline.f_vals[i-1][mt][mp]) / pow(Spline.r_vals[i] - Spline.r_vals[i-1], 2)
                            + (Spline.f_vals[i+1][mt][mp] - Spline.f_vals[i][mt][mp]) / pow(Spline.r_vals[i+1] - Spline.r_vals[i], 2) );
                
            new_c[i] = ci/(bi - new_c[i-1]*ai);
            new_d[i] = (di - new_d[i-1]*ai)/(bi - new_c[i-1]*ai);
        }
            
        ai = 1.0/(Spline.r_vals[Spline.length_r-1] - Spline.r_vals[Spline.length_r-2]);
        bi = 2.0/(Spline.r_vals[Spline.length_r-1] - Spline.r_vals[Spline.length_r-2]);
        di = 3.0 * (Spline.f_vals[Spline.length_r-1][mt][mp] - Spline.f_vals[Spline.length_r-2][mt][mp])
                        / pow(Spline.r_vals[Spline.length_r-1] - Spline.r_vals[Spline.length_r-2], 2);
            
        new_d[Spline.length_r-1] = (di - new_d[Spline.length_r - 2]*ai)/(bi - new_c[Spline.length_r - 2]*ai);
        
        Spline.f_slopes[Spline.length_r - 1][mt][mp] = new_d[Spline.length_r - 1];
        for(int i = Spline.length_r - 2; i >= 0; i--) Spline.f_slopes[i][mt][mp] = new_d[i] - new_c[i] * Spline.f_slopes[i+1][mt][mp];
    }}
    
    // Set df/dt and df/dp values to use for setting dfdt_slopes and dfdp_slopes
    double dfdt [Spline.length_r][Spline.length_t][Spline.length_p];
    double dfdp [Spline.length_r][Spline.length_t][Spline.length_p];
    
    for(int mt = 0; mt < Spline.length_t; mt++){
    for(int mp = 0; mp < Spline.length_p; mp++){
        int mt_up = min(mt + 1, Spline.length_t - 1);   int mt_dn = max(mt - 1, 0);
        int mp_up = min(mp + 1, Spline.length_p - 1);   int mp_dn = max(mp - 1, 0);
        for(int mr = 0; mr < Spline.length_r; mr++){
            dfdt[mr][mt][mp] = (Spline.f_vals[mr][mt_up][mp] - Spline.f_vals[mr][mt_dn][mp])/(Spline.t_vals[mt_up] - Spline.t_vals[mt_dn]);
            dfdp[mr][mt][mp] = (Spline.f_vals[mr][mt][mp_up] - Spline.f_vals[mr][mt][mp_dn])/(Spline.p_vals[mp_up] - Spline.p_vals[mp_dn]);
        }
    }}
    
    // Set dfdt_slopes
    for(int mt = 0; mt < Spline.length_t; mt++){
    for(int mp = 0; mp < Spline.length_p; mp++){
            
        bi = 2.0 / (Spline.r_vals[1] - Spline.r_vals[0]);
        ci = 1.0 / (Spline.r_vals[1] - Spline.r_vals[0]);
        di = 3.0 * (dfdt[1][mt][mp] - dfdt[0][mt][mp]) / pow(Spline.r_vals[1] - Spline.r_vals[0], 2);
            
        new_c[0] = ci/bi;
        new_d[0] = di/bi;
            
        for(int i = 1; i < Spline.length_r - 1; i++) {
            ai = 1.0/(Spline.r_vals[i] - Spline.r_vals[i-1]);
            bi = 2.0 * (1.0/(Spline.r_vals[i] - Spline.r_vals[i-1]) + 1.0/(Spline.r_vals[i+1] - Spline.r_vals[i]));
            ci = 1.0/(Spline.r_vals[i+1] - Spline.r_vals[i]);
            di = 3.0 * ((dfdt[i][mt][mp] - dfdt[i+1][mt][mp]) / pow(Spline.r_vals[i] - Spline.r_vals[i-1], 2)
                            + (dfdt[i+1][mt][mp] - dfdt[i][mt][mp]) / pow(Spline.r_vals[i+1] - Spline.r_vals[i], 2) );
                
            new_c[i] = ci/(bi - new_c[i-1]*ai);
            new_d[i] = (di - new_d[i-1]*ai)/(bi - new_c[i-1]*ai);
        }
            
        ai = 1.0/(Spline.r_vals[Spline.length_r-1] - Spline.r_vals[Spline.length_r-2]);
        bi = 2.0/(Spline.r_vals[Spline.length_r-1] - Spline.r_vals[Spline.length_r-2]);
        di = 3.0 * (dfdt[Spline.length_r-1][mt][mp] - dfdt[Spline.length_r-2][mt][mp])
                    / pow(Spline.r_vals[Spline.length_r-1] - Spline.r_vals[Spline.length_r-2], 2);
            
        new_d[Spline.length_r-1] = (di - new_d[Spline.length_r - 2]*ai)/(bi - new_c[Spline.length_r - 2]*ai);
        
        Spline.dfdt_slopes[Spline.length_r - 1][mt][mp] = new_d[Spline.length_r - 1];
        for(int i = Spline.length_r - 2; i >= 0; i--) Spline.dfdt_slopes[i][mt][mp] = new_d[i] - new_c[i] * Spline.dfdt_slopes[i+1][mt][mp];
    }}
    
    // Set dfdp_slopes
    for(int mt = 0; mt < Spline.length_t; mt++){
    for(int mp = 0; mp < Spline.length_p; mp++){
            
        bi = 2.0 / (Spline.r_vals[1] - Spline.r_vals[0]);
        ci = 1.0 / (Spline.r_vals[1] - Spline.r_vals[0]);
        di = 3.0 * (dfdp[1][mt][mp] - dfdp[0][mt][mp]) / pow(Spline.r_vals[1] - Spline.r_vals[0], 2);
            
        new_c[0] = ci/bi;
        new_d[0] = di/bi;
            
        for(int i = 1; i < Spline.length_r - 1; i++) {
            ai = 1.0/(Spline.r_vals[i] - Spline.r_vals[i-1]);
            bi = 2.0 * (1.0/(Spline.r_vals[i] - Spline.r_vals[i-1]) + 1.0/(Spline.r_vals[i+1] - Spline.r_vals[i]));
            ci = 1.0/(Spline.r_vals[i+1] - Spline.r_vals[i]);
            di = 3.0 * ((dfdp[i][mt][mp] - dfdp[i+1][mt][mp]) / pow(Spline.r_vals[i] - Spline.r_vals[i-1], 2)
                            + (dfdp[i+1][mt][mp] - dfdp[i][mt][mp]) / pow(Spline.r_vals[i+1] - Spline.r_vals[i], 2) );
            
            new_c[i] = ci/(bi - new_c[i-1]*ai);
            new_d[i] = (di - new_d[i-1]*ai)/(bi - new_c[i-1]*ai);
        }
            
        ai = 1.0/(Spline.r_vals[Spline.length_r-1] - Spline.r_vals[Spline.length_r-2]);
        bi = 2.0/(Spline.r_vals[Spline.length_r-1] - Spline.r_vals[Spline.length_r-2]);
        di = 3.0 * (dfdp[Spline.length_r-1][mt][mp] - dfdp[Spline.length_r-2][mt][mp])
                    / pow(Spline.r_vals[Spline.length_r-1] - Spline.r_vals[Spline.length_r-2], 2);
            
        new_d[Spline.length_r-1] = (di - new_d[Spline.length_r - 2]*ai)/(bi - new_c[Spline.length_r - 2]*ai);
            
        Spline.dfdp_slopes[Spline.length_r - 1][mt][mp] = new_d[Spline.length_r - 1];
        for(int i = Spline.length_r - 2; i >= 0; i--) Spline.dfdp_slopes[i][mt][mp] = new_d[i] - new_c[i] * Spline.dfdp_slopes[i+1][mt][mp];
    }}
}


//-----------------------------------------------//
//------Functions to Find the Segment Index------//
//-----and Evaluate the 1D Vertical Splines------//
//-----------------------------------------------//
int Find_Segment(double x, double* x_vals, int length, int & prev){
    int index = length + 1;
    bool done = false;
    
    if(x <= x_vals[length-1] && x >= x_vals[0]){
        // Check previous index and bounding segments
        if(x >= x_vals[prev] && x <= x_vals[prev+1]){
            done = true;
        }
        
        if(!done && prev+2 <= length-1){
            if(x >= x_vals[prev+1] && x <= x_vals[prev+2]){
                done = true;
                prev = prev + 1;
            }
        }
        if(!done && prev-1 >= 0){
            if(x >= x_vals[prev-1] && x <= x_vals[prev]){
                done = true;
                prev = prev - 1;
            }
        }
        
        if(!done){
            for(int i = 0; i < length; i++){
                if(x >= x_vals[i] && x <= x_vals[i+1]){
                    index = i;
                    break;
                }
                if (x >= x_vals[length - 2 - i] && x < x_vals[length - 1 - i]){
                    index = (length - 2) - i;
                    break;
                }
            }
            prev = index;
        }
        return prev;
    } else {
        cout << "Cannot interpolate outside of given bounds.  Value of " << setprecision(12) << x << " is outside of scope: " << setprecision(12) <<  x_vals[0] << " to " << setprecision(12) << x_vals[length-1] << "." << '\n';
    }
}


double Eval_Vert_Spline_f(double r, struct MultiDimSpline_3DGlobal Spline, int kr, int kt, int kp){
    double X = (r - Spline.r_vals[kr])/(Spline.r_vals[kr+1] - Spline.r_vals[kr]);
    double A = Spline.f_slopes[kr][kt][kp] * (Spline.r_vals[kr+1] - Spline.r_vals[kr]) - (Spline.f_vals[kr+1][kt][kp] - Spline.f_vals[kr][kt][kp]);
    double B = -Spline.f_slopes[kr+1][kt][kp] * (Spline.r_vals[kr+1] - Spline.r_vals[kr]) + (Spline.f_vals[kr+1][kt][kp] - Spline.f_vals[kr][kt][kp]);
    
    return (1.0 - X) * Spline.f_vals[kr][kt][kp] + X * Spline.f_vals[kr+1][kt][kp] + X * (1.0 - X) * (A * (1.0 - X ) + B * X);
}

double Eval_Vert_Spline_dfdr(double r, struct MultiDimSpline_3DGlobal Spline, int kr, int kt, int kp){
    double X = (r - Spline.r_vals[kr])/(Spline.r_vals[kr+1] - Spline.r_vals[kr]);
    double A = Spline.f_slopes[kr][kt][kp] * (Spline.r_vals[kr+1] - Spline.r_vals[kr]) - (Spline.f_vals[kr+1][kt][kp] - Spline.f_vals[kr][kt][kp]);
    double B = -Spline.f_slopes[kr+1][kt][kp] * (Spline.r_vals[kr+1] - Spline.r_vals[kr]) + (Spline.f_vals[kr+1][kt][kp] - Spline.f_vals[kr][kt][kp]);
    
    return (Spline.f_vals[kr+1][kt][kp] - Spline.f_vals[kr][kt][kp])/(Spline.r_vals[kr+1] - Spline.r_vals[kr])
                + (1.0 - 2.0 * X) * (A * (1.0 - X) + B * X)/(Spline.r_vals[kr+1] - Spline.r_vals[kr])
                    + X * (1.0 - X) * (B - A)/(Spline.r_vals[kr+1] - Spline.r_vals[kr]);
}

double Eval_Vert_Spline_dfdt(double r, struct MultiDimSpline_3DGlobal Spline, int kr, int kt, int kp){
    int kt_up = min(kt + 1, Spline.length_t - 1);
    int kt_dn = max(kt - 1, 0);
    
    double dfdt_kr = (Spline.f_vals[kr][kt_up][kp] - Spline.f_vals[kr][kt_dn][kp])/(Spline.t_vals[kt_up] - Spline.t_vals[kt_dn]);          // df/dt at kr
    double dfdt_krp1 = (Spline.f_vals[kr+1][kt_up][kp] - Spline.f_vals[kr+1][kt_dn][kp])/(Spline.t_vals[kt_up] - Spline.t_vals[kt_dn]);     // df/dt at kr + 1
    
    double X = (r - Spline.r_vals[kr])/(Spline.r_vals[kr+1] - Spline.r_vals[kr]);
    double A = Spline.dfdt_slopes[kr][kt][kp] * (Spline.r_vals[kr+1] - Spline.r_vals[kr]) - (dfdt_krp1 - dfdt_kr);
    double B = -Spline.dfdt_slopes[kr+1][kt][kp] * (Spline.r_vals[kr+1] - Spline.r_vals[kr]) + (dfdt_krp1 - dfdt_kr);
    
    return (1.0 - X) * dfdt_kr + X * dfdt_krp1 + X * (1.0 - X) * (A * (1.0 - X ) + B * X);
}

double Eval_Vert_Spline_dfdp(double r, struct MultiDimSpline_3DGlobal Spline, int kr, int kt, int kp){
    int kp_up = min(kp + 1, Spline.length_p - 1);
    int kp_dn = max(kp - 1, 0);
    
    double dfdp_kr = (Spline.f_vals[kr][kt][kp_up] - Spline.f_vals[kr][kt][kp_dn])/(Spline.p_vals[kp_up] - Spline.p_vals[kp_dn]);          // df/dp at kr
    double dfdp_krp1 = (Spline.f_vals[kr+1][kt][kp_up] - Spline.f_vals[kr+1][kt][kp_dn])/(Spline.p_vals[kp_up] - Spline.p_vals[kp_dn]);     // df/dp at kr + 1
    
    double X = (r - Spline.r_vals[kr])/(Spline.r_vals[kr+1] - Spline.r_vals[kr]);
    double A = Spline.dfdp_slopes[kr][kt][kp] * (Spline.r_vals[kr+1] - Spline.r_vals[kr]) - (dfdp_krp1 - dfdp_kr);
    double B = -Spline.dfdp_slopes[kr+1][kt][kp] * (Spline.r_vals[kr+1] - Spline.r_vals[kr]) + (dfdp_krp1 - dfdp_kr);
    
    return (1.0 - X) * dfdp_kr + X * dfdp_krp1 + X * (1.0 - X) * (A * (1.0 - X ) + B * X);
}

double Eval_Vert_Spline_ddfdrdr(double r, struct MultiDimSpline_3DGlobal Spline, int kr, int kt, int kp){
    double X = (r - Spline.r_vals[kr])/(Spline.r_vals[kr+1] - Spline.r_vals[kr]);
    double A = Spline.f_slopes[kr][kt][kp] * (Spline.r_vals[kr+1] - Spline.r_vals[kr]) - (Spline.f_vals[kr+1][kt][kp] - Spline.f_vals[kr][kt][kp]);
    double B = -Spline.f_slopes[kr+1][kt][kp] * (Spline.r_vals[kr+1] - Spline.r_vals[kr]) + (Spline.f_vals[kr+1][kt][kp] - Spline.f_vals[kr][kt][kp]);
    
    return 2.0 * (B - 2.0 * A + (A - B) * 3.0 * X)/pow(Spline.r_vals[kr+1] - Spline.r_vals[kr],2);
}

double Eval_Vert_Spline_ddfdrdt(double r, struct MultiDimSpline_3DGlobal Spline, int kr, int kt, int kp){
    int kt_up = min(kt + 1, Spline.length_t - 1);
    int kt_dn = max(kt - 1, 0);
    
    double dfdt_kr = (Spline.f_vals[kr][kt_up][kp] - Spline.f_vals[kr][kt_dn][kp])/(Spline.t_vals[kt_up] - Spline.t_vals[kt_dn]);          // df/dt at kr
    double dfdt_krp1 = (Spline.f_vals[kr+1][kt_up][kp] - Spline.f_vals[kr+1][kt_dn][kp])/(Spline.t_vals[kt_up] - Spline.t_vals[kt_dn]);     // df/dt at kr + 1
    
    double X = (r - Spline.r_vals[kr])/(Spline.r_vals[kr+1] - Spline.r_vals[kr]);
    double A = Spline.dfdt_slopes[kr][kt][kp] * (Spline.r_vals[kr+1] - Spline.r_vals[kr]) - (dfdt_krp1 - dfdt_kr);
    double B = -Spline.dfdt_slopes[kr+1][kt][kp] * (Spline.r_vals[kr+1] - Spline.r_vals[kr]) + (dfdt_krp1 - dfdt_kr);
    
    return (dfdt_krp1 - dfdt_krp1)/(Spline.r_vals[kr+1] - Spline.r_vals[kr])
                + (1.0 - 2.0 * X) * (A * (1.0 - X) + B * X)/(Spline.r_vals[kr+1] - Spline.r_vals[kr])
                    + X * (1.0 - X) * (B - A)/(Spline.r_vals[kr+1] - Spline.r_vals[kr]);
}

double Eval_Vert_Spline_ddfdrdp(double r, struct MultiDimSpline_3DGlobal Spline, int kr, int kt, int kp){
    int kp_up = min(kp + 1, Spline.length_p - 1);
    int kp_dn = max(kp - 1, 0);
    
    double dfdp_kr = (Spline.f_vals[kr][kt][kp_up] - Spline.f_vals[kr][kt][kp_dn])/(Spline.p_vals[kp_up] - Spline.p_vals[kp_dn]);          // df/dp at kr
    double dfdp_krp1 = (Spline.f_vals[kr+1][kt][kp_up] - Spline.f_vals[kr+1][kt][kp_dn])/(Spline.p_vals[kp_up] - Spline.p_vals[kp_dn]);     // df/dp at kr + 1
    
    double X = (r - Spline.r_vals[kr])/(Spline.r_vals[kr+1] - Spline.r_vals[kr]);
    double A = Spline.dfdp_slopes[kr][kt][kp] * (Spline.r_vals[kr+1] - Spline.r_vals[kr]) - (dfdp_krp1 - dfdp_kr);
    double B = -Spline.dfdp_slopes[kr+1][kt][kp] * (Spline.r_vals[kr+1] - Spline.r_vals[kr]) + (dfdp_krp1 - dfdp_kr);
    
    return (dfdp_krp1 - dfdp_krp1)/(Spline.r_vals[kr+1] - Spline.r_vals[kr])
                + (1.0 - 2.0 * X) * (A * (1.0 - X) + B * X)/(Spline.r_vals[kr+1] - Spline.r_vals[kr])
                    + X * (1.0 - X) * (B - A)/(Spline.r_vals[kr+1] - Spline.r_vals[kr]);
}

//----------------------------------------------------//
//-----------Evauation of the 2D Spline and-----------//
//-------its First and Second Order Derivatives-------//
//----------------------------------------------------//
double BiCubic_Deriv_dfdt(double r, struct MultiDimSpline_3DGlobal Spline, int kr, int kt, int kp){
    // Calculate discrete d/dt of f to obtain df/dt
    int kt_up = kt + 1;
    int kt_dn = kt - 1;
    
	if(kt_up > Spline.length_t - 1) kt_up = kt;           // Use backward t derivative if kt+1 is not defined
    if(kt_dn < 0)                   kt_dn = kt;           // Use forward t derivative if kt-1 is not defined
    
    return (Eval_Vert_Spline_f(r, Spline, kr, kt_up, kp) - Eval_Vert_Spline_f(r, Spline, kr, kt_dn, kp))
                /(Spline.t_vals[kt_up] - Spline.t_vals[kt_dn]);
}

double BiCubic_Deriv_dfdp(double r, struct MultiDimSpline_3DGlobal Spline, int kr, int kt, int kp){
    // Calculate discrete d/dp of f to obtain df/dp
    int kp_up = kp + 1;
    int kp_dn = kp - 1;
    
	if(kp_up > Spline.length_p - 1) kp_up = kp;           // Use backward p derivative if kp+1 is not defined
    if(kp_dn < 0)                   kp_dn = kp;           // Use forward p derivative if kp-1 is not defined
    
    return (Eval_Vert_Spline_f(r, Spline, kr, kt, kp_up) - Eval_Vert_Spline_f(r, Spline, kr, kt, kp_dn))
                /(Spline.p_vals[kp_up] - Spline.p_vals[kp_dn]);
    
}

double BiCubic_Deriv_ddfdtdt(double r, struct MultiDimSpline_3DGlobal Spline, int kr, int kt, int kp){
    // Calculate discrete d/dt of df/dt to obtain d^2f/dt^2
    int kt_up = kt + 1;
    int kt_dn = kt - 1;
    
	if(kt_up > Spline.length_t - 1) kt_up = kt;           // Use backward t derivative if kt+1 is not defined
    if(kt_dn < 0)                   kt_dn = kt;           // Use forward t derivative if kt-1 is not defined
    
    return (Eval_Vert_Spline_dfdt(r, Spline, kr, kt_up, kp) - Eval_Vert_Spline_dfdt(r, Spline, kr, kt_dn, kp))
                /(Spline.t_vals[kt_up] - Spline.t_vals[kt_dn]);
}

double BiCubic_Deriv_ddfdpdp(double r, struct MultiDimSpline_3DGlobal Spline, int kr, int kt, int kp){
    // Calculate discrete d/dp of df/dp to obtain d^2f/dp^2
    int kp_up = kp + 1;
    int kp_dn = kp - 1;
    
	if(kp_up > Spline.length_p - 1) kp_up = kp;           // Use backward p derivative if ky+1 is not defined
    if(kp_dn < 0)                   kp_dn = kp;           // Use forward p derivative if ky-1 is not defined
    
    return (Eval_Vert_Spline_dfdp(r, Spline, kr, kt, kp_up) - Eval_Vert_Spline_dfdp(r, Spline, kr, kt, kp_dn))
                /(Spline.p_vals[kp_up] - Spline.p_vals[kp_dn]);
    
}

double BiCubic_Deriv_ddfdtdp(double r, struct MultiDimSpline_3DGlobal Spline, int kr, int kt, int kp){
    // Calculate discrete d^2/dtdp of f to obtain d^2f/dtdp
    int kt_up = kt+1;    int kt_dn = kt-1;
    int kp_up = kp+1;    int kp_dn = kp-1;
    
	if(kt_up > Spline.length_t - 1) kt_up = kt;           // Use backward t derivative if kt+1 is not defined
    if(kt_dn < 0)                   kt_dn = kt;           // Use forward t derivative if kt-1 is not defined
    
	if(kp_up > Spline.length_p - 1) kp_up = kp;           // Use backward p derivative if kp+1 is not defined
    if(kp_dn < 0)                   kp_dn = kp;           // Use forward p derivative if kp-1 is not defined
    
    return (Eval_Vert_Spline_f(r, Spline, kr, kt_up, kp_up)
            - Eval_Vert_Spline_f(r, Spline, kr, kt_up, kp_dn)
                - Eval_Vert_Spline_f(r, Spline, kr, kt_dn, kp_up)
                    + Eval_Vert_Spline_f(r, Spline, kr, kt_dn, kp_dn))
                        /((Spline.t_vals[kt_up] - Spline.t_vals[kt_dn])*(Spline.p_vals[kp_up] - Spline.p_vals[kp_dn]));
}


double BiCubic_Deriv_dddfdrdrdt(double r, struct MultiDimSpline_3DGlobal Spline, int kr, int kt, int kp){
    // Calculate discrete d/dt of d^2f/dr^2 to obtain d^3f/dr^2dt
    int kt_up = kt + 1;
    int kt_dn = kt - 1;
    
	if(kt_up > Spline.length_t - 1) kt_up = kt;           // Use backward t derivative if kt+1 is not defined
    if(kt_dn < 0)                   kt_dn = kt;           // Use forward t derivative if kt-1 is not defined
    
    return (Eval_Vert_Spline_ddfdrdr(r, Spline, kr, kt_up, kp) - Eval_Vert_Spline_ddfdrdr(r, Spline, kr, kt_dn, kp))
                /(Spline.t_vals[kt_up] - Spline.t_vals[kt_dn]);
}

double BiCubic_Deriv_dddfdrdrdp(double r, struct MultiDimSpline_3DGlobal Spline, int kr, int kt, int kp){
    // Calculate discrete d/dp of d^2f/dr^2 to obtain d^3f/dr^2dp
    int kp_up = kp + 1;
    int kp_dn = kp - 1;
    
	if(kp_up > Spline.length_p - 1) kp_up = kp;           // Use backward t derivative if kt+1 is not defined
    if(kp_dn < 0)                   kp_dn = kp;           // Use forward t derivative if kt-1 is not defined
    
    return (Eval_Vert_Spline_ddfdrdr(r, Spline, kr, kt, kp_up) - Eval_Vert_Spline_ddfdrdr(r, Spline, kr, kt, kp_dn))
            /(Spline.p_vals[kp_up] - Spline.p_vals[kp_dn]);
}

double BiCubic_Deriv_dddfdrdtdp(double r, struct MultiDimSpline_3DGlobal Spline, int kr, int kt, int kp){
    // Calculate discrete d^2/dtdp of df/dr to obtain d^3f/drdtdp
    int kt_up = kt+1;    int kt_dn = kt-1;
    int kp_up = kp+1;    int kp_dn = kp-1;
    
	if(kt_up > Spline.length_t - 1) kt_up = kt;           // Use backward t derivative if kt+1 is not defined
    if(kt_dn < 0)                   kt_dn = kt;           // Use forward t derivative if kt-1 is not defined
    
	if(kp_up > Spline.length_p - 1) kp_up = kp;           // Use backward p derivative if kp+1 is not defined
    if(kp_dn < 0)                   kp_dn = kp;           // Use forward p derivative if kp-1 is not defined
    
    return (Eval_Vert_Spline_dfdr(r, Spline, kr, kt_up, kp_up)
            - Eval_Vert_Spline_dfdr(r, Spline, kr, kt_up, kp_dn)
                - Eval_Vert_Spline_dfdr(r, Spline, kr, kt_dn, kp_up)
                    + Eval_Vert_Spline_dfdr(r, Spline, kr, kt_dn, kp_dn))
                        /((Spline.t_vals[kt_up] - Spline.t_vals[kt_dn])*(Spline.p_vals[kp_up] - Spline.p_vals[kp_dn]));
}


double BiCubic_Deriv_dddfdrdpdp(double r, struct MultiDimSpline_3DGlobal Spline, int kr, int kt, int kp){
    // Calculate discrete d/dp of d^2f/drdr to obtain d^3f/dy^2dz
    int kt_up = kt + 1;
    int kt_dn = kt - 1;
    
	if(kt_up > Spline.length_t - 1) kt_up = kt;           // Use backward t derivative if kt+1 is not defined
    if(kt_dn < 0)                   kt_dn = kt;           // Use forward t derivative if kt-1 is not defined
    
    return (Eval_Vert_Spline_ddfdrdr(r, Spline, kr, kt_up, kp) - Eval_Vert_Spline_ddfdrdr(r, Spline, kr, kt_dn, kp))
                /(Spline.t_vals[kt_up] - Spline.t_vals[kt_dn]);
    
}


double BiCubic_Deriv_dddfdtdtdp(double r, struct MultiDimSpline_3DGlobal Spline, int kr, int kt, int kp){
    // Calculate discrete d^2/dtdp of df/dt to obtain d^3f/dt^2dp
    int kt_up = kt+1;    int kt_dn = kt-1;
    int kp_up = kp+1;    int kp_dn = kp-1;
    
	if(kt_up > Spline.length_t - 1) kt_up = kt;           // Use backward t derivative if kt+1 is not defined
    if(kt_dn < 0)                   kt_dn = kt;           // Use forward t derivative if kt-1 is not defined
    
	if(kp_up > Spline.length_p - 1) kp_up = kp;           // Use backward p derivative if kp+1 is not defined
    if(kp_dn < 0)                   kp_dn = kp;           // Use forward p derivative if kp-1 is not defined
    
    return (Eval_Vert_Spline_dfdt(r, Spline, kr, kt_up, kp_up)
            - Eval_Vert_Spline_dfdt(r, Spline, kr, kt_up, kp_dn)
                - Eval_Vert_Spline_dfdt(r, Spline, kr, kt_dn, kp_up)
                    + Eval_Vert_Spline_dfdt(r, Spline, kr, kt_dn, kp_dn))
                        /((Spline.t_vals[kt_up] - Spline.t_vals[kt_dn])*(Spline.p_vals[kp_up] - Spline.p_vals[kp_dn]));
}

double BiCubic_Deriv_dddfdtdpdp(double r, struct MultiDimSpline_3DGlobal Spline, int kr, int kt, int kp){
    // Calculate discrete d^2/dtdp of df/dp to obtain d^3f/dtdp^2
    int kt_up = kt+1;    int kt_dn = kt-1;
    int kp_up = kp+1;    int kp_dn = kp-1;
    
	if(kt_up > Spline.length_t - 1) kt_up = kt;           // Use backward t derivative if kt+1 is not defined
    if(kt_dn < 0)                   kt_dn = kt;           // Use forward t derivative if kt-1 is not defined
    
	if(kp_up > Spline.length_p - 1) kp_up = kp;           // Use backward p derivative if kp+1 is not defined
    if(kp_dn < 0)                   kp_dn = kp;           // Use forward p derivative if kp-1 is not defined
    
    return (Eval_Vert_Spline_dfdp(r, Spline, kr, kt_up, kp_up)
            - Eval_Vert_Spline_dfdp(r, Spline, kr, kt_up, kp_dn)
                - Eval_Vert_Spline_dfdp(r, Spline, kr, kt_dn, kp_up)
                    + Eval_Vert_Spline_dfdp(r, Spline, kr, kt_dn, kp_dn))
                        /((Spline.t_vals[kt_up] - Spline.t_vals[kt_dn])*(Spline.p_vals[kp_up] - Spline.p_vals[kp_dn]));
}

double BiCubic_Deriv_ddddfdrdrdtdp(double r, struct MultiDimSpline_3DGlobal Spline, int kr, int kt, int kp){
    // Calculate discrete d^2/dtdp of d^2f/dr^2 to obtain d^4f/dr^2dtdp
    int kt_up = kt+1;    int kt_dn = kt-1;
    int kp_up = kp+1;    int kp_dn = kp-1;
    
	if(kt_up > Spline.length_t - 1) kt_up = kt;           // Use backward t derivative if kt+1 is not defined
    if(kt_dn < 0)                   kt_dn = kt;           // Use forward t derivative if kt-1 is not defined
    
	if(kp_up > Spline.length_p - 1) kp_up = kp;           // Use backward p derivative if kp+1 is not defined
    if(kp_dn < 0)                   kp_dn = kp;           // Use forward p derivative if kp-1 is not defined
    
    return (Eval_Vert_Spline_ddfdrdr(r, Spline, kr, kt_up, kp_up)
            - Eval_Vert_Spline_ddfdrdr(r, Spline, kr, kt_up, kp_dn)
                - Eval_Vert_Spline_ddfdrdr(r, Spline, kr, kt_dn, kp_up)
                    + Eval_Vert_Spline_ddfdrdr(r, Spline, kr, kt_dn, kp_dn))
                        /((Spline.t_vals[kt_up] - Spline.t_vals[kt_dn])*(Spline.p_vals[kp_up] - Spline.p_vals[kp_dn]));
}



//----------------------------------------------------//
//-----------Evauation of the 2D Spline and-----------//
//-------its First and Second Order Derivatives-------//
//----------------------------------------------------//
double Eval_Spline_f(double r, double t, double p, struct MultiDimSpline_3DGlobal & Spline){
    double A_vec[16];
	double X_vec[16];
    
    // Locate kr, kt, kp such that (r[kr] <= r' <= r[kr+1]), (t[kt] <= t' <= t[kt+1]), and (p[kp] <= p' <= p[kp+1])
    int kr = Find_Segment(r, Spline.r_vals, Spline.length_r, Spline.accel[0]);
    int kt = Find_Segment(t, Spline.t_vals, Spline.length_t, Spline.accel[1]);
    int kp = Find_Segment(p, Spline.p_vals, Spline.length_p, Spline.accel[2]);
    
    // Set dt_scalar, dp_scalar, and t_scaled, p_scaled to equate to unit square
    double dt_scalar = Spline.t_vals[kt+1] - Spline.t_vals[kt];
    double dp_scalar = Spline.p_vals[kp+1] - Spline.p_vals[kp];
    
    double t_scaled = (t - Spline.t_vals[kt])/dt_scalar;
	double p_scaled = (p - Spline.p_vals[kp])/dp_scalar;
    
    X_vec[0] = Eval_Vert_Spline_f(r, Spline, kr, kt, kp);
    X_vec[1] = Eval_Vert_Spline_f(r, Spline, kr, kt+1, kp);
    X_vec[2] = Eval_Vert_Spline_f(r, Spline, kr, kt, kp+1);
    X_vec[3] = Eval_Vert_Spline_f(r, Spline, kr, kt+1, kp+1);
    
    X_vec[4] = BiCubic_Deriv_dfdt(r, Spline, kr, kt, kp)*dt_scalar;
    X_vec[5] = BiCubic_Deriv_dfdt(r, Spline, kr, kt+1, kp)*dt_scalar;
    X_vec[6] = BiCubic_Deriv_dfdt(r, Spline, kr, kt, kp+1)*dt_scalar;
    X_vec[7] = BiCubic_Deriv_dfdt(r, Spline, kr, kt+1, kp+1)*dt_scalar;
    
    X_vec[8] =  BiCubic_Deriv_dfdp(r, Spline, kr, kt, kp)*dp_scalar;
    X_vec[9] =  BiCubic_Deriv_dfdp(r, Spline, kr, kt+1, kp)*dp_scalar;
    X_vec[10] = BiCubic_Deriv_dfdp(r, Spline, kr, kt, kp+1)*dp_scalar;
    X_vec[11] = BiCubic_Deriv_dfdp(r, Spline, kr, kt+1, kp+1)*dp_scalar;
    
    X_vec[12]   = BiCubic_Deriv_ddfdtdp(r, Spline, kr, kt, kp)*dt_scalar*dp_scalar;
    X_vec[13]   = BiCubic_Deriv_ddfdtdp(r, Spline, kr, kt+1, kp)*dt_scalar*dp_scalar;
    X_vec[14]	= BiCubic_Deriv_ddfdtdp(r, Spline, kr, kt, kp+1)*dt_scalar*dp_scalar;
    X_vec[15]   = BiCubic_Deriv_ddfdtdp(r, Spline, kr, kt+1, kp+1)*dt_scalar*dp_scalar;
    
    // Set the coefficients by matrix multiplication
    for(int j = 0; j < 16; j++){
        A_vec[j] = 0;
        for(int k = 0; k < 16; k++){
            A_vec[j] += BiCubic_ConversionMatrix[j][k]*X_vec[k];
        }
    }
    
	double result = 0;
	for(int k1 = 0; k1 < 4; k1++){
    for(int k2 = 0; k2 < 4; k2++){
        result+=1.0*A_vec[k2*4 + k1]*pow(t_scaled,k1)*pow(p_scaled,k2);
    }}
	
    return result;
}

double Eval_Spline_df(double r, double t, double p, int index, struct MultiDimSpline_3DGlobal & Spline){
    double A_vec[16];
	double X_vec[16];
    
    // Locate kr, kt, kp such that (r[kr] <= r' <= r[kr+1]), (t[kt] <= t' <= t[kt+1]), and (p[kp] <= p' <= p[kp+1])
    int kr = Find_Segment(r, Spline.r_vals, Spline.length_r, Spline.accel[0]);
    int kt = Find_Segment(t, Spline.t_vals, Spline.length_t, Spline.accel[1]);
    int kp = Find_Segment(p, Spline.p_vals, Spline.length_p, Spline.accel[2]);
    
    // Set dt_scalar, dp_scalar, and t_scaled, p_scaled to equate to unit square
    double dt_scalar = Spline.t_vals[kt+1] - Spline.t_vals[kt];
    double dp_scalar = Spline.p_vals[kp+1] - Spline.p_vals[kp];
    
    double t_scaled = (t - Spline.t_vals[kt])/dt_scalar;
	double p_scaled = (p - Spline.p_vals[kp])/dp_scalar;
    
    if(index==0){
        X_vec[0] = Eval_Vert_Spline_dfdr(r, Spline, kr, kt, kp);
        X_vec[1] = Eval_Vert_Spline_dfdr(r, Spline, kr, kt+1, kp);
        X_vec[2] = Eval_Vert_Spline_dfdr(r, Spline, kr, kt, kp+1);
        X_vec[3] = Eval_Vert_Spline_dfdr(r, Spline, kr, kt+1, kp+1);
        
        X_vec[4] = Eval_Vert_Spline_ddfdrdt(r, Spline, kr, kt, kp)*dt_scalar;
        X_vec[5] = Eval_Vert_Spline_ddfdrdt(r, Spline, kr, kt+1, kp)*dt_scalar;
        X_vec[6] = Eval_Vert_Spline_ddfdrdt(r, Spline, kr, kt, kp+1)*dt_scalar;
        X_vec[7] = Eval_Vert_Spline_ddfdrdt(r, Spline, kr, kt+1, kp+1)*dt_scalar;
        
        X_vec[8] =  Eval_Vert_Spline_ddfdrdp(r, Spline, kr, kt, kp)*dp_scalar;
        X_vec[9] =  Eval_Vert_Spline_ddfdrdp(r, Spline, kr, kt+1, kp)*dp_scalar;
        X_vec[10] = Eval_Vert_Spline_ddfdrdp(r, Spline, kr, kt, kp+1)*dp_scalar;
        X_vec[11] = Eval_Vert_Spline_ddfdrdp(r, Spline, kr, kt+1, kp+1)*dp_scalar;
        
        X_vec[12]   = BiCubic_Deriv_dddfdrdtdp(r, Spline, kr, kt, kp)*dt_scalar*dp_scalar;
        X_vec[13]   = BiCubic_Deriv_dddfdrdtdp(r, Spline, kr, kt+1, kp)*dt_scalar*dp_scalar;
        X_vec[14]	= BiCubic_Deriv_dddfdrdtdp(r, Spline, kr, kt, kp+1)*dt_scalar*dp_scalar;
        X_vec[15]   = BiCubic_Deriv_dddfdrdtdp(r, Spline, kr, kt+1, kp+1)*dt_scalar*dp_scalar;

    } else if(index==1){
        X_vec[0] = BiCubic_Deriv_dfdt(r, Spline, kr, kt, kp);
        X_vec[1] = BiCubic_Deriv_dfdt(r, Spline, kr, kt+1, kp);
        X_vec[2] = BiCubic_Deriv_dfdt(r, Spline, kr, kt, kp+1);
        X_vec[3] = BiCubic_Deriv_dfdt(r, Spline, kr, kt+1, kp+1);
        
        X_vec[4] = BiCubic_Deriv_ddfdtdt(r, Spline, kr, kt, kp)*dt_scalar;
        X_vec[5] = BiCubic_Deriv_ddfdtdt(r, Spline, kr, kt+1, kp)*dt_scalar;
        X_vec[6] = BiCubic_Deriv_ddfdtdt(r, Spline, kr, kt, kp+1)*dt_scalar;
        X_vec[7] = BiCubic_Deriv_ddfdtdt(r, Spline, kr, kt+1, kp+1)*dt_scalar;
        
        X_vec[8] =  BiCubic_Deriv_ddfdtdp(r, Spline, kr, kt, kp)*dp_scalar;
        X_vec[9] =  BiCubic_Deriv_ddfdtdp(r, Spline, kr, kt+1, kp)*dp_scalar;
        X_vec[10] = BiCubic_Deriv_ddfdtdp(r, Spline, kr, kt, kp+1)*dp_scalar;
        X_vec[11] = BiCubic_Deriv_ddfdtdp(r, Spline, kr, kt+1, kp+1)*dp_scalar;
        
        X_vec[12]   = BiCubic_Deriv_dddfdtdtdp(r, Spline, kr, kt, kp)*dt_scalar*dp_scalar;
        X_vec[13]   = BiCubic_Deriv_dddfdtdtdp(r, Spline, kr, kt+1, kp)*dt_scalar*dp_scalar;
        X_vec[14]	= BiCubic_Deriv_dddfdtdtdp(r, Spline, kr, kt, kp+1)*dt_scalar*dp_scalar;
        X_vec[15]   = BiCubic_Deriv_dddfdtdtdp(r, Spline, kr, kt+1, kp+1)*dt_scalar*dp_scalar;
        
    } else if(index==2){
        X_vec[0] = BiCubic_Deriv_dfdp(r, Spline, kr, kt, kp);
        X_vec[1] = BiCubic_Deriv_dfdp(r, Spline, kr, kt+1, kp);
        X_vec[2] = BiCubic_Deriv_dfdp(r, Spline, kr, kt, kp+1);
        X_vec[3] = BiCubic_Deriv_dfdp(r, Spline, kr, kt+1, kp+1);
        
        X_vec[4] = BiCubic_Deriv_ddfdtdp(r, Spline, kr, kt, kp)*dt_scalar;
        X_vec[5] = BiCubic_Deriv_ddfdtdp(r, Spline, kr, kt+1, kp)*dt_scalar;
        X_vec[6] = BiCubic_Deriv_ddfdtdp(r, Spline, kr, kt, kp+1)*dt_scalar;
        X_vec[7] = BiCubic_Deriv_ddfdtdp(r, Spline, kr, kt+1, kp+1)*dt_scalar;
        
        X_vec[8] =  BiCubic_Deriv_ddfdpdp(r, Spline, kr, kt, kp)*dp_scalar;
        X_vec[9] =  BiCubic_Deriv_ddfdpdp(r, Spline, kr, kt+1, kp)*dp_scalar;
        X_vec[10] = BiCubic_Deriv_ddfdpdp(r, Spline, kr, kt, kp+1)*dp_scalar;
        X_vec[11] = BiCubic_Deriv_ddfdpdp(r, Spline, kr, kt+1, kp+1)*dp_scalar;
        
        X_vec[12]   = BiCubic_Deriv_dddfdtdpdp(r, Spline, kr, kt, kp)*dt_scalar*dp_scalar;
        X_vec[13]   = BiCubic_Deriv_dddfdtdpdp(r, Spline, kr, kt+1, kp)*dt_scalar*dp_scalar;
        X_vec[14]	= BiCubic_Deriv_dddfdtdpdp(r, Spline, kr, kt, kp+1)*dt_scalar*dp_scalar;
        X_vec[15]   = BiCubic_Deriv_dddfdtdpdp(r, Spline, kr, kt+1, kp+1)*dt_scalar*dp_scalar;
    }
    
    // Set the coefficients by matrix multiplication
    for(int j = 0; j < 16; j++){
        A_vec[j] = 0;
        for(int k = 0; k < 16; k++){
            A_vec[j] += BiCubic_ConversionMatrix[j][k]*X_vec[k];
        }
    }
    
	double result = 0;
    for(int k1 = 0; k1 < 4; k1++){
    for(int k2 = 0; k2 < 4; k2++){
        result+=1.0*A_vec[k2*4 + k1]*pow(t_scaled,k1)*pow(p_scaled,k2);
    }}
    
    return result;
}


double Eval_Spline_ddf(double r, double t, double p, int index1, int index2, struct MultiDimSpline_3DGlobal & Spline){
    double A_vec[16];
	double X_vec[16];
    
    // Locate kr, kt, kp such that (r[kr] <= r' <= r[kr+1]), (t[kt] <= t' <= t[kt+1]), and (p[kp] <= p' <= p[kp+1])
    int kr = Find_Segment(r, Spline.r_vals, Spline.length_r, Spline.accel[0]);
    int kt = Find_Segment(t, Spline.t_vals, Spline.length_t, Spline.accel[1]);
    int kp = Find_Segment(p, Spline.p_vals, Spline.length_p, Spline.accel[2]);
    
    // Set dt_scalar, dp_scalar, and t_scaled, p_scaled to equate to unit square
    double dt_scalar = Spline.t_vals[kt+1] - Spline.t_vals[kt];
    double dp_scalar = Spline.p_vals[kp+1] - Spline.p_vals[kp];
    
    double t_scaled = (t - Spline.t_vals[kt])/dt_scalar;
	double p_scaled = (p - Spline.p_vals[kp])/dp_scalar;
    
    if(index1==0 && index2==0){
        // To calculate d^2 f/ dr^2, bicubic spline the second r derivative of f
        X_vec[0] = Eval_Vert_Spline_ddfdrdr(r, Spline, kr, kt, kp);
        X_vec[1] = Eval_Vert_Spline_ddfdrdr(r, Spline, kr, kt+1, kp);
        X_vec[2] = Eval_Vert_Spline_ddfdrdr(r, Spline, kr, kt, kp+1);
        X_vec[3] = Eval_Vert_Spline_ddfdrdr(r, Spline, kr, kt+1, kp+1);
        
        X_vec[4] = BiCubic_Deriv_dddfdrdrdt(r, Spline, kr, kt, kp)*dt_scalar;
        X_vec[5] = BiCubic_Deriv_dddfdrdrdt(r, Spline, kr, kt+1, kp)*dt_scalar;
        X_vec[6] = BiCubic_Deriv_dddfdrdrdt(r, Spline, kr, kt, kp+1)*dt_scalar;
        X_vec[7] = BiCubic_Deriv_dddfdrdrdt(r, Spline, kr, kt+1, kp+1)*dt_scalar;
        
        X_vec[8] =  BiCubic_Deriv_dddfdrdrdp(r, Spline, kr, kt, kp)*dp_scalar;
        X_vec[9] =  BiCubic_Deriv_dddfdrdrdp(r, Spline, kr, kt+1, kp)*dp_scalar;
        X_vec[10] = BiCubic_Deriv_dddfdrdrdp(r, Spline, kr, kt, kp+1)*dp_scalar;
        X_vec[11] = BiCubic_Deriv_dddfdrdrdp(r, Spline, kr, kt+1, kp+1)*dp_scalar;
        
        X_vec[12]   = BiCubic_Deriv_ddddfdrdrdtdp(r, Spline, kr, kt, kp)*dt_scalar*dp_scalar;
        X_vec[13]   = BiCubic_Deriv_ddddfdrdrdtdp(r, Spline, kr, kt+1, kp)*dt_scalar*dp_scalar;
        X_vec[14]	= BiCubic_Deriv_ddddfdrdrdtdp(r, Spline, kr, kt, kp+1)*dt_scalar*dp_scalar;
        X_vec[15]   = BiCubic_Deriv_ddddfdrdrdtdp(r, Spline, kr, kt+1, kp+1)*dt_scalar*dp_scalar;
        
    } else if((index1==1 && index2==1) || (index1==1 && index2==2) || (index1==2 && index2==1)){
        // To calculate d^2 f/ dt^2 or d^2 f/ dtdp, bicubic spline df/dt and take a t or p derivative respectively (below)
        X_vec[0] = BiCubic_Deriv_dfdt(r, Spline, kr, kt, kp);
        X_vec[1] = BiCubic_Deriv_dfdt(r, Spline, kr, kt+1, kp);
        X_vec[2] = BiCubic_Deriv_dfdt(r, Spline, kr, kt, kp+1);
        X_vec[3] = BiCubic_Deriv_dfdt(r, Spline, kr, kt+1, kp+1);
        
        X_vec[4] = BiCubic_Deriv_ddfdtdt(r, Spline, kr, kt, kp)*dt_scalar;
        X_vec[5] = BiCubic_Deriv_ddfdtdt(r, Spline, kr, kt+1, kp)*dt_scalar;
        X_vec[6] = BiCubic_Deriv_ddfdtdt(r, Spline, kr, kt, kp+1)*dt_scalar;
        X_vec[7] = BiCubic_Deriv_ddfdtdt(r, Spline, kr, kt+1, kp+1)*dt_scalar;
        
        X_vec[8] =  BiCubic_Deriv_ddfdtdp(r, Spline, kr, kt, kp)*dp_scalar;
        X_vec[9] =  BiCubic_Deriv_ddfdtdp(r, Spline, kr, kt+1, kp)*dp_scalar;
        X_vec[10] = BiCubic_Deriv_ddfdtdp(r, Spline, kr, kt, kp+1)*dp_scalar;
        X_vec[11] = BiCubic_Deriv_ddfdtdp(r, Spline, kr, kt+1, kp+1)*dp_scalar;
        
        X_vec[12]   = BiCubic_Deriv_dddfdtdtdp(r, Spline, kr, kt, kp)*dt_scalar*dp_scalar;
        X_vec[13]   = BiCubic_Deriv_dddfdtdtdp(r, Spline, kr, kt+1, kp)*dt_scalar*dp_scalar;
        X_vec[14]	= BiCubic_Deriv_dddfdtdtdp(r, Spline, kr, kt, kp+1)*dt_scalar*dp_scalar;
        X_vec[15]   = BiCubic_Deriv_dddfdtdtdp(r, Spline, kr, kt+1, kp+1)*dt_scalar*dp_scalar;
        
    } else if(index1==2 && index2==2){
     // To calculate d^2 f/ dp^2, bicubic spline df/dp and take a p derivative (below)
        X_vec[0] = BiCubic_Deriv_dfdp(r, Spline, kr, kt, kp);
        X_vec[1] = BiCubic_Deriv_dfdp(r, Spline, kr, kt+1, kp);
        X_vec[2] = BiCubic_Deriv_dfdp(r, Spline, kr, kt, kp+1);
        X_vec[3] = BiCubic_Deriv_dfdp(r, Spline, kr, kt+1, kp+1);
        
        X_vec[4] = BiCubic_Deriv_ddfdtdp(r, Spline, kr, kt, kp)*dt_scalar;
        X_vec[5] = BiCubic_Deriv_ddfdtdp(r, Spline, kr, kt+1, kp)*dt_scalar;
        X_vec[6] = BiCubic_Deriv_ddfdtdp(r, Spline, kr, kt, kp+1)*dt_scalar;
        X_vec[7] = BiCubic_Deriv_ddfdtdp(r, Spline, kr, kt+1, kp+1)*dt_scalar;
        
        X_vec[8] =  BiCubic_Deriv_ddfdpdp(r, Spline, kr, kt, kp)*dp_scalar;
        X_vec[9] =  BiCubic_Deriv_ddfdpdp(r, Spline, kr, kt+1, kp)*dp_scalar;
        X_vec[10] = BiCubic_Deriv_ddfdpdp(r, Spline, kr, kt, kp+1)*dp_scalar;
        X_vec[11] = BiCubic_Deriv_ddfdpdp(r, Spline, kr, kt+1, kp+1)*dp_scalar;
        
        X_vec[12]   = BiCubic_Deriv_dddfdtdpdp(r, Spline, kr, kt, kp)*dt_scalar*dp_scalar;
        X_vec[13]   = BiCubic_Deriv_dddfdtdpdp(r, Spline, kr, kt+1, kp)*dt_scalar*dp_scalar;
        X_vec[14]	= BiCubic_Deriv_dddfdtdpdp(r, Spline, kr, kt, kp+1)*dt_scalar*dp_scalar;
        X_vec[15]   = BiCubic_Deriv_dddfdtdpdp(r, Spline, kr, kt+1, kp+1)*dt_scalar*dp_scalar;
     
    } else if((index1==0 && index2==1) || (index1==1 && index2==0) || (index1==0 && index2==2) || (index1==2 && index2==0)){
        // To calcualte d^2 f/drdt or d^2 f / drdp, bicubic spline df/dr and take a t or p derivative respectively (below)
        X_vec[0] = Eval_Vert_Spline_dfdr(r, Spline, kr, kt, kp);
        X_vec[1] = Eval_Vert_Spline_dfdr(r, Spline, kr, kt+1, kp);
        X_vec[2] = Eval_Vert_Spline_dfdr(r, Spline, kr, kt, kp+1);
        X_vec[3] = Eval_Vert_Spline_dfdr(r, Spline, kr, kt+1, kp+1);
        
        X_vec[4] = Eval_Vert_Spline_ddfdrdt(r, Spline, kr, kt, kp)*dt_scalar;
        X_vec[5] = Eval_Vert_Spline_ddfdrdt(r, Spline, kr, kt+1, kp)*dt_scalar;
        X_vec[6] = Eval_Vert_Spline_ddfdrdt(r, Spline, kr, kt, kp+1)*dt_scalar;
        X_vec[7] = Eval_Vert_Spline_ddfdrdt(r, Spline, kr, kt+1, kp+1)*dt_scalar;
        
        X_vec[8] =  Eval_Vert_Spline_ddfdrdp(r, Spline, kr, kt, kp)*dp_scalar;
        X_vec[9] =  Eval_Vert_Spline_ddfdrdp(r, Spline, kr, kt+1, kp)*dp_scalar;
        X_vec[10] = Eval_Vert_Spline_ddfdrdp(r, Spline, kr, kt, kp+1)*dp_scalar;
        X_vec[11] = Eval_Vert_Spline_ddfdrdp(r, Spline, kr, kt+1, kp+1)*dp_scalar;
        
        X_vec[12]   = BiCubic_Deriv_dddfdrdtdp(r, Spline, kr, kt, kp)*dt_scalar*dp_scalar;
        X_vec[13]   = BiCubic_Deriv_dddfdrdtdp(r, Spline, kr, kt+1, kp)*dt_scalar*dp_scalar;
        X_vec[14]	= BiCubic_Deriv_dddfdrdtdp(r, Spline, kr, kt, kp+1)*dt_scalar*dp_scalar;
        X_vec[15]   = BiCubic_Deriv_dddfdrdtdp(r, Spline, kr, kt+1, kp+1)*dt_scalar*dp_scalar;
        
    }
    
    // Set the coefficients by matrix multiplication
    for(int j = 0; j < 16; j++){
        A_vec[j] = 0;
        for(int k = 0; k < 16; k++){
            A_vec[j] += BiCubic_ConversionMatrix[j][k]*X_vec[k];
        }
    }
       
	double result = 0;
    if((index1==1 && index2==1) || (index1==0 && index2==1) || (index1==1 && index2==0)){
        // d^2 f / dt^2 and d^2 f/ drdt require taking a t derivative of the bicubic spline
        for(int k1 = 1; k1 < 4; k1++){
        for(int k2 = 0; k2 < 4; k2++){
            result+=1.0*k1*A_vec[k2*4 + k1]*pow(t_scaled,k1-1)*pow(p_scaled,k2)/dt_scalar;
        }}
        
    } else if((index1==2 && index2==2) || (index1==0 && index2==2) || (index1==2 && index2==0) || (index1==1 && index2==2) || (index1==2 && index2==1)){
        // d^2 f / dp^2, d^2 f / drdp, and d^2 f / dtdp require taking a p derivative
        for(int k1 = 0; k1 < 4; k1++){
        for(int k2 = 1; k2 < 4; k2++){
            result+=1.0*k2*A_vec[k2*4 + k1]*pow(t_scaled,k1)*pow(p_scaled,k2-1)/dp_scalar;
        }}
    } else {
        // The remaining result is obtained immediately without differentiating the bicubic spline
        for(int k1 = 0; k1 < 4; k1++){
        for(int k2 = 0; k2 < 4; k2++){
            result+=1.0*A_vec[k2*4 + k1]*pow(t_scaled,k1)*pow(p_scaled,k2);
        }}
    }
	
    return result;
}

void Eval_Spline_AllOrder1(double r, double t, double p, struct MultiDimSpline_3DGlobal & Spline, double & f, double & dfdr, double & dfdt, double & dfdp){
    double A_vec[16];
	double X_vec[16];
    
    // Locate kr, kt, kp such that (r[kr] <= r' <= r[kr+1]), (t[kt] <= t' <= t[kt+1]), and (p[kp] <= p' <= p[kp+1])
    double r_eval = min(r, r_max);  r_eval = max(r_eval, r_min);    // Check that r_min <= r_eval <= r_max
    double t_eval = min(t, t_max);  t_eval = max(t_eval, t_min);    // Check that t_min <= t_eval <= t_max
    double p_eval = min(p, p_max);  p_eval = max(p_eval, p_min);    // Check that p_min <= p_eval <= p_max
    
    int kr = Find_Segment(r_eval, Spline.r_vals, Spline.length_r, Spline.accel[0]);
    int kt = Find_Segment(t_eval, Spline.t_vals, Spline.length_t, Spline.accel[1]);
    int kp = Find_Segment(p_eval, Spline.p_vals, Spline.length_p, Spline.accel[2]);
    
    // Set dt_scalar, dp_scalar, and t_scaled, p_scaled to equate to unit square
    double dt_scalar = Spline.t_vals[kt+1] - Spline.t_vals[kt];
    double dp_scalar = Spline.p_vals[kp+1] - Spline.p_vals[kp];
    
    double t_scaled = (t_eval - Spline.t_vals[kt])/dt_scalar;
	double p_scaled = (p_eval - Spline.p_vals[kp])/dp_scalar;
    
    // Precompute finite difference df/dx, df/dy, and d^2f/dxdy which are used multiple times
    double Finite_dfdt_00 = BiCubic_Deriv_dfdt(r_eval, Spline, kr, kt, kp);
    double Finite_dfdt_10 = BiCubic_Deriv_dfdt(r_eval, Spline, kr, kt+1, kp);
    double Finite_dfdt_01 = BiCubic_Deriv_dfdt(r_eval, Spline, kr, kt, kp+1);
    double Finite_dfdt_11 = BiCubic_Deriv_dfdt(r_eval, Spline, kr, kt+1, kp+1);
    
    double Finite_dfdp_00 = BiCubic_Deriv_dfdp(r_eval, Spline, kr, kt, kp);
    double Finite_dfdp_10 = BiCubic_Deriv_dfdp(r_eval, Spline, kr, kt+1, kp);
    double Finite_dfdp_01 = BiCubic_Deriv_dfdp(r_eval, Spline, kr, kt, kp+1);
    double Finite_dfdp_11 = BiCubic_Deriv_dfdp(r_eval, Spline, kr, kt+1, kp+1);
    
    double Finite_ddfdtdp_00 = BiCubic_Deriv_ddfdtdp(r_eval, Spline, kr, kt, kp);
    double Finite_ddfdtdp_10 = BiCubic_Deriv_ddfdtdp(r_eval, Spline, kr, kt+1, kp);
    double Finite_ddfdtdp_01 = BiCubic_Deriv_ddfdtdp(r_eval, Spline, kr, kt, kp+1);
    double Finite_ddfdtdp_11 = BiCubic_Deriv_ddfdtdp(r_eval, Spline, kr, kt+1, kp+1);
    
    // Bicubic Spline f
        X_vec[0] = Eval_Vert_Spline_f(r_eval, Spline, kr, kt, kp);
        X_vec[1] = Eval_Vert_Spline_f(r_eval, Spline, kr, kt+1, kp);
        X_vec[2] = Eval_Vert_Spline_f(r_eval, Spline, kr, kt, kp+1);
        X_vec[3] = Eval_Vert_Spline_f(r_eval, Spline, kr, kt+1, kp+1);
    
        X_vec[4] = Finite_dfdt_00*dt_scalar;
        X_vec[5] = Finite_dfdt_10*dt_scalar;
        X_vec[6] = Finite_dfdt_01*dt_scalar;
        X_vec[7] = Finite_dfdt_11*dt_scalar;
    
        X_vec[8] =  Finite_dfdp_00*dp_scalar;
        X_vec[9] =  Finite_dfdp_10*dp_scalar;
        X_vec[10] = Finite_dfdp_01*dp_scalar;
        X_vec[11] = Finite_dfdp_11*dp_scalar;
    
        X_vec[12]   = Finite_ddfdtdp_00*dt_scalar*dp_scalar;
        X_vec[13]   = Finite_ddfdtdp_10*dt_scalar*dp_scalar;
        X_vec[14]	= Finite_ddfdtdp_01*dt_scalar*dp_scalar;
        X_vec[15]   = Finite_ddfdtdp_11*dt_scalar*dp_scalar;
    
        // Set A[j] by matrix multiplication
        for(int j = 0; j < 16; j++){
            A_vec[j] = 0;
            for(int k = 0; k < 16; k++){
                A_vec[j] += BiCubic_ConversionMatrix[j][k]*X_vec[k];
            }
        }

        f = 0;
        for(int k1 = 0; k1 < 4; k1++){
        for(int k2 = 0; k2 < 4; k2++){
            f+=1.0*A_vec[k1 + 4*k2]*pow(t_scaled,k1)*pow(p_scaled,k2);
        }}
    
    // Bicubic Spline df/dr
        X_vec[0] = Eval_Vert_Spline_dfdr(r_eval, Spline, kr, kt, kp);
        X_vec[1] = Eval_Vert_Spline_dfdr(r_eval, Spline, kr, kt+1, kp);
        X_vec[2] = Eval_Vert_Spline_dfdr(r_eval, Spline, kr, kt, kp+1);
        X_vec[3] = Eval_Vert_Spline_dfdr(r_eval, Spline, kr, kt+1, kp+1);
    
        X_vec[4] = Eval_Vert_Spline_ddfdrdt(r_eval, Spline, kr, kt, kp)*dt_scalar;
        X_vec[5] = Eval_Vert_Spline_ddfdrdt(r_eval, Spline, kr, kt+1, kp)*dt_scalar;
        X_vec[6] = Eval_Vert_Spline_ddfdrdt(r_eval, Spline, kr, kt, kp+1)*dt_scalar;
        X_vec[7] = Eval_Vert_Spline_ddfdrdt(r_eval, Spline, kr, kt+1, kp+1)*dt_scalar;
    
        X_vec[8] =  Eval_Vert_Spline_ddfdrdp(r_eval, Spline, kr, kt, kp)*dp_scalar;
        X_vec[9] =  Eval_Vert_Spline_ddfdrdp(r_eval, Spline, kr, kt+1, kp)*dp_scalar;
        X_vec[10] = Eval_Vert_Spline_ddfdrdp(r_eval, Spline, kr, kt, kp+1)*dp_scalar;
        X_vec[11] = Eval_Vert_Spline_ddfdrdp(r_eval, Spline, kr, kt+1, kp+1)*dp_scalar;
    
        X_vec[12]   = BiCubic_Deriv_dddfdrdtdp(r_eval, Spline, kr, kt, kp)*dt_scalar*dp_scalar;
        X_vec[13]   = BiCubic_Deriv_dddfdrdtdp(r_eval, Spline, kr, kt+1, kp)*dt_scalar*dp_scalar;
        X_vec[14]	= BiCubic_Deriv_dddfdrdtdp(r_eval, Spline, kr, kt, kp+1)*dt_scalar*dp_scalar;
        X_vec[15]   = BiCubic_Deriv_dddfdrdtdp(r_eval, Spline, kr, kt+1, kp+1)*dt_scalar*dp_scalar;
    
        // Set A[j] by matrix multiplication
        for(int j = 0; j < 16; j++){
            A_vec[j] = 0;
            for(int k = 0; k < 16; k++){
                A_vec[j] += BiCubic_ConversionMatrix[j][k]*X_vec[k];
            }
        }
    
        dfdr = 0;
        for(int k1 = 0; k1 < 4; k1++){
        for(int k2 = 0; k2 < 4; k2++){
            dfdr+=1.0*A_vec[k1 + 4*k2]*pow(t_scaled,k1)*pow(p_scaled,k2);
        }}
	
    // Bicubic Spline df/dt
        X_vec[0] = Finite_dfdt_00;
        X_vec[1] = Finite_dfdt_10;
        X_vec[2] = Finite_dfdt_01;
        X_vec[3] = Finite_dfdt_11;
    
        X_vec[4] = BiCubic_Deriv_ddfdtdt(r_eval, Spline, kr, kt, kp)*dt_scalar;
        X_vec[5] = BiCubic_Deriv_ddfdtdt(r_eval, Spline, kr, kt+1, kp)*dt_scalar;
        X_vec[6] = BiCubic_Deriv_ddfdtdt(r_eval, Spline, kr, kt, kp+1)*dt_scalar;
        X_vec[7] = BiCubic_Deriv_ddfdtdt(r_eval, Spline, kr, kt+1, kp+1)*dt_scalar;
    
        X_vec[8] =  Finite_ddfdtdp_00*dp_scalar;
        X_vec[9] =  Finite_ddfdtdp_10*dp_scalar;
        X_vec[10] = Finite_ddfdtdp_01*dp_scalar;
        X_vec[11] = Finite_ddfdtdp_11*dp_scalar;
    
        X_vec[12]   = BiCubic_Deriv_dddfdtdtdp(r_eval, Spline, kr, kt, kp)*dt_scalar*dp_scalar;
        X_vec[13]   = BiCubic_Deriv_dddfdtdtdp(r_eval, Spline, kr, kt+1, kp)*dt_scalar*dp_scalar;
        X_vec[14]	= BiCubic_Deriv_dddfdtdtdp(r_eval, Spline, kr, kt, kp+1)*dt_scalar*dp_scalar;
        X_vec[15]   = BiCubic_Deriv_dddfdtdtdp(r_eval, Spline, kr, kt+1, kp+1)*dt_scalar*dp_scalar;
    
        // Set A[j] by matrix multiplication
        for(int j = 0; j < 16; j++){
            A_vec[j] = 0;
            for(int k = 0; k < 16; k++){
                A_vec[j] += BiCubic_ConversionMatrix[j][k]*X_vec[k];
            }
        }
    
        dfdt = 0;
        for(int k1 = 0; k1 < 4; k1++){
        for(int k2 = 0; k2 < 4; k2++){
            dfdt+=1.0*A_vec[k1 + 4*k2]*pow(t_scaled,k1)*pow(p_scaled,k2);
        }}
    
    // Bicubic Spline df/dp
        X_vec[0] = Finite_dfdp_00;
        X_vec[1] = Finite_dfdp_10;
        X_vec[2] = Finite_dfdp_01;
        X_vec[3] = Finite_dfdp_11;
    
        X_vec[4] = Finite_ddfdtdp_00*dt_scalar;
        X_vec[5] = Finite_ddfdtdp_10*dt_scalar;
        X_vec[6] = Finite_ddfdtdp_01*dt_scalar;
        X_vec[7] = Finite_ddfdtdp_11*dt_scalar;
    
        X_vec[8] =  BiCubic_Deriv_ddfdpdp(r_eval, Spline, kr, kt, kp)*dp_scalar;
        X_vec[9] =  BiCubic_Deriv_ddfdpdp(r_eval, Spline, kr, kt+1, kp)*dp_scalar;
        X_vec[10] = BiCubic_Deriv_ddfdpdp(r_eval, Spline, kr, kt, kp+1)*dp_scalar;
        X_vec[11] = BiCubic_Deriv_ddfdpdp(r_eval, Spline, kr, kt+1, kp+1)*dp_scalar;
    
        X_vec[12]   = BiCubic_Deriv_dddfdtdpdp(r_eval, Spline, kr, kt, kp)*dt_scalar*dp_scalar;
        X_vec[13]   = BiCubic_Deriv_dddfdtdpdp(r_eval, Spline, kr, kt+1, kp)*dt_scalar*dp_scalar;
        X_vec[14]	= BiCubic_Deriv_dddfdtdpdp(r_eval, Spline, kr, kt, kp+1)*dt_scalar*dp_scalar;
        X_vec[15]   = BiCubic_Deriv_dddfdtdpdp(r_eval, Spline, kr, kt+1, kp+1)*dt_scalar*dp_scalar;
    
        // Set A[j] by matrix multiplication
        for(int j = 0; j < 16; j++){
            A_vec[j] = 0;
            for(int k = 0; k < 16; k++){
                A_vec[j] += BiCubic_ConversionMatrix[j][k]*X_vec[k];
            }
        }
    
        dfdp = 0;
        for(int k1 = 0; k1 < 4; k1++){
        for(int k2 = 0; k2 < 4; k2++){
            dfdp+=1.0*A_vec[k1 + 4*k2]*pow(t_scaled,k1)*pow(p_scaled,k2);
        }}
}

void Eval_Spline_AllOrder2(double r, double t, double p, struct MultiDimSpline_3DGlobal & Spline, double & f, double & dfdr, double & dfdt, double & dfdp, double & ddfdrdr, double & ddfdtdt, double & ddfdpdp, double & ddfdrdt, double & ddfdrdp, double & ddfdtdp){
    double A_vec[16];
	double X_vec[16];
    
    // Locate kr, kt, kp such that (r[kr] <= r' <= r[kr+1]), (t[kt] <= t' <= t[kt+1]), and (p[kp] <= p' <= p[kp+1])
    double r_eval = min(r, r_max);  r_eval = max(r_eval, r_min);    // Check that r_min <= r_eval <= r_max
    double t_eval = min(t, t_max);  t_eval = max(t_eval, t_min);    // Check that t_min <= t_eval <= t_max
    double p_eval = min(p, p_max);  p_eval = max(p_eval, p_min);    // Check that p_min <= p_eval <= p_max
    
    int kr = Find_Segment(r_eval, Spline.r_vals, Spline.length_r, Spline.accel[0]);
    int kt = Find_Segment(t_eval, Spline.t_vals, Spline.length_t, Spline.accel[1]);
    int kp = Find_Segment(p_eval, Spline.p_vals, Spline.length_p, Spline.accel[2]);
    
    // Set dt_scalar, dp_scalar, and t_scaled, p_scaled to equate to unit square
    double dt_scalar = Spline.t_vals[kt+1] - Spline.t_vals[kt];
    double dp_scalar = Spline.p_vals[kp+1] - Spline.p_vals[kp];
    
    double t_scaled = (t_eval - Spline.t_vals[kt])/dt_scalar;
	double p_scaled = (p_eval - Spline.p_vals[kp])/dp_scalar;
    
    // Precompute finite difference df/dx, df/dy, and d^2f/dxdy which are used multiple times
    double Finite_dfdt_00 = BiCubic_Deriv_dfdt(r_eval, Spline, kr, kt, kp);
    double Finite_dfdt_10 = BiCubic_Deriv_dfdt(r_eval, Spline, kr, kt+1, kp);
    double Finite_dfdt_01 = BiCubic_Deriv_dfdt(r_eval, Spline, kr, kt, kp+1);
    double Finite_dfdt_11 = BiCubic_Deriv_dfdt(r_eval, Spline, kr, kt+1, kp+1);
    
    double Finite_dfdp_00 = BiCubic_Deriv_dfdp(r_eval, Spline, kr, kt, kp);
    double Finite_dfdp_10 = BiCubic_Deriv_dfdp(r_eval, Spline, kr, kt+1, kp);
    double Finite_dfdp_01 = BiCubic_Deriv_dfdp(r_eval, Spline, kr, kt, kp+1);
    double Finite_dfdp_11 = BiCubic_Deriv_dfdp(r_eval, Spline, kr, kt+1, kp+1);
    
    double Finite_ddfdtdp_00 = BiCubic_Deriv_ddfdtdp(r_eval, Spline, kr, kt, kp);
    double Finite_ddfdtdp_10 = BiCubic_Deriv_ddfdtdp(r_eval, Spline, kr, kt+1, kp);
    double Finite_ddfdtdp_01 = BiCubic_Deriv_ddfdtdp(r_eval, Spline, kr, kt, kp+1);
    double Finite_ddfdtdp_11 = BiCubic_Deriv_ddfdtdp(r_eval, Spline, kr, kt+1, kp+1);
    
    // Bicubic Spline f
        X_vec[0] = Eval_Vert_Spline_f(r_eval, Spline, kr, kt, kp);
        X_vec[1] = Eval_Vert_Spline_f(r_eval, Spline, kr, kt+1, kp);
        X_vec[2] = Eval_Vert_Spline_f(r_eval, Spline, kr, kt, kp+1);
        X_vec[3] = Eval_Vert_Spline_f(r_eval, Spline, kr, kt+1, kp+1);
    
        X_vec[4] = Finite_dfdt_00*dt_scalar;
        X_vec[5] = Finite_dfdt_10*dt_scalar;
        X_vec[6] = Finite_dfdt_01*dt_scalar;
        X_vec[7] = Finite_dfdt_11*dt_scalar;
    
        X_vec[8] =  Finite_dfdp_00*dp_scalar;
        X_vec[9] =  Finite_dfdp_10*dp_scalar;
        X_vec[10] = Finite_dfdp_01*dp_scalar;
        X_vec[11] = Finite_dfdp_11*dp_scalar;
    
        X_vec[12]   = Finite_ddfdtdp_00*dt_scalar*dp_scalar;
        X_vec[13]   = Finite_ddfdtdp_10*dt_scalar*dp_scalar;
        X_vec[14]	= Finite_ddfdtdp_01*dt_scalar*dp_scalar;
        X_vec[15]   = Finite_ddfdtdp_11*dt_scalar*dp_scalar;
    
        for(int j = 0; j < 16; j++){
            A_vec[j] = 0;
            for(int k = 0; k < 16; k++){
                A_vec[j] += BiCubic_ConversionMatrix[j][k]*X_vec[k];
            }
        }
    
        f = 0;
        for(int k1 = 0; k1 < 4; k1++){
        for(int k2 = 0; k2 < 4; k2++){
            f+=1.0*A_vec[k1 + 4*k2]*pow(t_scaled,k1)*pow(p_scaled,k2);
        }}
    
    // Bicubic Spline df/dr, d^2f/drdt, and d^2f/drdp
        X_vec[0] = Eval_Vert_Spline_dfdr(r_eval, Spline, kr, kt, kp);
        X_vec[1] = Eval_Vert_Spline_dfdr(r_eval, Spline, kr, kt+1, kp);
        X_vec[2] = Eval_Vert_Spline_dfdr(r_eval, Spline, kr, kt, kp+1);
        X_vec[3] = Eval_Vert_Spline_dfdr(r_eval, Spline, kr, kt+1, kp+1);
    
        X_vec[4] = Eval_Vert_Spline_ddfdrdt(r_eval, Spline, kr, kt, kp)*dt_scalar;
        X_vec[5] = Eval_Vert_Spline_ddfdrdt(r_eval, Spline, kr, kt+1, kp)*dt_scalar;
        X_vec[6] = Eval_Vert_Spline_ddfdrdt(r_eval, Spline, kr, kt, kp+1)*dt_scalar;
        X_vec[7] = Eval_Vert_Spline_ddfdrdt(r_eval, Spline, kr, kt+1, kp+1)*dt_scalar;
    
        X_vec[8] =  Eval_Vert_Spline_ddfdrdp(r_eval, Spline, kr, kt, kp)*dp_scalar;
        X_vec[9] =  Eval_Vert_Spline_ddfdrdp(r_eval, Spline, kr, kt+1, kp)*dp_scalar;
        X_vec[10] = Eval_Vert_Spline_ddfdrdp(r_eval, Spline, kr, kt, kp+1)*dp_scalar;
        X_vec[11] = Eval_Vert_Spline_ddfdrdp(r_eval, Spline, kr, kt+1, kp+1)*dp_scalar;
    
        X_vec[12]   = BiCubic_Deriv_dddfdrdtdp(r_eval, Spline, kr, kt, kp)*dt_scalar*dp_scalar;
        X_vec[13]   = BiCubic_Deriv_dddfdrdtdp(r_eval, Spline, kr, kt+1, kp)*dt_scalar*dp_scalar;
        X_vec[14]	= BiCubic_Deriv_dddfdrdtdp(r_eval, Spline, kr, kt, kp+1)*dt_scalar*dp_scalar;
        X_vec[15]   = BiCubic_Deriv_dddfdrdtdp(r_eval, Spline, kr, kt+1, kp+1)*dt_scalar*dp_scalar;
    
        for(int j = 0; j < 16; j++){
            A_vec[j] = 0;
            for(int k = 0; k < 16; k++){
                A_vec[j] += BiCubic_ConversionMatrix[j][k]*X_vec[k];
            }
        }
    
        dfdr = 0;
        for(int k1 = 0; k1 < 4; k1++){
        for(int k2 = 0; k2 < 4; k2++){
            dfdr+=1.0*A_vec[k1 + 4*k2]*pow(t_scaled,k1)*pow(p_scaled,k2);
        }}
    
        ddfdrdt = 0;
        for(int k1 = 1; k1 < 4; k1++){
        for(int k2 = 0; k2 < 4; k2++){
            ddfdrdt+=1.0*k1*A_vec[k1 + 4*k2]*pow(t_scaled,k1-1)*pow(p_scaled,k2);
        }}
    
        ddfdrdp = 0;
        for(int k1 = 0; k1 < 4; k1++){
        for(int k2 = 1; k2 < 4; k2++){
            ddfdrdp+=1.0*k2*A_vec[k1 + 4*k2]*pow(t_scaled,k1)*pow(p_scaled,k2-1);
        }}
	
    // Bicubic Spline df/dt
        X_vec[0] = Finite_dfdt_00;
        X_vec[1] = Finite_dfdt_10;
        X_vec[2] = Finite_dfdt_01;
        X_vec[3] = Finite_dfdt_11;
    
        X_vec[4] = BiCubic_Deriv_ddfdtdt(r_eval, Spline, kr, kt, kp)*dt_scalar;
        X_vec[5] = BiCubic_Deriv_ddfdtdt(r_eval, Spline, kr, kt+1, kp)*dt_scalar;
        X_vec[6] = BiCubic_Deriv_ddfdtdt(r_eval, Spline, kr, kt, kp+1)*dt_scalar;
        X_vec[7] = BiCubic_Deriv_ddfdtdt(r_eval, Spline, kr, kt+1, kp+1)*dt_scalar;
    
        X_vec[8] =  Finite_ddfdtdp_00*dp_scalar;
        X_vec[9] =  Finite_ddfdtdp_10*dp_scalar;
        X_vec[10] = Finite_ddfdtdp_01*dp_scalar;
        X_vec[11] = Finite_ddfdtdp_11*dp_scalar;
    
        X_vec[12]   = BiCubic_Deriv_dddfdtdtdp(r_eval, Spline, kr, kt, kp)*dt_scalar*dp_scalar;
        X_vec[13]   = BiCubic_Deriv_dddfdtdtdp(r_eval, Spline, kr, kt+1, kp)*dt_scalar*dp_scalar;
        X_vec[14]	= BiCubic_Deriv_dddfdtdtdp(r_eval, Spline, kr, kt, kp+1)*dt_scalar*dp_scalar;
        X_vec[15]   = BiCubic_Deriv_dddfdtdtdp(r_eval, Spline, kr, kt+1, kp+1)*dt_scalar*dp_scalar;
    
        for(int j = 0; j < 16; j++){
            A_vec[j] = 0;
            for(int k = 0; k < 16; k++){
                A_vec[j] += BiCubic_ConversionMatrix[j][k]*X_vec[k];
            }
        }

        dfdt = 0;
        for(int k1 = 0; k1 < 4; k1++){
        for(int k2 = 0; k2 < 4; k2++){
            dfdt+=1.0*A_vec[k1 + 4*k2]*pow(t_scaled,k1)*pow(p_scaled,k2);
        }}
    
        ddfdtdt = 0;
        for(int k1 = 1; k1 < 4; k1++){
        for(int k2 = 0; k2 < 4; k2++){
            ddfdtdt+=1.0*k1*A_vec[k1 + 4*k2]*pow(t_scaled,k1-1)*pow(p_scaled,k2);
        }}
    
        ddfdtdp = 0;
        for(int k1 = 0; k1 < 4; k1++){
        for(int k2 = 1; k2 < 4; k2++){
            ddfdtdp+=1.0*k2*A_vec[k1 + 4*k2]*pow(t_scaled,k1)*pow(p_scaled,k2-1);
        }}
    
    // Bicubic Spline df/dp
        X_vec[0] = Finite_dfdp_00;
        X_vec[1] = Finite_dfdp_10;
        X_vec[2] = Finite_dfdp_01;
        X_vec[3] = Finite_dfdp_11;
    
        X_vec[4] = Finite_ddfdtdp_00*dt_scalar;
        X_vec[5] = Finite_ddfdtdp_10*dt_scalar;
        X_vec[6] = Finite_ddfdtdp_01*dt_scalar;
        X_vec[7] = Finite_ddfdtdp_11*dt_scalar;
    
        X_vec[8] =  BiCubic_Deriv_ddfdpdp(r_eval, Spline, kr, kt, kp)*dp_scalar;
        X_vec[9] =  BiCubic_Deriv_ddfdpdp(r_eval, Spline, kr, kt+1, kp)*dp_scalar;
        X_vec[10] = BiCubic_Deriv_ddfdpdp(r_eval, Spline, kr, kt, kp+1)*dp_scalar;
        X_vec[11] = BiCubic_Deriv_ddfdpdp(r_eval, Spline, kr, kt+1, kp+1)*dp_scalar;
    
        X_vec[12]   = BiCubic_Deriv_dddfdtdpdp(r_eval, Spline, kr, kt, kp)*dt_scalar*dp_scalar;
        X_vec[13]   = BiCubic_Deriv_dddfdtdpdp(r_eval, Spline, kr, kt+1, kp)*dt_scalar*dp_scalar;
        X_vec[14]	= BiCubic_Deriv_dddfdtdpdp(r_eval, Spline, kr, kt, kp+1)*dt_scalar*dp_scalar;
        X_vec[15]   = BiCubic_Deriv_dddfdtdpdp(r_eval, Spline, kr, kt+1, kp+1)*dt_scalar*dp_scalar;
    
        for(int j = 0; j < 16; j++){
            A_vec[j] = 0;
            for(int k = 0; k < 16; k++){
                A_vec[j] += BiCubic_ConversionMatrix[j][k]*X_vec[k];
            }
        }
    
        dfdp = 0;
        for(int k1 = 0; k1 < 4; k1++){
        for(int k2 = 0; k2 < 4; k2++){
            dfdp+=1.0*A_vec[k1 + 4*k2]*pow(t_scaled,k1)*pow(p_scaled,k2);
        }}
    
        ddfdpdp = 0;
        for(int k1 = 0; k1 < 4; k1++){
        for(int k2 = 1; k2 < 4; k2++){
            ddfdpdp+=1.0*k2*A_vec[k1 + 4*k2]*pow(t_scaled,k1)*pow(p_scaled,k2-1);
        }}
    
    // Bicubic spline d^2f/dr^2
        X_vec[0] = Eval_Vert_Spline_ddfdrdr(r_eval, Spline, kr, kt, kp);
        X_vec[1] = Eval_Vert_Spline_ddfdrdr(r_eval, Spline, kr, kt+1, kp);
        X_vec[2] = Eval_Vert_Spline_ddfdrdr(r_eval, Spline, kr, kt, kp+1);
        X_vec[3] = Eval_Vert_Spline_ddfdrdr(r_eval, Spline, kr, kt+1, kp+1);
    
        X_vec[4] = BiCubic_Deriv_dddfdrdrdt(r_eval, Spline, kr, kt, kp)*dt_scalar;
        X_vec[5] = BiCubic_Deriv_dddfdrdrdt(r_eval, Spline, kr, kt+1, kp)*dt_scalar;
        X_vec[6] = BiCubic_Deriv_dddfdrdrdt(r_eval, Spline, kr, kt, kp+1)*dt_scalar;
        X_vec[7] = BiCubic_Deriv_dddfdrdrdt(r_eval, Spline, kr, kt+1, kp+1)*dt_scalar;
    
        X_vec[8] =  BiCubic_Deriv_dddfdrdrdp(r_eval, Spline, kr, kt, kp)*dp_scalar;
        X_vec[9] =  BiCubic_Deriv_dddfdrdrdp(r_eval, Spline, kr, kt+1, kp)*dp_scalar;
        X_vec[10] = BiCubic_Deriv_dddfdrdrdp(r_eval, Spline, kr, kt, kp+1)*dp_scalar;
        X_vec[11] = BiCubic_Deriv_dddfdrdrdp(r_eval, Spline, kr, kt+1, kp+1)*dp_scalar;
    
        X_vec[12]   = BiCubic_Deriv_ddddfdrdrdtdp(r_eval, Spline, kr, kt, kp)*dt_scalar*dp_scalar;
        X_vec[13]   = BiCubic_Deriv_ddddfdrdrdtdp(r_eval, Spline, kr, kt+1, kp)*dt_scalar*dp_scalar;
        X_vec[14]	= BiCubic_Deriv_ddddfdrdrdtdp(r_eval, Spline, kr, kt, kp+1)*dt_scalar*dp_scalar;
        X_vec[15]   = BiCubic_Deriv_ddddfdrdrdtdp(r_eval, Spline, kr, kt+1, kp+1)*dt_scalar*dp_scalar;
    
        // Set A[j] by matrix multiplication
        for(int j = 0; j < 16; j++){
            A_vec[j] = 0;
            for(int k = 0; k < 16; k++){
                A_vec[j] += BiCubic_ConversionMatrix[j][k]*X_vec[k];
            }
        }
    
        // Evaluate the bicubic spline
        ddfdrdr = 0;
        for(int k1 = 0; k1 < 4; k1++){
        for(int k2 = 0; k2 < 4; k2++){
            ddfdrdr+=1.0*A_vec[k1 + 4*k2]*pow(t_scaled,k1)*pow(p_scaled,k2);
        }}
}


//----------------------------------------------------//
//-------------Combined Function to Input-------------//
//------G2S Files and Generate the Interpolation------//
//----------------------------------------------------//
struct MultiDimSpline_3DGlobal Temp_Spline;
struct MultiDimSpline_3DGlobal Windu_Spline;
struct MultiDimSpline_3DGlobal Windv_Spline;
struct MultiDimSpline_3DGlobal Density_Spline;

void Spline_Multi_G2S(char* profile_prefix, char* lat_file, char* lon_file, char* option){
    SetUp_G2S_Arrays(profile_prefix, lat_file, lon_file);
    Load_G2S_Multi(profile_prefix, lat_file, lon_file, option);
    
    Temp_Spline.length_r = r_cnt;           Windu_Spline.length_r = r_cnt;          Windv_Spline.length_r = r_cnt;          Density_Spline.length_r = r_cnt;
    Temp_Spline.length_t = t_cnt;           Windu_Spline.length_t = t_cnt;          Windv_Spline.length_t = t_cnt;          Density_Spline.length_t = t_cnt;
    Temp_Spline.length_p = p_cnt;           Windu_Spline.length_p = p_cnt;          Windv_Spline.length_p = p_cnt;          Density_Spline.length_p = p_cnt;
    Temp_Spline.r_vals = r_vals;            Windu_Spline.r_vals = r_vals;           Windv_Spline.r_vals = r_vals;           Density_Spline.r_vals = r_vals;
    Temp_Spline.t_vals = t_vals;            Windu_Spline.t_vals = t_vals;           Windv_Spline.t_vals = t_vals;           Density_Spline.t_vals = t_vals;
    Temp_Spline.p_vals = p_vals;            Windu_Spline.p_vals = p_vals;           Windv_Spline.p_vals = p_vals;           Density_Spline.p_vals = p_vals;
    Temp_Spline.f_vals = T_vals;            Windu_Spline.f_vals = u_vals;           Windv_Spline.f_vals = v_vals;           Density_Spline.f_vals = rho_vals;
    Temp_Spline.f_slopes = T_slopes;        Windu_Spline.f_slopes = u_slopes;       Windv_Spline.f_slopes = v_slopes;       Density_Spline.f_slopes = rho_slopes;
    Temp_Spline.dfdt_slopes = T_slopes_dt;  Windu_Spline.dfdt_slopes = u_slopes_dt; Windv_Spline.dfdt_slopes = v_slopes_dt; Density_Spline.dfdt_slopes = rho_slopes_dt;
    Temp_Spline.dfdp_slopes = T_slopes_dp;  Windu_Spline.dfdp_slopes = u_slopes_dp; Windv_Spline.dfdp_slopes = v_slopes_dp; Density_Spline.dfdp_slopes = rho_slopes_dp;
    GeoAc_SetPropRegion();
    
    BuildSlopesArrays(Temp_Spline);     BuildSlopesArrays(Windu_Spline);
    BuildSlopesArrays(Density_Spline);  BuildSlopesArrays(Windv_Spline);
    
    Set_Slopes_Multi(Temp_Spline);      Set_Slopes_Multi(Windu_Spline);
    Set_Slopes_Multi(Density_Spline);   Set_Slopes_Multi(Windv_Spline);
}

void ClearAll(){
    ClearSlopesArrays(Temp_Spline);     ClearSlopesArrays(Windu_Spline);
    ClearSlopesArrays(Density_Spline);  ClearSlopesArrays(Windv_Spline);
    Clear_G2S_Arrays();
}

//-------------------------------------//
//---------Atmospheric Density---------//
//-------------------------------------//
double rho(double r, double theta, double phi){
    double r_eval = min(r, r_max);      r_eval = max(r_eval, r_min);    // Check that r_min <= r_eval <= r_max
    double t_eval = min(theta, t_max);  t_eval = max(t_eval, t_min);    // Check that t_min <= t_eval <= t_max
    double p_eval = min(phi, p_max);    p_eval = max(p_eval, p_min);    // Check that p_min <= p_eval <= p_max
    
    return Eval_Spline_f(r_eval,t_eval,p_eval,Density_Spline);
}

//-------------------------------------------//
//---------Thermodynamic Sound Speed---------//
//------------and its Derivatives------------//
//-------------------------------------------//

double gamR = 0.00040187; // gamma * R in km^2/s^2 * 1/K, c(x,y,z) = sqrt(gamma*r*T(x,y,z))

double c(double r, double theta, double phi){
    double r_eval = min(r, r_max);      r_eval = max(r_eval, r_min);    // Check that r_min <= r_eval <= r_max
    double t_eval = min(theta, t_max);  t_eval = max(t_eval, t_min);    // Check that t_min <= t_eval <= t_max
    double p_eval = min(phi, p_max);    p_eval = max(p_eval, p_min);    // Check that p_min <= p_eval <= p_max
    
    return sqrt(gamR * Eval_Spline_f(r_eval,t_eval,p_eval,Temp_Spline));
}

double c_diff(double r, double theta, double phi, int n){
    double r_eval = min(r, r_max);      r_eval = max(r_eval, r_min);    // Check that r_min <= r_eval <= r_max
    double t_eval = min(theta, t_max);  t_eval = max(t_eval, t_min);    // Check that t_min <= t_eval <= t_max
    double p_eval = min(phi, p_max);    p_eval = max(p_eval, p_min);    // Check that p_min <= p_eval <= p_max
    
    return gamR / (2.0 *c(r,theta,phi)) * Eval_Spline_df(r_eval,t_eval,p_eval,n,Temp_Spline);
}

double c_ddiff(double r, double theta, double phi, int n1, int n2){
    double r_eval = min(r, r_max);      r_eval = max(r_eval, r_min);    // Check that r_min <= r_eval <= r_max
    double t_eval = min(theta, t_max);  t_eval = max(t_eval, t_min);    // Check that t_min <= t_eval <= t_max
    double p_eval = min(phi, p_max);    p_eval = max(p_eval, p_min);    // Check that p_min <= p_eval <= p_max
    double SndSpd = c(r,theta,phi);
    
    return gamR / (2.0 * SndSpd) * Eval_Spline_ddf(r_eval,t_eval,p_eval,n1,n2,Temp_Spline)
             - pow(gamR,2)/(4.0 * pow(SndSpd,3)) * Eval_Spline_df(r_eval,t_eval,p_eval,n1,Temp_Spline)*Eval_Spline_df(r_eval,t_eval,p_eval,n2,Temp_Spline);
}


//------------------------------------------//
//---------East-West Wind Component---------//
//------------and its Derivatives-----------//
//------------------------------------------//
double u(double r, double theta, double phi){
    double r_eval = min(r, r_max);      r_eval = max(r_eval, r_min);    // Check that r_min <= r_eval <= r_max
    double t_eval = min(theta, t_max);  t_eval = max(t_eval, t_min);    // Check that t_min <= t_eval <= t_max
    double p_eval = min(phi, p_max);    p_eval = max(p_eval, p_min);    // Check that p_min <= p_eval <= p_max

    return Eval_Spline_f(r_eval, t_eval, p_eval, Windu_Spline);
}


double u_diff(double r, double theta, double phi, int n){
    double r_eval = min(r, r_max);      r_eval = max(r_eval, r_min);    // Check that r_min <= r_eval <= r_max
    double t_eval = min(theta, t_max);  t_eval = max(t_eval, t_min);    // Check that t_min <= t_eval <= t_max
    double p_eval = min(phi, p_max);    p_eval = max(p_eval, p_min);    // Check that p_min <= p_eval <= p_max
    
    return Eval_Spline_df(r_eval, t_eval, p_eval, n, Windu_Spline);

}

double u_ddiff(double r, double theta, double phi, int n1, int n2){
    double r_eval = min(r, r_max);      r_eval = max(r_eval, r_min);    // Check that r_min <= r_eval <= r_max
    double t_eval = min(theta, t_max);  t_eval = max(t_eval, t_min);    // Check that t_min <= t_eval <= t_max
    double p_eval = min(phi, p_max);    p_eval = max(p_eval, p_min);    // Check that p_min <= p_eval <= p_max
    
    return Eval_Spline_ddf(r_eval, t_eval, p_eval, n1, n2, Windu_Spline);
}


//--------------------------------------------//
//---------North-South Wind Component---------//
//-------------and its Derivatives------------//
//--------------------------------------------//
double v(double r, double theta, double phi){
    double r_eval = min(r, r_max);      r_eval = max(r_eval, r_min);    // Check that r_min <= r_eval <= r_max
    double t_eval = min(theta, t_max);  t_eval = max(t_eval, t_min);    // Check that t_min <= t_eval <= t_max
    double p_eval = min(phi, p_max);    p_eval = max(p_eval, p_min);    // Check that p_min <= p_eval <= p_max
    
    return Eval_Spline_f(r_eval, t_eval, p_eval, Windv_Spline);
}


double v_diff(double r, double theta, double phi, int n){
    double r_eval = min(r, r_max);      r_eval = max(r_eval, r_min);    // Check that r_min <= r_eval <= r_max
    double t_eval = min(theta, t_max);  t_eval = max(t_eval, t_min);    // Check that t_min <= t_eval <= t_max
    double p_eval = min(phi, p_max);    p_eval = max(p_eval, p_min);    // Check that p_min <= p_eval <= p_max
    
    return Eval_Spline_df(r_eval, t_eval, p_eval, n, Windv_Spline);
}

double v_ddiff(double r, double theta, double phi, int n1, int n2){
    double r_eval = min(r, r_max);      r_eval = max(r_eval, r_min);    // Check that r_min <= r_eval <= r_max
    double t_eval = min(theta, t_max);  t_eval = max(t_eval, t_min);    // Check that t_min <= t_eval <= t_max
    double p_eval = min(phi, p_max);    p_eval = max(p_eval, p_min);    // Check that p_min <= p_eval <= p_max
    
    return Eval_Spline_ddf(r_eval, t_eval, p_eval, n1, n2, Windv_Spline);
}


//-----------------------------------------//
//---------Vertical Wind Component---------//
//-----------and its Derivatives-----------//
//-----------------------------------------//
double w(double r, double theta, double phi){                         return 0.0;}
double w_diff(double r, double theta, double phi, int n){             return 0.0;}
double w_ddiff(double r, double theta, double phi, int n1, int n2){   return 0.0;}

#endif /* G2S_GLIOBAL_RNGDEP_ATMOSPHERE_CPP_ */
