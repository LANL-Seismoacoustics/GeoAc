# ifndef G2S_SPLINE1D_CPP_
# define G2S_SPLINE1D_CPP_

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <sstream>
#include <math.h>

#include "G2S_Spline1D.h"
#include "Atmo_State.h"
#include "GeoAc.Parameters.h"

using namespace std;
//-----------------------------------------------//
//---------Define the Propagation Region---------//
//-----------------------------------------------//
double z_min, z_max;

void GeoAc_SetPropRegion(){
    z_min = Windu_Spline.x_vals[0];
    z_max = Windu_Spline.x_vals[Windu_Spline.length-1];
    
    GeoAc_vert_limit  =	z_max;
    GeoAc_range_limit = 10000.0;
}

//-----------------------------------------------//
//---------Topographical Ground Function---------//
//-----------------------------------------------//
double z_grnd = 0.0;

double GroundTopography(double x, double y){
    return z_grnd;
}

//----------------------------------------//
//------Parameters for Interpolation------//
//----------------------------------------//
int z_cnt;          // Number of vertical points
int accel;          // Acceleration index
double* z_vals;     // z_k elements (z_cnt length)

double* T_vals;    // Temperature at (r_i, z_k) (r_cnt x z_cnt)
double* u_vals;    // E-W winds at (r_i, z_k) (r_cnt x z_cnt)
double* v_vals;    // N-S winds at (r_i, z_k) (r_cnt x z_cnt)
double* rho_vals;  // Density at (r_i, z_k) (r_cnt x z_cnt)

double* T_slopes;    // Slopes to compute interpoalted temperature at each horizontal node
double* u_slopes;    // Slopes to compute interpolated E-W winds at each horizontal node
double* v_slopes;    // Slopes to compute interpolated N-S winds at each horizontal node
double* rho_slopes;  // Slopes to compute interpolated density at each horizontal node

//----------------------------------------//
//----------File IO Manipulation----------//
//----------------------------------------//
int file_length(string file_name){
	ifstream file_in;	file_in.open(file_name.c_str() );
	if(!file_in.is_open()){
		cout << "Error opening file, check file name" << '\n';
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
		cout << "Error opening file, check file name" << '\n';
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
void SetUp_G2S_Arrays(char* file_name){
    z_cnt = file_length(file_name);
    
    z_vals = new double [z_cnt];
    T_vals = new double [z_cnt];   rho_vals = new double [z_cnt];
    u_vals = new double [z_cnt];   v_vals = new double [z_cnt];
}

void Load_G2S(char* file_name, char* option){
    ifstream file_in; double temp;
    file_in.open(file_name);
    
    if (strncmp(option, "zTuvdp",6) == 0){
        for (int nz = 0; nz < z_cnt; nz++){
            file_in >> z_vals[nz];     // Extract z_i value
            file_in >> T_vals[nz];     // Extract T(z_i)
            file_in >> u_vals[nz];     // Extract u(z_i)
            file_in >> v_vals[nz];     // Extract v(z_i)
            file_in >> rho_vals[nz];   // Extract rho(z_i)
            file_in >> temp;           // Extract p(z_i) but don't store it
            
            // Convert winds m/s -> km/s and scale near the ground to guarantee u(z_g), v(z_g) = 0
            u_vals[nz] *= (2.0 / (1.0 + exp(-(z_vals[nz] - z_grnd)/0.2)) - 1.0) / 1000.0;
            v_vals[nz] *= (2.0 / (1.0 + exp(-(z_vals[nz] - z_grnd)/0.2)) - 1.0) / 1000.0;
        }
    } else if (strncmp(option, "zuvwTdp",7) == 0){
        for (int nz = 0; nz < z_cnt; nz++){
            file_in >> z_vals[nz];     // Extract z_i value
            file_in >> u_vals[nz];     // Extract u(z_i)
            file_in >> v_vals[nz];     // Extract v(z_i)
            file_in >> temp;           // Extract w(z_i) but don't store it
            file_in >> T_vals[nz];     // Extract T(z_i)
            file_in >> rho_vals[nz];   // Extract rho(z_i)
            file_in >> temp;           // Extract p(z_i) but don't store it
            
            // Convert winds m/s -> km/s and scale near the ground to guarantee u(z_g), v(z_g) = 0
            u_vals[nz] *= (2.0 / (1.0 + exp(-(z_vals[nz] - z_grnd)/0.2)) - 1.0) / 1000.0;
            v_vals[nz] *= (2.0 / (1.0 + exp(-(z_vals[nz] - z_grnd)/0.2)) - 1.0) / 1000.0;
        }
    } else { cout << "Unrecognized profile option: " << option << ".  Valid options are: zTuvdp and zuvwTdp" << '\n';}
    file_in.close();
}

void Clear_G2S_Arrays(){
    delete z_vals;
    delete T_vals;      delete rho_vals;
    delete u_vals;      delete v_vals;
}

//----------------------------------------------//
//------------Functions to Build and------------//
//-------Clear the Input and Slope Arrays-------//
//----------------------------------------------//
void BuildSlopeArray(struct NaturalCubicSpline_1D & Spline){        Spline.slopes = new double [Spline.length];}
void ClearSlopeArray(struct NaturalCubicSpline_1D & Spline){        delete Spline.slopes;}

//--------------------------------------//
//-----------Functions to Set-----------//
//-------the Interpolation Slopes-------//
//--------------------------------------//
void Set_Slopes(struct NaturalCubicSpline_1D & Spline){
    double ai, bi, ci, di;
    double new_c[Spline.length-1];
    double new_d[Spline.length];
    
    bi = 2.0 / (Spline.x_vals[1] - Spline.x_vals[0]);
    ci = 1.0 / (Spline.x_vals[1] - Spline.x_vals[0]);
    di = 3.0 * (Spline.f_vals[1] - Spline.f_vals[0]) / pow(Spline.x_vals[1] - Spline.x_vals[0], 2);
    
    new_c[0] = ci/bi;
    new_d[0] = di/bi;
    
    for(int i = 1; i < Spline.length - 1; i++) {
        ai = 1.0/(Spline.x_vals[i] - Spline.x_vals[i-1]);
        bi = 2.0 * (1.0/(Spline.x_vals[i] - Spline.x_vals[i-1]) + 1.0/(Spline.x_vals[i+1] - Spline.x_vals[i]));
        ci = 1.0/(Spline.x_vals[i+1] - Spline.x_vals[i]);
        di = 3.0 * ((Spline.f_vals[i] - Spline.f_vals[i-1]) / pow(Spline.x_vals[i] - Spline.x_vals[i-1], 2)
                    + (Spline.f_vals[i+1] - Spline.f_vals[i]) / pow(Spline.x_vals[i+1] - Spline.x_vals[i], 2) );
        
        new_c[i] = ci/(bi - new_c[i-1]*ai);
        new_d[i] = (di - new_d[i-1]*ai)/(bi - new_c[i-1]*ai);
    }
    
    ai = 1.0/(Spline.x_vals[Spline.length-1] - Spline.x_vals[Spline.length-2]);
    bi = 2.0/(Spline.x_vals[Spline.length-1] - Spline.x_vals[Spline.length-2]);
    di = 3.0 * (Spline.f_vals[Spline.length-1] - Spline.f_vals[Spline.length-2])
            / pow(Spline.x_vals[Spline.length-1] - Spline.x_vals[Spline.length-2], 2);
    
    new_d[Spline.length-1] = (di - new_d[Spline.length - 2]*ai)/(bi - new_c[Spline.length - 2]*ai);
    
    Spline.slopes[Spline.length - 1] = new_d[Spline.length - 1];
    for(int i = Spline.length - 2; i > -1; i--){
        Spline.slopes[i] = new_d[i] - new_c[i] * Spline.slopes[i+1];
    }

}

//-----------------------------------------------//
//------Functions to Find the Segment Index------//
//-----and Evaluate the 1D Vertical Splines------//
//-----------------------------------------------//
int Find_Segment(double x, double* x_vals, int length, int & prev){
    int index = length + 1;
    bool done = false;
    
    if(x > x_vals[length-1] || x < x_vals[0]){
        cout << "Cannot interpolate outside of given bounds.  x = " << x << " is invalid." << '\n';
    } else {
        
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
    }
}

double Eval_Spline_f(double x, struct NaturalCubicSpline_1D & Spline){
    int k = Find_Segment(x, Spline.x_vals, Spline.length, Spline.accel);
    
    if(k < Spline.length){
        double X = (x - Spline.x_vals[k])/(Spline.x_vals[k+1] - Spline.x_vals[k]);
        double A = Spline.slopes[k] * (Spline.x_vals[k+1] - Spline.x_vals[k]) - (Spline.f_vals[k+1] - Spline.f_vals[k]);
        double B = -Spline.slopes[k+1] * (Spline.x_vals[k+1] - Spline.x_vals[k]) + (Spline.f_vals[k+1] - Spline.f_vals[k]);
        
        return (1.0 - X) * Spline.f_vals[k] + X * Spline.f_vals[k+1] + X * (1.0 - X) * (A * (1.0 - X ) + B * X);
    } else { return 0.0;}
}

double Eval_Spline_df(double x, struct NaturalCubicSpline_1D & Spline){
    int k = Find_Segment(x, Spline.x_vals, Spline.length, Spline.accel);
    
    if(k < Spline.length){
        double X = (x - Spline.x_vals[k])/(Spline.x_vals[k+1] - Spline.x_vals[k]);
        double A = Spline.slopes[k] * (Spline.x_vals[k+1] - Spline.x_vals[k]) - (Spline.f_vals[k+1] - Spline.f_vals[k]);
        double B = -Spline.slopes[k+1] * (Spline.x_vals[k+1] - Spline.x_vals[k]) + (Spline.f_vals[k+1] - Spline.f_vals[k]);
        
        return (Spline.f_vals[k+1] - Spline.f_vals[k])/(Spline.x_vals[k+1] - Spline.x_vals[k])
        + (1.0 - 2.0 * X) * (A * (1.0 - X) + B * X)/(Spline.x_vals[k+1] - Spline.x_vals[k])
        + X * (1.0 - X) * (B - A)/(Spline.x_vals[k+1] - Spline.x_vals[k]);
    } else { return 0.0;}
}

double Eval_Spline_ddf(double x, struct NaturalCubicSpline_1D & Spline){
    int k = Find_Segment(x, Spline.x_vals, Spline.length, Spline.accel);
    
    if(k < Spline.length){
        double X = (x - Spline.x_vals[k])/(Spline.x_vals[k+1] - Spline.x_vals[k]);
        double A = Spline.slopes[k] * (Spline.x_vals[k+1] - Spline.x_vals[k]) - (Spline.f_vals[k+1] - Spline.f_vals[k]);
        double B = -Spline.slopes[k+1] * (Spline.x_vals[k+1] - Spline.x_vals[k]) + (Spline.f_vals[k+1] - Spline.f_vals[k]);
        
        return 2.0 * (B - 2.0 * A + (A - B) * 3.0 * X)/pow(Spline.x_vals[k+1] - Spline.x_vals[k],2);
    } else { return 0.0;}
}


//----------------------------------------------------//
//-------------Combined Function to Input-------------//
//------G2S Files and Generate the Interpolation------//
//----------------------------------------------------//
struct NaturalCubicSpline_1D Temp_Spline;
struct NaturalCubicSpline_1D Windu_Spline;
struct NaturalCubicSpline_1D Windv_Spline;
struct NaturalCubicSpline_1D Density_Spline;

void Spline_Single_G2S(char* file_name, char* format){
    SetUp_G2S_Arrays(file_name);
    
    Temp_Spline.length = z_cnt;     Windu_Spline.length = z_cnt;    Windv_Spline.length = z_cnt;    Density_Spline.length = z_cnt;
    Temp_Spline.accel = accel;      Windu_Spline.accel = accel;     Windv_Spline.accel = accel;     Density_Spline.accel = accel;
    Temp_Spline.x_vals = z_vals;    Windu_Spline.x_vals = z_vals;   Windv_Spline.x_vals = z_vals;   Density_Spline.x_vals = z_vals;
    Temp_Spline.f_vals = T_vals;    Windu_Spline.f_vals = u_vals;   Windv_Spline.f_vals = v_vals;   Density_Spline.f_vals = rho_vals;
    Temp_Spline.slopes = T_slopes;  Windu_Spline.slopes = u_slopes; Windv_Spline.slopes = v_slopes; Density_Spline.slopes = rho_slopes;
    
    BuildSlopeArray(Temp_Spline);       BuildSlopeArray(Windu_Spline);
    BuildSlopeArray(Density_Spline);    BuildSlopeArray(Windv_Spline);

    Load_G2S(file_name, format);
    GeoAc_SetPropRegion();
    
    Set_Slopes(Temp_Spline);            Set_Slopes(Windu_Spline);
    Set_Slopes(Density_Spline);         Set_Slopes(Windv_Spline);
}

void ClearAll(){
    ClearSlopeArray(Temp_Spline);       ClearSlopeArray(Windu_Spline);
    ClearSlopeArray(Density_Spline);    ClearSlopeArray(Windv_Spline);    
    Clear_G2S_Arrays();
}


//-------------------------------------//
//---------Atmospheric Density---------//
//-------------------------------------//
double rho(double x, double y, double z){
    double z_eval = min(z, z_max);  z_eval = max(z_eval, z_min);    // Check that z_min <= z_eval <= z_max
    
    return Eval_Spline_f(z_eval,Density_Spline);
}

//-------------------------------------------//
//---------Thermodynamic Sound Speed---------//
//------------and its Derivatives------------//
//-------------------------------------------//

double gamR = 0.00040187; // gamma * R in km^2/s^2 * 1/K, c(x,y,z) = sqrt(gamma*r*T(x,y,z))

double c(double x, double y, double z){
    double z_eval = min(z, z_max);  z_eval = max(z_eval, z_min);    // Check that z_min <= z_eval <= z_max
    
    return sqrt(gamR * Eval_Spline_f(z_eval,Temp_Spline));
}

double c_diff(double x, double y, double z, int n){
    double z_eval = min(z, z_max);  z_eval = max(z_eval, z_min);    // Check that z_min <= z_eval <= z_max
    
    if(n==2){   return gamR / (2.0 * c(x,y,z)) * Eval_Spline_df(z_eval,Temp_Spline);}
    else {      return 0.0;}
}

double c_ddiff(double x, double y, double z, int n1, int n2){
    double z_eval = min(z, z_max);  z_eval = max(z_eval, z_min);    // Check that z_min <= z_eval <= z_max
    
    if(n1==2 && n2==2){
        return gamR / (2.0 * c(x,y,z)) * Eval_Spline_ddf(z_eval,Temp_Spline)
                - pow(gamR,2)/(4.0 * pow(c(x,y,z),3)) * pow(Eval_Spline_df(z_eval,Temp_Spline),2);
    } else {
        return 0.0;
    }
}


//------------------------------------------//
//---------East-West Wind Component---------//
//------------and its Derivatives-----------//
//------------------------------------------//
double u(double x, double y, double z){
    double z_eval = min(z, z_max);  z_eval = max(z_eval, z_min);    // Check that z_min <= z_eval <= z_max
    
    return Eval_Spline_f(z_eval, Windu_Spline);
}


double u_diff(double x, double y, double z, int n){
    double z_eval = min(z, z_max);  z_eval = max(z_eval, z_min);    // Check that z_min <= z_eval <= z_max
    
    if(n==2)  return Eval_Spline_df(z_eval, Windu_Spline);
    else        return 0.0;
}

double u_ddiff(double x, double y, double z, int n1, int n2){
    double z_eval = min(z, z_max);  z_eval = max(z_eval, z_min);    // Check that z_min <= z_eval <= z_max
    
    if(n1==2&&n2==2)    return Eval_Spline_ddf(z_eval, Windu_Spline);
    else                return 0.0;
}


//--------------------------------------------//
//---------North-South Wind Component---------//
//-------------and its Derivatives------------//
//--------------------------------------------//
double v(double x, double y, double z){
    double z_eval = min(z, z_max);  z_eval = max(z_eval, z_min);    // Check that z_min <= z_eval <= z_max
    
    return Eval_Spline_f(z_eval, Windv_Spline);
}

double v_diff(double x, double y, double z, int n){
    double z_eval = min(z, z_max);  z_eval = max(z_eval, z_min);    // Check that z_min <= z_eval <= z_max
    
    if(n==2)   return Eval_Spline_df(z_eval, Windv_Spline);
    else       return 0.0;
}

double v_ddiff(double x, double y, double z, int n1, int n2){
    double z_eval = min(z, z_max);  z_eval = max(z_eval, z_min);    // Check that z_min <= z_eval <= z_max
    
    if(n1==2&&n2==2)    return Eval_Spline_ddf(z_eval, Windv_Spline);
    else                return 0.0;
}


//-----------------------------------------//
//---------Vertical Wind Component---------//
//-----------and its Derivatives-----------//
//-----------------------------------------//
double w(double x, double y, double z){                         return 0.0;}
double w_diff(double x, double y, double z, int n){             return 0.0;}
double w_ddiff(double x, double y, double z, int n1, int n2){   return 0.0;}

#endif /* G2S_SPLINE1D_CPP_ */
