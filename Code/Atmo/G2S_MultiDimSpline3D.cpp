# ifndef G2S_RANGEDEPENDENT_ATMOSPHERE_CPP_
# define G2S_RANGEDEPENDENT_ATMOSPHERE_CPP_

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <sstream>
#include <math.h>

#include "G2S_MultiDimSpline3D.h"
#include "Atmo_State.h"
#include "GeoAc.Parameters.h"

using namespace std;

//-----------------------------------------------//
//---------Define the Propagation Region---------//
//-----------------------------------------------//
double x_min, x_max;
double y_min, y_max;
double z_min, z_max;

void GeoAc_SetPropRegion(){
    x_min = Windu_Spline.x_vals[0], x_max = Windu_Spline.x_vals[Windu_Spline.length_x-1];
    y_min = Windu_Spline.y_vals[0], y_max = Windu_Spline.y_vals[Windu_Spline.length_y-1];
    z_min = Windu_Spline.z_vals[0], z_max = Windu_Spline.z_vals[Windu_Spline.length_z-1];

    GeoAc_vert_limit  =	z_max;
    GeoAc_x_min_limit = x_min;      GeoAc_x_max_limit = x_max;
    GeoAc_y_min_limit = y_min;      GeoAc_y_max_limit = y_max;
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
int x_cnt;  // Number of x data points (assumed constant for all y and z values)
int y_cnt;  // Number of y data points (assumed constant for all x and z values)
int z_cnt;  // Number of z data points (assumed constant for all x and y values)

double* x_vals;     // x_i elements (x_cnt length)
double* y_vals;     // y_j elements (y_cnt length)
double* z_vals;     // z_k elements (z_cnt length)

double*** T_vals;           // Temperature at (x_i, y_j, z_k) (x_cnt x y_cnt x z_cnt)
double*** T_slopes;         // Slopes for vertical splines of temperature at each x[], y[] node
double*** T_slopes_dx;      // Slopes for d/dx of vertical splines of temperature at each x[], y[] node
double*** T_slopes_dy;      // Slopes for d/dy of vertical splines of temperature at each x[], y[] node

double*** u_vals;           // E-W winds at (x_i, y_j, z_k) (x_cnt x y_cnt x z_cnt)
double*** u_slopes;         // Slopes for vertical splines of temperature at each x[], y[] node
double*** u_slopes_dx;      // Slopes for d/dx of vertical splines of temperature at each x[], y[] node
double*** u_slopes_dy;      // Slopes for d/dy of vertical splines of temperature at each x[], y[] node

double*** v_vals;           // N-S winds at (x_i, y_j, z_k) (x_cnt x y_cnt x z_cnt)
double*** v_slopes;         // Slopes for vertical splines of temperature at each x[], y[] node
double*** v_slopes_dx;      // Slopes for d/dx of vertical splines of temperature at each x[], y[] node
double*** v_slopes_dy;      // Slopes for d/dy of vertical splines of temperature at each x[], y[] node

double*** rho_vals;         // Density at (x_i, y_j, z_k) (x_cnt x y_cnt x z_cnt)
double*** rho_slopes;         // Slopes for vertical splines of temperature at each x[], y[] node
double*** rho_slopes_dx;      // Slopes for d/dx of vertical splines of temperature at each x[], y[] node
double*** rho_slopes_dy;      // Slopes for d/dy of vertical splines of temperature at each x[], y[] node


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
void SetUp_G2S_Arrays(char* file_prefix, char* locx_file, char* locy_file){
    char output_buffer [50];
    
    x_cnt = file_length(locx_file);
    y_cnt = file_length(locy_file);
    
    // Open the first file and determine its length
    sprintf(output_buffer, "%s%i.met", file_prefix, 0);
    if(Check_G2S_Format(output_buffer)) z_cnt = file_length(output_buffer);
    
    // Set the sizes of the various arrays
    x_vals = new double [x_cnt];
    y_vals = new double [y_cnt];
    z_vals = new double [z_cnt];
    
    T_vals = new double** [x_cnt];   rho_vals = new double** [x_cnt];
    u_vals = new double** [x_cnt];   v_vals = new double** [x_cnt];
    
    for(int nx = 0; nx < x_cnt; nx++){
        T_vals[nx] = new double* [y_cnt];    rho_vals[nx] = new double* [y_cnt];
        u_vals[nx] = new double* [y_cnt];    v_vals[nx] = new double* [y_cnt];
        for(int ny = 0; ny < y_cnt; ny++){
            T_vals[nx][ny] = new double [z_cnt];    rho_vals[nx][ny] = new double [z_cnt];
            u_vals[nx][ny] = new double [z_cnt];    v_vals[nx][ny] = new double [z_cnt];
        }
    }
}

void Load_G2S_Multi(char* file_prefix, char* locx_file, char* locy_file, char* option){
    ifstream file_in; double temp;
    char output_buffer [50];
    
    file_in.open(locx_file);
    for (int i = 0; i < x_cnt; i++){file_in >> x_vals[i];}
    file_in.close();
    
    file_in.open(locy_file);
    for (int i = 0; i < y_cnt; i++){file_in >> y_vals[i];}
    file_in.close();
    
    // Copy the G2S data into the various arrays
    for(int ny = 0; ny < y_cnt; ny++){
    for(int nx = 0; nx < x_cnt; nx++){
        sprintf(output_buffer, "%s%i.met", file_prefix, nx + ny*x_cnt);
        ifstream file_in; file_in.open(output_buffer);
            
        if (strncmp(option, "zTuvdp",6) == 0){
            for (int nz = 0; nz < z_cnt; nz++){
                file_in >> z_vals[nz];           // Extract z_i value
                file_in >> T_vals[nx][ny][nz];   // Extract T(z_i)
                file_in >> u_vals[nx][ny][nz];   // Extract u(z_i)
                file_in >> v_vals[nx][ny][nz];   // Extract v(z_i)
                file_in >> rho_vals[nx][ny][nz]; // Extract rho(z_i)
                file_in >> temp;                 // Extract p(z_i) but don't store it
                
                // Convert winds m/s -> km/s and scale near the ground to guarantee u(z_g), v(z_g) = 0 (currently set at z_g = 0)
                u_vals[nx][ny][nz] *= (2.0 / (1.0 + exp(-(z_vals[nz] - z_grnd)/0.05)) - 1.0) / 1000.0;
                v_vals[nx][ny][nz] *= (2.0 / (1.0 + exp(-(z_vals[nz] - z_grnd)/0.05)) - 1.0) / 1000.0;
            }
        } else if (strncmp(option, "zuvwTdp",7) == 0){
            for (int nz = 0; nz < z_cnt; nz++){
                file_in >> z_vals[nz];           // Extract z_i value
                file_in >> u_vals[nx][ny][nz];   // Extract u(z_i)
                file_in >> v_vals[nx][ny][nz];   // Extract v(z_i)
                file_in >> temp;                 // Extract w(z_i) but don't store it
                file_in >> T_vals[nx][ny][nz];   // Extract T(z_i)
                file_in >> rho_vals[nx][ny][nz]; // Extract rho(z_i)
                file_in >> temp;                 // Extract p(z_i) but don't store it
                    
                // Convert winds m/s -> km/s and scale near the ground to guarantee u(z_g), v(z_g) = 0 (currently set at z_g = 0)
                u_vals[nx][ny][nz] *= (2.0 / (1.0 + exp(-(z_vals[nz] - z_grnd)/0.05)) - 1.0) / 1000.0;
                v_vals[nx][ny][nz] *= (2.0 / (1.0 + exp(-(z_vals[nz] - z_grnd)/0.05)) - 1.0) / 1000.0;
            }
        } else {
            cout << "Unrecognized profile option: " << option << ".  Valid options are: zTuvdp and zuvwTdp" << '\n';
        }
        file_in.close();
    }}
}

void Clear_G2S_Arrays(){
    // Clear the arrays
    delete x_vals;
    delete y_vals;
    delete z_vals;
    
    for(int nx = 0; nx < x_cnt; nx++){
        for(int ny = 0; ny < y_cnt; ny++){
            delete T_vals[nx][ny];  delete rho_vals[nx][ny];
            delete u_vals[nx][ny];  delete v_vals[nx][ny];
        }
        delete T_vals[nx];  delete rho_vals[nx];
        delete u_vals[nx];  delete v_vals[nx];
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


void BuildInputArrays(struct MultiDimSpline_3D & Spline){
    Spline.x_vals = new double [Spline.length_x];
    Spline.y_vals = new double [Spline.length_y];
    Spline.z_vals = new double [Spline.length_z];
    
    Spline.f_vals = new double** [Spline.length_x];
    
    for(int i = 0; i < Spline.length_x; i++){
        Spline.f_vals[i] = new double* [Spline.length_y];
        for(int j = 0; j < Spline.length_y; j++){
            Spline.f_vals[i][j] = new double [Spline.length_z];
        }
    }
}
void BuildSlopesArrays(struct MultiDimSpline_3D & Spline){
    Spline.f_slopes = new double** [Spline.length_x];
    Spline.dfdx_slopes = new double** [Spline.length_x];
    Spline.dfdy_slopes = new double** [Spline.length_x];
    
    for(int i = 0; i < Spline.length_x; i++){
        Spline.f_slopes[i] = new double* [Spline.length_y];
        Spline.dfdx_slopes[i] = new double* [Spline.length_y];
        Spline.dfdy_slopes[i] = new double* [Spline.length_y];
        for(int j = 0; j < Spline.length_y; j++){
            Spline.f_slopes[i][j] = new double [Spline.length_z];
            Spline.dfdx_slopes[i][j] = new double [Spline.length_z];
            Spline.dfdy_slopes[i][j] = new double [Spline.length_z];
        }
    }
}

void ClearInputArrays(struct MultiDimSpline_3D & Spline){
    delete Spline.x_vals;
    delete Spline.y_vals;
    delete Spline.z_vals;
    
    for(int i = 0; i < Spline.length_x; i++){
        for(int j = 0; j < Spline.length_y; j++){
            delete Spline.f_vals[i][j];
        }
        delete Spline.f_vals[i];
    }
    delete Spline.f_vals;
    
}

void ClearSlopesArrays(struct MultiDimSpline_3D & Spline){
    for(int i = 0; i < Spline.length_x; i++){
        for(int j = 0; j < Spline.length_y; j++){
            delete Spline.f_slopes[i][j];
            delete Spline.dfdx_slopes[i][j];
            delete Spline.dfdy_slopes[i][j];
        }
        delete Spline.f_slopes[i];
        delete Spline.dfdx_slopes[i];
        delete Spline.dfdy_slopes[i];
    }
    delete Spline.f_slopes;
    delete Spline.dfdx_slopes;
    delete Spline.dfdy_slopes;
    
}


//--------------------------------------//
//-----------Functions to Set-----------//
//-------the Interpolation Slopes-------//
//--------------------------------------//
void Set_Slopes_Multi(struct MultiDimSpline_3D & Spline){
    double ai, bi, ci, di;
    
    double new_c[Spline.length_z-1];
    double new_d[Spline.length_z];
    
    // Set f_slopes
    for(int mx = 0; mx < Spline.length_x; mx++){
        for(int my = 0; my < Spline.length_y; my++){
            
            bi = 2.0 / (Spline.z_vals[1] - Spline.z_vals[0]);
            ci = 1.0 / (Spline.z_vals[1] - Spline.z_vals[0]);
            di = 3.0 * (Spline.f_vals[mx][my][1] - Spline.f_vals[mx][my][0]) / pow(Spline.z_vals[1] - Spline.z_vals[0], 2);
            
            new_c[0] = ci/bi;
            new_d[0] = di/bi;
            
            for(int i = 1; i < Spline.length_z - 1; i++) {
                ai = 1.0/(Spline.z_vals[i] - Spline.z_vals[i-1]);
                bi = 2.0 * (1.0/(Spline.z_vals[i] - Spline.z_vals[i-1]) + 1.0/(Spline.z_vals[i+1] - Spline.z_vals[i]));
                ci = 1.0/(Spline.z_vals[i+1] - Spline.z_vals[i]);
                di = 3.0 * ((Spline.f_vals[mx][my][i] - Spline.f_vals[mx][my][i-1]) / pow(Spline.z_vals[i] - Spline.z_vals[i-1], 2)
                            + (Spline.f_vals[mx][my][i+1] - Spline.f_vals[mx][my][i]) / pow(Spline.z_vals[i+1] - Spline.z_vals[i], 2) );
                
                new_c[i] = ci/(bi - new_c[i-1]*ai);
                new_d[i] = (di - new_d[i-1]*ai)/(bi - new_c[i-1]*ai);
            }
            
            ai = 1.0/(Spline.z_vals[Spline.length_z-1] - Spline.z_vals[Spline.length_z-2]);
            bi = 2.0/(Spline.z_vals[Spline.length_z-1] - Spline.z_vals[Spline.length_z-2]);
            di = 3.0 * (Spline.f_vals[mx][my][Spline.length_z-1] - Spline.f_vals[mx][my][Spline.length_z-2])
            / pow(Spline.z_vals[Spline.length_z-1] - Spline.z_vals[Spline.length_z-2], 2);
            
            new_d[Spline.length_z-1] = (di - new_d[Spline.length_z - 2]*ai)/(bi - new_c[Spline.length_z - 2]*ai);
            
            Spline.f_slopes[mx][my][Spline.length_z - 1] = new_d[Spline.length_z - 1];
            for(int i = Spline.length_z - 2; i >= 0; i--) Spline.f_slopes[mx][my][i] = new_d[i] - new_c[i] * Spline.f_slopes[mx][my][i+1];
        }}
    
    // Set df/dx and df/dy values to use for setting dfdx_slopes and dfdy_slopes
    double dfdx [Spline.length_x][Spline.length_y][Spline.length_z];
    double dfdy [Spline.length_x][Spline.length_y][Spline.length_z];
    
    for(int mx = 0; mx < Spline.length_x; mx++){
        for(int my = 0; my < Spline.length_y; my++){
            for(int mz = 0; mz < Spline.length_z; mz++){
                int mx_up = min(mx + 1, Spline.length_x - 1);   int mx_dn = max(mx - 1, 0);
                int my_up = min(my + 1, Spline.length_y - 1);   int my_dn = max(my - 1, 0);
                
                dfdx[mx][my][mz] = (Spline.f_vals[mx_up][my][mz] - Spline.f_vals[mx_dn][my][mz])/(Spline.x_vals[mx_up] - Spline.x_vals[mx_dn]);
                dfdy[mx][my][mz] = (Spline.f_vals[mx][my_up][mz] - Spline.f_vals[mx][my_dn][mz])/(Spline.y_vals[my_up] - Spline.y_vals[my_dn]);
            }}}
    
    
    // Set dfdx_slopes
    for(int mx = 0; mx < Spline.length_x; mx++){
        for(int my = 0; my < Spline.length_y; my++){
            
            bi = 2.0 / (Spline.z_vals[1] - Spline.z_vals[0]);
            ci = 1.0 / (Spline.z_vals[1] - Spline.z_vals[0]);
            di = 3.0 * (dfdx[mx][my][1] - dfdx[mx][my][0]) / pow(Spline.z_vals[1] - Spline.z_vals[0], 2);
            
            new_c[0] = ci/bi;
            new_d[0] = di/bi;
            
            for(int i = 1; i < Spline.length_z - 1; i++) {
                ai = 1.0/(Spline.z_vals[i] - Spline.z_vals[i-1]);
                bi = 2.0 * (1.0/(Spline.z_vals[i] - Spline.z_vals[i-1]) + 1.0/(Spline.z_vals[i+1] - Spline.z_vals[i]));
                ci = 1.0/(Spline.z_vals[i+1] - Spline.z_vals[i]);
                di = 3.0 * ((dfdx[mx][my][i] - dfdx[mx][my][i-1]) / pow(Spline.z_vals[i] - Spline.z_vals[i-1], 2)
                            + (dfdx[mx][my][i+1] - dfdx[mx][my][i]) / pow(Spline.z_vals[i+1] - Spline.z_vals[i], 2) );
                
                new_c[i] = ci/(bi - new_c[i-1]*ai);
                new_d[i] = (di - new_d[i-1]*ai)/(bi - new_c[i-1]*ai);
            }
            
            ai = 1.0/(Spline.z_vals[Spline.length_z-1] - Spline.z_vals[Spline.length_z-2]);
            bi = 2.0/(Spline.z_vals[Spline.length_z-1] - Spline.z_vals[Spline.length_z-2]);
            di = 3.0 * (dfdx[mx][my][Spline.length_z-1] - dfdx[mx][my][Spline.length_z-2])
            / pow(Spline.z_vals[Spline.length_z-1] - Spline.z_vals[Spline.length_z-2], 2);
            
            new_d[Spline.length_z-1] = (di - new_d[Spline.length_z - 2]*ai)/(bi - new_c[Spline.length_z - 2]*ai);
            
            Spline.dfdx_slopes[mx][my][Spline.length_z - 1] = new_d[Spline.length_z - 1];
            for(int i = Spline.length_z - 2; i >= 0; i--) Spline.dfdx_slopes[mx][my][i] = new_d[i] - new_c[i] * Spline.dfdx_slopes[mx][my][i+1];
        }}
    
    // Set dfdy_slopes
    for(int mx = 0; mx < Spline.length_x; mx++){
        for(int my = 0; my < Spline.length_y; my++){
            
            bi = 2.0 / (Spline.z_vals[1] - Spline.z_vals[0]);
            ci = 1.0 / (Spline.z_vals[1] - Spline.z_vals[0]);
            di = 3.0 * (dfdy[mx][my][1] - dfdy[mx][my][0]) / pow(Spline.z_vals[1] - Spline.z_vals[0], 2);
            
            new_c[0] = ci/bi;
            new_d[0] = di/bi;
            
            for(int i = 1; i < Spline.length_z - 1; i++) {
                ai = 1.0/(Spline.z_vals[i] - Spline.z_vals[i-1]);
                bi = 2.0 * (1.0/(Spline.z_vals[i] - Spline.z_vals[i-1]) + 1.0/(Spline.z_vals[i+1] - Spline.z_vals[i]));
                ci = 1.0/(Spline.z_vals[i+1] - Spline.z_vals[i]);
                di = 3.0 * ((dfdy[mx][my][i] - dfdy[mx][my][i-1]) / pow(Spline.z_vals[i] - Spline.z_vals[i-1], 2)
                            + (dfdy[mx][my][i+1] - dfdy[mx][my][i]) / pow(Spline.z_vals[i+1] - Spline.z_vals[i], 2) );
                
                new_c[i] = ci/(bi - new_c[i-1]*ai);
                new_d[i] = (di - new_d[i-1]*ai)/(bi - new_c[i-1]*ai);
            }
            
            ai = 1.0/(Spline.z_vals[Spline.length_z-1] - Spline.z_vals[Spline.length_z-2]);
            bi = 2.0/(Spline.z_vals[Spline.length_z-1] - Spline.z_vals[Spline.length_z-2]);
            di = 3.0 * (dfdy[mx][my][Spline.length_z-1] - dfdy[mx][my][Spline.length_z-2])
            / pow(Spline.z_vals[Spline.length_z-1] - Spline.z_vals[Spline.length_z-2], 2);
            
            new_d[Spline.length_z-1] = (di - new_d[Spline.length_z - 2]*ai)/(bi - new_c[Spline.length_z - 2]*ai);
            
            Spline.dfdy_slopes[mx][my][Spline.length_z - 1] = new_d[Spline.length_z - 1];
            for(int i = Spline.length_z - 2; i >= 0; i--) Spline.dfdy_slopes[mx][my][i] = new_d[i] - new_c[i] * Spline.dfdy_slopes[mx][my][i+1];
        }}
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


double Eval_Vert_Spline_f(double z, struct MultiDimSpline_3D Spline, int kx, int ky, int kz){
    double X = (z - Spline.z_vals[kz])/(Spline.z_vals[kz+1] - Spline.z_vals[kz]);
    double A = Spline.f_slopes[kx][ky][kz] * (Spline.z_vals[kz+1] - Spline.z_vals[kz]) - (Spline.f_vals[kx][ky][kz+1] - Spline.f_vals[kx][ky][kz]);
    double B = -Spline.f_slopes[kx][ky][kz+1] * (Spline.z_vals[kz+1] - Spline.z_vals[kz]) + (Spline.f_vals[kx][ky][kz+1] - Spline.f_vals[kx][ky][kz]);
    
    return (1.0 - X) * Spline.f_vals[kx][ky][kz] + X * Spline.f_vals[kx][ky][kz+1] + X * (1.0 - X) * (A * (1.0 - X ) + B * X);
}

double Eval_Vert_Spline_dfdx(double z, struct MultiDimSpline_3D Spline, int kx, int ky, int kz){
    int kx_up = min(kx + 1, Spline.length_x - 1);
    int kx_dn = max(kx - 1, 0);
    
    double dfdx_kz = (Spline.f_vals[kx_up][ky][kz] - Spline.f_vals[kx_dn][ky][kz])/(Spline.x_vals[kx_up] - Spline.x_vals[kx_dn]);           // df/dx at kz
    double dfdx_kzp1 = (Spline.f_vals[kx_up][ky][kz+1] - Spline.f_vals[kx_dn][ky][kz+1])/(Spline.x_vals[kx_up] - Spline.x_vals[kx_dn]);     // df/dx at kz + 1
    
    double X = (z - Spline.z_vals[kz])/(Spline.z_vals[kz+1] - Spline.z_vals[kz]);
    double A = Spline.dfdx_slopes[kx][ky][kz] * (Spline.z_vals[kz+1] - Spline.z_vals[kz]) - (dfdx_kzp1 - dfdx_kz);
    double B = -Spline.dfdx_slopes[kx][ky][kz+1] * (Spline.z_vals[kz+1] - Spline.z_vals[kz]) + (dfdx_kzp1 - dfdx_kz);
    
    return (1.0 - X) * dfdx_kz + X * dfdx_kzp1 + X * (1.0 - X) * (A * (1.0 - X ) + B * X);
}

double Eval_Vert_Spline_dfdy(double z, struct MultiDimSpline_3D Spline, int kx, int ky, int kz){
    int ky_up = min(ky + 1, Spline.length_y - 1);
    int ky_dn = max(ky - 1, 0);
    
    double dfdy_kz = (Spline.f_vals[kx][ky_up][kz] - Spline.f_vals[kx][ky_dn][kz])/(Spline.y_vals[ky_up] - Spline.y_vals[ky_dn]);       // df/dy at kz
    double dfdy_kzp1 = (Spline.f_vals[kx][ky_up][kz+1] - Spline.f_vals[kx][ky_dn][kz+1])/(Spline.y_vals[ky_up] - Spline.y_vals[ky_dn]); // df/dy at kz + 1
    
    double X = (z - Spline.z_vals[kz])/(Spline.z_vals[kz+1] - Spline.z_vals[kz]);
    double A = Spline.dfdy_slopes[kx][ky][kz] * (Spline.z_vals[kz+1] - Spline.z_vals[kz]) - (dfdy_kzp1 - dfdy_kz);
    double B = -Spline.dfdy_slopes[kx][ky][kz+1] * (Spline.z_vals[kz+1] - Spline.z_vals[kz]) + (dfdy_kzp1 - dfdy_kz);
    
    return (1.0 - X) * dfdy_kz + X * dfdy_kzp1 + X * (1.0 - X) * (A * (1.0 - X ) + B * X);
}

double Eval_Vert_Spline_dfdz(double z, struct MultiDimSpline_3D Spline, int kx, int ky, int kz){
    double X = (z - Spline.z_vals[kz])/(Spline.z_vals[kz+1] - Spline.z_vals[kz]);
    double A = Spline.f_slopes[kx][ky][kz] * (Spline.z_vals[kz+1] - Spline.z_vals[kz]) - (Spline.f_vals[kx][ky][kz+1] - Spline.f_vals[kx][ky][kz]);
    double B = -Spline.f_slopes[kx][ky][kz+1] * (Spline.z_vals[kz+1] - Spline.z_vals[kz]) + (Spline.f_vals[kx][ky][kz+1] - Spline.f_vals[kx][ky][kz]);
    
    return (Spline.f_vals[kx][ky][kz+1] - Spline.f_vals[kx][ky][kz])/(Spline.z_vals[kz+1] - Spline.z_vals[kz])
    + (1.0 - 2.0 * X) * (A * (1.0 - X) + B * X)/(Spline.z_vals[kz+1] - Spline.z_vals[kz])
    + X * (1.0 - X) * (B - A)/(Spline.z_vals[kz+1] - Spline.z_vals[kz]);
}

double Eval_Vert_Spline_ddfdxdz(double z, struct MultiDimSpline_3D Spline, int kx, int ky, int kz){
    int kx_up = min(kx + 1, Spline.length_x - 1);
    int kx_dn = max(kx - 1, 0);
    
    double dfdx_kz = (Spline.f_vals[kx_up][ky][kz] - Spline.f_vals[kx_dn][ky][kz])/(Spline.x_vals[kx_up] - Spline.x_vals[kx_dn]);           // df/dx at kz
    double dfdx_kzp1 = (Spline.f_vals[kx_up][ky][kz+1] - Spline.f_vals[kx_dn][ky][kz+1])/(Spline.x_vals[kx_up] - Spline.x_vals[kx_dn]);     // df/dx at kz + 1
    
    double X = (z - Spline.z_vals[kz])/(Spline.z_vals[kz+1] - Spline.z_vals[kz]);
    double A = Spline.dfdx_slopes[kx][ky][kz] * (Spline.z_vals[kz+1] - Spline.z_vals[kz]) - (dfdx_kzp1 - dfdx_kz);
    double B = -Spline.dfdx_slopes[kx][ky][kz+1] * (Spline.z_vals[kz+1] - Spline.z_vals[kz]) + (dfdx_kzp1 - dfdx_kz);
    
    return (dfdx_kzp1 - dfdx_kz)/(Spline.z_vals[kz+1] - Spline.z_vals[kz])
    + (1.0 - 2.0 * X) * (A * (1.0 - X) + B * X)/(Spline.z_vals[kz+1] - Spline.z_vals[kz])
    + X * (1.0 - X) * (B - A)/(Spline.z_vals[kz+1] - Spline.z_vals[kz]);
    
}

double Eval_Vert_Spline_ddfdydz(double z, struct MultiDimSpline_3D Spline, int kx, int ky, int kz){
    int ky_up = min(ky + 1, Spline.length_y - 1);
    int ky_dn = max(ky - 1, 0);
    
    double dfdy_kz = (Spline.f_vals[kx][ky_up][kz] - Spline.f_vals[kx][ky_dn][kz])/(Spline.y_vals[ky_up] - Spline.y_vals[ky_dn]);       // df/dy at kz
    double dfdy_kzp1 = (Spline.f_vals[kx][ky_up][kz+1] - Spline.f_vals[kx][ky_dn][kz+1])/(Spline.y_vals[ky_up] - Spline.y_vals[ky_dn]); // df/dy at kz + 1
    
    double X = (z - Spline.z_vals[kz])/(Spline.z_vals[kz+1] - Spline.z_vals[kz]);
    double A = Spline.dfdy_slopes[kx][ky][kz] * (Spline.z_vals[kz+1] - Spline.z_vals[kz]) - (dfdy_kzp1 - dfdy_kz);
    double B = -Spline.dfdy_slopes[kx][ky][kz+1] * (Spline.z_vals[kz+1] - Spline.z_vals[kz]) + (dfdy_kzp1 - dfdy_kz);
    
    return (dfdy_kzp1 - dfdy_kz)/(Spline.z_vals[kz+1] - Spline.z_vals[kz])
    + (1.0 - 2.0 * X) * (A * (1.0 - X) + B * X)/(Spline.z_vals[kz+1] - Spline.z_vals[kz])
    + X * (1.0 - X) * (B - A)/(Spline.z_vals[kz+1] - Spline.z_vals[kz]);
    
}

double Eval_Vert_Spline_ddfdzdz(double z, struct MultiDimSpline_3D Spline, int kx, int ky, int kz){
    double X = (z - Spline.z_vals[kz])/(Spline.z_vals[kz+1] - Spline.z_vals[kz]);
    double A = Spline.f_slopes[kx][ky][kz] * (Spline.z_vals[kz+1] - Spline.z_vals[kz]) - (Spline.f_vals[kx][ky][kz+1] - Spline.f_vals[kx][ky][kz]);
    double B = -Spline.f_slopes[kx][ky][kz+1] * (Spline.z_vals[kz+1] - Spline.z_vals[kz]) + (Spline.f_vals[kx][ky][kz+1] - Spline.f_vals[kx][ky][kz]);
    
    return 2.0 * (B - 2.0 * A + (A - B) * 3.0 * X)/pow(Spline.z_vals[kz+1] - Spline.z_vals[kz],2);
}

//----------------------------------------------------//
//-----------Evauation of the 2D Spline and-----------//
//-------its First and Second Order Derivatives-------//
//----------------------------------------------------//
double BiCubic_Deriv_dfdx(double z, struct MultiDimSpline_3D Spline, int kx, int ky, int kz){
    // Calculate discrete d/dx of f to obtain df/dx
    int kx_up = kx + 1;
    int kx_dn = kx - 1;
    
	if(kx_up > Spline.length_x - 1) kx_up = kx;           // Use backward x derivative if kx+1 is not defined
    if(kx_dn < 0)                   kx_dn = kx;           // Use forward x derivative if kx-1 is not defined
    
    return (Eval_Vert_Spline_f(z, Spline, kx_up, ky, kz) - Eval_Vert_Spline_f(z, Spline, kx_dn, ky, kz))
                /(Spline.x_vals[kx_up] - Spline.x_vals[kx_dn]);
}

double BiCubic_Deriv_dfdy(double z, struct MultiDimSpline_3D Spline, int kx, int ky, int kz){
    // Calculate discrete d/dy of f to obtain df/dy
    int ky_up = ky + 1;
    int ky_dn = ky - 1;
    
	if(ky_up > Spline.length_y - 1) ky_up = ky;           // Use backward y derivative if ky+1 is not defined
    if(ky_dn < 0)                   ky_dn = ky;           // Use forward y derivative if ky-1 is not defined
    
    return (Eval_Vert_Spline_f(z, Spline, kx, ky_up, kz) - Eval_Vert_Spline_f(z, Spline, kx, ky_dn, kz))
                /(Spline.y_vals[ky_up] - Spline.y_vals[ky_dn]);
    
}

double BiCubic_Deriv_ddfdxdx(double z, struct MultiDimSpline_3D Spline, int kx, int ky, int kz){
    // Calculate discrete d/dx of df/dx to obtain d^2f/dx^2
    int kx_up = kx + 1;
    int kx_dn = kx - 1;
    
	if(kx_up > Spline.length_x - 1) kx_up = kx;           // Use backward x derivative if kx+1 is not defined
    if(kx_dn < 0)                   kx_dn = kx;           // Use forward x derivative if kx-1 is not defined
    
    return (Eval_Vert_Spline_dfdx(z, Spline, kx_up, ky, kz) - Eval_Vert_Spline_dfdx(z, Spline, kx_dn, ky, kz))
                /(Spline.x_vals[kx_up] - Spline.x_vals[kx_dn]);
}

double BiCubic_Deriv_ddfdydy(double z, struct MultiDimSpline_3D Spline, int kx, int ky, int kz){
    // Calculate discrete d/dy of df/dy to obtain d^2f/dy^2
    int ky_up = ky + 1;
    int ky_dn = ky - 1;
    
	if(ky_up > Spline.length_y - 1) ky_up = ky;           // Use backward y derivative if ky+1 is not defined
    if(ky_dn < 0)                   ky_dn = ky;           // Use forward y derivative if ky-1 is not defined
    
    return (Eval_Vert_Spline_dfdy(z, Spline, kx, ky_up, kz) - Eval_Vert_Spline_dfdy(z, Spline, kx, ky_dn, kz))
                /(Spline.y_vals[ky_up] - Spline.y_vals[ky_dn]);
    
}

double BiCubic_Deriv_ddfdxdy(double z, struct MultiDimSpline_3D Spline, int kx, int ky, int kz){
    // Calculate discrete d^2/dxdy of f to obtain d^2f/dxdy
    int kx_up = kx+1;    int kx_dn = kx-1;
    int ky_up = ky+1;    int ky_dn = ky-1;
    
	if(kx_up > Spline.length_x - 1) kx_up = kx;           // Use backward x derivative if kx+1 is not defined
    if(kx_dn < 0)                   kx_dn = kx;           // Use forward x derivative if kx-1 is not defined
    
	if(ky_up > Spline.length_y - 1) ky_up = ky;           // Use backward y derivative if ky+1 is not defined
    if(ky_dn < 0)                   ky_dn = ky;           // Use forward y derivative if ky-1 is not defined
    
    return (Eval_Vert_Spline_f(z, Spline, kx_up, ky_up, kz)
            - Eval_Vert_Spline_f(z, Spline, kx_up, ky_dn, kz)
                - Eval_Vert_Spline_f(z, Spline, kx_dn, ky_up, kz)
                    + Eval_Vert_Spline_f(z, Spline, kx_dn, ky_dn, kz))
                        /((Spline.x_vals[kx_up] - Spline.x_vals[kx_dn])*(Spline.y_vals[ky_up] - Spline.y_vals[ky_dn]));
    
}

double BiCubic_Deriv_dddfdxdxdy(double z, struct MultiDimSpline_3D Spline, int kx, int ky, int kz){
    // Calculate discrete d^2/dxdy of df/dx to obtain d^3f/dx^2dy
    int kx_up = kx+1;    int kx_dn = kx-1;
    int ky_up = ky+1;    int ky_dn = ky-1;
    
	if(kx_up > Spline.length_x - 1) kx_up = kx;           // Use backward x derivative if kx+1 is not defined
    if(kx_dn < 0)                   kx_dn = kx;           // Use forward x derivative if kx-1 is not defined
    
	if(ky_up > Spline.length_y - 1) ky_up = ky;           // Use backward y derivative if ky+1 is not defined
    if(ky_dn < 0)                   ky_dn = ky;           // Use forward y derivative if ky-1 is not defined
    
    return (Eval_Vert_Spline_dfdx(z, Spline, kx_up, ky_up, kz)
            - Eval_Vert_Spline_dfdx(z, Spline, kx_up, ky_dn, kz)
                - Eval_Vert_Spline_dfdx(z, Spline, kx_dn, ky_up, kz)
                    + Eval_Vert_Spline_dfdx(z, Spline, kx_dn, ky_dn, kz))
                        /((Spline.x_vals[kx_up] - Spline.x_vals[kx_dn])*(Spline.y_vals[ky_up] - Spline.y_vals[ky_dn]));
}

double BiCubic_Deriv_dddfdxdxdz(double z, struct MultiDimSpline_3D Spline, int kx, int ky, int kz){
    // Calculate discrete d/dx of d^2f/dxdz to obtain d^3f/dx^2dz
    int kx_up = kx+1;    int kx_dn = kx-1;
    int ky_up = ky+1;    int ky_dn = ky-1;
    
	if(kx_up > Spline.length_x - 1) kx_up = kx;           // Use backward x derivative if kx+1 is not defined
    if(kx_dn < 0)                   kx_dn = kx;           // Use forward x derivative if kx-1 is not defined
    
    return (Eval_Vert_Spline_ddfdxdz(z, Spline, kx_up, ky, kz) - Eval_Vert_Spline_ddfdxdz(z, Spline, kx_dn, ky, kz))
                /(Spline.x_vals[kx_up] - Spline.x_vals[kx_dn]);
}
double BiCubic_Deriv_dddfdxdydy(double z, struct MultiDimSpline_3D Spline, int kx, int ky, int kz){
    // Calculate discrete d^2/dxdy of df/dy to obtain d^3f/dxdy^2
    int kx_up = kx+1;    int kx_dn = kx-1;
    int ky_up = ky+1;    int ky_dn = ky-1;
    
	if(kx_up > Spline.length_x - 1) kx_up = kx;           // Use backward x derivative if kx+1 is not defined
    if(kx_dn < 0)                   kx_dn = kx;           // Use forward x derivative if kx-1 is not defined
    
	if(ky_up > Spline.length_y - 1) ky_up = ky;           // Use backward y derivative if ky+1 is not defined
    if(ky_dn < 0)                   ky_dn = ky;           // Use forward y derivative if ky-1 is not defined
    
    return (Eval_Vert_Spline_dfdy(z, Spline, kx_up, ky_up, kz)
            - Eval_Vert_Spline_dfdy(z, Spline, kx_up, ky_dn, kz)
                - Eval_Vert_Spline_dfdy(z, Spline, kx_dn, ky_up, kz)
                    + Eval_Vert_Spline_dfdy(z, Spline, kx_dn, ky_dn, kz))
                        /((Spline.x_vals[kx_up] - Spline.x_vals[kx_dn])*(Spline.y_vals[ky_up] - Spline.y_vals[ky_dn]));
}

double BiCubic_Deriv_dddfdxdydz(double z, struct MultiDimSpline_3D Spline, int kx, int ky, int kz){
    // Calculate discrete d^2/dxdy of df/dz to obtain d^3f/dxdydz
    int kx_up = kx+1;    int kx_dn = kx-1;
    int ky_up = ky+1;    int ky_dn = ky-1;
    
	if(kx_up > Spline.length_x - 1) kx_up = kx;           // Use backward x derivative if kx+1 is not defined
    if(kx_dn < 0)                   kx_dn = kx;           // Use forward x derivative if kx-1 is not defined
    
	if(ky_up > Spline.length_y - 1) ky_up = ky;           // Use backward y derivative if ky+1 is not defined
    if(ky_dn < 0)                   ky_dn = ky;           // Use forward y derivative if ky-1 is not defined
    
    return (Eval_Vert_Spline_dfdz(z, Spline, kx_up, ky_up, kz)
            - Eval_Vert_Spline_dfdz(z, Spline, kx_up, ky_dn, kz)
                - Eval_Vert_Spline_dfdz(z, Spline, kx_dn, ky_up, kz)
                    + Eval_Vert_Spline_dfdz(z, Spline, kx_dn, ky_dn, kz))
                        /((Spline.x_vals[kx_up] - Spline.x_vals[kx_dn])*(Spline.y_vals[ky_up] - Spline.y_vals[ky_dn]));
}


double BiCubic_Deriv_dddfdxdzdz(double z, struct MultiDimSpline_3D Spline, int kx, int ky, int kz){
    // Calculate discrete d/dx of d^2f/dz^2 to obtain d^3f/dxdz^2
    int kx_up = kx + 1;
    int kx_dn = kx - 1;
    
	if(kx_up > Spline.length_x - 1) kx_up = kx;           // Use backward x derivative if kx+1 is not defined
    if(kx_dn < 0)                   kx_dn = kx;           // Use forward x derivative if kx-1 is not defined
    
    return (Eval_Vert_Spline_ddfdzdz(z, Spline, kx_up, ky, kz) - Eval_Vert_Spline_ddfdzdz(z, Spline, kx_dn, ky, kz))
                /(Spline.x_vals[kx_up] - Spline.x_vals[kx_dn]);
}

double BiCubic_Deriv_dddfdydydz(double z, struct MultiDimSpline_3D Spline, int kx, int ky, int kz){
    // Calculate discrete d/dy of d^2f/dydz to obtain d^3f/dy^2dz
    int ky_up = ky + 1;
    int ky_dn = ky - 1;
    
	if(ky_up > Spline.length_y - 1) ky_up = ky;           // Use backward y derivative if ky+1 is not defined
    if(ky_dn < 0)                   ky_dn = ky;           // Use forward y derivative if ky-1 is not defined
    
    return (Eval_Vert_Spline_ddfdydz(z, Spline, kx, ky_up, kz) - Eval_Vert_Spline_ddfdydz(z, Spline, kx, ky_dn, kz))
                /(Spline.y_vals[ky_up] - Spline.y_vals[ky_dn]);
    
}

double BiCubic_Deriv_dddfdydzdz(double z, struct MultiDimSpline_3D Spline, int kx, int ky, int kz){
    // Calculate discrete d/dy of d^2f/dz^2 to obtain d^3f/dydz^2
    int ky_up = ky + 1;
    int ky_dn = ky - 1;
    
	if(ky_up > Spline.length_y - 1) ky_up = ky;           // Use backward y derivative if ky+1 is not defined
    if(ky_dn < 0)                   ky_dn = ky;           // Use forward y derivative if ky-1 is not defined
    
    return (Eval_Vert_Spline_ddfdzdz(z, Spline, kx, ky_up, kz) - Eval_Vert_Spline_ddfdzdz(z, Spline, kx, ky_dn, kz))
                /(Spline.y_vals[ky_up] - Spline.y_vals[ky_dn]);
    
}

double BiCubic_Deriv_ddddfdxdxdydz(double z, struct MultiDimSpline_3D Spline, int kx, int ky, int kz){
    // Calculate discrete d^2/dxdy of d^2f/dxdz to obtain d^4f/dx^2dydz
    
    int kx_up = kx+1;    int kx_dn = kx-1;
    int ky_up = ky+1;    int ky_dn = ky-1;
    
	if(kx_up > Spline.length_x - 1) kx_up = kx;           // Use backward x derivative if kx+1 is not defined
    if(kx_dn < 0)                   kx_dn = kx;           // Use forward x derivative if kx-1 is not defined
    
	if(ky_up > Spline.length_y - 1) ky_up = ky;           // Use backward y derivative if ky+1 is not defined
    if(ky_dn < 0)                   ky_dn = ky;           // Use forward y derivative if ky-1 is not defined
    
    return (Eval_Vert_Spline_ddfdxdz(z, Spline, kx_up, ky_up, kz)
            - Eval_Vert_Spline_ddfdxdz(z, Spline, kx_up, ky_dn, kz)
                - Eval_Vert_Spline_ddfdxdz(z, Spline, kx_dn, ky_up, kz)
                    + Eval_Vert_Spline_ddfdxdz(z, Spline, kx_dn, ky_dn, kz))
                        /((Spline.x_vals[kx_up] - Spline.x_vals[kx_dn])*(Spline.y_vals[ky_up] - Spline.y_vals[ky_dn]));
}


double BiCubic_Deriv_ddddfdxdydydz(double z, struct MultiDimSpline_3D Spline, int kx, int ky, int kz){
    // Calculate discrete d^2/dxdy of d^2f/dydz to obtain d^4f/dxdy^2dz
    
    int kx_up = kx+1;    int kx_dn = kx-1;
    int ky_up = ky+1;    int ky_dn = ky-1;
    
	if(kx_up > Spline.length_x - 1) kx_up = kx;           // Use backward x derivative if kx+1 is not defined
    if(kx_dn < 0)                   kx_dn = kx;           // Use forward x derivative if kx-1 is not defined
    
	if(ky_up > Spline.length_y - 1) ky_up = ky;           // Use backward y derivative if ky+1 is not defined
    if(ky_dn < 0)                   ky_dn = ky;           // Use forward y derivative if ky-1 is not defined
    
    return (Eval_Vert_Spline_ddfdydz(z, Spline, kx_up, ky_up, kz)
            - Eval_Vert_Spline_ddfdydz(z, Spline, kx_up, ky_dn, kz)
                - Eval_Vert_Spline_ddfdydz(z, Spline, kx_dn, ky_up, kz)
                    + Eval_Vert_Spline_ddfdydz(z, Spline, kx_dn, ky_dn, kz))
                        /((Spline.x_vals[kx_up] - Spline.x_vals[kx_dn])*(Spline.y_vals[ky_up] - Spline.y_vals[ky_dn]));
}




double BiCubic_Deriv_ddddfdxdydzdz(double z, struct MultiDimSpline_3D Spline, int kx, int ky, int kz){
    // Calculate discrete d^2/dxdy of d^2f/dz^2 to obtain d^4f/dxdydz^2
    
    int kx_up = kx+1;    int kx_dn = kx-1;
    int ky_up = ky+1;    int ky_dn = ky-1;
    
	if(kx_up > Spline.length_x - 1) kx_up = kx;           // Use backward x derivative if kx+1 is not defined
    if(kx_dn < 0)                   kx_dn = kx;           // Use forward x derivative if kx-1 is not defined
    
	if(ky_up > Spline.length_y - 1) ky_up = ky;           // Use backward y derivative if ky+1 is not defined
    if(ky_dn < 0)                   ky_dn = ky;           // Use forward y derivative if ky-1 is not defined
    
    return (Eval_Vert_Spline_ddfdzdz(z, Spline, kx_up, ky_up, kz)
            - Eval_Vert_Spline_ddfdzdz(z, Spline, kx_up, ky_dn, kz)
                - Eval_Vert_Spline_ddfdzdz(z, Spline, kx_dn, ky_up, kz)
                    + Eval_Vert_Spline_ddfdzdz(z, Spline, kx_dn, ky_dn, kz))
                        /((Spline.x_vals[kx_up] - Spline.x_vals[kx_dn])*(Spline.y_vals[ky_up] - Spline.y_vals[ky_dn]));
}

//----------------------------------------------------//
//-----------Evauation of the 2D Spline and-----------//
//-------its First and Second Order Derivatives-------//
//----------------------------------------------------//
double Eval_Spline_f(double x, double y, double z, struct MultiDimSpline_3D & Spline){
    double A_vec[16];
	double X_vec[16];
    double BicCoeff[4][4];
    
    // Locate kx, ky, kz such that (x[kx] <= x' <= x[kx+1]), (y[ky] <= y' <= y[ky+1]), and (z[kz] <= z' <= z[kz+1])
    int kx = Find_Segment(x, Spline.x_vals, Spline.length_x, Spline.accel[0]);
    int ky = Find_Segment(y, Spline.y_vals, Spline.length_y, Spline.accel[1]);
    int kz = Find_Segment(z, Spline.z_vals, Spline.length_z, Spline.accel[2]);
    
    double dx_scalar = Spline.x_vals[kx+1] - Spline.x_vals[kx]; // Set dr scaling to produce unit square relations
    double dy_scalar = Spline.y_vals[ky+1] - Spline.y_vals[ky]; // Set dz scaling to produce unit square relations
    
    X_vec[0] = Eval_Vert_Spline_f(z, Spline, kx, ky, kz);
    X_vec[1] = Eval_Vert_Spline_f(z, Spline, kx+1, ky, kz);
    X_vec[2] = Eval_Vert_Spline_f(z, Spline, kx, ky+1, kz);
    X_vec[3] = Eval_Vert_Spline_f(z, Spline, kx+1, ky+1, kz);
    
    X_vec[4] = BiCubic_Deriv_dfdx(z, Spline, kx, ky,kz)*dx_scalar;
    X_vec[5] = BiCubic_Deriv_dfdx(z, Spline, kx+1, ky,kz)*dx_scalar;
    X_vec[6] = BiCubic_Deriv_dfdx(z, Spline, kx, ky+1,kz)*dx_scalar;
    X_vec[7] = BiCubic_Deriv_dfdx(z, Spline, kx+1, ky+1,kz)*dx_scalar;
    
    X_vec[8] =  BiCubic_Deriv_dfdy(z, Spline, kx, ky,kz)*dx_scalar;
    X_vec[9] =  BiCubic_Deriv_dfdy(z, Spline, kx+1, ky,kz)*dx_scalar;
    X_vec[10] = BiCubic_Deriv_dfdy(z, Spline, kx, ky+1,kz)*dx_scalar;
    X_vec[11] = BiCubic_Deriv_dfdy(z, Spline, kx+1, ky+1,kz)*dx_scalar;
    
    X_vec[12]   = BiCubic_Deriv_ddfdxdy(z, Spline, kx, ky, kz)*dx_scalar*dy_scalar;
    X_vec[13]   = BiCubic_Deriv_ddfdxdy(z, Spline, kx+1, ky, kz)*dx_scalar*dy_scalar;
    X_vec[14]	= BiCubic_Deriv_ddfdxdy(z, Spline, kx, ky+1, kz)*dx_scalar*dy_scalar;
    X_vec[15]   = BiCubic_Deriv_ddfdxdy(z, Spline, kx+1, ky+1, kz)*dx_scalar*dy_scalar;
    
    // Set the coefficients by matrix multiplication
    for(int j = 0; j < 16; j++){
        A_vec[j] = 0;
        for(int k = 0; k < 16; k++){
            A_vec[j] += BiCubic_ConversionMatrix[j][k]*X_vec[k];
        }
    }
    
    // Copy the a_vec elements into the coefficients array
    for(int k1=0; k1 < 4; k1++){
    for(int k2=0; k2 < 4; k2++){
        BicCoeff[k1][k2] = A_vec[k2*4 + k1];
    }}
    
    double x_scaled = (x - Spline.x_vals[kx])/(Spline.x_vals[kx+1] - Spline.x_vals[kx]);
	double y_scaled = (y - Spline.y_vals[ky])/(Spline.y_vals[ky+1] - Spline.y_vals[ky]);
    
	double result = 0;
	for(int k1 = 0; k1 < 4; k1++){
    for(int k2 = 0; k2 < 4; k2++){
        result+=1.0*BicCoeff[k1][k2]*pow(x_scaled,k1)*pow(y_scaled,k2);
    }}
	
    return result;
}

double Eval_Spline_df(double x, double y, double z, int index, struct MultiDimSpline_3D & Spline){
    double A_vec[16];
	double X_vec[16];
    double BicCoeff[4][4];
    
    // Locate kx, ky, kz such that (x[kx] <= x' <= x[kx+1]), (y[ky] <= y' <= y[ky+1]), and (z[kz] <= z' <= z[kz+1])
    int kx = Find_Segment(x, Spline.x_vals, Spline.length_x, Spline.accel[0]);
    int ky = Find_Segment(y, Spline.y_vals, Spline.length_y, Spline.accel[1]);
    int kz = Find_Segment(z, Spline.z_vals, Spline.length_z, Spline.accel[2]);
    
    double dx_scalar = Spline.x_vals[kx+1] - Spline.x_vals[kx]; // Set dx scaling to produce unit square relations
    double dy_scalar = Spline.y_vals[ky+1] - Spline.y_vals[ky]; // Set dy scaling to produce unit square relations
    
    if(index==0){
        X_vec[0] = BiCubic_Deriv_dfdx(z, Spline, kx, ky, kz);
        X_vec[1] = BiCubic_Deriv_dfdx(z, Spline, kx+1, ky, kz);
        X_vec[2] = BiCubic_Deriv_dfdx(z, Spline, kx, ky+1, kz);
        X_vec[3] = BiCubic_Deriv_dfdx(z, Spline, kx+1, ky+1, kz);
        
        X_vec[4] = BiCubic_Deriv_ddfdxdx(z, Spline, kx, ky,kz)*dx_scalar;
        X_vec[5] = BiCubic_Deriv_ddfdxdx(z, Spline, kx+1, ky,kz)*dx_scalar;
        X_vec[6] = BiCubic_Deriv_ddfdxdx(z, Spline, kx, ky+1,kz)*dx_scalar;
        X_vec[7] = BiCubic_Deriv_ddfdxdx(z, Spline, kx+1, ky+1,kz)*dx_scalar;
        
        X_vec[8] =  BiCubic_Deriv_ddfdxdy(z, Spline, kx, ky,kz)*dx_scalar;
        X_vec[9] =  BiCubic_Deriv_ddfdxdy(z, Spline, kx+1, ky,kz)*dx_scalar;
        X_vec[10] = BiCubic_Deriv_ddfdxdy(z, Spline, kx, ky+1,kz)*dx_scalar;
        X_vec[11] = BiCubic_Deriv_ddfdxdy(z, Spline, kx+1, ky+1,kz)*dx_scalar;
        
        X_vec[12]   = BiCubic_Deriv_dddfdxdxdy(z, Spline, kx, ky, kz)*dx_scalar*dy_scalar;
        X_vec[13]   = BiCubic_Deriv_dddfdxdxdy(z, Spline, kx+1, ky, kz)*dx_scalar*dy_scalar;
        X_vec[14]	= BiCubic_Deriv_dddfdxdxdy(z, Spline, kx, ky+1, kz)*dx_scalar*dy_scalar;
        X_vec[15]   = BiCubic_Deriv_dddfdxdxdy(z, Spline, kx+1, ky+1, kz)*dx_scalar*dy_scalar;
        
    } else if(index==1){
        X_vec[0] = BiCubic_Deriv_dfdy(z, Spline, kx, ky, kz);
        X_vec[1] = BiCubic_Deriv_dfdy(z, Spline, kx+1, ky, kz);
        X_vec[2] = BiCubic_Deriv_dfdy(z, Spline, kx, ky+1, kz);
        X_vec[3] = BiCubic_Deriv_dfdy(z, Spline, kx+1, ky+1, kz);
        
        X_vec[4] = BiCubic_Deriv_ddfdxdy(z, Spline, kx, ky,kz)*dx_scalar;
        X_vec[5] = BiCubic_Deriv_ddfdxdy(z, Spline, kx+1, ky,kz)*dx_scalar;
        X_vec[6] = BiCubic_Deriv_ddfdxdy(z, Spline, kx, ky+1,kz)*dx_scalar;
        X_vec[7] = BiCubic_Deriv_ddfdxdy(z, Spline, kx+1, ky+1,kz)*dx_scalar;
        
        X_vec[8] =  BiCubic_Deriv_ddfdydy(z, Spline, kx, ky,kz)*dx_scalar;
        X_vec[9] =  BiCubic_Deriv_ddfdydy(z, Spline, kx+1, ky,kz)*dx_scalar;
        X_vec[10] = BiCubic_Deriv_ddfdydy(z, Spline, kx, ky+1,kz)*dx_scalar;
        X_vec[11] = BiCubic_Deriv_ddfdydy(z, Spline, kx+1, ky+1,kz)*dx_scalar;
        
        X_vec[12]   = BiCubic_Deriv_dddfdxdydy(z, Spline, kx, ky, kz)*dx_scalar*dy_scalar;
        X_vec[13]   = BiCubic_Deriv_dddfdxdydy(z, Spline, kx+1, ky, kz)*dx_scalar*dy_scalar;
        X_vec[14]	= BiCubic_Deriv_dddfdxdydy(z, Spline, kx, ky+1, kz)*dx_scalar*dy_scalar;
        X_vec[15]   = BiCubic_Deriv_dddfdxdydy(z, Spline, kx+1, ky+1, kz)*dx_scalar*dy_scalar;
        
    } else if(index==2){
        X_vec[0] = Eval_Vert_Spline_dfdz(z, Spline, kx, ky, kz);
        X_vec[1] = Eval_Vert_Spline_dfdz(z, Spline, kx+1, ky, kz);
        X_vec[2] = Eval_Vert_Spline_dfdz(z, Spline, kx, ky+1, kz);
        X_vec[3] = Eval_Vert_Spline_dfdz(z, Spline, kx+1, ky+1, kz);
        
        X_vec[4] = Eval_Vert_Spline_ddfdxdz(z, Spline, kx, ky,kz)*dx_scalar;
        X_vec[5] = Eval_Vert_Spline_ddfdxdz(z, Spline, kx+1, ky,kz)*dx_scalar;
        X_vec[6] = Eval_Vert_Spline_ddfdxdz(z, Spline, kx, ky+1,kz)*dx_scalar;
        X_vec[7] = Eval_Vert_Spline_ddfdxdz(z, Spline, kx+1, ky+1,kz)*dx_scalar;
        
        X_vec[8] =  Eval_Vert_Spline_ddfdydz(z, Spline, kx, ky,kz)*dx_scalar;
        X_vec[9] =  Eval_Vert_Spline_ddfdydz(z, Spline, kx+1, ky,kz)*dx_scalar;
        X_vec[10] = Eval_Vert_Spline_ddfdydz(z, Spline, kx, ky+1,kz)*dx_scalar;
        X_vec[11] = Eval_Vert_Spline_ddfdydz(z, Spline, kx+1, ky+1,kz)*dx_scalar;
        
        X_vec[12]   = BiCubic_Deriv_dddfdxdydz(z, Spline, kx, ky, kz)*dx_scalar*dy_scalar;
        X_vec[13]   = BiCubic_Deriv_dddfdxdydz(z, Spline, kx+1, ky, kz)*dx_scalar*dy_scalar;
        X_vec[14]	= BiCubic_Deriv_dddfdxdydz(z, Spline, kx, ky+1, kz)*dx_scalar*dy_scalar;
        X_vec[15]   = BiCubic_Deriv_dddfdxdydz(z, Spline, kx+1, ky+1, kz)*dx_scalar*dy_scalar;
    }
    
    // Set the coefficients by matrix multiplication
    for(int j = 0; j < 16; j++){
        A_vec[j] = 0;
        for(int k = 0; k < 16; k++){
            A_vec[j] += BiCubic_ConversionMatrix[j][k]*X_vec[k];
        }
    }
    
    // Copy the a_vec elements into the coefficients array
    for(int k1=0; k1 < 4; k1++){
    for(int k2=0; k2 < 4; k2++){
        BicCoeff[k1][k2] = A_vec[k2*4 + k1];
    }}
    
    double x_scaled = (x - Spline.x_vals[kx])/(Spline.x_vals[kx+1] - Spline.x_vals[kx]);
	double y_scaled = (y - Spline.y_vals[ky])/(Spline.y_vals[ky+1] - Spline.y_vals[ky]);
    
	double result = 0;
    for(int k1 = 0; k1 < 4; k1++){
    for(int k2 = 0; k2 < 4; k2++){
        result+=1.0*BicCoeff[k1][k2]*pow(x_scaled,k1)*pow(y_scaled,k2);
    }}
    
    return result;
}


double Eval_Spline_ddf(double x, double y, double z, int index1, int index2, struct MultiDimSpline_3D & Spline){
    double A_vec[16];
	double X_vec[16];
    double BicCoeff[4][4];
    
    // Locate kx, ky, kz such that (x[kx] <= x' <= x[kx+1]), (y[ky] <= y' <= y[ky+1]), and (z[kz] <= z' <= z[kz+1])
    int kx = Find_Segment(x, Spline.x_vals, Spline.length_x, Spline.accel[0]);
    int ky = Find_Segment(y, Spline.y_vals, Spline.length_y, Spline.accel[1]);
    int kz = Find_Segment(z, Spline.z_vals, Spline.length_z, Spline.accel[2]);
    
    double dx_scalar = Spline.x_vals[kx+1] - Spline.x_vals[kx]; // Set dr scaling to produce unit square relations
    double dy_scalar = Spline.y_vals[ky+1] - Spline.y_vals[ky]; // Set dz scaling to produce unit square relations
    
    if(index1==0 && index2==0){
        // To calculate d^2 f/ dx^2, bicubic spline df/dx and take an x derivative (below)
        X_vec[0] = Eval_Vert_Spline_dfdx(z, Spline, kx, ky, kz);
        X_vec[1] = Eval_Vert_Spline_dfdx(z, Spline, kx+1, ky, kz);
        X_vec[2] = Eval_Vert_Spline_dfdx(z, Spline, kx, ky+1, kz);
        X_vec[3] = Eval_Vert_Spline_dfdx(z, Spline, kx+1, ky+1, kz);
        
        X_vec[4] = BiCubic_Deriv_ddfdxdx(z, Spline, kx, ky,kz)*dx_scalar;
        X_vec[5] = BiCubic_Deriv_ddfdxdx(z, Spline, kx+1, ky,kz)*dx_scalar;
        X_vec[6] = BiCubic_Deriv_ddfdxdx(z, Spline, kx, ky+1,kz)*dx_scalar;
        X_vec[7] = BiCubic_Deriv_ddfdxdx(z, Spline, kx+1, ky+1,kz)*dx_scalar;
        
        X_vec[8] =  BiCubic_Deriv_ddfdxdy(z, Spline, kx, ky,kz)*dx_scalar;
        X_vec[9] =  BiCubic_Deriv_ddfdxdy(z, Spline, kx+1, ky,kz)*dx_scalar;
        X_vec[10] = BiCubic_Deriv_ddfdxdy(z, Spline, kx, ky+1,kz)*dx_scalar;
        X_vec[11] = BiCubic_Deriv_ddfdxdy(z, Spline, kx+1, ky+1,kz)*dx_scalar;
        
        X_vec[12]   = BiCubic_Deriv_dddfdxdxdy(z, Spline, kx, ky, kz)*dx_scalar*dy_scalar;
        X_vec[13]   = BiCubic_Deriv_dddfdxdxdy(z, Spline, kx+1, ky, kz)*dx_scalar*dy_scalar;
        X_vec[14]	= BiCubic_Deriv_dddfdxdxdy(z, Spline, kx, ky+1, kz)*dx_scalar*dy_scalar;
        X_vec[15]   = BiCubic_Deriv_dddfdxdxdy(z, Spline, kx+1, ky+1, kz)*dx_scalar*dy_scalar;
        
    } else if(index1==1 && index2==1){
        // To calculate d^2 f/ dy^2, bicubic spline df/dy and take a y derivative (below)
        X_vec[0] = Eval_Vert_Spline_dfdy(z, Spline, kx, ky, kz);
        X_vec[1] = Eval_Vert_Spline_dfdy(z, Spline, kx+1, ky, kz);
        X_vec[2] = Eval_Vert_Spline_dfdy(z, Spline, kx, ky+1, kz);
        X_vec[3] = Eval_Vert_Spline_dfdy(z, Spline, kx+1, ky+1, kz);
        
        X_vec[4] = BiCubic_Deriv_ddfdxdy(z, Spline, kx, ky,kz)*dx_scalar;
        X_vec[5] = BiCubic_Deriv_ddfdxdy(z, Spline, kx+1, ky,kz)*dx_scalar;
        X_vec[6] = BiCubic_Deriv_ddfdxdy(z, Spline, kx, ky+1,kz)*dx_scalar;
        X_vec[7] = BiCubic_Deriv_ddfdxdy(z, Spline, kx+1, ky+1,kz)*dx_scalar;
        
        X_vec[8] =  BiCubic_Deriv_ddfdydy(z, Spline, kx, ky,kz)*dx_scalar;
        X_vec[9] =  BiCubic_Deriv_ddfdydy(z, Spline, kx+1, ky,kz)*dx_scalar;
        X_vec[10] = BiCubic_Deriv_ddfdydy(z, Spline, kx, ky+1,kz)*dx_scalar;
        X_vec[11] = BiCubic_Deriv_ddfdydy(z, Spline, kx+1, ky+1,kz)*dx_scalar;
        
        X_vec[12]   = BiCubic_Deriv_dddfdxdydy(z, Spline, kx, ky, kz)*dx_scalar*dy_scalar;
        X_vec[13]   = BiCubic_Deriv_dddfdxdydy(z, Spline, kx+1, ky, kz)*dx_scalar*dy_scalar;
        X_vec[14]	= BiCubic_Deriv_dddfdxdydy(z, Spline, kx, ky+1, kz)*dx_scalar*dy_scalar;
        X_vec[15]   = BiCubic_Deriv_dddfdxdydy(z, Spline, kx+1, ky+1, kz)*dx_scalar*dy_scalar;
        
    } else if(index1==2 && index2==2){
        // To calculate d^2 f/ dz^2, bicubic spline the second z derivative of f
        X_vec[0] = Eval_Vert_Spline_ddfdzdz(z, Spline, kx, ky, kz);
        X_vec[1] = Eval_Vert_Spline_ddfdzdz(z, Spline, kx+1, ky, kz);
        X_vec[2] = Eval_Vert_Spline_ddfdzdz(z, Spline, kx, ky+1, kz);
        X_vec[3] = Eval_Vert_Spline_ddfdzdz(z, Spline, kx+1, ky+1, kz);
        
        X_vec[4] = BiCubic_Deriv_dddfdxdzdz(z, Spline, kx, ky,kz)*dx_scalar;
        X_vec[5] = BiCubic_Deriv_dddfdxdzdz(z, Spline, kx+1, ky,kz)*dx_scalar;
        X_vec[6] = BiCubic_Deriv_dddfdxdzdz(z, Spline, kx, ky+1,kz)*dx_scalar;
        X_vec[7] = BiCubic_Deriv_dddfdxdzdz(z, Spline, kx+1, ky+1,kz)*dx_scalar;
        
        X_vec[8] =  BiCubic_Deriv_dddfdydzdz(z, Spline, kx, ky,kz)*dx_scalar;
        X_vec[9] =  BiCubic_Deriv_dddfdydzdz(z, Spline, kx+1, ky,kz)*dx_scalar;
        X_vec[10] = BiCubic_Deriv_dddfdydzdz(z, Spline, kx, ky+1,kz)*dx_scalar;
        X_vec[11] = BiCubic_Deriv_dddfdydzdz(z, Spline, kx+1, ky+1,kz)*dx_scalar;
        
        X_vec[12]   = BiCubic_Deriv_ddddfdxdydzdz(z, Spline, kx, ky, kz)*dx_scalar*dy_scalar;
        X_vec[13]   = BiCubic_Deriv_ddddfdxdydzdz(z, Spline, kx+1, ky, kz)*dx_scalar*dy_scalar;
        X_vec[14]	= BiCubic_Deriv_ddddfdxdydzdz(z, Spline, kx, ky+1, kz)*dx_scalar*dy_scalar;
        X_vec[15]   = BiCubic_Deriv_ddddfdxdydzdz(z, Spline, kx+1, ky+1, kz)*dx_scalar*dy_scalar;
        
    } else if((index1==0 && index2==1) || (index1==1 && index2==0)){
        // To calculate d^2 f/ dxdy, bicubic spline df/dx and take a y derivative (below)
        X_vec[0] = Eval_Vert_Spline_dfdx(z, Spline, kx, ky, kz);
        X_vec[1] = Eval_Vert_Spline_dfdx(z, Spline, kx+1, ky, kz);
        X_vec[2] = Eval_Vert_Spline_dfdx(z, Spline, kx, ky+1, kz);
        X_vec[3] = Eval_Vert_Spline_dfdx(z, Spline, kx+1, ky+1, kz);
        
        X_vec[4] = BiCubic_Deriv_ddfdxdx(z, Spline, kx, ky,kz)*dx_scalar;
        X_vec[5] = BiCubic_Deriv_ddfdxdx(z, Spline, kx+1, ky,kz)*dx_scalar;
        X_vec[6] = BiCubic_Deriv_ddfdxdx(z, Spline, kx, ky+1,kz)*dx_scalar;
        X_vec[7] = BiCubic_Deriv_ddfdxdx(z, Spline, kx+1, ky+1,kz)*dx_scalar;
        
        X_vec[8] =  BiCubic_Deriv_ddfdxdy(z, Spline, kx, ky,kz)*dx_scalar;
        X_vec[9] =  BiCubic_Deriv_ddfdxdy(z, Spline, kx+1, ky,kz)*dx_scalar;
        X_vec[10] = BiCubic_Deriv_ddfdxdy(z, Spline, kx, ky+1,kz)*dx_scalar;
        X_vec[11] = BiCubic_Deriv_ddfdxdy(z, Spline, kx+1, ky+1,kz)*dx_scalar;
        
        X_vec[12]   = BiCubic_Deriv_dddfdxdxdy(z, Spline, kx, ky, kz)*dx_scalar*dy_scalar;
        X_vec[13]   = BiCubic_Deriv_dddfdxdxdy(z, Spline, kx+1, ky, kz)*dx_scalar*dy_scalar;
        X_vec[14]	= BiCubic_Deriv_dddfdxdxdy(z, Spline, kx, ky+1, kz)*dx_scalar*dy_scalar;
        X_vec[15]   = BiCubic_Deriv_dddfdxdxdy(z, Spline, kx+1, ky+1, kz)*dx_scalar*dy_scalar;
        
    } else if((index1==0 && index2==2) || (index1==2 && index2==0)){
        // To calcualte d^2 f/dxdz, bicubic spline the z derivative of df/dx
        X_vec[0] = Eval_Vert_Spline_ddfdxdz(z, Spline, kx, ky, kz);
        X_vec[1] = Eval_Vert_Spline_ddfdxdz(z, Spline, kx+1, ky, kz);
        X_vec[2] = Eval_Vert_Spline_ddfdxdz(z, Spline, kx, ky+1, kz);
        X_vec[3] = Eval_Vert_Spline_ddfdxdz(z, Spline, kx+1, ky+1, kz);
        
        X_vec[4] = BiCubic_Deriv_dddfdxdxdz(z, Spline, kx, ky,kz)*dx_scalar;
        X_vec[5] = BiCubic_Deriv_dddfdxdxdz(z, Spline, kx+1, ky,kz)*dx_scalar;
        X_vec[6] = BiCubic_Deriv_dddfdxdxdz(z, Spline, kx, ky+1,kz)*dx_scalar;
        X_vec[7] = BiCubic_Deriv_dddfdxdxdz(z, Spline, kx+1, ky+1,kz)*dx_scalar;
        
        X_vec[8] =  BiCubic_Deriv_dddfdxdydz(z, Spline, kx, ky,kz)*dx_scalar;
        X_vec[9] =  BiCubic_Deriv_dddfdxdydz(z, Spline, kx+1, ky,kz)*dx_scalar;
        X_vec[10] = BiCubic_Deriv_dddfdxdydz(z, Spline, kx, ky+1,kz)*dx_scalar;
        X_vec[11] = BiCubic_Deriv_dddfdxdydz(z, Spline, kx+1, ky+1,kz)*dx_scalar;
        
        X_vec[12]   = BiCubic_Deriv_ddddfdxdxdydz(z, Spline, kx, ky, kz)*dx_scalar*dy_scalar;
        X_vec[13]   = BiCubic_Deriv_ddddfdxdxdydz(z, Spline, kx+1, ky, kz)*dx_scalar*dy_scalar;
        X_vec[14]	= BiCubic_Deriv_ddddfdxdxdydz(z, Spline, kx, ky+1, kz)*dx_scalar*dy_scalar;
        X_vec[15]   = BiCubic_Deriv_ddddfdxdxdydz(z, Spline, kx+1, ky+1, kz)*dx_scalar*dy_scalar;
        
    } else if((index1==1 && index2==2) || (index1==2 && index2==1)){
        // To calculate d^2 f/dydz, bicubic spline the z derivative of df/dy
        X_vec[0] = Eval_Vert_Spline_ddfdydz(z, Spline, kx, ky, kz);
        X_vec[1] = Eval_Vert_Spline_ddfdydz(z, Spline, kx+1, ky, kz);
        X_vec[2] = Eval_Vert_Spline_ddfdydz(z, Spline, kx, ky+1, kz);
        X_vec[3] = Eval_Vert_Spline_ddfdydz(z, Spline, kx+1, ky+1, kz);
        
        X_vec[4] = BiCubic_Deriv_dddfdxdydz(z, Spline, kx, ky,kz)*dx_scalar;
        X_vec[5] = BiCubic_Deriv_dddfdxdydz(z, Spline, kx+1, ky,kz)*dx_scalar;
        X_vec[6] = BiCubic_Deriv_dddfdxdydz(z, Spline, kx, ky+1,kz)*dx_scalar;
        X_vec[7] = BiCubic_Deriv_dddfdxdydz(z, Spline, kx+1, ky+1,kz)*dx_scalar;
        
        X_vec[8] =  BiCubic_Deriv_dddfdydydz(z, Spline, kx, ky,kz)*dx_scalar;
        X_vec[9] =  BiCubic_Deriv_dddfdydydz(z, Spline, kx+1, ky,kz)*dx_scalar;
        X_vec[10] = BiCubic_Deriv_dddfdydydz(z, Spline, kx, ky+1,kz)*dx_scalar;
        X_vec[11] = BiCubic_Deriv_dddfdydydz(z, Spline, kx+1, ky+1,kz)*dx_scalar;
        
        X_vec[12]   = BiCubic_Deriv_ddddfdxdydydz(z, Spline, kx, ky, kz)*dx_scalar*dy_scalar;
        X_vec[13]   = BiCubic_Deriv_ddddfdxdydydz(z, Spline, kx+1, ky, kz)*dx_scalar*dy_scalar;
        X_vec[14]	= BiCubic_Deriv_ddddfdxdydydz(z, Spline, kx, ky+1, kz)*dx_scalar*dy_scalar;
        X_vec[15]   = BiCubic_Deriv_ddddfdxdydydz(z, Spline, kx+1, ky+1, kz)*dx_scalar*dy_scalar;
    }
    
    // Set the coefficients by matrix multiplication
    for(int j = 0; j < 16; j++){
        A_vec[j] = 0;
        for(int k = 0; k < 16; k++){
            A_vec[j] += BiCubic_ConversionMatrix[j][k]*X_vec[k];
        }
    }
    
    // Copy the a_vec elements into the coefficients array
    for(int k1=0; k1 < 4; k1++){
    for(int k2=0; k2 < 4; k2++){
        BicCoeff[k1][k2] = A_vec[k2*4 + k1];
    }}
    
    double x_scaled = (x - Spline.x_vals[kx])/(Spline.x_vals[kx+1] - Spline.x_vals[kx]);
	double y_scaled = (y - Spline.y_vals[ky])/(Spline.y_vals[ky+1] - Spline.y_vals[ky]);
    
	double result = 0;
    if(index1==0 && index2==0){
        // d^2 f / dx^2 requires taking an x derivative df/dx
        for(int k1 = 1; k1 < 4; k1++){
        for(int k2 = 0; k2 < 4; k2++){
            result+=1.0*k1*BicCoeff[k1][k2]*pow(x_scaled,k1-1)*pow(y_scaled,k2)/(Spline.x_vals[kx+1] - Spline.x_vals[kx]);
        }}
    } else if((index1==1 && index2==1) || (index1==0 && index2==1) || (index1==1 && index2==0)){
        // d^2 f / dy^2 and d^2 f / dxdy require taking a y derivative of df/dy and df/dx respectively
        for(int k1 = 0; k1 < 4; k1++){
        for(int k2 = 1; k2 < 4; k2++){
            result+=1.0*k2*BicCoeff[k1][k2]*pow(x_scaled,k1)*pow(y_scaled,k2-1)/(Spline.y_vals[ky+1] - Spline.y_vals[ky]);
        }}
    } else {
        // All other results are obtained immediately without differentiating the bicubic spline
        for(int k1 = 0; k1 < 4; k1++){
        for(int k2 = 0; k2 < 4; k2++){
            result+=1.0*BicCoeff[k1][k2]*pow(x_scaled,k1)*pow(y_scaled,k2);
        }}
    }
	
    return result;
}

void Eval_Spline_AllOrder1(double x, double y, double z, struct MultiDimSpline_3D & Spline, double & f, double & dfdx, double & dfdy, double & dfdz){
    double A_vec[16];
	double X_vec[16];
    
    // Locate kx, ky, kz such that (x[kx] <= x' <= x[kx+1]), (y[ky] <= y' <= y[ky+1]), and (z[kz] <= z' <= z[kz+1])
    double x_eval = min(x, x_max);  x_eval = max(x_eval, x_min);    // Check that x_min <= x_eval <= x_max
    double y_eval = min(y, y_max);  y_eval = max(y_eval, y_min);    // Check that y_min <= y_eval <= y_max
    double z_eval = min(z, z_max);  z_eval = max(z_eval, z_min);    // Check that z_min <= z_eval <= z_max
    
    int kx = Find_Segment(x_eval, Spline.x_vals, Spline.length_x, Spline.accel[0]);
    int ky = Find_Segment(y_eval, Spline.y_vals, Spline.length_y, Spline.accel[1]);
    int kz = Find_Segment(z_eval, Spline.z_vals, Spline.length_z, Spline.accel[2]);
    
    // Compute the dx and dy scaling so that the grid square can be treated as a unit grid then the scaled versions of x and y within the unit grid
    double dx_scalar = Spline.x_vals[kx+1] - Spline.x_vals[kx];
    double dy_scalar = Spline.y_vals[ky+1] - Spline.y_vals[ky];
    
    double x_scaled = (x_eval - Spline.x_vals[kx])/(Spline.x_vals[kx+1] - Spline.x_vals[kx]);
	double y_scaled = (y_eval - Spline.y_vals[ky])/(Spline.y_vals[ky+1] - Spline.y_vals[ky]);
    
    // Precompute finite difference df/dx, df/dy, and d^2f/dxdy which are used multiple times
    double Finite_dfdx_00 = BiCubic_Deriv_dfdx(z_eval, Spline, kx, ky, kz);
    double Finite_dfdx_10 = BiCubic_Deriv_dfdx(z_eval, Spline, kx+1, ky, kz);
    double Finite_dfdx_01 = BiCubic_Deriv_dfdx(z_eval, Spline, kx, ky+1, kz);
    double Finite_dfdx_11 = BiCubic_Deriv_dfdx(z_eval, Spline, kx+1, ky+1, kz);
    
    double Finite_dfdy_00 = BiCubic_Deriv_dfdy(z_eval, Spline, kx, ky, kz);
    double Finite_dfdy_10 = BiCubic_Deriv_dfdy(z_eval, Spline, kx+1, ky, kz);
    double Finite_dfdy_01 = BiCubic_Deriv_dfdy(z_eval, Spline, kx, ky+1, kz);
    double Finite_dfdy_11 = BiCubic_Deriv_dfdy(z_eval, Spline, kx+1, ky+1, kz);
    
    double Finite_ddfdxdy_00 = BiCubic_Deriv_ddfdxdy(z_eval, Spline, kx, ky, kz);
    double Finite_ddfdxdy_10 = BiCubic_Deriv_ddfdxdy(z_eval, Spline, kx+1, ky, kz);
    double Finite_ddfdxdy_01 = BiCubic_Deriv_ddfdxdy(z_eval, Spline, kx, ky+1, kz);
    double Finite_ddfdxdy_11 = BiCubic_Deriv_ddfdxdy(z_eval, Spline, kx+1, ky+1, kz);
    
    // Bicubic Spline f and evaluate it
        X_vec[0] = Eval_Vert_Spline_f(z_eval, Spline, kx, ky, kz);
        X_vec[1] = Eval_Vert_Spline_f(z_eval, Spline, kx+1, ky, kz);
        X_vec[2] = Eval_Vert_Spline_f(z_eval, Spline, kx, ky+1, kz);
        X_vec[3] = Eval_Vert_Spline_f(z_eval, Spline, kx+1, ky+1, kz);
    
        X_vec[4] = Finite_dfdx_00*dx_scalar;
        X_vec[5] = Finite_dfdx_10*dx_scalar;
        X_vec[6] = Finite_dfdx_01*dx_scalar;
        X_vec[7] = Finite_dfdx_11*dx_scalar;
    
        X_vec[8] =  Finite_dfdy_00*dy_scalar;
        X_vec[9] =  Finite_dfdy_10*dy_scalar;
        X_vec[10] = Finite_dfdy_01*dy_scalar;
        X_vec[11] = Finite_dfdy_11*dy_scalar;
    
        X_vec[12]   = Finite_ddfdxdy_00*dx_scalar*dy_scalar;
        X_vec[13]   = Finite_ddfdxdy_10*dx_scalar*dy_scalar;
        X_vec[14]	= Finite_ddfdxdy_01*dx_scalar*dy_scalar;
        X_vec[15]   = Finite_ddfdxdy_11*dx_scalar*dy_scalar;
    
        // Set A[j] by matrix multiplication
        for(int j = 0; j < 16; j++){
            A_vec[j] = 0;
            for(int k = 0; k < 16; k++){
                A_vec[j] += BiCubic_ConversionMatrix[j][k]*X_vec[k];
            }
        }
    
        // Evaluate f(x,y,z) = Sum{k1,k2 = 0...3} c_{k1,k2}(z) x^k1 y^k2
        // where c_{k1,k2}(z) = a[k1 + 4*k2](z)
        f = 0;
        for(int k1 = 0; k1 < 4; k1++){
        for(int k2 = 0; k2 < 4; k2++){
            f+=1.0*A_vec[k1 + 4*k2]*pow(x_scaled,k1)*pow(y_scaled,k2);
        }}
	
    // Bicubic Spline df/dx
        X_vec[0] = Finite_dfdx_00;
        X_vec[1] = Finite_dfdx_10;
        X_vec[2] = Finite_dfdx_01;
        X_vec[3] = Finite_dfdx_11;
    
        X_vec[4] = BiCubic_Deriv_ddfdxdx(z_eval, Spline, kx, ky,kz)*dx_scalar;
        X_vec[5] = BiCubic_Deriv_ddfdxdx(z_eval, Spline, kx+1, ky,kz)*dx_scalar;
        X_vec[6] = BiCubic_Deriv_ddfdxdx(z_eval, Spline, kx, ky+1,kz)*dx_scalar;
        X_vec[7] = BiCubic_Deriv_ddfdxdx(z_eval, Spline, kx+1, ky+1,kz)*dx_scalar;
    
        X_vec[8]   =  Finite_ddfdxdy_00*dy_scalar;
        X_vec[9]   =  Finite_ddfdxdy_10*dy_scalar;
        X_vec[10]	= Finite_ddfdxdy_01*dy_scalar;
        X_vec[11]   = Finite_ddfdxdy_11*dy_scalar;
    
        X_vec[12]   = BiCubic_Deriv_dddfdxdxdy(z_eval, Spline, kx, ky, kz)*dx_scalar*dy_scalar;
        X_vec[13]   = BiCubic_Deriv_dddfdxdxdy(z_eval, Spline, kx+1, ky, kz)*dx_scalar*dy_scalar;
        X_vec[14]	= BiCubic_Deriv_dddfdxdxdy(z_eval, Spline, kx, ky+1, kz)*dx_scalar*dy_scalar;
        X_vec[15]   = BiCubic_Deriv_dddfdxdxdy(z_eval, Spline, kx+1, ky+1, kz)*dx_scalar*dy_scalar;
    
        // Set A[j] by matrix multiplication
        for(int j = 0; j < 16; j++){
            A_vec[j] = 0;
            for(int k = 0; k < 16; k++){
                A_vec[j] += BiCubic_ConversionMatrix[j][k]*X_vec[k];
            }
        }
    
        // Evaluate approriate derivatives of df/dx(x,y) = Sum{k1,k2 = 0...3} c_{k1,k2} x^k1 y^k2
        // where c_{k1,k2} = a[k1 + 4*k2]
        dfdx = 0;
        for(int k1 = 0; k1 < 4; k1++){
        for(int k2 = 0; k2 < 4; k2++){
            dfdx+=1.0*A_vec[k1 + 4*k2]*pow(x_scaled,k1)*pow(y_scaled,k2);
        }}
    
    // Bicubic Spline df/dy
        X_vec[0] = Finite_dfdy_00;
        X_vec[1] = Finite_dfdy_10;
        X_vec[2] = Finite_dfdy_01;
        X_vec[3] = Finite_dfdy_11;
    
        X_vec[4]   = Finite_ddfdxdy_00*dx_scalar;
        X_vec[5]   = Finite_ddfdxdy_10*dx_scalar;
        X_vec[6]   = Finite_ddfdxdy_01*dx_scalar;
        X_vec[7]   = Finite_ddfdxdy_11*dx_scalar;
    
        X_vec[8] =  BiCubic_Deriv_ddfdydy(z_eval, Spline, kx, ky,kz)*dy_scalar;
        X_vec[9] =  BiCubic_Deriv_ddfdydy(z_eval, Spline, kx+1, ky,kz)*dy_scalar;
        X_vec[10] = BiCubic_Deriv_ddfdydy(z_eval, Spline, kx, ky+1,kz)*dy_scalar;
        X_vec[11] = BiCubic_Deriv_ddfdydy(z_eval, Spline, kx+1, ky+1,kz)*dy_scalar;
    
        X_vec[12]   = BiCubic_Deriv_dddfdxdydy(z_eval, Spline, kx, ky, kz)*dx_scalar*dy_scalar;
        X_vec[13]   = BiCubic_Deriv_dddfdxdydy(z_eval, Spline, kx+1, ky, kz)*dx_scalar*dy_scalar;
        X_vec[14]	= BiCubic_Deriv_dddfdxdydy(z_eval, Spline, kx, ky+1, kz)*dx_scalar*dy_scalar;
        X_vec[15]   = BiCubic_Deriv_dddfdxdydy(z_eval, Spline, kx+1, ky+1, kz)*dx_scalar*dy_scalar;
    
        // Set A[j] by matrix multiplication
        for(int j = 0; j < 16; j++){
            A_vec[j] = 0;
            for(int k = 0; k < 16; k++){
                A_vec[j] += BiCubic_ConversionMatrix[j][k]*X_vec[k];
            }
        }
    
        // Evaluate approriate derivatives of df/dx(x,y) = Sum{k1,k2 = 0...3} c_{k1,k2} x^k1 y^k2
        // where c_{k1,k2} = a[k1 + 4*k2]
        dfdy = 0;
        for(int k1 = 0; k1 < 4; k1++){
        for(int k2 = 0; k2 < 4; k2++){
            dfdy+=1.0*A_vec[k1 + 4*k2]*pow(x_scaled,k1)*pow(y_scaled,k2);
        }}
    
    // Bicubic Spline df/dz and evaluate df/dz, d^2f/dxdz, and d^2f/dydz from the Bicubic Spline
        X_vec[0] = Eval_Vert_Spline_dfdz(z_eval, Spline, kx, ky, kz);
        X_vec[1] = Eval_Vert_Spline_dfdz(z_eval, Spline, kx+1, ky, kz);
        X_vec[2] = Eval_Vert_Spline_dfdz(z_eval, Spline, kx, ky+1, kz);
        X_vec[3] = Eval_Vert_Spline_dfdz(z_eval, Spline, kx+1, ky+1, kz);
    
        X_vec[4] = Eval_Vert_Spline_ddfdxdz(z_eval, Spline, kx, ky,kz)*dx_scalar;
        X_vec[5] = Eval_Vert_Spline_ddfdxdz(z_eval, Spline, kx+1, ky,kz)*dx_scalar;
        X_vec[6] = Eval_Vert_Spline_ddfdxdz(z_eval, Spline, kx, ky+1,kz)*dx_scalar;
        X_vec[7] = Eval_Vert_Spline_ddfdxdz(z_eval, Spline, kx+1, ky+1,kz)*dx_scalar;
    
        X_vec[8] =  Eval_Vert_Spline_ddfdydz(z_eval, Spline, kx, ky,kz)*dy_scalar;
        X_vec[9] =  Eval_Vert_Spline_ddfdydz(z_eval, Spline, kx+1, ky,kz)*dy_scalar;
        X_vec[10] = Eval_Vert_Spline_ddfdydz(z_eval, Spline, kx, ky+1,kz)*dy_scalar;
        X_vec[11] = Eval_Vert_Spline_ddfdydz(z_eval, Spline, kx+1, ky+1,kz)*dy_scalar;
    
        X_vec[12]   = BiCubic_Deriv_dddfdxdydz(z_eval, Spline, kx, ky, kz)*dx_scalar*dy_scalar;
        X_vec[13]   = BiCubic_Deriv_dddfdxdydz(z_eval, Spline, kx+1, ky, kz)*dx_scalar*dy_scalar;
        X_vec[14]	= BiCubic_Deriv_dddfdxdydz(z_eval, Spline, kx, ky+1, kz)*dx_scalar*dy_scalar;
        X_vec[15]   = BiCubic_Deriv_dddfdxdydz(z_eval, Spline, kx+1, ky+1, kz)*dx_scalar*dy_scalar;
    
        // Set A[j] by matrix multiplication
        for(int j = 0; j < 16; j++){
            A_vec[j] = 0;
            for(int k = 0; k < 16; k++){
                A_vec[j] += BiCubic_ConversionMatrix[j][k]*X_vec[k];
            }
        }
    
        // Evaluate approriate derivatives of df/dx(x,y) = Sum{k1,k2 = 0...3} c_{k1,k2} x^k1 y^k2
        // where c_{k1,k2} = a[k1 + 4*k2]
        dfdz = 0;
        for(int k1 = 0; k1 < 4; k1++){
        for(int k2 = 0; k2 < 4; k2++){
            dfdz+=1.0*A_vec[k1 + 4*k2]*pow(x_scaled,k1)*pow(y_scaled,k2);
        }}
}

void Eval_Spline_AllOrder2(double x, double y, double z, struct MultiDimSpline_3D & Spline, double & f, double & dfdx, double & dfdy, double & dfdz, double & ddfdxdx, double & ddfdydy, double & ddfdzdz, double & ddfdxdy, double & ddfdxdz, double & ddfdydz){
    double A_vec[16];
	double X_vec[16];
    
    // Locate kx, ky, kz such that (x[kx] <= x' <= x[kx+1]), (y[ky] <= y' <= y[ky+1]), and (z[kz] <= z' <= z[kz+1])
    double x_eval = min(x, x_max);  x_eval = max(x_eval, x_min);    // Check that x_min <= x_eval <= x_max
    double y_eval = min(y, y_max);  y_eval = max(y_eval, y_min);    // Check that y_min <= y_eval <= y_max
    double z_eval = min(z, z_max);  z_eval = max(z_eval, z_min);    // Check that z_min <= z_eval <= z_max
    
    int kx = Find_Segment(x_eval, Spline.x_vals, Spline.length_x, Spline.accel[0]);
    int ky = Find_Segment(y_eval, Spline.y_vals, Spline.length_y, Spline.accel[1]);
    int kz = Find_Segment(z_eval, Spline.z_vals, Spline.length_z, Spline.accel[2]);
    
    // Compute the dx and dy scaling so that the grid square can be treated as a unit grid then the scaled versions of x and y within the unit grid
    double dx_scalar = Spline.x_vals[kx+1] - Spline.x_vals[kx];
    double dy_scalar = Spline.y_vals[ky+1] - Spline.y_vals[ky];
    
    double x_scaled = (x_eval - Spline.x_vals[kx])/dx_scalar;
	double y_scaled = (y_eval - Spline.y_vals[ky])/dy_scalar;
    
    // Precompute finite difference df/dx, df/dy, and d^2f/dxdy which are used multiple times
    double Finite_dfdx_00 = BiCubic_Deriv_dfdx(z_eval, Spline, kx, ky, kz);
    double Finite_dfdx_10 = BiCubic_Deriv_dfdx(z_eval, Spline, kx+1, ky, kz);
    double Finite_dfdx_01 = BiCubic_Deriv_dfdx(z_eval, Spline, kx, ky+1, kz);
    double Finite_dfdx_11 = BiCubic_Deriv_dfdx(z_eval, Spline, kx+1, ky+1, kz);

    double Finite_dfdy_00 = BiCubic_Deriv_dfdy(z_eval, Spline, kx, ky, kz);
    double Finite_dfdy_10 = BiCubic_Deriv_dfdy(z_eval, Spline, kx+1, ky, kz);
    double Finite_dfdy_01 = BiCubic_Deriv_dfdy(z_eval, Spline, kx, ky+1, kz);
    double Finite_dfdy_11 = BiCubic_Deriv_dfdy(z_eval, Spline, kx+1, ky+1, kz);

    double Finite_ddfdxdy_00 = BiCubic_Deriv_ddfdxdy(z_eval, Spline, kx, ky, kz);
    double Finite_ddfdxdy_10 = BiCubic_Deriv_ddfdxdy(z_eval, Spline, kx+1, ky, kz);
    double Finite_ddfdxdy_01 = BiCubic_Deriv_ddfdxdy(z_eval, Spline, kx, ky+1, kz);
    double Finite_ddfdxdy_11 = BiCubic_Deriv_ddfdxdy(z_eval, Spline, kx+1, ky+1, kz);

    // Bicubic Spline f and evaluate it
        X_vec[0] = Eval_Vert_Spline_f(z_eval, Spline, kx, ky, kz);
        X_vec[1] = Eval_Vert_Spline_f(z_eval, Spline, kx+1, ky, kz);
        X_vec[2] = Eval_Vert_Spline_f(z_eval, Spline, kx, ky+1, kz);
        X_vec[3] = Eval_Vert_Spline_f(z_eval, Spline, kx+1, ky+1, kz);
    
        X_vec[4] = Finite_dfdx_00*dx_scalar;
        X_vec[5] = Finite_dfdx_10*dx_scalar;
        X_vec[6] = Finite_dfdx_01*dx_scalar;
        X_vec[7] = Finite_dfdx_11*dx_scalar;
    
        X_vec[8] =  Finite_dfdy_00*dy_scalar;
        X_vec[9] =  Finite_dfdy_10*dy_scalar;
        X_vec[10] = Finite_dfdy_01*dy_scalar;
        X_vec[11] = Finite_dfdy_11*dy_scalar;
    
        X_vec[12]   = Finite_ddfdxdy_00*dx_scalar*dy_scalar;
        X_vec[13]   = Finite_ddfdxdy_10*dx_scalar*dy_scalar;
        X_vec[14]	= Finite_ddfdxdy_01*dx_scalar*dy_scalar;
        X_vec[15]   = Finite_ddfdxdy_11*dx_scalar*dy_scalar;
    
        // Set A[j] by matrix multiplication
        for(int j = 0; j < 16; j++){
            A_vec[j] = 0;
            for(int k = 0; k < 16; k++){
                A_vec[j] += BiCubic_ConversionMatrix[j][k]*X_vec[k];
            }
        }
    
        // Evaluate f(x,y,z) = Sum{k1,k2 = 0...3} c_{k1,k2}(z) x^k1 y^k2
        // where c_{k1,k2}(z) = a[k1 + 4*k2](z)
        f = 0;
        for(int k1 = 0; k1 < 4; k1++){
        for(int k2 = 0; k2 < 4; k2++){
            f+=1.0*A_vec[k1 + 4*k2]*pow(x_scaled,k1)*pow(y_scaled,k2);
        }}
	
    // Bicubic Spline df/dx and evaluate df/dx, d^2f/dx^2, and d^2f/dxdy
        X_vec[0] = Finite_dfdx_00;
        X_vec[1] = Finite_dfdx_10;
        X_vec[2] = Finite_dfdx_01;
        X_vec[3] = Finite_dfdx_11;
    
        X_vec[4] = BiCubic_Deriv_ddfdxdx(z_eval, Spline, kx, ky,kz)*dx_scalar;
        X_vec[5] = BiCubic_Deriv_ddfdxdx(z_eval, Spline, kx+1, ky,kz)*dx_scalar;
        X_vec[6] = BiCubic_Deriv_ddfdxdx(z_eval, Spline, kx, ky+1,kz)*dx_scalar;
        X_vec[7] = BiCubic_Deriv_ddfdxdx(z_eval, Spline, kx+1, ky+1,kz)*dx_scalar;
    
        X_vec[8]   =  Finite_ddfdxdy_00*dy_scalar;
        X_vec[9]   =  Finite_ddfdxdy_10*dy_scalar;
        X_vec[10]	= Finite_ddfdxdy_01*dy_scalar;
        X_vec[11]   = Finite_ddfdxdy_11*dy_scalar;
    
        X_vec[12]   = BiCubic_Deriv_dddfdxdxdy(z_eval, Spline, kx, ky, kz)*dx_scalar*dy_scalar;
        X_vec[13]   = BiCubic_Deriv_dddfdxdxdy(z_eval, Spline, kx+1, ky, kz)*dx_scalar*dy_scalar;
        X_vec[14]	= BiCubic_Deriv_dddfdxdxdy(z_eval, Spline, kx, ky+1, kz)*dx_scalar*dy_scalar;
        X_vec[15]   = BiCubic_Deriv_dddfdxdxdy(z_eval, Spline, kx+1, ky+1, kz)*dx_scalar*dy_scalar;

    
        // Set A[j] by matrix multiplication
        for(int j = 0; j < 16; j++){
            A_vec[j] = 0;
            for(int k = 0; k < 16; k++){
                A_vec[j] += BiCubic_ConversionMatrix[j][k]*X_vec[k];
            }
        }
    
        // Evaluate approriate derivatives of df/dx(x,y) = Sum{k1,k2 = 0...3} c_{k1,k2} x^k1 y^k2
        // where c_{k1,k2} = a[k1 + 4*k2]
        dfdx = 0;
        for(int k1 = 0; k1 < 4; k1++){
        for(int k2 = 0; k2 < 4; k2++){
            dfdx+=1.0*A_vec[k1 + 4*k2]*pow(x_scaled,k1)*pow(y_scaled,k2);
        }}
    
        ddfdxdx = 0;
        for(int k1 = 1; k1 < 4; k1++){
        for(int k2 = 0; k2 < 4; k2++){
            ddfdxdx+=1.0*k1*A_vec[k1 + 4*k2]*pow(x_scaled,k1-1)*pow(y_scaled,k2)/dx_scalar;
        }}
    
        ddfdxdy = 0;
        for(int k1 = 0; k1 < 4; k1++){
        for(int k2 = 1; k2 < 4; k2++){
            ddfdxdy+=1.0*k2*A_vec[k1 + 4*k2]*pow(x_scaled,k1)*pow(y_scaled,k2-1)/dy_scalar;
        }}
    
    // Bicubic Spline df/dy and evaluate df/dy and d^2f/dy^2 from the Bicubic Spline
        X_vec[0] = Finite_dfdy_00;
        X_vec[1] = Finite_dfdy_10;
        X_vec[2] = Finite_dfdy_01;
        X_vec[3] = Finite_dfdy_11;
    
        X_vec[4]   = Finite_ddfdxdy_00*dx_scalar;
        X_vec[5]   = Finite_ddfdxdy_10*dx_scalar;
        X_vec[6]   = Finite_ddfdxdy_01*dx_scalar;
        X_vec[7]   = Finite_ddfdxdy_11*dx_scalar;
    
        X_vec[8] =  BiCubic_Deriv_ddfdydy(z_eval, Spline, kx, ky,kz)*dy_scalar;
        X_vec[9] =  BiCubic_Deriv_ddfdydy(z_eval, Spline, kx+1, ky,kz)*dy_scalar;
        X_vec[10] = BiCubic_Deriv_ddfdydy(z_eval, Spline, kx, ky+1,kz)*dy_scalar;
        X_vec[11] = BiCubic_Deriv_ddfdydy(z_eval, Spline, kx+1, ky+1,kz)*dy_scalar;
    
        X_vec[12]   = BiCubic_Deriv_dddfdxdydy(z_eval, Spline, kx, ky, kz)*dx_scalar*dy_scalar;
        X_vec[13]   = BiCubic_Deriv_dddfdxdydy(z_eval, Spline, kx+1, ky, kz)*dx_scalar*dy_scalar;
        X_vec[14]	= BiCubic_Deriv_dddfdxdydy(z_eval, Spline, kx, ky+1, kz)*dx_scalar*dy_scalar;
        X_vec[15]   = BiCubic_Deriv_dddfdxdydy(z_eval, Spline, kx+1, ky+1, kz)*dx_scalar*dy_scalar;
    
        // Set A[j] by matrix multiplication
        for(int j = 0; j < 16; j++){
            A_vec[j] = 0;
            for(int k = 0; k < 16; k++){
                A_vec[j] += BiCubic_ConversionMatrix[j][k]*X_vec[k];
            }
        }
    
        // Evaluate approriate derivatives of df/dx(x,y) = Sum{k1,k2 = 0...3} c_{k1,k2} x^k1 y^k2
        // where c_{k1,k2} = a[k1 + 4*k2]
        dfdy = 0;
        for(int k1 = 0; k1 < 4; k1++){
        for(int k2 = 0; k2 < 4; k2++){
            dfdy+=1.0*A_vec[k1 + 4*k2]*pow(x_scaled,k1)*pow(y_scaled,k2);
        }}
    
        ddfdydy = 0;
        for(int k1 = 0; k1 < 4; k1++){
        for(int k2 = 1; k2 < 4; k2++){
            ddfdydy+=1.0*k2*A_vec[k1 + 4*k2]*pow(x_scaled,k1)*pow(y_scaled,k2-1)/dy_scalar;
        }}
    
    
    // Bicubic Spline df/dz and evaluate df/dz, d^2f/dxdz, and d^2f/dydz from the Bicubic Spline
        X_vec[0] = Eval_Vert_Spline_dfdz(z_eval, Spline, kx, ky, kz);
        X_vec[1] = Eval_Vert_Spline_dfdz(z_eval, Spline, kx+1, ky, kz);
        X_vec[2] = Eval_Vert_Spline_dfdz(z_eval, Spline, kx, ky+1, kz);
        X_vec[3] = Eval_Vert_Spline_dfdz(z_eval, Spline, kx+1, ky+1, kz);
    
        X_vec[4] = Eval_Vert_Spline_ddfdxdz(z_eval, Spline, kx, ky,kz)*dx_scalar;
        X_vec[5] = Eval_Vert_Spline_ddfdxdz(z_eval, Spline, kx+1, ky,kz)*dx_scalar;
        X_vec[6] = Eval_Vert_Spline_ddfdxdz(z_eval, Spline, kx, ky+1,kz)*dx_scalar;
        X_vec[7] = Eval_Vert_Spline_ddfdxdz(z_eval, Spline, kx+1, ky+1,kz)*dx_scalar;
    
        X_vec[8] =  Eval_Vert_Spline_ddfdydz(z_eval, Spline, kx, ky,kz)*dy_scalar;
        X_vec[9] =  Eval_Vert_Spline_ddfdydz(z_eval, Spline, kx+1, ky,kz)*dy_scalar;
        X_vec[10] = Eval_Vert_Spline_ddfdydz(z_eval, Spline, kx, ky+1,kz)*dy_scalar;
        X_vec[11] = Eval_Vert_Spline_ddfdydz(z_eval, Spline, kx+1, ky+1,kz)*dy_scalar;
    
        X_vec[12]   = BiCubic_Deriv_dddfdxdydz(z_eval, Spline, kx, ky, kz)*dx_scalar*dy_scalar;
        X_vec[13]   = BiCubic_Deriv_dddfdxdydz(z_eval, Spline, kx+1, ky, kz)*dx_scalar*dy_scalar;
        X_vec[14]	= BiCubic_Deriv_dddfdxdydz(z_eval, Spline, kx, ky+1, kz)*dx_scalar*dy_scalar;
        X_vec[15]   = BiCubic_Deriv_dddfdxdydz(z_eval, Spline, kx+1, ky+1, kz)*dx_scalar*dy_scalar;
    
        // Set A[j] by matrix multiplication
        for(int j = 0; j < 16; j++){
            A_vec[j] = 0;
            for(int k = 0; k < 16; k++){
                A_vec[j] += BiCubic_ConversionMatrix[j][k]*X_vec[k];
            }
        }
    
        // Evaluate approriate derivatives of df/dx(x,y) = Sum{k1,k2 = 0...3} c_{k1,k2} x^k1 y^k2
        // where c_{k1,k2} = a[k1 + 4*k2]
        dfdz = 0;
        for(int k1 = 0; k1 < 4; k1++){
        for(int k2 = 0; k2 < 4; k2++){
            dfdz+=1.0*A_vec[k1 + 4*k2]*pow(x_scaled,k1)*pow(y_scaled,k2);
        }}
    
        ddfdxdz = 0;
        for(int k1 = 1; k1 < 4; k1++){
        for(int k2 = 0; k2 < 4; k2++){
            ddfdxdz+=1.0*k1*A_vec[k1 + 4*k2]*pow(x_scaled,k1-1)*pow(y_scaled,k2)/dx_scalar;
        }}
    
        ddfdydz = 0;
        for(int k1 = 0; k1 < 4; k1++){
        for(int k2 = 1; k2 < 4; k2++){
            ddfdydz+=1.0*k2*A_vec[k1 + 4*k2]*pow(x_scaled,k1)*pow(y_scaled,k2-1)/dy_scalar;
        }}   
    
    // Bicubic Spline d^2f/dz^2
        X_vec[0] = Eval_Vert_Spline_ddfdzdz(z_eval, Spline, kx, ky, kz);
        X_vec[1] = Eval_Vert_Spline_ddfdzdz(z_eval, Spline, kx+1, ky, kz);
        X_vec[2] = Eval_Vert_Spline_ddfdzdz(z_eval, Spline, kx, ky+1, kz);
        X_vec[3] = Eval_Vert_Spline_ddfdzdz(z_eval, Spline, kx+1, ky+1, kz);
    
        X_vec[4] = BiCubic_Deriv_dddfdxdzdz(z_eval, Spline, kx, ky,kz)*dx_scalar;
        X_vec[5] = BiCubic_Deriv_dddfdxdzdz(z_eval, Spline, kx+1, ky,kz)*dx_scalar;
        X_vec[6] = BiCubic_Deriv_dddfdxdzdz(z_eval, Spline, kx, ky+1,kz)*dx_scalar;
        X_vec[7] = BiCubic_Deriv_dddfdxdzdz(z_eval, Spline, kx+1, ky+1,kz)*dx_scalar;
    
        X_vec[8] =  BiCubic_Deriv_dddfdydzdz(z_eval, Spline, kx, ky,kz)*dx_scalar;
        X_vec[9] =  BiCubic_Deriv_dddfdydzdz(z_eval, Spline, kx+1, ky,kz)*dx_scalar;
        X_vec[10] = BiCubic_Deriv_dddfdydzdz(z_eval, Spline, kx, ky+1,kz)*dx_scalar;
        X_vec[11] = BiCubic_Deriv_dddfdydzdz(z_eval, Spline, kx+1, ky+1,kz)*dx_scalar;
    
        X_vec[12]   = BiCubic_Deriv_ddddfdxdydzdz(z_eval, Spline, kx, ky, kz)*dx_scalar*dy_scalar;
        X_vec[13]   = BiCubic_Deriv_ddddfdxdydzdz(z_eval, Spline, kx+1, ky, kz)*dx_scalar*dy_scalar;
        X_vec[14]	= BiCubic_Deriv_ddddfdxdydzdz(z_eval, Spline, kx, ky+1, kz)*dx_scalar*dy_scalar;
        X_vec[15]   = BiCubic_Deriv_ddddfdxdydzdz(z_eval, Spline, kx+1, ky+1, kz)*dx_scalar*dy_scalar;
    
        // Set A[j] by matrix multiplication
        for(int j = 0; j < 16; j++){
            A_vec[j] = 0;
            for(int k = 0; k < 16; k++){
                A_vec[j] += BiCubic_ConversionMatrix[j][k]*X_vec[k];
            }
        }
    
        // Evaluate d^2f/dz^2(x,y) = Sum{k1,k2 = 0...3} c_{k1,k2} x^k1 y^k2
        // where c_{k1,k2} = a[k1 + 4*k2]
        ddfdzdz = 0;
        for(int k1 = 0; k1 < 4; k1++){
        for(int k2 = 0; k2 < 4; k2++){
            ddfdzdz+=1.0*A_vec[k1 + 4*k2]*pow(x_scaled,k1)*pow(y_scaled,k2);
        }}
}


//----------------------------------------------------//
//-------------Combined Function to Input-------------//
//------G2S Files and Generate the Interpolation------//
//----------------------------------------------------//
struct MultiDimSpline_3D Temp_Spline;       struct MultiDimSpline_3D Windu_Spline;
struct MultiDimSpline_3D Density_Spline;    struct MultiDimSpline_3D Windv_Spline;

void Spline_Multi_G2S(char* profile_prefix, char* locx, char* locy, char* option){
    SetUp_G2S_Arrays(profile_prefix, locx, locy);
    Load_G2S_Multi(profile_prefix, locx, locy, option);
    
    Temp_Spline.length_x = x_cnt;           Windu_Spline.length_x = x_cnt;          Windv_Spline.length_x = x_cnt;          Density_Spline.length_x = x_cnt;
    Temp_Spline.length_y = y_cnt;           Windu_Spline.length_y = y_cnt;          Windv_Spline.length_y = y_cnt;          Density_Spline.length_y = y_cnt;
    Temp_Spline.length_z = z_cnt;           Windu_Spline.length_z = z_cnt;          Windv_Spline.length_z = z_cnt;          Density_Spline.length_z = z_cnt;
    Temp_Spline.x_vals = x_vals;            Windu_Spline.x_vals = x_vals;           Windv_Spline.x_vals = x_vals;           Density_Spline.x_vals = x_vals;
    Temp_Spline.y_vals = y_vals;            Windu_Spline.y_vals = y_vals;           Windv_Spline.y_vals = y_vals;           Density_Spline.y_vals = y_vals;
    Temp_Spline.z_vals = z_vals;            Windu_Spline.z_vals = z_vals;           Windv_Spline.z_vals = z_vals;           Density_Spline.z_vals = z_vals;
    Temp_Spline.f_vals = T_vals;            Windu_Spline.f_vals = u_vals;           Windv_Spline.f_vals = v_vals;           Density_Spline.f_vals = rho_vals;
    Temp_Spline.f_slopes = T_slopes;        Windu_Spline.f_slopes = u_slopes;       Windv_Spline.f_slopes = v_slopes;       Density_Spline.f_slopes = rho_slopes;
    Temp_Spline.dfdx_slopes = T_slopes_dx;  Windu_Spline.dfdx_slopes = u_slopes_dx; Windv_Spline.dfdx_slopes = v_slopes_dx; Density_Spline.dfdx_slopes = rho_slopes_dx;
    Temp_Spline.dfdy_slopes = T_slopes_dy;  Windu_Spline.dfdy_slopes = u_slopes_dy; Windv_Spline.dfdy_slopes = v_slopes_dy; Density_Spline.dfdy_slopes = rho_slopes_dy;
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
double rho(double x, double y, double z){
    double x_eval = min(x, x_max);  x_eval = max(x_eval, x_min);    // Check that x_min <= x_eval <= x_max
    double y_eval = min(y, y_max);  y_eval = max(y_eval, y_min);    // Check that y_min <= y_eval <= y_max
    double z_eval = min(z, z_max);  z_eval = max(z_eval, z_min);    // Check that z_min <= z_eval <= z_max
    
    return Eval_Spline_f(x_eval,y_eval,z_eval,Density_Spline);
}

//-------------------------------------------//
//---------Thermodynamic Sound Speed---------//
//------------and its Derivatives------------//
//-------------------------------------------//

double gamR = 0.00040187; // gamma * R in km^2/s^2 * 1/K, c(x,y,z) = sqrt(gamma*r*T(x,y,z))

double c(double x, double y, double z){
    double x_eval = min(x, x_max);  x_eval = max(x_eval, x_min);    // Check that x_min <= x_eval <= x_max
    double y_eval = min(y, y_max);  y_eval = max(y_eval, y_min);    // Check that y_min <= y_eval <= y_max
    double z_eval = min(z, z_max);  z_eval = max(z_eval, z_min);    // Check that z_min <= z_eval <= z_max
    
    return sqrt(gamR * Eval_Spline_f(x_eval,y_eval,z_eval,Temp_Spline));
}

double c_diff(double x, double y, double z, int n){
    double x_eval = min(x, x_max);  x_eval = max(x_eval, x_min);    // Check that x_min <= x_eval <= x_max
    double y_eval = min(y, y_max);  y_eval = max(y_eval, y_min);    // Check that y_min <= y_eval <= y_max
    double z_eval = min(z, z_max);  z_eval = max(z_eval, z_min);    // Check that z_min <= z_eval <= z_max
    
    return gamR / (2.0 * c(x,y,z)) * Eval_Spline_df(x_eval,y_eval,z_eval,n,Temp_Spline);
}

double c_ddiff(double x, double y, double z, int n1, int n2){
    double x_eval = min(x, x_max);  x_eval = max(x_eval, x_min);    // Check that x_min <= x_eval <= x_max
    double y_eval = min(y, y_max);  y_eval = max(y_eval, y_min);    // Check that y_min <= y_eval <= y_max
    double z_eval = min(z, z_max);  z_eval = max(z_eval, z_min);    // Check that z_min <= z_eval <= z_max
    double SndSpd = c(x,y,z);
    
    return gamR / (2.0 * SndSpd) * Eval_Spline_ddf(x_eval,y_eval,z_eval,n1,n2,Temp_Spline)
             - pow(gamR,2)/(4.0 * pow(SndSpd,3)) * Eval_Spline_df(x_eval,y_eval,z_eval,n1,Temp_Spline)*Eval_Spline_df(x_eval,y_eval,z_eval,n2,Temp_Spline);
}


//------------------------------------------//
//---------East-West Wind Component---------//
//------------and its Derivatives-----------//
//------------------------------------------//
double u(double x, double y, double z){
    double x_eval = min(x, x_max);  x_eval = max(x_eval, x_min);    // Check that x_min <= x_eval <= x_max
    double y_eval = min(y, y_max);  y_eval = max(y_eval, y_min);    // Check that y_min <= y_eval <= y_max
    double z_eval = min(z, z_max);  z_eval = max(z_eval, z_min);    // Check that z_min <= z_eval <= z_max

    return Eval_Spline_f(x_eval, y_eval, z_eval, Windu_Spline);
}


double u_diff(double x, double y, double z, int n){
    double x_eval = min(x, x_max);  x_eval = max(x_eval, x_min);    // Check that x_min <= x_eval <= x_max
    double y_eval = min(y, y_max);  y_eval = max(y_eval, y_min);    // Check that y_min <= y_eval <= y_max
    double z_eval = min(z, z_max);  z_eval = max(z_eval, z_min);    // Check that z_min <= z_eval <= z_max
    
    return Eval_Spline_df(x_eval, y_eval, z_eval, n, Windu_Spline);

}

double u_ddiff(double x, double y, double z, int n1, int n2){
    double x_eval = min(x, x_max);  x_eval = max(x_eval, x_min);    // Check that x_min <= x_eval <= x_max
    double y_eval = min(y, y_max);  y_eval = max(y_eval, y_min);    // Check that y_min <= y_eval <= y_max
    double z_eval = min(z, z_max);  z_eval = max(z_eval, z_min);    // Check that z_min <= z_eval <= z_max
    
    return Eval_Spline_ddf(x_eval, y_eval, z_eval, n1, n2, Windu_Spline);
}


//--------------------------------------------//
//---------North-South Wind Component---------//
//-------------and its Derivatives------------//
//--------------------------------------------//
double v(double x, double y, double z){
    double x_eval = min(x, x_max);  x_eval = max(x_eval, x_min);    // Check that x_min <= x_eval <= x_max
    double y_eval = min(y, y_max);  y_eval = max(y_eval, y_min);    // Check that y_min <= y_eval <= y_max
    double z_eval = min(z, z_max);  z_eval = max(z_eval, z_min);    // Check that z_min <= z_eval <= z_max
    
    return Eval_Spline_f(x_eval, y_eval, z_eval, Windv_Spline);
}


double v_diff(double x, double y, double z, int n){
    double x_eval = min(x, x_max);  x_eval = max(x_eval, x_min);    // Check that x_min <= x_eval <= x_max
    double y_eval = min(y, y_max);  y_eval = max(y_eval, y_min);    // Check that y_min <= y_eval <= y_max
    double z_eval = min(z, z_max);  z_eval = max(z_eval, z_min);    // Check that z_min <= z_eval <= z_max
    
    return Eval_Spline_df(x_eval, y_eval, z_eval, n, Windv_Spline);
    
}

double v_ddiff(double x, double y, double z, int n1, int n2){
    double x_eval = min(x, x_max);  x_eval = max(x_eval, x_min);    // Check that x_min <= x_eval <= x_max
    double y_eval = min(y, y_max);  y_eval = max(y_eval, y_min);    // Check that y_min <= y_eval <= y_max
    double z_eval = min(z, z_max);  z_eval = max(z_eval, z_min);    // Check that z_min <= z_eval <= z_max
    
    return Eval_Spline_ddf(x_eval, y_eval, z_eval, n1, n2, Windv_Spline);
}


//-----------------------------------------//
//---------Vertical Wind Component---------//
//-----------and its Derivatives-----------//
//-----------------------------------------//
double w(double x, double y, double z){                         return 0.0;}
double w_diff(double x, double y, double z, int n){             return 0.0;}
double w_ddiff(double x, double y, double z, int n1, int n2){   return 0.0;}



#endif /* G2S_RANGEDEPENDENT_ATMOSPHERE_CPP_ */
