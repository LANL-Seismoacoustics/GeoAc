# ifndef G2S_GLIOBAL_RNGDEP_ATMOSPHERE_H_
# define G2S_GLIOBAL_RNGDEP_ATMOSPHERE_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>

using namespace std;

//-----------------------------------------------//
//---------Define the Propagation Region---------//
//-----------------------------------------------//
extern double r_min, r_max;
extern double t_min, t_max;
extern double p_min, p_max;

//----------------------------------------//
//------Parameters for Interpolation------//
//----------------------------------------//
extern int r_cnt;       // Number of altitude(r - R_Ground) data points
extern int t_cnt;       // Number of latitude (theta -> t) data points
extern int p_cnt;       // Number of longitude (phi -> p) data points

extern int r_accel;     // Index for r[k] < r' < r[k+1]
extern int t_accel;     // Index for t[k] < t' < t[k+1]
extern int p_accel;     // Index for p[k] < p' < p[k+1]

extern double* r_vals;  // r_i elements (r_cnt length)
extern double* t_vals;  // t_j elements (t_cnt length)
extern double* p_vals;  // p_k elements (p_cnt length)

extern double*** T_vals;        // Temperature at (r_i, t_j, p_k) (r_cnt x t_cnt x p_cnt)
extern double*** T_slopes;      // Slopes for vertical splines of temperature at each t[], p[] node
extern double*** T_slopes_dt;   // Slopes for d/dp of vertical splines of temperature at each t[], p[] node
extern double*** T_slopes_dp;   // Slopes for d/dt of vertical splines of temperature at each t[], p[] node

extern double*** u_vals;        // E-W winds at (r_i, t_j, p_k) (r_cnt x t_cnt x p_cnt)
extern double*** u_slopes;      // Slopes for vertical splines of temperature at each t[], p[] node
extern double*** u_slopes_dt;   // Slopes for d/dt of vertical splines of temperature at each t[], p[] node
extern double*** u_slopes_dp;   // Slopes for d/dp of vertical splines of temperature at each t[], p[] node

extern double*** v_vals;        // N-S winds at (r_i, t_j, p_k) (r_cnt x t_cnt x p_cnt)
extern double*** v_slopes;      // Slopes for vertical splines of temperature at each t[], p[] node
extern double*** v_slopes_dt;   // Slopes for d/dt of vertical splines of temperature at each t[], p[] node
extern double*** v_slopes_dp;   // Slopes for d/dp of vertical splines of temperature at each t[], p[] node

extern double*** rho_vals;      // Density at (r_i, t_j, p_k) (r_cnt x t_cnt x p_cnt)
extern double*** rho_slopes;    // Slopes for vertical splines of temperature at each t[], p[] node
extern double*** rho_slopes_dt; // Slopes for d/dt of vertical splines of temperature at each t[], p[] node
extern double*** rho_slopes_dp; // Slopes for d/dp of vertical splines of temperature at each t[], p[] node

//----------------------------------------//
//----------File IO Manipulation----------//
//----------------------------------------//
int file_length(string);            // Function to determine length (number of rows) in input file
int file_width(string);             // Function to determine width (number of columns) in input file
bool Check_G2S_Format(string);      // Function to check that there are columns for z, u, v, T, p, rho
                                    // Currently only check number of columns; might program to check for "realistic" values

//----------------------------------------//
//---------G2S Array Manipulation---------//
//----------------------------------------//
void SetUp_G2S_Arrays(char*, char*, char*);         // Take the file prefix in char* and determine r_cnt, z_cnt; Set arrays to appropriate sizes
void Load_G2S_Multi(char*, char*, char*, char*);    // Take the file prefix in char* and r values in file name (string) and load all data into the arrays
void Clear_G2S_Arrays();                            // Clear and delete the G2S arrays

//-------------------------------------//
//----------Three-Dimensional----------//
//-------Interpolation Structure-------//
//-------------------------------------//
extern double BiCubic_ConversionMatrix [16][16];

struct MultiDimSpline_3DGlobal{
	int length_r;           // Number of r nodes
    int length_t;           // Number of theta (t) nodes; corresponds to latitudes
    int length_p;           // Number of phi (p) nodes; corresponds to longitudes
    
    int accel[3];          // Indices of current r, t, p point
    
	double* r_vals;         // 1D array of r values
    double* t_vals;         // 1D array of t values
    double* p_vals;         // 1D array of p values
	
    double*** f_vals;       // 3D array of f(r,t,p) values, f_vals[i][j][k] = f(r[i], t[j], p[k])
	double*** f_slopes;     // Slopes used to generate natural cubic spline solution at each t and p = constant node, describes f(r, t[i], p[i])
    double*** dfdt_slopes;  // Slopes used to generate df/dt at each node point, describes df/dt @ r, t[i], p[j]
    double*** dfdp_slopes;  // Slopes used to generate df/dp at each node point, describes df/dp @ r, t[i], p[j]
};

//----------------------------------------------//
//------------Functions to Build and------------//
//-------Clear the Input and Slope Arrays-------//
//----------------------------------------------//
void BuildInputArrays(struct MultiDimSpline_3DGlobal &);
void BuildSlopesArrays(struct MultiDimSpline_3DGlobal &);

void ClearInputArrays(struct MultiDimSpline_3DGlobal &);
void ClearSlopesArrays(struct MultiDimSpline_3DGlobal &);

//--------------------------------------//
//-----------Functions to Set-----------//
//-------the Interpolation Slopes-------//
//--------------------------------------//
void Set_Slopes_Multi(struct MultiDimSpline_3DGlobal &);      // Set the slopes at all constant r values in the f_vals array

//-----------------------------------------------//
//------Functions to Find the Segment Index------//
//-----and Evaluate the 1D Vertical Splines------//
//-----------------------------------------------//
int Find_Segment(double, double*, int, int&);                                             // Find k such that x[k] <= x <= x]k+1]

double Eval_Vert_Spline_f(double, struct MultiDimSpline_3DGlobal, int, int, int);         // Evaluate the vertical spline of f at t_i, p_j

double Eval_Vert_Spline_dfdr(double, struct MultiDimSpline_3DGlobal, int, int, int);      // Evaluate the vertical spline of df/dr at t_i, p_j

double Eval_Vert_Spline_ddfdrdr(double, struct MultiDimSpline_3DGlobal, int, int, int);   // Evaluate the vertical spline of d^2f/dr^2 at t_i, p_j
double Eval_Vert_Spline_ddfdrdt(double, struct MultiDimSpline_3DGlobal, int, int, int);   // Evaluate the vertical spline of d^2f/drdt at t_i, p_j
double Eval_Vert_Spline_ddfdrdp(double, struct MultiDimSpline_3DGlobal, int, int, int);   // Evaluate the vertical spline of d^2f/drdp at t_i, p_j

double BiCubic_Deriv_dfdt(double, struct MultiDimSpline_3DGlobal, int, int, int);         // Evaluate df/dt using discrete derivatives
double BiCubic_Deriv_dfdp(double, struct MultiDimSpline_3DGlobal, int, int, int);         // Evaluate df/dp using discrete derivatives

double BiCubic_Deriv_ddfdtdt(double, struct MultiDimSpline_3DGlobal, int, int, int);      // Evaluate d^2f/dt^2 using discrete derivatives
double BiCubic_Deriv_ddfdpdp(double, struct MultiDimSpline_3DGlobal, int, int, int);      // Evaluate d^2f/dp^2 using discrete derivatives
double BiCubic_Deriv_ddfdtdp(double, struct MultiDimSpline_3DGlobal, int, int, int);      // Evaluate d^2f/dtdp using discrete derivatives

double BiCubic_Deriv_dddfdrdrdt(double, struct MultiDimSpline_3DGlobal, int, int, int);   // Evaluate d^3f/dr^2dt using discrete derivatives
double BiCubic_Deriv_dddfdrdrdp(double, struct MultiDimSpline_3DGlobal, int, int, int);   // Evaluate d^3f/dr^2dp using discrete derivatives
double BiCubic_Deriv_dddfdrdtdp(double, struct MultiDimSpline_3DGlobal, int, int, int);   // Evaluate d^3f/drdtdp using discrete derivatives
double BiCubic_Deriv_dddfdtdtdp(double, struct MultiDimSpline_3DGlobal, int, int, int);   // Evaluate d^3f/dt^2dp using discrete derivatives
double BiCubic_Deriv_dddfdtdpdp(double, struct MultiDimSpline_3DGlobal, int, int, int);   // Evaluate d^3f/dtdp^2 using discrete derivatives

double BiCubic_Deriv_ddddfdrdrdtdp(double, struct MultiDimSpline_3DGlobal, int, int, int);    // Evaluate d^4f/dx^2dydz using discrete derivatives

//----------------------------------------------------//
//-----------Evauation of the 2D Spline and-----------//
//-------its First and Second Order Derivatives-------//
//----------------------------------------------------//
double Eval_Spline_f(double, double, double, struct MultiDimSpline_3DGlobal &);                           // Evaluate the spline at x, y, z
double Eval_Spline_df(double, double, double, int, struct MultiDimSpline_3DGlobal &);                     // Evaluate df/dx_i at x, y, z
double Eval_Spline_ddf(double, double, double, int, int, struct MultiDimSpline_3DGlobal &);               // Evaluate d^2 f/dx_i dx_j at x, y, z
void Eval_Spline_AllOrder1(double, double, double, struct MultiDimSpline_3DGlobal &, double &, double &, double &, double &);                  // Evaluate f, df/dx, df/dy, and df/dz
void Eval_Spline_AllOrder2(double, double, double, struct MultiDimSpline_3DGlobal &, double &, double &, double &, double &, double &, double &, double &, double &, double &, double &);
        // Evaluate f, df/dx, df/dy, df/dz, d^2f/dx^2, d^2fdy^2, d^2f/dz^2, d^2f/dxdy, d^2f/dxdz, and d^2f/dydz



//----------------------------------------------------//
//-------------Combined Function to Input-------------//
//------G2S Files and Generate the Interpolation------//
//----------------------------------------------------//
extern struct MultiDimSpline_3DGlobal Temp_Spline;      // Temperature interpolation
extern struct MultiDimSpline_3DGlobal Windu_Spline;     // E-W wind component interpolation
extern struct MultiDimSpline_3DGlobal Windv_Spline;     // N-S wind component interpolation
extern struct MultiDimSpline_3DGlobal Density_Spline;   // Density interpolation

void Spline_Multi_G2S(char*, char*, char*, char*);      // Input the profile prefix and x_vals/y_vals information to create splines of 3D G2S data
void ClearAll();                                        // Function to clear the G2S arrays and slopes from interpolation


#endif /* G2S_GLIOBAL_RNGDEP_ATMOSPHERE_H_ */
