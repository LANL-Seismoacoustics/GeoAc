# ifndef G2S_MULTIDIMSPLINE3D_H_
# define G2S_MULTIDIMSPLINE3D_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>

using namespace std;

//----------------------------------------//
//------Parameters for Interpolation------//
//----------------------------------------//
extern int x_cnt;  // Number of x data points (assumed constant for all y and z values)
extern int y_cnt;  // Number of y data points (assumed constant for all x and z values)
extern int z_cnt;  // Number of z data points (assumed constant for all x and y values)

extern int x_accel;         // Index for x[k] < x' < x[k+1]
extern int y_accel;         // Index for y[k] < y' < y[k+1]
extern int z_accel;         // Index for z[k] < z' < z[k+1]

extern double* x_vals;     // x_i elements (x_cnt length)
extern double* y_vals;     // y_j elements (y_cnt length)
extern double* z_vals;     // z_k elements (z_cnt length)

extern double*** T_vals;           // Temperature at (x_i, y_j, z_k) (x_cnt x y_cnt x z_cnt)
extern double*** T_slopes;         // Slopes for vertical splines of temperature at each x[], y[] node
extern double*** T_slopes_dx;      // Slopes for d/dx of vertical splines of temperature at each x[], y[] node
extern double*** T_slopes_dy;      // Slopes for d/dy of vertical splines of temperature at each x[], y[] node

extern double*** u_vals;           // E-W winds at (x_i, y_j, z_k) (x_cnt x y_cnt x z_cnt)
extern double*** u_slopes;         // Slopes for vertical splines of temperature at each x[], y[] node
extern double*** u_slopes_dx;      // Slopes for d/dx of vertical splines of temperature at each x[], y[] node
extern double*** u_slopes_dy;      // Slopes for d/dy of vertical splines of temperature at each x[], y[] node

extern double*** v_vals;           // N-S winds at (x_i, y_j, z_k) (x_cnt x y_cnt x z_cnt)
extern double*** v_slopes;         // Slopes for vertical splines of temperature at each x[], y[] node
extern double*** v_slopes_dx;      // Slopes for d/dx of vertical splines of temperature at each x[], y[] node
extern double*** v_slopes_dy;      // Slopes for d/dy of vertical splines of temperature at each x[], y[] node

extern double*** rho_vals;         // Density at (x_i, y_j, z_k) (x_cnt x y_cnt x z_cnt)
extern double*** rho_slopes;         // Slopes for vertical splines of temperature at each x[], y[] node
extern double*** rho_slopes_dx;      // Slopes for d/dx of vertical splines of temperature at each x[], y[] node
extern double*** rho_slopes_dy;      // Slopes for d/dy of vertical splines of temperature at each x[], y[] node

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
void SetUp_G2S_Arrays(char*, char*, char*);           // Take the file prefix in char* and determine r_cnt, z_cnt; Set arrays to appropriate sizes
void Load_G2S_Multi(char*, char*, char*);     // Take the file prefix in char* and r values in file name (string) and load all data into the arrays
void Clear_G2S_Arrays();                // Clear and delete the G2S arrays

//-------------------------------------//
//----------Three-Dimensional----------//
//-------Interpolation Structure-------//
//-------------------------------------//
extern double BiCubic_ConversionMatrix [16][16];

struct MultiDimSpline_3D{
	int length_x;           // Number of x nodes
    int length_y;           // Number of y nodes
    int length_z;           // Number of z nodes
    
    int accel[3];          // Indices of x,y,z point
    
	double* x_vals;         // 1D array of x values
    double* y_vals;         // 1D array of y values
    double* z_vals;         // 1D array of z values
	
    double*** f_vals;       // 3D array of f(x,y,z) values, f_vals[i][j][k] = f(x[i], y[j], z[k])
	double*** f_slopes;     // Slopes used to generate natural cubic spline solution at each x and y = constant node, describes f(x[i], y[i], z)
    double*** dfdx_slopes;  // Slopes used to generate df/dx at each node point, describes df/dx @ x[i], y[j], z
    double*** dfdy_slopes;  // Slopes used to generate df/dy at each node point, describes df/dx @ x[i], y[j], z
};

//----------------------------------------------//
//------------Functions to Build and------------//
//-------Clear the Input and Slope Arrays-------//
//----------------------------------------------//
void BuildInputArrays(struct MultiDimSpline_3D &);
void BuildSlopesArrays(struct MultiDimSpline_3D &);

void ClearInputArrays(struct MultiDimSpline_3D &);
void ClearSlopesArrays(struct MultiDimSpline_3D &);

//--------------------------------------//
//-----------Functions to Set-----------//
//-------the Interpolation Slopes-------//
//--------------------------------------//
void Set_Slopes_Multi(struct MultiDimSpline_3D &);      // Set the slopes at all constant r values in the f_vals array

//-----------------------------------------------//
//------Functions to Find the Segment Index------//
//-----and Evaluate the 1D Vertical Splines------//
//-----------------------------------------------//
int Find_Segment(double, double*, int, int&);                                             // Find k such that x[k] <= x <= x]k+1]

double Eval_Vert_Spline_f(double, struct MultiDimSpline_3D, int, int, int);         // Evaluate the vertical spline of f at x_i, y_j
double Eval_Vert_Spline_dfdx(double, struct MultiDimSpline_3D, int, int, int);      // Evaluate the vertical spline of df/dx at x_i, y_j
double Eval_Vert_Spline_dfdy(double, struct MultiDimSpline_3D, int, int, int);      // Evaluate the vertical spline of df/dy at x_i, y_j
double Eval_Vert_Spline_dfdz(double, struct MultiDimSpline_3D, int, int, int);      // Evaluate the vertical spline of df/dz at x_i, y_j
double Eval_Vert_Spline_ddfdxdz(double, struct MultiDimSpline_3D, int, int, int);   // Evaluate the vertical spline of d^2f/dxdz at x_i, y_j
double Eval_Vert_Spline_ddfdydz(double, struct MultiDimSpline_3D, int, int, int);   // Evaluate the vertical spline of d^2f/dydz at x_i, y_j
double Eval_Vert_Spline_ddfdzdz(double, struct MultiDimSpline_3D, int, int, int);   // Evaluate the vertical spline of d^2f/dz^2 at x_i, y_j

double BiCubic_Deriv_dfdx(double, struct MultiDimSpline_3D, int, int, int);         // Evaluate df/dx using discrete derivatives
double BiCubic_Deriv_dfdy(double, struct MultiDimSpline_3D, int, int, int);         // Evaluate df/dy using discrete derivatives

double BiCubic_Deriv_ddfdxdx(double, struct MultiDimSpline_3D, int, int, int);      // Evaluate d^2f/dx^2 using discrete derivatives
double BiCubic_Deriv_ddfdydy(double, struct MultiDimSpline_3D, int, int, int);      // Evaluate d^2f/dy^2 using discrete derivatives
double BiCubic_Deriv_ddfdxdy(double, struct MultiDimSpline_3D, int, int, int);      // Evaluate d^2f/dxdy using discrete derivatives

double BiCubic_Deriv_dddfdxdxdy(double, struct MultiDimSpline_3D, int, int, int);   // Evaluate d^3f/dx^2dy using discrete derivatives
double BiCubic_Deriv_dddfdxdxdz(double, struct MultiDimSpline_3D, int, int, int);   // Evaluate d^3f/dx^2dz using discrete derivatives
double BiCubic_Deriv_dddfdxdydy(double, struct MultiDimSpline_3D, int, int, int);   // Evaluate d^3f/dxdy^2 using discrete derivatives
double BiCubic_Deriv_dddfdxdydz(double, struct MultiDimSpline_3D, int, int, int);   // Evaluate d^3f/dxdydz using discrete derivatives
double BiCubic_Deriv_dddfdxdzdz(double, struct MultiDimSpline_3D, int, int, int);   // Evaluate d^3f/dxdz^2 using discrete derivatives
double BiCubic_Deriv_dddfdydydz(double, struct MultiDimSpline_3D, int, int, int);   // Evaluate d^3f/dy^2dz using discrete derivatives
double BiCubic_Deriv_dddfdydzdz(double, struct MultiDimSpline_3D, int, int, int);   // Evaluate d^3f/dydz^2 using discrete derivatives

double BiCubic_Deriv_ddddfdxdxdydz(double, struct MultiDimSpline_3D, int, int, int);    // Evaluate d^4f/dx^2dydz using discrete derivatives
double BiCubic_Deriv_ddddfdxdydydz(double, struct MultiDimSpline_3D, int, int, int);    // Evaluate d^4f/dxdy^2dz using discrete derivatives
double BiCubic_Deriv_ddddfdxdydzdz(double, struct MultiDimSpline_3D, int, int, int);    // Evaluate d^4f/dxdydz^2 using discrete derivatives

//----------------------------------------------------//
//-----------Evauation of the 2D Spline and-----------//
//-------its First and Second Order Derivatives-------//
//----------------------------------------------------//
double Eval_Spline_f(double, double, double, struct MultiDimSpline_3D &);                           // Evaluate the spline at x, y, z
double Eval_Spline_df(double, double, double, int, struct MultiDimSpline_3D &);                     // Evaluate df/dx_i at x, y, z
double Eval_Spline_ddf(double, double, double, int, int, struct MultiDimSpline_3D &);               // Evaluate d^2 f/dx_i dx_j at x, y, z
void Eval_Spline_AllOrder1(double, double, double, struct MultiDimSpline_3D &, double &, double &, double &, double &);                  // Evaluate f, df/dx, df/dy, and df/dz
void Eval_Spline_AllOrder2(double, double, double, struct MultiDimSpline_3D &, double &, double &, double &, double &, double &, double &, double &, double &, double &, double &);
        // Evaluate f, df/dx, df/dy, df/dz, d^2f/dx^2, d^2fdy^2, d^2f/dz^2, d^2f/dxdy, d^2f/dxdz, and d^2f/dydz


//----------------------------------------------------//
//-------------Combined Function to Input-------------//
//------G2S Files and Generate the Interpolation------//
//----------------------------------------------------//
extern struct MultiDimSpline_3D Temp_Spline;    // Temperature interpolation
extern struct MultiDimSpline_3D Windu_Spline;   // E-W wind component interpolation
extern struct MultiDimSpline_3D Windv_Spline;   // N-S wind component interpolation
extern struct MultiDimSpline_3D Density_Spline; // Density interpolation

void Spline_Multi_G2S(char*, char*, char*);     // Input the profile prefix and x_vals/y_vals information to create splines of 3D G2S data
void ClearAll();                                // Function to clear the G2S arrays and slopes from interpolation

extern double gamR;
#endif /* G2S_MULTIDIMSPLINE3D_H_ */