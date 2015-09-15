# ifndef G2S_SPLINE1D_H_
# define G2S_SPLINE1D_H_

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
extern int z_cnt;           // Number of vertical points
extern double* z_vals;      // z_k elements

extern double* T_vals;      // Temperature at z_k
extern double* u_vals;      // E-W winds at z_k
extern double* v_vals;      // N-S winds at z_k
extern double* rho_vals;    // Density at z_k

extern double* T_slopes;    // Slopes to compute interpoalted temperature
extern double* u_slopes;    // Slopes to compute interpolated E-W winds
extern double* v_slopes;    // Slopes to compute interpolated N-S winds
extern double* rho_slopes;  // Slopes to compute interpolated density

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
void SetUp_G2S_Arrays(char*);           // Take the file prefix in char* and determine r_cnt, z_cnt; Set arrays to appropriate sizes
void Load_G2S_Single(char*, char*);     // Take the file prefix in char* and r values in file name (string) and load all data into the arrays
void Clear_G2S_Arrays();                // Clear and delete the G2S arrays

//-------------------------------------//
//-----------One-Dimensional-----------//
//-------Interpolation Structure-------//
//-------------------------------------//
struct NaturalCubicSpline_1D{
	int length;			// Length of input files (x and f(x))
    int accel;          // Index of previous table look up; used to increase spline speed
	double* x_vals;     // 1D array of x values
	double* f_vals;     // 1D array of f(x) values, f_vals[i] = f(x[i])
	double* slopes;     // Slopes used to generate natural cubic spline solution
};

//----------------------------------------------//
//------------Functions to Build and------------//
//-------Clear the Input and Slope Arrays-------//
//----------------------------------------------//
void BuildSlopeArrays(struct NaturalCubicSpline_1D &);
void ClearSlopeArrays(struct NaturalCubicSpline_1D &);

//--------------------------------------//
//-----------Functions to Set-----------//
//-------the Interpolation Slopes-------//
//--------------------------------------//
void Set_Slopes(struct NaturalCubicSpline_1D &);                        // Set the slopes at all constant r values in the f_vals array

//-----------------------------------------------//
//------Functions to Find the Segment Index------//
//-----and Evaluate the 1D Vertical Splines------//
//-----------------------------------------------//
int Find_Segment(double, double*, int, int &);                   // Find k such that x[k] <= x <= x[k+1]
double Eval_Spline_f(double, struct NaturalCubicSpline_1D &);    // Evaluate f @ x
double Eval_Spline_df(double, struct NaturalCubicSpline_1D &);   // Evaluate df/dx
double Eval_Spline_ddf(double, struct NaturalCubicSpline_1D &);  // Evaluate d^2f/dx^2

//----------------------------------------------------//
//-------------Combined Function to Input-------------//
//------G2S Files and Generate the Interpolation------//
//----------------------------------------------------//
extern struct NaturalCubicSpline_1D Temp_Spline;
extern struct NaturalCubicSpline_1D Windu_Spline;
extern struct NaturalCubicSpline_1D Windv_Spline;
extern struct NaturalCubicSpline_1D Density_Spline;

void Spline_Single_G2S(char*, char*);   // Input the profile prefix and r_vals information to create splines of 2D G2S data
void ClearAll();                        // Function to clear the G2S arrays and slopes from interpolation


#endif /* G2S_SPLINE2D_H_ */