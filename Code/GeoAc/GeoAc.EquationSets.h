#ifndef GEOAC_EQSETS_H_
#define GEOAC_EQSETS_H_

using namespace std;

void    GeoAc_SetSystem();                                              // Function to set dimensions and stratified

void    GeoAc_SetInitialConditions(double**&);                          // Function sets initial conditions for a source at the origin with launch angles GeoAc_theta, GeoAc_phi
void    GeoAc_SetInitialConditions(double**&, double, double);          // Function sets initial conditions for a source at (r_0, z_0) with launch angles GeoAc_theta, GeoAc_phi
void    GeoAc_SetInitialConditions(double**&, double, double, double);  // Function sets initial conditions for a source at (x_0, y_0, z_0) with launch angles GeoAc_theta, GeoAc_phi

void    GeoAc_ApproximateIntercept(double**, int, double*&);            // Function uses linear interpolation to estimate ground intercept values
void    GeoAc_SetReflectionConditions(double**&, int);                  // Function sets reflection conditions

void    GeoAc_UpdateSources(double, double*);                           // Function to update source equation parameters which are used multiple times
double  GeoAc_Set_ds(double*);                                          // Function to modify the solver step size
double  GeoAc_EvalSrcEq(double, double*, int);                          // Function to evaluate the source equations

double  GeoAc_EvalHamiltonian(double, double*);                         // Function to evaluate the hamiltonian (eikonal) at s, current_values
double  GeoAc_EvalHamiltonian(double**, int);                           // Function to evaluate the hamiltonian at index of solution**

bool    GeoAc_BreakCheck(double **, int);                               // Function to check for ray leaving propagation region
bool    GeoAc_GroundCheck(double **, int);                              // Function to check for ray returning to ground

double  GeoAc_TravelTime(double **, int);                               // Function to calculate travel time from ray origin to s = ds * index
void    GeoAc_TravelTimeSegment(double&, double**, int, int);           // Function increment travel time from ds*k_1 to d2*k_2
double  GeoAc_SB_Atten(double **, int, double);                         // Function to calculate atmospheric attenuation through a ray path
void    GeoAc_SB_AttenSegment(double&, double **, int , int, double);   // Function to increment atmospheric attenuation through a ray path
double  GeoAc_Jacobian(double **, int);                                 // Function to calculate the Jacobian determinant
double  GeoAc_Amplitude(double **, int);                                // Function to calculate the transport equation coefficient
int     GeoAc_CausticCnt(double **, int, int);                          // Function to count caustics from ray origin to s = ds * index

#endif /* GEOAC_EQSETS_H_ */
