#ifndef GEOAC_EIGENRAY_CPP_
#define GEOAC_EIGENRAY_CPP_

#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "Atmo_State.h"
#include "GeoAc.Parameters.h"
#include "GeoAc.EquationSets.h"
#include "GeoAc.Solver.h"
#include "GeoAc.Interface.h"

using namespace std;

// Define global parameters used in calculating eigenrays (modifiable via command line)
bool verbose_output = false;
int eigenray_count = 0;
double d_theta_big = 0.25;
double d_theta_small = 0.002;
ofstream results;

double Modify_d_theta(double dr, double dr_dtheta){
    double width = 2.0*pow(dr_dtheta,2);
    
    return d_theta_big - (d_theta_big - d_theta_small) * exp(-dr*dr/width);
}

bool GeoAc_EstimateEigenray(double Source_Loc[3], double Receiver_Loc[2], double theta_min, double theta_max, double & theta_estimate, double & phi_estimate, double & theta_next, int bounces, double azimuth_error_limit){

    double r_rcvr = sqrt(pow(Receiver_Loc[0] - Source_Loc[0], 2) + pow(Receiver_Loc[1] - Source_Loc[1], 2));
    double phi = 180.0/3.14159 * atan2(Receiver_Loc[1] - Source_Loc[1], Receiver_Loc[0] - Source_Loc[0]);
    
    if(verbose_output){
        cout << '\t' << "Estimating eigenray angles for source-receiver separated by " << r_rcvr << " km, and azimuth " << 90.0 - phi;
        cout << " degrees from N.  Inclination limits: [" << theta_min << ", " << theta_max << "]." << '\n';
    }
    
    int	iterations = 0, k, GeoAc_length = GeoAc_ray_limit * int(1.0/(GeoAc_ds_min*10));
    theta_estimate = theta_max;
    
    double** solution;
    GeoAc_ConfigureCalcAmp(false);                      // Only calculate the ray path geometry to accelerate the search.
    GeoAc_BuildSolutionArray(solution, GeoAc_length);   // Build the array for the ray path.
    
    double r, r_prev, d_theta = d_theta_big, d_phi = 10.0;
    bool BreakCheck, success, theta_max_reached;
    
    success = 0;
    theta_max_reached = false;
    while(fabs(d_phi) > azimuth_error_limit && iterations < 5){
        GeoAc_phi =     phi*Pi/180.0;
        r = r_rcvr;
        r_prev = r_rcvr;

        // Cycle through inclinations until theta < theta_max
        for(double theta = theta_min; theta <= theta_max; theta+=d_theta){
            if(theta + d_theta >= theta_max)  theta_max_reached = true;
            
            GeoAc_theta =	theta*Pi/180.0;
            GeoAc_SetInitialConditions(solution, Source_Loc[0], Source_Loc[1], Source_Loc[2]);
        
            k = GeoAc_Propagate_RK4(solution, BreakCheck);
            if(!BreakCheck){
                for(int n_bnc = 1; n_bnc <= bounces; n_bnc++){
                    GeoAc_SetReflectionConditions(solution,k);
                    k = GeoAc_Propagate_RK4(solution, BreakCheck);
                    if(BreakCheck) break;
                }
            }

            if(verbose_output){
                cout << '\t' << '\t' << "Ray launched at " << theta << " degrees arrives at range " << sqrt(pow(solution[k][0] - Source_Loc[0],2) + pow(solution[k][1] - Source_Loc[1],2));
                cout << " km after " << bounces << " reflections." << '\t' << "Exact arrival at " << solution[k][0] << " km East, " << solution[k][1] << " km North" << '\n';
            }
        
            if(BreakCheck){
                r = r_rcvr;
                r_prev = r_rcvr;
                success = false;
            } else { r = sqrt(pow(solution[k][0] - Source_Loc[0],2) + pow(solution[k][1] - Source_Loc[1],2));}
        
            if((r - r_rcvr)*(r_prev - r_rcvr) < 0.0){
                if(iterations==0) theta_next = theta;
                
                // Inclination angles intersect the desired range, check azimuth deviation
                d_phi = (atan2(Receiver_Loc[1] - Source_Loc[1], Receiver_Loc[0] - Source_Loc[0])
                         - atan2(solution[k][1] - Source_Loc[1], solution[k][0] - Source_Loc[0]))*180.0/Pi;
                
                while(d_phi > 180.0)  d_phi-=360.0;
                while(d_phi < -180.0) d_phi+=360.0;
                
                if(fabs(d_phi) < azimuth_error_limit){
                    if(verbose_output) cout << '\t' << '\t' << "Azimuth deviation = " << d_phi << ".  Less than " << azimuth_error_limit << " degrees.  Estimates acceptable." << '\n' << '\n';
                    theta_estimate = theta - d_theta;
                    phi_estimate = phi;
                    return true;
                } else {
                    if(verbose_output) cout << '\t' << '\t' << "Azimuth deviation = " << d_phi << ".  Greater than " << azimuth_error_limit << " degrees.  Compensating and searching inclinations again." << '\n' << '\n';
                    phi += d_phi*0.9;
                    theta_min=max(theta - 7.5, theta_min);
                }
                break;
            }
            if(iterations >= 3){ d_theta = Modify_d_theta(r - r_rcvr, (r - r_prev)/(2.0 * d_theta));}
            r_prev = r;
        }
        if(theta_max_reached){
            theta_next = theta_max;
            break;
        }
        iterations++;
        if(iterations >= 1 && iterations < 3){ d_theta=d_theta_big/2.0;}

    }
    GeoAc_DeleteSolutionArray(solution, GeoAc_length);
    
    if(verbose_output) cout << '\t' << '\t' << "Reached maximum inclination angle or iteration limit." << '\n';
    return false;
}

void GeoAc_3DEigenray_LM(double Source_Loc[3], double Receiver_Loc[2], double & theta, double & phi, double freq, int bnc_cnt, int iterate_limit, char title[]){
	bool BreakCheck;
    char output_buffer [60];
    ofstream raypath;
    double D, attenuation, back_az, back_az_dev, arrival_incl, dr, dr_prev = 10000.0, travel_time;
    double M_Comps[3], nu0[3], M, nu0_xy[2];

    if(GeoAc_AtmoStrat){
        M_Comps[0] = u(Source_Loc[0], Source_Loc[1], Source_Loc[2]) / c(Source_Loc[0], Source_Loc[1], Source_Loc[2]);
        M_Comps[1] = v(Source_Loc[0], Source_Loc[1], Source_Loc[2]) / c(Source_Loc[0], Source_Loc[1], Source_Loc[2]);
        M_Comps[2] = w(Source_Loc[0], Source_Loc[1], Source_Loc[2]) / c(Source_Loc[0], Source_Loc[1], Source_Loc[2]);
    }
    
	double tolerance = 0.1;						// Absolute distance between arrivals and receiver below which to stop searching (in km)
    
	double theta_lim_step = 0.2;				// Define limiting step size and step scalar
	double phi_lim_step = 0.2;
    double step_scalar = 1.0;
    
	long double x, y, dx, dy, dx_dt, dy_dt, dx_dp, dy_dp;
	long double det, dt, dp;
    
    int	k, GeoAc_length = GeoAc_ray_limit * int(1.0/(GeoAc_ds_min*10));
    double** solution;
    
    GeoAc_ConfigureCalcAmp(true);
    GeoAc_BuildSolutionArray(solution, GeoAc_length);
    
    if(verbose_output) cout << '\t' << '\t' << "Searching for exact eigenray using auxiliary parameters." << '\n';
	for(int n = 0; n <= iterate_limit; n++){
        if(n == iterate_limit){
            if(verbose_output){cout << '\t' <<'\t' << '\t' << "Search for exact eigenray maxed out iterations.  No eigneray idenfied." << '\n' << '\n';}
            break;
        }
        
		// Initialize the solution array and calculate the ray path
        GeoAc_theta =	theta*Pi/180.0;
		GeoAc_phi = 	phi*Pi/180.0;
        
        if(GeoAc_AtmoStrat){
            nu0[0] = cos(GeoAc_theta) * cos(GeoAc_phi);
            nu0[1] = cos(GeoAc_theta) * sin(GeoAc_phi);
            nu0[2] = sin(GeoAc_theta);
        
            M = 1.0 + (nu0[0] * M_Comps[0] + nu0[1] * M_Comps[1] + nu0[2] * M_Comps[2]);
            nu0_xy[0] = nu0[0] / M;
            nu0_xy[1] = nu0[1] / M;
        }
        
        GeoAc_SetInitialConditions(solution, Source_Loc[0], Source_Loc[1], Source_Loc[2]);
        if(verbose_output) cout << '\t' << '\t' << "Plotting ray path with theta = " << theta << ", phi = " << 90.0 - phi;

        k = GeoAc_Propagate_RK4(solution, BreakCheck);
        if(BreakCheck){
            if(verbose_output) cout << '\t' << "Ray path left propagation region." << '\n';
            break;
        }
        for(int n_bnc = 1; n_bnc <= bnc_cnt; n_bnc++){
            GeoAc_SetReflectionConditions(solution,k);
            k = GeoAc_Propagate_RK4(solution, BreakCheck);
            if(BreakCheck){
                if(verbose_output) cout << '\t' << "Ray path left propagation region." << '\n';
                break;
            }
        }
        if(BreakCheck) break;
        
		// Determine arrival location and check if it's within the defined tolerance
        x = solution[k][0];     dx = Receiver_Loc[0] - x;
        y = solution[k][1];     dy = Receiver_Loc[1] - y;
        dr = sqrt(dx*dx+dy*dy);
        if(verbose_output) cout << '\t' << '\t' << "Arrival after " << bnc_cnt << " reflections at (" << x << ", " << y << "), distance to receiver = " << dr << " km." << '\n';

        if(dr < tolerance){
            sprintf(output_buffer, "%s_Eigenray-%i.dat", title, eigenray_count);
            raypath.open(output_buffer);
            
            raypath << "# x [km]";
            raypath << '\t' << "y [km]";
            raypath << '\t' << "z [km]";
            raypath << '\t' << "Geo. Atten. [dB]";
            raypath << '\t' << "Atmo. Atten. [dB]";
            raypath << '\t' << "Travel Time [s]" << '\n';
            
            attenuation = 0.0;
            travel_time = 0.0;
            
            GeoAc_SetInitialConditions(solution, Source_Loc[0], Source_Loc[1], Source_Loc[2]);
            k = GeoAc_Propagate_RK4(solution, BreakCheck);
            
            for(int m = 1; m < k ; m++){
                GeoAc_TravelTimeSegment(travel_time, solution, m-1,m);
                GeoAc_SB_AttenSegment(attenuation, solution, m-1, m, freq);
                
                if(m % 25 == 0){
                    raypath << solution[m][0];
                    raypath << '\t' << solution[m][1];
                    raypath << '\t' << max(solution[m][2],0.0);
                    raypath << '\t' << 20.0*log10(GeoAc_Amplitude(solution,m));
                    raypath << '\t' << -attenuation;
                    raypath << '\t' << travel_time << '\n';

                }
            }
            for(int n_bnc = 1; n_bnc <= bnc_cnt; n_bnc++){
                GeoAc_SetReflectionConditions(solution,k);
                k = GeoAc_Propagate_RK4(solution, BreakCheck);

                for(int m = 1; m < k; m++){
                    GeoAc_TravelTimeSegment(travel_time, solution, m-1,m);
                    GeoAc_SB_AttenSegment(attenuation, solution, m-1, m, freq);
                    
                    if(m % 25 == 0){
                        raypath << solution[m][0];
                        raypath << '\t' << solution[m][1];
                        raypath << '\t' << solution[m][2];
                        raypath << '\t' << 20.0*log10(GeoAc_Amplitude(solution,m));
                        raypath << '\t' << attenuation;
                        raypath << '\t' << travel_time << '\n';
                    }
                }
                
            }
            raypath.close();
            
            if(GeoAc_AtmoStrat){
                back_az = (90.0 - GeoAc_phi * 180.0/Pi) + 180.0;
                arrival_incl = - asin(c(solution[k][0], solution[k][1], z_grnd) / c(Source_Loc[0], Source_Loc[1], Source_Loc[2]) * solution[k][3]) * 180.0 / Pi;
            } else {
                back_az = 90.0 - atan2(-solution[k][4], -solution[k][3]) * 180.0 / Pi;
                arrival_incl = - asin(c(solution[k][0], solution[k][1], z_grnd) / c(Source_Loc[0], Source_Loc[1], Source_Loc[2]) * solution[k][5]) * 180.0 / Pi;
            }
            
            back_az_dev = back_az - (90.0 - atan2(Source_Loc[1] - Receiver_Loc[1], Source_Loc[0] - Receiver_Loc[0]) * 180.0/Pi);
            while(back_az > 180.0)      back_az-=360.0;
            while(back_az < -180.0)     back_az+=360.0;
            while(back_az_dev > 180.0)  back_az_dev-=360.0;
            while(back_az_dev < -180.0) back_az_dev+=360.0;
            
            if(!verbose_output) cout << '\t' << "Eigenray identified:" << '\t' << "theta, phi = " << setprecision(8) << theta << ", " << 90.0 - phi << " degrees." << '\n';
            if(verbose_output){
                cout << '\t' << '\t' << "Eigenray Identified:" << '\n';
                cout << '\t' << '\t' << '\t' << "theta, phi = " << setprecision(8) << theta << ", " << 90.0 - phi << " degrees." << '\n';
                cout << '\t' << '\t' << '\t' << "Travel Time = " << travel_time << " seconds." << '\n';
                cout << '\t' << '\t' << '\t' <<  "Celerity = " << sqrt(pow(solution[k][0] - Source_Loc[0],2) + pow(solution[k][1] - Source_Loc[1],2))/travel_time << " km/s." << '\n';
                cout << '\t' << '\t' << '\t' <<  "Amplitude (geometric) = " << 20.0*log10(GeoAc_Amplitude(solution,k)) << " dB." << '\n';
                cout << '\t' << '\t' << '\t' <<  "Atmospheric Attenuation = " << -attenuation << " dB." << '\n';
                cout << '\t' << '\t' << '\t' << "Arrival inclination = " << arrival_incl << " degrees." << '\n';
                cout << '\t' << '\t' << '\t' <<  "Azimuth to source = " << 90.0 - atan2(Source_Loc[1] - Receiver_Loc[1], Source_Loc[0] - Receiver_Loc[0]) * 180.0/Pi << '\n';
                cout << '\t' << '\t' << '\t' <<  "Back Azimuth of arrival = " << back_az << '\n';
                cout << '\t' << '\t' << '\t' <<  "Azimuth Deviation = " << back_az_dev  << " degrees." << '\n' << '\n';
            }
            
            results << "Eigenray-" << eigenray_count << ".  " << bnc_cnt << " bounce(s)." << '\n';
            results << '\t' << "theta, phi = " << setprecision(8) << theta << ", " << 90.0 - phi << " degrees." << '\n';
            results << '\t' << "Travel Time = " << travel_time << " seconds." << '\n';
            results << '\t' << "Celerity = " << sqrt(pow(solution[k][0] - Source_Loc[0],2) + pow(solution[k][1] - Source_Loc[1],2))/travel_time << " km/s." << '\n';
            results << '\t' << "Amplitude (geometric) = " << 20.0*log10(GeoAc_Amplitude(solution,k)) << " dB." << '\n';
            results << '\t' << "Atmospheric attenuation = " << -attenuation << " dB." << '\n';
            results << '\t' << "Arrival inclination = " << arrival_incl << " degrees." << '\n';
            results << '\t' << "Azimuth to source = " << 90.0 - atan2(Source_Loc[1] - Receiver_Loc[1], Source_Loc[0] - Receiver_Loc[0]) * 180.0/Pi << '\n';
            results << '\t' << "Back azimuth of arrival = " << back_az << '\n';
            results << '\t' << "Azimuth deviation = " << back_az_dev << " degrees." << '\n' << '\n';
            
            eigenray_count++;
            
            break;
        } else if(n > 0 && dr > dr_prev){
            // If the range to the receiver has increased, undo the previous changes to theta and phi,
            // half the step scalar and repeat the step using the new scaled increments
            theta-=dt*step_scalar;
            phi-=dp*step_scalar;
            
            step_scalar/=2.0;
            
            if(sqrt(dt*dt+dp*dp) * step_scalar < 1.0e-12){
                if (verbose_output) cout << '\t' << '\t' <<  '\t' << "Step size too small, psuedo-critical ray path likely." << '\n' << '\n';
                break;
            }
        } else {
            step_scalar = min(1.0, step_scalar * 1.25);
            
            // Calculate the transformation matrix to obtain dt and dp from dx and dy
            if(GeoAc_AtmoStrat){
                dx_dt = solution[k][4] - nu0_xy[0] / solution[k][3] * solution[k][6];
                dy_dt = solution[k][5] - nu0_xy[1] / solution[k][3] * solution[k][6];
                dx_dp = solution[k][8] - nu0_xy[0] / solution[k][3] * solution[k][10];
                dy_dp = solution[k][9] - nu0_xy[1] / solution[k][3] * solution[k][10];
            } else {
                dx_dt = solution[k][6] - solution[k][3] / solution[k][5] * solution[k][8];
                dy_dt = solution[k][7] - solution[k][4] / solution[k][5] * solution[k][8];
                dx_dp = solution[k][12] - solution[k][3] / solution[k][5] * solution[k][14];
                dy_dp = solution[k][13] - solution[k][4] / solution[k][5] * solution[k][14];
            }
            det = dx_dt * dy_dp - dx_dp * dy_dt;
            
            dt = 1.0/det * (dy_dp * dx - dx_dp * dy) * 180.0/Pi;
            dp = 1.0/det * (dx_dt * dy - dy_dt * dx) * 180.0/Pi;
                       
            // Correct the steps if they are too large or have gone outside the acceptable range
            if(dt > theta_lim_step)		dt =  theta_lim_step;   if(dp > phi_lim_step)		dp =  phi_lim_step;
            if(dt < -theta_lim_step)	dt = -theta_lim_step;   if(dp < -phi_lim_step)		dp = -phi_lim_step;
            
            // Update the angles, copy the current dr into dr_prev
            theta += dt * step_scalar;
            phi += dp * step_scalar;
            dr_prev = dr;
        }
        // Clear the solution array and prepare to trace the new ray path
        GeoAc_ClearSolutionArray(solution,k);
	}
    GeoAc_DeleteSolutionArray(solution, GeoAc_length);
}

#endif /* GEOAC_EIGENRAY_CPP_ */
