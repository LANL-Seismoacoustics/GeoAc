#ifndef GEOAC_EIGENRAY_H_
#define GEOAC_EIGENRAY_H_

#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "GeoAc.Parameters.h"
#include "Atmo_State.h"
#include "GeoAc.EquationSets.h"
#include "GeoAc.Solver.h"
#include "GeoAc.Interface.h"


using namespace std;

// Define global parameters used in calculating eigenrays
bool verbose_output = false;
int eigenray_count = 0;
double d_theta_big = 0.25;
double d_theta_small = 0.002;
ofstream results;

double Calc_Bearing(double lat1, double long1, double lat2, double long2){
    double term1 = sin((long2 - long1)*Pi/180.0);
    double term2 = cos(lat1*Pi/180.0)*tan(lat2*Pi/180.0) - sin(lat1*Pi/180.0)*cos((long2-long1)*Pi/180.0);
    
    return atan2(term1, term2)*180.0/Pi;
}

double Calc_GC_Distance(double lat1, double long1, double lat2, double long2){
    double term1 = pow(sin((lat2 - lat1)*Pi/180.0/2.0),2);
    double term2 = cos(lat1*Pi/180.0) * cos(lat2*Pi/180.0) * pow(sin((long2 - long1)*Pi/180.0/2.0),2);
    
    return 2.0 * r_earth * asin(sqrt(term1+term2));
}

double Modify_d_theta(double dr, double dr_dtheta){
    double width = 1.0/2.0 * pow(dr_dtheta,2);
    
    return d_theta_big - (d_theta_big - d_theta_small) * exp(-dr*dr/width);
}


bool GeoAc_EstimateEigenray(double Source_Loc[3], double Receiver_Loc[2], double theta_min, double theta_max, double & theta_estimate, double & phi_estimate, double & theta_next, int bounces, double azimuth_error_limit){
    
    double GC_r_rcvr = Calc_GC_Distance(Source_Loc[0], Source_Loc[1], Receiver_Loc[0], Receiver_Loc[1]);
    double phi = Calc_Bearing(Source_Loc[0], Source_Loc[1], Receiver_Loc[0], Receiver_Loc[1]);
    
    if(verbose_output){
        cout << '\t' << "Estimating eigenray angles for source-receiver separated by great circle distance " << GC_r_rcvr << " km, and azimuth " << phi;
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
        GeoAc_phi =     (90.0 - phi)*Pi/180.0;
        r = GC_r_rcvr;
        r_prev = GC_r_rcvr;
        
        // Cycle through inclinations until theta < theta_max
        for(double theta = theta_min; theta < theta_max; theta+=d_theta){
            if(theta+d_theta >= theta_max) theta_max_reached = true;
            
            GeoAc_theta =	theta*Pi/180.0;
            GeoAc_SetInitialConditions(solution, Source_Loc[2], Source_Loc[0]*Pi/180.0, Source_Loc[1]*Pi/180.0);
            
            k = GeoAc_Propagate_RK4(solution, BreakCheck);
            if(!BreakCheck){
                for(int n_bnc = 1; n_bnc <= bounces; n_bnc++){
                    GeoAc_SetReflectionConditions(solution,k);
                    k = GeoAc_Propagate_RK4(solution, BreakCheck);
                    if(BreakCheck) break;
                }
            }
            
            if(BreakCheck){
                r = GC_r_rcvr;
                r_prev = GC_r_rcvr;
                success = 0;
            } else {
                r = Calc_GC_Distance(Source_Loc[0], Source_Loc[1],  solution[k][1]*180.0/Pi, solution[k][2]*180.0/Pi);
            }
            
            if(verbose_output){
                cout << '\t' << '\t' << "Ray launched at inclination=" << GeoAc_theta * 180.0/Pi << " degrees arrives at range " << r;
                cout << " km after " << bounces << " bounces.  Exact arrival at " << solution[k][1]*180.0/Pi << " degrees N latitude, " << solution[k][2]*180.0/Pi << " degrees E longitude" << '\n';
            }
            
            if((r - GC_r_rcvr)*(r_prev - GC_r_rcvr) < 0.0){
                if(iterations==0) theta_next = theta;
                
                d_phi  = Calc_Bearing(Source_Loc[0], Source_Loc[1], Receiver_Loc[0], Receiver_Loc[1]);
                d_phi -= Calc_Bearing(Source_Loc[0], Source_Loc[1], solution[k][1]*180.0/Pi, solution[k][2]*180.0/Pi);
                while(d_phi > 180.0){  d_phi-=360.0;}
                while(d_phi < -180.0){ d_phi+=360.0;}
                
                if(fabs(d_phi) < azimuth_error_limit){
                    if(verbose_output) cout << '\t' << '\t' << "Azimuth deviation less than " << azimuth_error_limit << " degrees.  Estimates acceptable." << '\n' << '\n';
                    theta_estimate = theta - d_theta;
                    phi_estimate = 90.0 - phi;
                    return true;
                } else {
                    if(verbose_output) cout << '\t' << '\t' << "Azimuth deviation greater than " << azimuth_error_limit << " degrees.  Compensating and searching inclinations again." << '\n' << '\n';
                    phi += d_phi*0.9;
                    theta_min=max(theta - 7.5, theta_min);
                }
                break;
            }
            if(iterations >= 3){ d_theta = Modify_d_theta(r - GC_r_rcvr, (r - r_prev)/(2.0*d_theta));}
            r_prev = r;
        }
        if(theta_max_reached){
            theta_next = theta_max;
            break;
        }
        iterations++;
        if(iterations >= 1 && iterations < 3){ d_theta = d_theta_big/2.0;}
    }
    GeoAc_DeleteSolutionArray(solution, GeoAc_length);
    
    if(verbose_output) cout << '\t' << '\t' << "Reached maximum inclination angle or iteration limit." << '\n' << '\n';
    return false;
}


void GeoAc_3DEigenray_LM(double Source_Loc[3], double Receiver_Loc[2], double & lt, double & lp, double freq, int bnc_cnt, int iterate_limit, char title[]){
	bool BreakCheck;
    char output_buffer [60];
    ofstream raypath;
    double D, attenuation, travel_time, back_az_dev, arrival_incl;
    double dr, dr_prev = 10000.0;
    
	double tolerance = 0.1;					// Absolute distance in km between arrivals and receiver below which to stop searching
    
	double lt_lim_step = 0.2;				// Define limiting step size and step scalar
	double lp_lim_step = 0.2;
    double step_scalar = 1.0;
    
	long double lat, lon, d_lat, d_lon, d_lat_dlt, d_lon_dlt, d_lat_dlp, d_lon_dlp, det, dlt, dlp;
    
    int	iterations = 0, k, GeoAc_length = GeoAc_ray_limit * int(1.0/(GeoAc_ds_min*10));
    double** solution;
    
    GeoAc_ConfigureCalcAmp(true);
    GeoAc_BuildSolutionArray(solution, GeoAc_length);

    if(verbose_output) cout << '\t' << '\t' << "Searching for exact eigenray using auxiliary parameters." << '\n';
	for(int n = 0; n <= iterate_limit; n++){
        if(n == iterate_limit){
            if(verbose_output){cout << '\t' <<'\t' << '\t' << "Search for exact eigenray maxed out iterations.  No eigneray idenfied." << '\n';}
            break;
        }
        
        // Initialize the solution array and calculate the ray path
        GeoAc_theta =	lt*Pi/180.0;
		GeoAc_phi = 	lp*Pi/180.0;
		GeoAc_SetInitialConditions(solution, Source_Loc[2], Source_Loc[0]*Pi/180.0, Source_Loc[1]*Pi/180.0);
        if(verbose_output) cout << '\t' << '\t' << "Plotting ray path with theta = " << lt << ", phi = " << 90.0 - lp;
		
        k = GeoAc_Propagate_RK4(solution, BreakCheck);
        if(BreakCheck) break;
        for(int n_bnc = 1; n_bnc <= bnc_cnt; n_bnc++){
            GeoAc_SetReflectionConditions(solution,k);
            k = GeoAc_Propagate_RK4(solution, BreakCheck);
            if(BreakCheck) break;
        }
        if(BreakCheck) break;
        
		// Determine arrival location and check if it's within the defined tolerance
        lat = solution[k][1];   lon = solution[k][2];
        dr = Calc_GC_Distance(lat*180.0/Pi, lon*180.0/Pi, Receiver_Loc[0], Receiver_Loc[1]);
        if(verbose_output) cout << '\t' << '\t' << "Arrival at (" << setprecision(8) << lat*180.0/Pi << ", " << lon*180.0/Pi << "), distance to receiver = " << dr << " km." << '\n';

        if(dr < tolerance){
            sprintf(output_buffer, "%s_Eigenray-%i.dat", title, eigenray_count);
            raypath.open(output_buffer);

            raypath << "# z [km]";
            raypath << '\t' << "lat [deg]";
            raypath << '\t' << "lon [deg]";
            raypath << '\t' << "Geo. Atten. [dB]";
            raypath << '\t' << "Atmo. Atten. [dB]";
            raypath << '\t' << "Travel Time [s]";
            raypath << '\n';
                
            attenuation = 0.0;
            travel_time = 0.0;
            
            GeoAc_SetInitialConditions(solution, Source_Loc[2], Source_Loc[0]*Pi/180.0, Source_Loc[1]*Pi/180.0);
            k = GeoAc_Propagate_RK4(solution, BreakCheck);

            for(int m=1;m<k;m++){
                GeoAc_TravelTimeSegment(travel_time, solution, m-1,m);
                GeoAc_SB_AttenSegment(attenuation, solution, m-1, m, freq);
                    
                if(m % 25 == 0){
                    raypath << solution[m][0] - r_earth;
                    raypath << '\t' << setprecision(8) << solution[m][1] * 180.0/Pi;
                    raypath << '\t' << setprecision(8) << solution[m][2] * 180.0/Pi;
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
                        raypath << solution[m][0] - r_earth;
                        raypath << '\t' << solution[m][1]*180.0/3.14159;
                        raypath << '\t' << solution[m][2]*180.0/3.14159;
                        raypath << '\t' << 20.0*log10(GeoAc_Amplitude(solution,m));
                        raypath << '\t' << attenuation;
                        raypath << '\t' << travel_time << '\n';
                    }
                }
            }
            raypath.close();
                        
            arrival_incl = - asin(c(solution[k][0], solution[k][1], solution[k][2]) / c(r_earth + Source_Loc[2], Source_Loc[0] * Pi / 180.0, Source_Loc[1] * Pi / 180.0) * solution[k][3]) * 180.0 / Pi;
            back_az_dev = (90.0 - atan2(-solution[k][4], -solution[k][5])*180.0/Pi) - Calc_Bearing(Receiver_Loc[0], Receiver_Loc[1], Source_Loc[0], Source_Loc[1]);
            if(back_az_dev > 180.0)  back_az_dev-=360.0;
            if(back_az_dev < -180.0) back_az_dev+=360.0;
                
            
            if(!verbose_output){cout << '\t' << "Eigenray identified:" << '\t' << "theta, phi = " << setprecision(8) << lt << ", " << 90.0 - lp << " degrees." << '\n';
            } else {
                cout << '\t' << '\t' << "Eigenray-" << eigenray_count << ".  " << bnc_cnt << " bounce(s)." << '\n';
                cout << '\t' << '\t' << '\t' << "theta, phi = " << setprecision(8) << lt << ", " << 90.0 - lp << " degrees." << '\n';
                cout << '\t' << '\t' << '\t' << "Travel Time = " << travel_time << " seconds." << '\n';
                cout << '\t' << '\t' << '\t' << "Celerity = " << Calc_GC_Distance(Source_Loc[0], Source_Loc[1],Receiver_Loc[0],Receiver_Loc[1])/travel_time << " km/s." << '\n';
                cout << '\t' << '\t' << '\t' << "Amplitude = " << 20.0*log10(GeoAc_Amplitude(solution,k)) << " dB." << '\n';
                cout << '\t' << '\t' << '\t' << "Atmospheric Attenuation = " << -attenuation << " dB." << '\n';
                cout << '\t' << '\t' << '\t' << "Arrival inclination = " << arrival_incl << " degrees." << '\n';
                cout << '\t' << '\t' << '\t' << "Bearing to source = " << Calc_Bearing(Receiver_Loc[0], Receiver_Loc[1], Source_Loc[0], Source_Loc[1]) << " degrees." << '\n';
                cout << '\t' << '\t' << '\t' << "Back azimuth of arrival = " << 90.0 - atan2(-solution[k][4], -solution[k][5]) * 180.0 / Pi << " degrees." << '\n';
                cout << '\t' << '\t' << '\t' << "Azimuth Deviation = " << back_az_dev  << " degrees." << '\n' << '\n';
            }
            
            results << "Eigenray-" << eigenray_count << ".  " << bnc_cnt << " bounce(s)." << '\n';
            results << '\t' << "theta, phi = " << setprecision(8) << lt << ", " << 90.0 - lp << " degrees." << '\n';
            results << '\t' << "Travel Time = " << travel_time << " seconds." << '\n';
            results << '\t' << "Celerity = " << Calc_GC_Distance(Source_Loc[0], Source_Loc[1],Receiver_Loc[0],Receiver_Loc[1])/travel_time << " km/s." << '\n';
            results << '\t' << "Amplitude (geometric) = " << 20.0*log10(GeoAc_Amplitude(solution,k)) << " dB." << '\n';
            results << '\t' << "Atmospheric attenuation = " << -attenuation << " dB." << '\n';
            results << '\t' << "Arrival inclination = " << arrival_incl << " degrees." << '\n';
            results << '\t' << "Bearing to source = " << Calc_Bearing(Receiver_Loc[0], Receiver_Loc[1], Source_Loc[0], Source_Loc[1]) << " degrees." << '\n';
            results << '\t' << "Back azimuth of arrival = " << 90.0 - atan2(-solution[k][4], -solution[k][5]) * 180.0 / Pi << " degrees." << '\n';
            results << '\t' << "Azimuth deviation = " << back_az_dev  << " degrees." << '\n' << '\n';
            
            eigenray_count++;
            
            break;
        } else if(n > 0 && dr > dr_prev){
            // If the range to the receiver has increased, undo the previous changes to theta and phi,
            // half the step scalar and repeat the step using the new scaled increments
            lt -= dlt * step_scalar;
            lp -= dlp * step_scalar;
            
            step_scalar/=2.0;
            
            // In the case that the resulting step size is near machine precision, exit search
            if(sqrt(dlt*dlt + dlp*dlp) * step_scalar < 1.0e-12){
                if (verbose_output) cout << '\t' << '\t' <<  '\t' << "Step size too small, psuedo-critical ray path likely." << '\n' << '\n';
                break;
            }
            
        } else {
            step_scalar = min(1.0, step_scalar * 1.25);
            
            // Calculate the transformation matrix to obtain dlt and dlp from d_lat and d_lon
            d_lat = Receiver_Loc[0] * Pi/180.0 - lat;
            d_lon = Receiver_Loc[1] * Pi/180.0 - lon;
            
            d_lat_dlt = solution[k][7]  - 1.0 / (r_earth + z_grnd) * solution[k][4] / solution[k][3] * solution[k][6];
            d_lat_dlp = solution[k][13] - 1.0 / (r_earth + z_grnd) * solution[k][4] / solution[k][3] * solution[k][12];
            d_lon_dlt = solution[k][8]  - 1.0 / ((r_earth + z_grnd) * cos(lat)) * solution[k][5] / solution[k][3] * solution[k][6];
            d_lon_dlp = solution[k][14] - 1.0 / ((r_earth + z_grnd) * cos(lat)) * solution[k][5] / solution[k][3] * solution[k][12];
             
            det = d_lat_dlt * d_lon_dlp - d_lat_dlp * d_lon_dlt;
            dlt = (d_lon_dlp * d_lat - d_lat_dlp * d_lon) / det * 180.0/Pi;
            dlp = (-d_lon_dlt * d_lat + d_lat_dlt * d_lon) / det * 180.0/Pi;
            
            // if(verbose_output)  cout << '\t' << '\t' << '\t' << '\t' << "Calculated dlt = " << dlt << ", dlp = " << dlp << "." << '\n';
            
            // Correct the steps if they are too large or have gone outside the acceptable range
            if(dlt >  lt_lim_step) dlt =  lt_lim_step;
            if(dlt < -lt_lim_step) dlt = -lt_lim_step;
            if(dlp >  lp_lim_step) dlp =  lp_lim_step;
            if(dlp < -lp_lim_step) dlp = -lp_lim_step;
           
            // Update the angles, copy the current dr into dr_prev
            lt += dlt * step_scalar;
            lp += dlp * step_scalar;
            dr_prev = dr;
        }
        // Clear the solution array and prepare to trace the new ray path
        GeoAc_ClearSolutionArray(solution,k);
	}
    GeoAc_DeleteSolutionArray(solution, GeoAc_length);
}
#endif /* GEOAC_EIGENRAY_H_ */
