#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>

#include "Atmo/G2S_GlobalMultiDimSpline3D.h"
#include "Atmo/Atmo_State.h"

#include "GeoAc/GeoAc.Parameters.h"
#include "GeoAc/GeoAc.EquationSets.h"
#include "GeoAc/GeoAc.Solver.h"
#include "GeoAc/GeoAc.Interface.h"
#include "GeoAc/GeoAc.Eigenray.h"

using namespace std;


void GeoAcGlobal_RngDep_Usage(){
    cout << '\n';
    cout << '\t' << "#################################################" << '\n';
    cout << '\t' << "####           GeoAcGlobal.RngDep            ####" << '\n';
    cout << '\t' << "####      Three-Dimensional Ray Tracing      ####" << '\n';
    cout << '\t' << "####       Using Spherical Coordinates       ####" << '\n';
    cout << '\t' << "#### Through an Inhomogeneous, Moving Medium ####" << '\n';
    cout << '\t' << "#################################################" << '\n' << '\n';
    
    cout << '\n';
    cout << "Usage: GeoAcGlobal.RngDep [option] profile_prefix nodes-lat.loc nodes-lon.loc [parameters]" << '\n';
    cout << '\t' << '\t' << "Enter only 1 option." << '\n';
    cout << '\t' << '\t' << "Each profile_prefix##.met is expected to contain columns describing {z[km]  T[K]  u(zonal wind)[m/s]  v(meridional wind)[m/s]  density[g/cm^3]  p[mbar]} " << '\n';
    cout << '\t' << '\t' << "The files nodes-lat.loc and nodes-lon.loc are expected to contain columns describing the latitude and longitude locations of the .met files." << '\n';
    cout << '\t' << '\t' << '\t' << "These can be generated using the included Python script for a single .bin file." << '\n';
    cout << '\t' << '\t' << '\t' << "Profiles are ordered such that profile_prefixN.met describes the atmosphere at lat = lat_i, lon = lon_j, N = i + n_lat * j" << '\n';
    cout << '\t' << '\t' << '\t' << "Profile format can be modified, see manual document for details." << '\n';
    cout << '\t' << '\t' << "Parameter calls are expected using the format: parameter_name=value." << '\n' << '\n';
    
    cout << "Options and parameters are:" << '\n';
    cout << '\t' << "-prop (generate ray paths for propagations at multiple azimuth and inclination angles)" << '\n';
    cout << '\t' << '\t' << "Parameter"  << '\t' << "Units" << '\t' << '\t' << "Default Value" << '\n';
    cout << '\t' << '\t' << "---------------------------------------------" << '\n';
    cout << '\t' << '\t' << "theta_min"         << '\t' << "degrees"            << '\t' << '\t' << "0.5" << '\n';
    cout << '\t' << '\t' << "theta_max"         << '\t' << "degrees"            << '\t' << '\t' << "45.0" << '\n';
    cout << '\t' << '\t' << "theta_step"        << '\t' << "degrees"            << '\t' << '\t' << "0.5"  << '\n';
    cout << '\t' << '\t' << "phi_min"           << '\t' << '\t' << "degrees"    << '\t' << '\t' << "-90.0" << '\n';
    cout << '\t' << '\t' << "phi_max"           << '\t' << '\t' << "degrees"    << '\t' << '\t' << "-90.0" << '\n';
    cout << '\t' << '\t' << "phi_step"          << '\t' << "degrees"            << '\t' << '\t' << "1.0"  << '\n';
    cout << '\t' << '\t' << "azimuth"           << '\t' << '\t' << "See Manual" << '\t' << "-90.0" << '\n';
    cout << '\t' << '\t' << "bounces"           << '\t' << '\t' << "integer"    << '\t' << '\t' << "2" << '\n';
    cout << '\t' << '\t' << "lat_src"           << '\t' << '\t' << "degrees"    << '\t' << '\t' << "30.0" << '\n';
    cout << '\t' << '\t' << "lon_src"           << '\t' << '\t' << "degrees"    << '\t' << '\t' << "0.0" << '\n';
    cout << '\t' << '\t' << "z_src"             << '\t' << '\t' << "km"         << '\t' << '\t' << "0.0" << '\n';
    cout << '\t' << '\t' << "WriteAtmo"         << '\t' << "True/False"         << '\t' << "False" << '\n';
    cout << '\t' << '\t' << "WriteRays"         << '\t' << "True/False"         << '\t' << "True" << '\n' << '\n';
    
    cout << '\t' << "-interactive (set a fixed source location and select individual ray paths to generate)" << '\n';
    cout << '\t' << '\t' << "Parameter"  << '\t' << "Units" << '\t' << '\t' << "Default Value" << '\n';
    cout << '\t' << '\t' << "---------------------------------------------" << '\n';
    cout << '\t' << '\t' << "lat_src"           << '\t' << '\t' << "degrees"    << '\t' << '\t' << "30.0" << '\n';
    cout << '\t' << '\t' << "lon_src"           << '\t' << '\t' << "degrees"    << '\t' << '\t' << "0.0" << '\n';
    cout << '\t' << '\t' << "z_src"             << '\t' << '\t' << "km"         << '\t' << '\t' << "0.0" << '\n' << '\n';
    
    cout << '\t' << "-eig_search (Search for all eigenrays connecting a source at (lat_src, lon_src, z_src) to a receiver " << '\n';
    cout << '\t' << '\t' << '\t' << "at (lat_rcvr, lon_rcvr, z_grnd) which have inclinations and ground reflections within specified limits)" << '\n';
    cout << '\t' << '\t' << "Parameter"  << '\t' << "Units" << '\t' << '\t' << "Default Value" << '\n';
    cout << '\t' << '\t' << "---------------------------------------------" << '\n';
    cout << '\t' << '\t' << "theta_min"         << '\t' << "degrees"            << '\t' << '\t' << "1.0" << '\n';
    cout << '\t' << '\t' << "theta_max"         << '\t' << "degrees"            << '\t' << '\t' << "45.0" << '\n';
    cout << '\t' << '\t' << "bnc_min"           << '\t' << '\t' << "integer"    << '\t' << '\t' << "0" << '\n';
    cout << '\t' << '\t' << "bnc_max"           << '\t' << '\t' << "integer"    << '\t' << '\t' << "0" << '\n';
    cout << '\t' << '\t' << "bounces"           << '\t' << '\t' << "See Manual" << '\t' << "0" << '\n';
    cout << '\t' << '\t' << "lat_src"           << '\t' << '\t' << "degrees"    << '\t' << '\t' << "30.0" << '\n';
    cout << '\t' << '\t' << "lon_src"           << '\t' << '\t' << "degrees"    << '\t' << '\t' << "0.0" << '\n';
    cout << '\t' << '\t' << "z_src"             << '\t' << '\t' << "km"         << '\t' << '\t' << "0.0" << '\n';
    cout << '\t' << '\t' << "lat_rcvr"          << '\t' << "degrees"            << '\t' << '\t' << "30" << '\n';
    cout << '\t' << '\t' << "lon_rcvr"          << '\t' << "degrees"            << '\t' << '\t' << "-2.5" << '\n';
    cout << '\t' << '\t' << "Verbose"           << '\t' << '\t' << "True/False" << '\t' << "False" << '\n';
    cout << '\t' << '\t' << "iterations"        << '\t' << "integer"            << '\t' << '\t' << "25" << '\n';
    cout << '\t' << '\t' << "azimuth_err_lim"   << '\t' << "degrees"            << '\t' << '\t' << "2.0" << '\n' << '\n';
    
    cout << '\t' << "-eig_direct (Search for a single eigenray connecting a source at (lat_src, lon_src, z_src) to a receiver " << '\n';
    cout << '\t' << '\t' << '\t' << "at (lat_rcvr, lon_rcvr, z_grnd) near an estiamted azimuth and inclination assuming a specific number of ground reflections)" << '\n';
    cout << '\t' << '\t' << "Parameter"  << '\t' << "Units" << '\t' << '\t' << "Default Value" << '\n';
    cout << '\t' << '\t' << "---------------------------------------------" << '\n';
    cout << '\t' << '\t' << "theta_est"         << '\t' << "degrees"            << '\t' << '\t' << "15.0" << '\n';
    cout << '\t' << '\t' << "phi_est"           << '\t' << '\t' << "degrees"    << '\t' << '\t' << "Great Circle Bearing From Source to Receiver" << '\n';
    cout << '\t' << '\t' << "bounces"           << '\t' << '\t' << "integer"    << '\t' << '\t' << "2" << '\n';
    cout << '\t' << '\t' << "lat_src"           << '\t' << '\t' << "degrees"    << '\t' << '\t' << "30.0" << '\n';
    cout << '\t' << '\t' << "lon_src"           << '\t' << '\t' << "degrees"    << '\t' << '\t' << "0.0" << '\n';
    cout << '\t' << '\t' << "z_src"             << '\t' << '\t' << "km"         << '\t' << '\t' << "0.0" << '\n';
    cout << '\t' << '\t' << "lat_rcvr"          << '\t' << "degrees"            << '\t' << '\t' << "30" << '\n';
    cout << '\t' << '\t' << "lon_rcvr"          << '\t' << "degrees"            << '\t' << '\t' << "10.0" << '\n';
    cout << '\t' << '\t' << "Verbose"           << '\t' << '\t' << "True/False" << '\t' << "False" << '\n';
    cout << '\t' << '\t' << "iterations"        << '\t' << "integer"            << '\t' << '\t' << "25" << '\n' << '\n';
    
    cout << '\t' << "Additional Parameters"  << '\t' << "Units/Options" << '\t' << "Default Value" << '\n';
    cout << '\t' << "---------------------------------------------" << '\n';
    cout << '\t' << "freq"              << '\t' << '\t' << '\t' << "Hz"         << '\t' << '\t' << "0.1" << '\n';
    cout << '\t' << "abs_coeff"         << '\t' << '\t' << "scalar"             << '\t' << '\t' << "0.3" << '\n';
    cout << '\t' << "z_grnd"            << '\t' << '\t' << '\t' << "km"         << '\t' << '\t' << "0.0" << '\n';
    cout << '\t' << "profile_format"    << '\t' << '\t' << "See Manual"         << '\t' << "zTuvdp" << '\n';
    cout << '\t' << "WriteCaustics*"    << '\t' << '\t' << "True/False"         << '\t' << "False" << '\n';
    cout << '\t' << "CalcAmp*"          << '\t' << '\t' << "True/False"         << '\t' << "True" << '\n';
    cout << '\t' << "alt_max"           << '\t' << '\t' << '\t' << "km"         << '\t' << "Interpolation Limits" << '\n';
    cout << '\t' << "lat_min"           << '\t' << '\t' << '\t' << "deg"        << '\t' << "Interpolation Limits" << '\n';
    cout << '\t' << "lat_max"           << '\t' << '\t' << '\t' << "deg"        << '\t' << "Interpolation Limits" << '\n';
    cout << '\t' << "lon_min"           << '\t' << '\t' << '\t' << "deg"        << '\t' << "Interpolation Limits" << '\n';
    cout << '\t' << "lon_max"           << '\t' << '\t' << '\t' << "deg"        << '\t' << "Interpolation Limits" << '\n';
    cout << '\t' << "* denotes parameters not available in eigenray methods." << '\n' << '\n';
    
    cout << "Examples (Note - Profiles are not provided, these examples will not run without user specified profiles):" << '\n';
    cout << '\t' << "./GeoAcGlobal.RngDep -prop Profiles/Profile Profiles/loc_lat.dat Profiles/loc_lon.dat theta_step=2.0 bounces=2 azimuth=-45.0" << '\n';
    cout << '\t' << "./GeoAcGlobal.RngDep -eig_search Profiles/Profile Profiles/loc_lat.dat Profiles/loc_lon.dat lat_src=45.0 lon_src=-110.0 lat_rcvr=42.5 y_rcvr=-112.5 verbose=True" << '\n' << '\n';
    
}

void GeoAcGlobal_RngDep_RunProp(char* inputs[], int count){
    double theta_min = 0.5, theta_max=45.0, theta_step=0.5;
    double phi_min=-90.0, phi_max=-90.0, phi_step=1.0;
    int bounces=2;
    bool CalcAmp=true, WriteAtmo=false, WriteRays=true, WriteCaustics=false;
    double freq=0.1;
    char* ProfileFormat = "zTuvdp";
    char input_check;
    z_grnd=0.0;
    
    for(int i = 5; i < count; i++) if (strncmp(inputs[i], "profile_format=",15) == 0){    ProfileFormat = inputs[i]+15;}
    Spline_Multi_G2S(inputs[2], inputs[3], inputs[4], ProfileFormat);
    
    double  lat_src=(t_vals[0] + t_vals[t_cnt-1])/2.0 * 180.0/Pi,
            lon_src=(p_vals[0] + p_vals[p_cnt-1])/2.0 * 180.0/Pi,
            z_src=0.0;
    
    for(int i = 5; i < count; i++){
        if (strncmp(inputs[i], "theta_min=",10) == 0){              theta_min = atof(inputs[i]+10);}
        else if (strncmp(inputs[i], "theta_max=",10) == 0){         theta_max = atof(inputs[i]+10);}
        else if (strncmp(inputs[i], "theta_step=",11) == 0){        theta_step = atof(inputs[i]+11);}
        else if (strncmp(inputs[i], "phi_min=",8) == 0){            phi_min = atof(inputs[i]+8);}
        else if (strncmp(inputs[i], "phi_max=",8) == 0){            phi_max = atof(inputs[i]+8);}
        else if (strncmp(inputs[i], "phi_step=",9) == 0){           phi_step = atof(inputs[i]+9);}
        else if (strncmp(inputs[i], "azimuth=",8) == 0){
            phi_min = atof(inputs[i]+8);
            phi_max = atof(inputs[i]+8);
            phi_step = 1.0;
        }
        else if (strncmp(inputs[i], "bounces=",8) == 0){            bounces = atoi(inputs[i]+8);}
        else if (strncmp(inputs[i], "lat_src=",8) == 0){            lat_src = atof(inputs[i]+8);}
        else if (strncmp(inputs[i], "lon_src=",8) == 0){            lon_src = atof(inputs[i]+8);}
        else if (strncmp(inputs[i], "z_src=",6) == 0){              z_src = atof(inputs[i]+6);}
        else if (strncmp(inputs[i], "WriteAtmo=",10) == 0){         WriteAtmo = string2bool(inputs[i]+10);}
        else if (strncmp(inputs[i], "WriteRays=",10) == 0){         WriteRays = string2bool(inputs[i]+10);}

        else if (strncmp(inputs[i], "freq=",5) == 0){               freq = atof(inputs[i]+5);}
        else if (strncmp(inputs[i], "abs_coeff=",10) == 0){         tweak_abs = max(0.0, atof(inputs[i]+10));}
        else if (strncmp(inputs[i], "z_grnd=",7) == 0){             z_grnd = atof(inputs[i]+7);}
        else if (strncmp(inputs[i], "profile_format=",15) == 0){    ProfileFormat = inputs[i]+15;}
        else if (strncmp(inputs[i], "WriteCaustics=",14) == 0){     WriteCaustics = string2bool(inputs[i]+14);}
        else if (strncmp(inputs[i], "CalcAmp=",8) == 0){            CalcAmp = string2bool(inputs[i]+8);}
        
        else if (strncmp(inputs[i], "alt_max=",8) == 0){            GeoAc_vert_limit = atof(inputs[i]+8);}
        else if (strncmp(inputs[i], "lat_min=",8) == 0){            GeoAc_lat_min_limit = atof(inputs[i]+8);}
        else if (strncmp(inputs[i], "lat_max=",8) == 0){            GeoAc_lat_max_limit = atof(inputs[i]+8);}
        else if (strncmp(inputs[i], "lon_min=",8) == 0){            GeoAc_lon_min_limit = atof(inputs[i]+8);}
        else if (strncmp(inputs[i], "lon_max=",8) == 0){            GeoAc_lon_max_limit = atof(inputs[i]+8);}
        
        else{
            cout << "***WARNING*** Unrecognized parameter entry: " << inputs[i] << '\n';
            cout << "Continue? (y/n):"; cin >> input_check;
            if(input_check!='y' && input_check!='Y') return;
        }
    }
    z_src = max(z_src,z_grnd);
    if(WriteCaustics) CalcAmp=true;
    GeoAc_ConfigureCalcAmp(CalcAmp);
    
    char file_title[50];
    for(int m = 0; m < 50; m++){
        if(inputs[2][m]=='.'){
            file_title[m]='\0';
            break;
        }
        file_title[m]=inputs[2][m];
    }
    
    // Define variables used for analysis
	double D, D_prev, travel_time_sum, attenuation, r_max;
	int k, length = GeoAc_ray_limit * int(1.0/(GeoAc_ds_min*10));
	bool BreakCheck;
	char output_buffer [60];
	
    // Write the profile to file if neceessary
    if(WriteAtmo) GeoAc_WriteProfile("atmo.dat", lat_src*Pi/180.0, lon_src*Pi/180.0, 90.0 - phi_min);
        
	ofstream results;
    ofstream raypath;
    ofstream caustics;
    
    sprintf(output_buffer, "%s_results.dat", file_title);
    results.open(output_buffer);
    results << "# theta [deg]";
    results << '\t' << "phi [deg]";
    results << '\t' << "Bounces";
    results << '\t' << "lat_0 [deg]";
    results << '\t' << "lon_0 [deg]";
    results << '\t' << "Travel Time [s]";
    results << '\t' << "Celerity [km/s]";
    results << '\t' << "Turning Height [km]";
    results << '\t' << "Inclination [deg]";
    results << '\t' << "Back Azimuth [deg]";
    results << '\t' << "Geo. Atten. [dB]";
    results << '\t' << "Atmo. Atten. [dB]";
    results << '\n';
    
	if(WriteRays){
        sprintf(output_buffer, "%s_raypaths.dat", file_title);
        raypath.open(output_buffer);
        
        raypath << "# z [km]";
        raypath << '\t' << "Lat [deg]";
        raypath << '\t' << "Long [deg]";
        raypath << '\t' << "Geo. Atten. [dB]";
        raypath << '\t' << "Atmo. Atten. [dB]";
        raypath << '\t' << "Travel Time [s]";
        raypath << '\n';
        
    }
    if(WriteCaustics){
        for (int bnc = 0; bnc <= bounces; bnc++){
            // Clear caustics files
            sprintf(output_buffer, "%s_caustics-path%i.dat", file_title, bnc);
            caustics.open(output_buffer);
            
            caustics << "# z [km]";
            caustics << '\t' << "Lat [deg]";
            caustics << '\t' << "Long [deg]";
            caustics << '\t' << "Travel Time [s]";
            caustics << '\n';
            
            caustics.close();
        }
    }
    
	double** solution;
	GeoAc_BuildSolutionArray(solution,length);
    
    for(double phi = phi_min;       phi <= phi_max;     phi+=phi_step){
        for(double theta = theta_min;   theta <= theta_max; theta+=theta_step){
            cout << "Plotting ray path w/ theta = " << theta << ", phi = " << phi << '\n';
            GeoAc_theta = theta*Pi/180.0;
            GeoAc_phi = Pi/2.0 - phi*Pi/180.0;
            
            GeoAc_SetInitialConditions(solution, z_src, lat_src*Pi/180.0, lon_src*Pi/180.0);
            travel_time_sum = 0.0;
            attenuation = 0.0;
            
            for(int bnc_cnt = 0; bnc_cnt <= bounces; bnc_cnt++){
                
                k = GeoAc_Propagate_RK4(solution, BreakCheck);
                
                if(WriteRays || WriteCaustics){
                    if(WriteCaustics){
                        sprintf(output_buffer, "%s_caustics-path%i.dat", file_title, bnc_cnt);
                        caustics.open(output_buffer,fstream::app);
                        
                        D_prev = GeoAc_Jacobian(solution,1);
                    }
                    
                    for(int m = 1; m < k ; m++){     // write profiles to data files and vector arrays
                        if(WriteCaustics) D = GeoAc_Jacobian(solution,m);
                        GeoAc_TravelTimeSegment(travel_time_sum, solution, m-1,m);
                        GeoAc_SB_AttenSegment(attenuation, solution, m-1, m, freq);
                        
                        if(WriteRays && m % 25 == 0){
                            raypath << solution[m][0] - r_earth;
                            raypath << '\t' << setprecision(8) << solution[m][1] * 180.0/Pi;
                            raypath << '\t' << setprecision(8) << solution[m][2] * 180.0/Pi;
                            if(CalcAmp){    raypath << '\t' << 20.0*log10(GeoAc_Amplitude(solution,m));}
                            else{           raypath << '\t' << 0.0;}
                            raypath << '\t' << -attenuation;
                            raypath << '\t' << travel_time_sum << '\n';
                        }
                        
                        if(WriteCaustics && D*D_prev < 0.0){
                            caustics << solution[m][0] - r_earth;
                            caustics << '\t' << setprecision(8) << solution[m][1] * 180.0/Pi;
                            caustics << '\t' << setprecision(8) << solution[m][2] * 180.0/Pi;
                            caustics << '\t' << travel_time_sum << '\n';
                        }
                        if(WriteCaustics) D_prev = D;
                    }
                    if(WriteCaustics) caustics.close();
                } else {
                    travel_time_sum+= GeoAc_TravelTime(solution, k);
                    attenuation+= GeoAc_SB_Atten(solution,k,freq);
                }
                
                if(BreakCheck) break;
                
                r_max = 0.0;
                for(int m = 0; m < k ; m++) r_max = max (r_max, solution[m][0] - r_earth);
                
                double GC_Dist1 = pow(sin((solution[k][1] - lat_src*Pi/180.0)/2.0),2);
                double GC_Dist2 = cos(lat_src*Pi/180.0) * cos(solution[k][1]) * pow(sin((solution[k][2] - lon_src*Pi/180.0)/2.0),2);
                
                double inclination = asin(c(solution[k][0], solution[k][1], solution[k][2]) / c(r_earth + z_src, lat_src*Pi/180.0, lon_src*Pi/180.0) * solution[k][3]) * 180.0 / Pi;
                double back_az = 90.0 - atan2(-solution[k][4], -solution[k][5]) * 180.0 / Pi;
                if(back_az < -180.0) back_az +=360.0;
                if(back_az >  180.0) back_az -=360.0;
                
                results << theta;
                results << '\t' << phi;
                results << '\t' << bnc_cnt;
                results << '\t' << setprecision(8) << solution[k][1] * 180.0/Pi;
                results << '\t' << setprecision(8) << solution[k][2] * 180.0/Pi;
                results << '\t' << travel_time_sum;
                results << '\t' << 2.0 * r_earth * asin(sqrt(GC_Dist1+GC_Dist2)) / travel_time_sum;
                results << '\t' << r_max;
                results << '\t' << inclination;
                results << '\t' << back_az;
                if(CalcAmp){    results << '\t' << 20.0*log10(GeoAc_Amplitude(solution,k));}
                else{           results << '\t' << 0.0;}
                results << '\t' << -attenuation;
                results << '\n';
                
                GeoAc_SetReflectionConditions(solution,k);
            }
            if(WriteRays){raypath << '\n';}
            GeoAc_ClearSolutionArray(solution,k);
        }
        results << '\n';
    }
	
	raypath.close();
	results.close();
	caustics.close();
	GeoAc_DeleteSolutionArray(solution, length);
    
}

void GeoAcGlobal_RngDep_RunInteractive(char* inputs[], int count){
    
    char* ProfileFormat = "zTuvdp";
    char input_check;
    z_grnd=0.0;
    for(int i = 5; i < count; i++){
        if (strncmp(inputs[i], "profile_format=",15) == 0){     ProfileFormat = inputs[i]+15;}
        else if (strncmp(inputs[i], "z_grnd=",7) == 0){         z_grnd = atof(inputs[i]+7);}
    }
    
    Spline_Multi_G2S(inputs[2], inputs[3], inputs[4], ProfileFormat);
    GeoAc_SetPropRegion();
    GeoAc_ConfigureCalcAmp(true);
    
    double  lat_src=(t_vals[0] + t_vals[t_cnt-1])/2.0 * 180.0/Pi,
    lon_src=(p_vals[0] + p_vals[p_cnt-1])/2.0 * 180.0/Pi,
    z_src=0.0;
    double freq=0.1, D, D_prev;
    bool CalcAmp=true, WriteCaustics=false;
    
    for(int i = 5; i < count; i++){
        if (strncmp(inputs[i], "lat_src=",8) == 0){                 lat_src = atof(inputs[i]+8);}
        else if (strncmp(inputs[i], "lon_src=",8) == 0){            lon_src = atof(inputs[i]+8);}
        else if (strncmp(inputs[i], "z_src=",6) == 0){              z_src = atof(inputs[i]+6);}
        
        else if (strncmp(inputs[i], "freq=",5) == 0){               freq = atof(inputs[i]+5);}
        else if (strncmp(inputs[i], "abs_coeff=",10) == 0){         tweak_abs = max(0.0, atof(inputs[i]+10));}
        else if (strncmp(inputs[i], "z_grnd=",7) == 0){             z_grnd = atof(inputs[i]+7);}
        else if (strncmp(inputs[i], "profile_format=",15) == 0){    ProfileFormat = inputs[i]+15;}
        else if (strncmp(inputs[i], "WriteCaustics=",14) == 0){     WriteCaustics = string2bool(inputs[i]+14);}
        else if (strncmp(inputs[i], "CalcAmp=",8) == 0){            CalcAmp = string2bool(inputs[i]+8);}
        
        else if (strncmp(inputs[i], "alt_max=",8) == 0){            GeoAc_vert_limit = atof(inputs[i]+8);}
        else if (strncmp(inputs[i], "lat_min=",8) == 0){            GeoAc_lat_min_limit = atof(inputs[i]+8);}
        else if (strncmp(inputs[i], "lat_max=",8) == 0){            GeoAc_lat_max_limit = atof(inputs[i]+8);}
        else if (strncmp(inputs[i], "lon_min=",8) == 0){            GeoAc_lon_min_limit = atof(inputs[i]+8);}
        else if (strncmp(inputs[i], "lon_max=",8) == 0){            GeoAc_lon_max_limit = atof(inputs[i]+8);}
        
        else{
            cout << "***WARNING*** Unrecognized parameter entry: " << inputs[i] << '\n';
            cout << "Continue? (y/n):"; cin >> input_check;
            if(input_check!='y' && input_check!='Y') return;
        }
    }
    z_src=max(z_src,z_grnd);
    if(WriteCaustics) CalcAmp=true;
    
    char file_title[50];
    for(int m = 0; m < 50; m++){
        if(inputs[2][m]=='.'){
            file_title[m]='\0';
            break;
        }
        file_title[m]=inputs[2][m];
    }
    double theta, phi;
    int bnc_cnt;
    char keepgoing;
    
    // Define variables used for analysis
    double travel_time_sum, attenuation, r_max;
    int k, length = GeoAc_ray_limit * int(1.0/(GeoAc_ds_min*10));
    bool BreakCheck;
    
    ofstream raypath;
    ofstream caustics;
    
    double** solution;
    GeoAc_BuildSolutionArray(solution,length);
    
    keepgoing='y';
    while(keepgoing=='y' || keepgoing=='Y'){
        raypath.open("raypath.dat");
        
        raypath << "# z [km]";
        raypath << '\t' << "Lat [deg]";
        raypath << '\t' << "Long [deg]";
        raypath << '\t' << "Geo. Atten. [dB]";
        raypath << '\t' << "Atmo. Atten. [dB]";
        raypath << '\t' << "Travel Time [s]";
        raypath << '\n';
        
        if(WriteCaustics){
            caustics.open("caustics.dat");
            
            caustics << "# x [km]";
            caustics << '\t' << "y [km]";
            caustics << '\t' << "z [km]";
            caustics << '\t' << "Travel Time [s]";
            caustics << '\n';
        }
        
        cout << '\t' << "Enter inclination angle [degrees]: ";  cin >> theta;
        cout << '\t' << "Enter azimuth angle [degrees]: ";      cin >> phi;
        cout << '\t' << "Enter number of bounces: ";            cin >> bnc_cnt;
        
        cout << '\n';
        cout << '\t' << "Plotting ray path w/ theta = " << theta << ", phi = " << phi << '\n';
        GeoAc_theta = theta*Pi/180.0;
        GeoAc_phi = Pi/2.0 - phi*Pi/180.0;
        
        GeoAc_SetInitialConditions(solution, z_src, lat_src*Pi/180.0, lon_src*Pi/180.0);
        travel_time_sum = 0.0;
        attenuation = 0.0;
        r_max = 0.0;
		
        for(int bnc = 0; bnc <= bnc_cnt; bnc++){
            k = GeoAc_Propagate_RK4(solution, BreakCheck);
            if(WriteCaustics){ D_prev = GeoAc_Jacobian(solution,1);}
            
            for(int m = 1; m < k ; m++){     // write profiles to data files and vector arrays
                GeoAc_TravelTimeSegment(travel_time_sum, solution, m-1,m);
                GeoAc_SB_AttenSegment(attenuation, solution, m-1, m, freq);
                r_max = max(r_max, solution[m][0]-r_earth);
                if(WriteCaustics){ D = GeoAc_Jacobian(solution,m);}
                
                if(m % 10 == 0){
                    raypath << solution[m][0] - r_earth;
                    raypath << '\t' << setprecision(8) << solution[m][1] * 180.0/Pi;
                    raypath << '\t' << setprecision(8) << solution[m][2] * 180.0/Pi;
                    if(CalcAmp){    raypath << '\t' << 20.0*log10(GeoAc_Amplitude(solution,m));}
                    else{           raypath << '\t' << 0.0;}
                    raypath << '\t' << -attenuation;
                    raypath << '\t' << travel_time_sum;
                    raypath << '\n';
                }
                if(WriteCaustics && D*D_prev < 0.0){
                    caustics << solution[m][0];
                    caustics << '\t' << solution[m][1];
                    caustics << '\t' << max(solution[m][2],0.0);
                    caustics <<'\t' << travel_time_sum << '\n';
                    caustics << '\n';
                }
                if(WriteCaustics){D_prev = D;}
            }
            if(BreakCheck) break;
            GeoAc_SetReflectionConditions(solution,k);
        }
        
        double GC_Dist1 = pow(sin((solution[k][1] - lat_src*Pi/180.0)/2.0),2);
        double GC_Dist2 = cos(lat_src*Pi/180.0) * cos(solution[k][1]) * pow(sin((solution[k][2] - lon_src*Pi/180.0)/2.0),2);
        double range = 2.0 * r_earth * asin(sqrt(GC_Dist1+GC_Dist2));
        
        double arg1 = sin((solution[k][2] - lon_src*Pi/180.0));
        double arg2 = cos(lat_src*Pi/180.0)*tan(solution[k][1]) - sin(lat_src*Pi/180.0)*cos((solution[k][2]-lon_src*Pi/180.0));
        double bearing = atan2(arg1, arg2)*180.0/Pi;
        
        double  back_az = 90.0 - atan2(-solution[k][4], -solution[k][5]) * 180.0 / Pi;
        if(back_az < -180.0) back_az +=360.0;
        if(back_az >  180.0) back_az -=360.0;
        
        if(!BreakCheck){
            cout << '\t' << '\t' << "Ray path arrived at " << solution[k][1]*180.0/Pi << " degrees latitude, " << solution[k][2]*180.0/Pi << " degrees longitude." << '\n';
            cout << '\t' << '\t' << "Arrival range = " << range << " km at azimuth " << bearing << " degrees from N." << '\n';
            if(CalcAmp) cout << '\t' << '\t' << "Geometric attenuation = " << 20.0*log10(GeoAc_Amplitude(solution,k)) << " dB." << '\n';
            cout << '\t' << '\t' << "Atmospheric attenuation = " << -attenuation << " dB." << '\n';
            cout << '\t' << '\t' << "Arrival celerity = " << range/travel_time_sum << "km/sec." << '\n';
            cout << '\t' << '\t' << "Back azimuth of the arrival = " << back_az  << " degrees from N. " << '\n' << '\n';
        } else {
            cout << '\t' << '\t' << "Ray path does not return to the ground " << '\n' << '\n';
        }
        
        raypath.close();
        if(WriteCaustics) caustics.close();
            
        cout << "Continue plotting other ray paths? (y/n): ";
        cin >> keepgoing;
    }
    GeoAc_DeleteSolutionArray(solution, length);

    
}

void GeoAcGlobal_RngDep_RunEigSearch(char* inputs[], int count){
    double theta_min = 0.5, theta_max = 45.0;
    int bnc_min = 0, bnc_max = 0;
    int iterations=25;
    double azimuth_err_lim=2.0;
    verbose_output=false;
    double freq = 0.1;
    char* ProfileFormat = "zTuvdp";
    char input_check;
    z_grnd=0.0;

    for(int i = 5; i < count; i++){
        if (strncmp(inputs[i], "profile_format=",15) == 0){     ProfileFormat = inputs[i]+15;}
        else if (strncmp(inputs[i], "z_grnd=",7) == 0){         z_grnd = atof(inputs[i]+7);}
    }
    
    Spline_Multi_G2S(inputs[2], inputs[3], inputs[4], ProfileFormat);
    GeoAc_SetPropRegion();
    
    double Source_Loc [3]   = { (t_vals[0] + t_vals[t_cnt-1])/2.0 * 180.0/Pi,
                                (p_vals[0] + p_vals[p_cnt-1])/2.0 * 180.0/Pi,
                                0.0};

    double Receiver_Loc [2] = { (t_vals[0] + t_vals[t_cnt-1]) * 0.75 * 180.0/Pi,
                                (p_vals[0] + p_vals[p_cnt-1]) * 0.75 * 180.0/Pi};
    
    for(int i = 5; i < count; i++){
        if (strncmp(inputs[i], "theta_min=",10) == 0){              theta_min = atof(inputs[i]+10);}
        else if (strncmp(inputs[i], "theta_max=",10) == 0){         theta_max = atof(inputs[i]+10);}
        else if (strncmp(inputs[i], "bnc_min=",8) == 0){            bnc_min = atoi(inputs[i]+8);}
        else if (strncmp(inputs[i], "bnc_max=",8) == 0){            bnc_max = atoi(inputs[i]+8);}
        else if (strncmp(inputs[i], "bounces=",8) == 0){
            bnc_min = atoi(inputs[i]+8);
            bnc_max = atoi(inputs[i]+8);
        }
        else if (strncmp(inputs[i], "lat_src=",8) == 0){            Source_Loc[0] = atof(inputs[i]+8);}
        else if (strncmp(inputs[i], "lon_src=",8) == 0){            Source_Loc[1] = atof(inputs[i]+8);}
        else if (strncmp(inputs[i], "z_src=",6) == 0){              Source_Loc[2] = atof(inputs[i]+6);}
        else if (strncmp(inputs[i], "lat_rcvr=",9) == 0){           Receiver_Loc[0] = atof(inputs[i]+9);}
        else if (strncmp(inputs[i], "lon_rcvr=",9) == 0){           Receiver_Loc[1] = atof(inputs[i]+9);}
        else if (strncmp(inputs[i], "Verbose=",8) == 0){            verbose_output = string2bool(inputs[i]+8);}
        else if (strncmp(inputs[i], "verbose=",8) == 0){            verbose_output = string2bool(inputs[i]+8);}
        else if (strncmp(inputs[i], "azimuth_err_lim=",16) == 0){   azimuth_err_lim=atof(inputs[i]+16);}
        else if (strncmp(inputs[i], "iterations=",11) == 0){        iterations=atof(inputs[i]+11);}
        
        else if (strncmp(inputs[i], "freq=",5) == 0){               freq = atof(inputs[i]+5);}
        else if (strncmp(inputs[i], "abs_coeff=",10) == 0){         tweak_abs = max(0.0, atof(inputs[i]+10));}
        else if (strncmp(inputs[i], "z_grnd=",7) == 0){             z_grnd = atof(inputs[i]+7);}
        else if (strncmp(inputs[i], "profile_format=",15) == 0){    ProfileFormat = inputs[i]+15;}
        
        else if (strncmp(inputs[i], "alt_max=",8) == 0){            GeoAc_vert_limit = atof(inputs[i]+8);}
        else if (strncmp(inputs[i], "lat_min=",8) == 0){            GeoAc_lat_min_limit = atof(inputs[i]+8);}
        else if (strncmp(inputs[i], "lat_max=",8) == 0){            GeoAc_lat_max_limit = atof(inputs[i]+8);}
        else if (strncmp(inputs[i], "lon_min=",8) == 0){            GeoAc_lon_min_limit = atof(inputs[i]+8);}
        else if (strncmp(inputs[i], "lon_max=",8) == 0){            GeoAc_lon_max_limit = atof(inputs[i]+8);}
        
        else{
            cout << "***WARNING*** Unrecognized parameter entry: " << inputs[i] << '\n';
            cout << "Continue? (y/n):"; cin >> input_check;
            if(input_check!='y' && input_check!='Y') return;
        }
    }
    Source_Loc[2]=max(Source_Loc[2],z_grnd);
    
    // Extract the file name from the input and use it to distinguish the resulting output
    char output_buffer[60];
    char file_title[50];
    for(int m = 0; m < 50; m++){
        if(inputs[2][m]=='.'){
            file_title[m]='\0';
            break;
        }
        file_title[m]=inputs[2][m];
    }
    
    double theta_start, theta_next, theta_est, phi_est;
    bool estimate_success;
    
    sprintf(output_buffer, "%s_results.dat", file_title);
    results.open(output_buffer);
    
    results << "GeoAcGlobal.RngDep - Eigenray Run Summary:" << '\n';
    results << '\t' << "Profile used: " << inputs[2] << '\n';
    results << '\t' << "Source Location (lat, lon, elev) : (" << Source_Loc[0] << ", " << Source_Loc[1] << ", " <<  Source_Loc[2] << ")." << '\n';
    results << '\t' << "Receiver Location (lat, lon, elev) : (" << Receiver_Loc[0] << ", " << Receiver_Loc[1] << ", " << z_grnd << ")." << '\n';
    results << '\t' << "Inclination range (degrees): " << theta_min << " - " << theta_max << "." << '\n';
    results << '\t' << "Ground reflection (bounce) limits: " << bnc_min << " - " << bnc_max << "." << '\n' << '\n';
    
    for(int n_bnc = bnc_min; n_bnc <= bnc_max; n_bnc++){
        cout << "Searching for " << n_bnc << " bounce eigenrays." << '\n';
        
        theta_start = theta_min;
        while(theta_start < theta_max){
            estimate_success = GeoAc_EstimateEigenray(Source_Loc, Receiver_Loc, theta_start, theta_max, theta_est, phi_est, theta_next, n_bnc, azimuth_err_lim);
            if(estimate_success)  GeoAc_3DEigenray_LM(Source_Loc, Receiver_Loc, theta_est, phi_est, freq, n_bnc, iterations, file_title);
            
            theta_start = theta_next;
        }
    }
    
    results.close();
    cout << "Identified " << eigenray_count << " eigenray(s)." << '\n';
    
}

void GeoAcGlobal_RngDep_RunEigDirect(char* inputs[], int count){
    
    double theta_est = 10.0, phi_est=45.0;
    int bounces = 0, iterations=25;
    verbose_output=false;
    double freq = 0.1;
    char* ProfileFormat = "zTuvdp";
    char input_check;
    z_grnd=0.0;
    
    for(int i = 5; i < count; i++){
        if (strncmp(inputs[i], "profile_format=",15) == 0){     ProfileFormat = inputs[i]+15;}
        else if (strncmp(inputs[i], "z_grnd=",7) == 0){         z_grnd = atof(inputs[i]+7);}
    }
    Spline_Multi_G2S(inputs[2], inputs[3], inputs[4], ProfileFormat);
    GeoAc_SetPropRegion();
    
    double Source_Loc [3]   = { (t_vals[0] + t_vals[t_cnt-1])/2.0 * 180.0/Pi,
                                (p_vals[0] + p_vals[p_cnt-1])/2.0 * 180.0/Pi,
                                0.0};
    
    double Receiver_Loc [2] = { (t_vals[0] + t_vals[t_cnt-1]) * 0.75 * 180.0/Pi,
                                (p_vals[0] + p_vals[p_cnt-1]) * 0.75 * 180.0/Pi};

    
    for(int i = 5; i < count; i++){
        if (strncmp(inputs[i], "theta_est=",10) == 0){              theta_est = atof(inputs[i]+10);}
        else if (strncmp(inputs[i], "bounces=",8) == 0){            bounces = atoi(inputs[i]+8);}
        else if (strncmp(inputs[i], "lat_src=",8) == 0){            Source_Loc[0] = atof(inputs[i]+8);}
        else if (strncmp(inputs[i], "lon_src=",8) == 0){            Source_Loc[1] = atof(inputs[i]+8);}
        else if (strncmp(inputs[i], "z_src=",6) == 0){              Source_Loc[2] = atof(inputs[i]+6);}
        else if (strncmp(inputs[i], "lat_rcvr=",9) == 0){           Receiver_Loc[0] = atof(inputs[i]+9);}
        else if (strncmp(inputs[i], "lon_rcvr=",9) == 0){           Receiver_Loc[1] = atof(inputs[i]+9);}
        else if (strncmp(inputs[i], "Verbose=",8) == 0){            verbose_output = string2bool(inputs[i]+8);}
        else if (strncmp(inputs[i], "verbose=",8) == 0){            verbose_output = string2bool(inputs[i]+8);}
        else if (strncmp(inputs[i], "iterations=",11) == 0){        iterations=atof(inputs[i]+11);}

        else if (strncmp(inputs[i], "freq=",5) == 0){               freq = atof(inputs[i]+5);}
        else if (strncmp(inputs[i], "abs_coeff=",10) == 0){         tweak_abs = max(0.0, atof(inputs[i]+10));}
        else if (strncmp(inputs[i], "z_grnd=",7) == 0){             z_grnd = atof(inputs[i]+7);}
        else if (strncmp(inputs[i], "profile_format=",15) == 0){    ProfileFormat = inputs[i]+15;}
        
        else if (strncmp(inputs[i], "alt_max=",8) == 0){            GeoAc_vert_limit = atof(inputs[i]+8);}
        else if (strncmp(inputs[i], "lat_min=",8) == 0){            GeoAc_lat_min_limit = atof(inputs[i]+8);}
        else if (strncmp(inputs[i], "lat_max=",8) == 0){            GeoAc_lat_max_limit = atof(inputs[i]+8);}
        else if (strncmp(inputs[i], "lon_min=",8) == 0){            GeoAc_lon_min_limit = atof(inputs[i]+8);}
        else if (strncmp(inputs[i], "lon_max=",8) == 0){            GeoAc_lon_max_limit = atof(inputs[i]+8);}
        
        else{
            cout << "***WARNING*** Unrecognized parameter entry: " << inputs[i] << '\n';
            cout << "Continue? (y/n):"; cin >> input_check;
            if(input_check!='y' && input_check!='Y') return;
        }
    }
    Source_Loc[2] = max(Source_Loc[2], z_grnd);
    
    // Set the phi_start angle by the great circle bearing to the receiver and then change it if it's been input
    phi_est=90.0 - Calc_Bearing(Source_Loc[0], Source_Loc[1], Receiver_Loc[0], Receiver_Loc[1]);
    for(int i = 3; i < count; i++){
        if (strncmp(inputs[i], "phi_est=",8) == 0){                 phi_est = 90.0 - atof(inputs[i]+8);}
    }
    
    // Extract the file name from the input and use it to distinguish the resulting output
    char file_title[50];
    for(int m = 0; m < 50; m++){
        if(inputs[2][m]=='.'){
            file_title[m]='\0';
            break;
        }
        file_title[m]=inputs[2][m];
    }
    
    GeoAc_3DEigenray_LM(Source_Loc, Receiver_Loc, theta_est, phi_est, freq, bounces, iterations, file_title);
}



int main(int argc, char* argv[]){
    if(argc < 5){
        GeoAcGlobal_RngDep_Usage();
        return 0;
    }
    
    if (strncmp(argv[1], "-prop",5) == 0){
        GeoAcGlobal_RngDep_RunProp(argv, argc);
    
    } else if (strncmp(argv[1], "-interactive",12) == 0){
        GeoAcGlobal_RngDep_RunInteractive(argv, argc);
        
    } else if (strncmp(argv[1], "-eig_search",11) == 0){
        GeoAcGlobal_RngDep_RunEigSearch(argv, argc);
        
    } else if (strncmp(argv[1], "-eig_direct",11) == 0){
        GeoAcGlobal_RngDep_RunEigDirect(argv, argc);
        
    } else {
        cout << "Unrecognized option." << '\n';
        return 0;
    }
    
    
    
    return 0;
}
