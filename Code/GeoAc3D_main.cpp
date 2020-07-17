#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include <math.h>

#include "GeoAc/GeoAc.Parameters.h"
#include "Atmo/Atmo_State.h"
#include "Atmo/G2S_Spline1D.h"
#include "GeoAc/GeoAc.EquationSets.h"
#include "GeoAc/GeoAc.Solver.h"
#include "GeoAc/GeoAc.Interface.h"
#include "GeoAc/GeoAc.Eigenray.h"

using namespace std;

void GeoAc3D_Usage(){
    cout << '\n';
    cout << "NOTE: GeoAc methods are no longer updated and have been replaced by the infraGA/GeoAc" << '\n';
	cout << "methods available at https://github.com/LANL-Seismoacoustics/infraGA" << '\n' << '\n';

    cout << '\t' << "#############################################" << '\n';
    cout << '\t' << "####               GeoAc3D               ####" << '\n';
    cout << '\t' << "####    Three-Dimensional Ray Tracing    ####" << '\n';
    cout << '\t' << "#### Through a Stratified, Moving Medium ####" << '\n';
    cout << '\t' << "#############################################" << '\n' << '\n';
    
    cout << '\n';
    cout << "Usage: GeoAc3D [option] profile.met [parameters]" << '\n';
    cout << '\t' << '\t' << "Enter only 1 option." << '\n';
    cout << '\t' << '\t' << "Profile.met is expected to contain columns describing {z[km]  T[K]  u(zonal wind)[m/s]  v(meridional wind)[m/s]  density[g/cm^3]  p[mbar]} " << '\n';
    cout << '\t' << '\t' << '\t' << "Profile format can be modified, see manual document for details." << '\n';
    cout << '\t' << '\t' << "Parameter calls are expected using the format: parameter_name=value." << '\n' << '\n';
    
    cout << "Options and parameters are:" << '\n';
    cout << '\t' << "-prop (generate ray paths for propagations at multiple azimuth and inclination angles)" << '\n';
    cout << '\t' << '\t' << "Parameter"  << '\t' << "Units/Options" << '\t' << "Default Value" << '\n';
    cout << '\t' << '\t' << "---------------------------------------------" << '\n';
    cout << '\t' << '\t' << "theta_min"         << '\t' << "degrees"            << '\t' << '\t' << "0.5" << '\n';
    cout << '\t' << '\t' << "theta_max"         << '\t' << "degrees"            << '\t' << '\t' << "45.0" << '\n';
    cout << '\t' << '\t' << "theta_step"        << '\t' << "degrees"            << '\t' << '\t' << "0.5"  << '\n';
    cout << '\t' << '\t' << "phi_min"           << '\t' << '\t' << "degrees"    << '\t' << '\t' << "-90.0" << '\n';
    cout << '\t' << '\t' << "phi_max"           << '\t' << '\t' << "degrees"    << '\t' << '\t' << "-90.0" << '\n';
    cout << '\t' << '\t' << "phi_step"          << '\t' << "degrees"            << '\t' << '\t' << "1.0"  << '\n';
    cout << '\t' << '\t' << "azimuth"           << '\t' << '\t' << "See Manual" << '\t' << "-90.0" << '\n';
    cout << '\t' << '\t' << "bounces"           << '\t' << '\t' << "integer"    << '\t' << '\t' << "2" << '\n';
    cout << '\t' << '\t' << "z_src"             << '\t' << '\t' << "km"         << '\t' << '\t' << "0.0" << '\n';
    cout << '\t' << '\t' << "WriteAtmo"         << '\t' << "True/False"         << '\t' << "False" << '\n';
    cout << '\t' << '\t' << "WriteRays"         << '\t' << "True/False"         << '\t' << "True" << '\n' << '\n';
    
    cout << '\t' << "-interactive (set a fixed source location and select individual ray paths to generate)" << '\n';
    cout << '\t' << '\t' << "Parameter"  << '\t' << "Units/Options" << '\t' << "Default Value" << '\n';
    cout << '\t' << '\t' << "---------------------------------------------" << '\n';
    cout << '\t' << '\t' << "z_src"             << '\t' << '\t' << "km" << '\t' << '\t' << "0.0" << '\n' << '\n';
    
    cout << '\t' << "-eig_search (Search for all eigenrays connecting a source at (0.0, 0.0, z_src) to a receiver " << '\n';
    cout << '\t' << '\t' << '\t' << "at (x_rcvr, y_rcvr, z_grnd) which have inclinations and ground reflections within specified limits)" << '\n';
    cout << '\t' << '\t' << "Parameter"  << '\t' << "Units/Options" << '\t' << "Default Value" << '\n';
    cout << '\t' << '\t' << "---------------------------------------------" << '\n';
    cout << '\t' << '\t' << "theta_min"         << '\t' << "degrees"            << '\t' << '\t' << "0.5" << '\n';
    cout << '\t' << '\t' << "theta_max"         << '\t' << "degrees"            << '\t' << '\t' << "45.0" << '\n';
    cout << '\t' << '\t' << "bnc_min"           << '\t' << '\t' << "integer"    << '\t' << '\t' << "0" << '\n';
    cout << '\t' << '\t' << "bnc_max"           << '\t' << '\t' << "integer"    << '\t' << '\t' << "0" << '\n';
    cout << '\t' << '\t' << "bounces"           << '\t' << '\t' << "See Manual" << '\t' << "0" << '\n';
    cout << '\t' << '\t' << "z_src"             << '\t' << '\t' << "km"         << '\t' << '\t' << "0.0" << '\n';
    cout << '\t' << '\t' << "x_rcvr"            << '\t' << '\t' << "km"         << '\t' << '\t' << "250.0" << '\n';
    cout << '\t' << '\t' << "y_rcvr"            << '\t' << '\t' << "km"         << '\t' << '\t' << "0.0" << '\n';
    cout << '\t' << '\t' << "Verbose"           << '\t' << '\t' << "True/False" << '\t' << "False" << '\n';
    cout << '\t' << '\t' << "iterations"        << '\t' << "integer"            << '\t' << '\t' << "25" << '\n';
    cout << '\t' << '\t' << "azimuth_err_lim"   << '\t' << "degrees"            << '\t' << '\t' << "2.0" << '\n' << '\n';
  
    cout << '\t' << "-eig_direct (Search for a single eigenray connecting a source at (0.0, 0.0, z_src) to a receiver at (x_rcvr, y_rcvr, z_grnd)" << '\n';
    cout << '\t' << '\t' << '\t' << "near an estimated azimuth and inclination pair assuming a specific number of ground reflections)" << '\n';
    cout << '\t' << '\t' << "Parameter"  << '\t' << "Units/Options" << '\t' << "Default Value" << '\n';
    cout << '\t' << '\t' << "---------------------------------------------" << '\n';
    cout << '\t' << '\t' << "theta_est"         << '\t' << "degrees"            << '\t' << '\t' << "15.0" << '\n';
    cout << '\t' << '\t' << "phi_est"           << '\t' <<'\t' << "degrees"     << '\t' << '\t' << "atan2(y_rvcr, x_rcvr)" << '\n';
    cout << '\t' << '\t' << "bounces"           << '\t' << '\t' << "integer"    << '\t' << '\t' << "0" << '\n';
    cout << '\t' << '\t' << "z_src"             << '\t' << '\t' << "km"         << '\t' << '\t' << "0.0" << '\n';
    cout << '\t' << '\t' << "x_rcvr"            << '\t' << '\t' << "km"         << '\t' << '\t' << "250.0" << '\n';
    cout << '\t' << '\t' << "y_rcvr"            << '\t' << '\t' << "km"         << '\t' << '\t' << "0.0" << '\n';
    cout << '\t' << '\t' << "Verbose"           << '\t' << '\t' << "True/False" << '\t' << "False" << '\n';
    cout << '\t' << '\t' << "iterations"        << '\t' << "integer"            << '\t' << '\t' << "25" << '\n' << '\n';
    
    cout << '\t' << "Additional Parameters"  << '\t' << "Units/Options" << '\t' << "Default Value" << '\n';
    cout << '\t' << "---------------------------------------------" << '\n';
    cout << '\t' << "freq"              << '\t' << '\t' << '\t' << "Hz"         << '\t' << '\t' << "0.1" << '\n';
    cout << '\t' << "abs_coeff"         << '\t' << '\t' << "scalar"             << '\t' << '\t' << "0.3" << '\n';
    cout << '\t' << "z_grnd"            << '\t' << '\t' << '\t' << "km"         << '\t' << '\t' << "0.0" << '\n';
    cout << '\t' << "profile_format"    << '\t' << '\t' << "See Manual"         << '\t' << "zTuvdp" << '\n';
    cout << '\t' << "WriteCaustics*"    << '\t' << '\t' << "True/False"         << '\t' << "False" << '\n';
    cout << '\t' << "CalcAmp*"          << '\t' << '\t' << "True/False" << '\t' << "True" << '\n';
    cout << '\t' << "alt_max"           << '\t' << '\t' << '\t' << "km"         << '\t' << '\t' << "Interpolation Max" << '\n';
    cout << '\t' << "rng_max"           << '\t' << '\t' << '\t' << "km"         << '\t' << '\t' << "10,000" << '\n' << '\n';
    cout << '\t' << "* denotes parameters not available in eigenray methods." << '\n' << '\n';

    cout << "Examples:" << '\n';
    cout << '\t' << "./GeoAc3D -prop ToyAtmo.met theta_step=1.0 bounces=2 azimuth=0.0" << '\n';
    cout << '\t' << "./GeoAc3D -prop ToyAtmo.met bounces=0 phi_min=-120.0 phi_max=-60.0 phi_step=3.0 WriteRays=False" << '\n';
    cout << '\t' << "./GeoAc3D -interactive ToyAtmo.met z_grnd=1.5 freq=0.2" << '\n';
    cout << '\t' << "./GeoAc3D -eig_search ToyAtmo.met x_rcvr=300.0 y_rcvr=150.0 verbose=True" << '\n';
    cout << '\t' << "./GeoAc3D -eig_direct ToyAtmo.met x_rcvr=-240.0 y_rcvr=75.0 theta_est=8.0" << '\n' << '\n';
}

void GeoAc3D_RunProp(char* inputs[], int count){
    double theta_min = 0.5, theta_max=45.0, theta_step=0.5;
    double phi_min=-90.0, phi_max=-90.0, phi_step=1.0;
    int bounces=2;
    double x_src=0.0, y_src=0.0, z_src=0.0;
    bool CalcAmp=true, WriteAtmo=false, WriteRays=true, WriteCaustics=false;
    double freq=0.1;
    char* ProfileFormat = "zTuvdp";
    char input_check;
    z_grnd=0.0;
    tweak_abs = 0.3;

    for(int i = 3; i < count; i++) if (strncmp(inputs[i], "profile_format=",15) == 0){ ProfileFormat = inputs[i]+15;}
    Spline_Single_G2S(inputs[2], ProfileFormat);
    
    for(int i = 3; i < count; i++){
        if (strncmp(inputs[i], "theta_min=",10) == 0){          theta_min = atof(inputs[i]+10);}
        else if (strncmp(inputs[i], "theta_max=",10) == 0){     theta_max = atof(inputs[i]+10);}
        else if (strncmp(inputs[i], "theta_step=",11) == 0){    theta_step = atof(inputs[i]+11);}
        else if (strncmp(inputs[i], "phi_min=",8) == 0){        phi_min = atof(inputs[i]+8);}
        else if (strncmp(inputs[i], "phi_max=",8) == 0){        phi_max = atof(inputs[i]+8);}
        else if (strncmp(inputs[i], "phi_step=",9) == 0){       phi_step = atof(inputs[i]+9);}
        else if (strncmp(inputs[i], "azimuth=",8) == 0){
            phi_min = atof(inputs[i]+8);
            phi_max = atof(inputs[i]+8);
            phi_step = 1.0;
        }
        else if (strncmp(inputs[i], "bounces=",8) == 0){        bounces = atoi(inputs[i]+8);}
        else if (strncmp(inputs[i], "z_src=",6) == 0){          z_src = atof(inputs[i]+6);}
        else if (strncmp(inputs[i], "WriteAtmo=",10) == 0){     WriteAtmo = string2bool(inputs[i]+10);}
        else if (strncmp(inputs[i], "WriteRays=",10) == 0){     WriteRays = string2bool(inputs[i]+10);}
        
        else if (strncmp(inputs[i], "freq=",5) == 0){               freq = atof(inputs[i]+5);}
        else if (strncmp(inputs[i], "abs_coeff=",10) == 0){         tweak_abs = max(0.0, atof(inputs[i]+10));}
        else if (strncmp(inputs[i], "z_grnd=",7) == 0){             z_grnd = atof(inputs[i]+7);}
        else if (strncmp(inputs[i], "profile_format=",15) == 0){    ProfileFormat = inputs[i]+15;}
        else if (strncmp(inputs[i], "WriteCaustics=",14) == 0){     WriteCaustics = string2bool(inputs[i]+14);}
        else if (strncmp(inputs[i], "CalcAmp=",8) == 0){            CalcAmp = string2bool(inputs[i]+8);}
        else if (strncmp(inputs[i], "alt_max=",8) == 0){            GeoAc_vert_limit = atof(inputs[i]+8);}
        else if (strncmp(inputs[i], "rng_max=",8) == 0){            GeoAc_range_limit = atof(inputs[i]+8);}
        else{
            cout << "***WARNING*** Unrecognized parameter entry: " << inputs[i] << '\n';
            cout << "Continue? (y/n):"; cin >> input_check;
            if(input_check!='y' && input_check!='Y') return;
        }
    }   
    z_src=max(z_grnd,z_src);
    if(WriteCaustics) CalcAmp=true;
    GeoAc_ConfigureCalcAmp(CalcAmp);
    
    // Extract the file name from the input and use it to distinguish the resulting output
    char file_title[50];
    for(int m = 0; m < 50; m++){
        if(inputs[2][m]=='.'){
            file_title[m]='\0';
            break;
        }
        file_title[m]=inputs[2][m];
    }
    
    // Define variables used for analysis
	double D, D_prev, travel_time_sum, attenuation, z_max, inclination, back_az;
	int k, length = GeoAc_ray_limit * int(1.0/(GeoAc_ds_min*10));
	bool BreakCheck;
	char output_buffer [60];
    
    // Write the profile to file if neceessary
    if(WriteAtmo) GeoAc_WriteProfile("atmo.dat", 90.0 - phi_min);
    
	ofstream results;
    ofstream raypath;
    ofstream caustics;
    
    sprintf(output_buffer, "%s_results.dat", file_title);
    results.open(output_buffer);
    results << "# theta [deg]";
    results << '\t' << "phi [deg]";
    results << '\t' << "n_b";
    results << '\t' << "x_0 [km]";
    results << '\t' << "y_0 [km]";
    results << '\t' << "Travel Time [s]";
    results << '\t' << "Turning Height [km]";
    results << '\t' << "Inclination [deg]";
    results << '\t' << "Back Azimuth [deg]";
    results << '\t' << "Geo. Atten. [dB]";
    results << '\t' << "Atmo. Atten. [dB]";
    results << '\n';
    
	if(WriteRays){
        sprintf(output_buffer, "%s_raypaths.dat", file_title);
        raypath.open(output_buffer);
        
        raypath << "# x [km]";
        raypath << '\t' << "y [km]";
        raypath << '\t' << "z [km]";
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
            
            caustics << "# x [km]";
            caustics << '\t' << "y [km]";
            caustics << '\t' << "z [km]";
            caustics << '\t' << "Travel Time [s]";
            caustics << '\n';
            
            caustics.close();
        }
    }
    
	double** solution;
	GeoAc_BuildSolutionArray(solution,length);
    
    for(double phi = phi_min;       phi <= phi_max;     phi+=phi_step){
    for(double theta = theta_min;   theta <= theta_max; theta+=theta_step){
        GeoAc_theta = theta*Pi/180.0;
        GeoAc_phi = Pi/2.0 - phi*Pi/180.0;
            
        GeoAc_SetInitialConditions(solution, x_src, y_src, z_src);
        travel_time_sum = 0.0;
        attenuation = 0.0;
        z_max = 0.0;
        
        cout << "Plotting ray path w/ theta = " << theta << ", phi = " << phi << '\n';
        for(int bnc_cnt = 0; bnc_cnt <= bounces; bnc_cnt++){
            k = GeoAc_Propagate_RK4(solution, BreakCheck);
                
            if(WriteRays || WriteCaustics){
                if(WriteCaustics){
                    sprintf(output_buffer, "%s_caustics-path%i.dat", file_title, bnc_cnt);
                    caustics.open(output_buffer,fstream::app);

                    D_prev = GeoAc_Jacobian(solution,1);
                }
                
                
                for(int m = 1; m < k ; m++){     // write profiles to data files and vector arrays
                    if(WriteCaustics){ D = GeoAc_Jacobian(solution,m);}
                    GeoAc_TravelTimeSegment(travel_time_sum, solution, m-1,m);
                    GeoAc_SB_AttenSegment(attenuation, solution, m-1, m, freq);
                        
                    if(WriteRays && m % 25 == 0){
                        raypath << solution[m][0];
                        raypath << '\t' << solution[m][1];
                        raypath << '\t' << max(solution[m][2],0.0);
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
                    }
                    if(WriteCaustics){D_prev = D;}
                }
                if(WriteCaustics) caustics.close();
            } else {
                travel_time_sum+= GeoAc_TravelTime(solution, k);
                attenuation+= GeoAc_SB_Atten(solution,k,freq);
            }
            if(BreakCheck) break;
                
            for(int m = 0; m < k ; m++){ z_max = max (z_max, solution[m][2]);}
                
            back_az = phi + 180.0;
            inclination = - asin(c(solution[k][0], solution[k][1], z_grnd) / c(x_src, y_src, z_src) * solution[k][3]) * 180.0 / Pi;
            while(back_az > 180.0)  back_az-=360.0;
            while(back_az < -180.0) back_az+=360.0;
                
            results << theta;
            results << '\t' << phi;
            results << '\t' << bnc_cnt;
            results << '\t' << solution[k][0];
            results << '\t' << solution[k][1];
            results << '\t' << travel_time_sum;
            results << '\t' << z_max;
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

void GeoAc3D_RunInteractive(char* inputs[], int count){
    double x_src=0.0, y_src=0.0, z_src=0.0;
    double freq=0.1, D, D_prev;
    bool CalcAmp=true, WriteCaustics=false;
    char* ProfileFormat = "zTuvdp";
    char input_check;
    z_grnd=0.0;
    tweak_abs = 0.3;
    
    for(int i = 3; i < count; i++) if (strncmp(inputs[i], "profile_format=",15) == 0){ ProfileFormat = inputs[i]+15;}
    Spline_Single_G2S(inputs[2], ProfileFormat);
    
    for(int i = 3; i < count; i++){
        if (strncmp(inputs[i], "z_src=",6) == 0){              z_src = atof(inputs[i]+6);}

        else if (strncmp(inputs[i], "freq=",5) == 0){               freq = atof(inputs[i]+5);}
        else if (strncmp(inputs[i], "abs_coeff=",10) == 0){         tweak_abs = max(0.0, atof(inputs[i]+10));}
        else if (strncmp(inputs[i], "z_grnd=",7) == 0){             z_grnd = atof(inputs[i]+7);}
        else if (strncmp(inputs[i], "profile_format=",15) == 0){    ProfileFormat = inputs[i]+15;}
        else if (strncmp(inputs[i], "WriteCaustics=",14) == 0){     WriteCaustics = string2bool(inputs[i]+14);}
        else if (strncmp(inputs[i], "CalcAmp=",8) == 0){            CalcAmp = string2bool(inputs[i]+8);}
        else if (strncmp(inputs[i], "alt_max=",8) == 0){            GeoAc_vert_limit = atof(inputs[i]+8);}
        else if (strncmp(inputs[i], "rng_max=",8) == 0){            GeoAc_range_limit = atof(inputs[i]+8);}
        else{
            cout << "***WARNING*** Unrecognized parameter entry: " << inputs[i] << '\n';
            cout << "Continue? (y/n):"; cin >> input_check;
            if(input_check!='y' && input_check!='Y') return;
        }
    }
    z_src=max(z_grnd,z_src);
    if(WriteCaustics) CalcAmp=true;
    GeoAc_ConfigureCalcAmp(CalcAmp);
    
    // Define variables used for analysis
    double theta, phi;
    int bounces;
    char keepgoing;
    
    double travel_time_sum, attenuation, z_max;
    int k, length = GeoAc_ray_limit * int(1.0/(GeoAc_ds_min*10));
    bool BreakCheck;
    
    ofstream raypath;
    ofstream caustics;
    
    double** solution;
    GeoAc_BuildSolutionArray(solution,length);
    
    keepgoing='y';
    while(keepgoing=='y' || keepgoing=='Y'){
        raypath.open("raypath.dat");
        
        raypath << "# x [km]";
        raypath << '\t' << "y [km]";
        raypath << '\t' << "z [km]";
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
        cout << '\t' << "Enter number of bounces: ";            cin >> bounces;
        
        cout << '\n';
        cout << '\t' << "Plotting ray path w/ theta = " << theta << ", phi = " << phi << '\n';
        GeoAc_theta = theta*Pi/180.0;
        GeoAc_phi = Pi/2.0 - phi*Pi/180.0;
        
        GeoAc_SetInitialConditions(solution, x_src, y_src, z_src);
        travel_time_sum = 0.0;
        attenuation = 0.0;
        z_max = 0.0;
		
        for(int bnc_cnt = 0; bnc_cnt <= bounces; bnc_cnt++){
            k = GeoAc_Propagate_RK4(solution, BreakCheck);
            if(WriteCaustics){ D_prev = GeoAc_Jacobian(solution,1);}

            for(int m = 1; m < k ; m++){     // write profiles to data files and vector arrays
                GeoAc_TravelTimeSegment(travel_time_sum, solution, m-1,m);
                GeoAc_SB_AttenSegment(attenuation, solution, m-1, m, freq);
                z_max = max (z_max, solution[m][2]);
                if(WriteCaustics){ D = GeoAc_Jacobian(solution,m);}
                
                if(m % 10 == 0){
                    raypath << solution[m][0];
                    raypath << '\t' << solution[m][1];
                    raypath << '\t' << max(solution[m][2],0.0);
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
                }
                if(WriteCaustics){D_prev = D;}
            }
            if(BreakCheck) break;
            GeoAc_SetReflectionConditions(solution,k);
        }
        
        double range = sqrt(pow(solution[k][0],2) + pow(solution[k][1],2));
        
        double back_az = phi - 180.0;
        while(back_az < -180.0) back_az +=360.0;
        while(back_az >  180.0) back_az -=360.0;
        
        if(!BreakCheck){
            cout << '\t' << '\t' << "Ray path arrived at " << solution[k][0] << " km E-W, " << solution[k][1] << " km N-S." << '\n';
            cout << '\t' << '\t' << "Arrival range = " << range << " km at azimuth " << atan2(solution[k][1], solution[k][0])*180.0/Pi << " degrees from N." << '\n';
            if(CalcAmp){ cout << '\t' << '\t' << "Geometric attenuation = " << 20.0*log10(GeoAc_Amplitude(solution,k)) << " dB." << '\n';}
            cout << '\t' << '\t' << "Atmospheric attenuation = " << -attenuation << " dB." << '\n';
            cout << '\t' << '\t' << "Arrival celerity = " << range/travel_time_sum << " km/sec." << '\n';
            cout << '\t' << '\t' << "Turning height of the ray = " << z_max << " km." << '\n';
            cout << '\t' << '\t' << "Back azimuth of the arrival = " << back_az  << " degrees (relative to N). " << '\n' << '\n';
        } else {
            cout << '\t' << '\t' << "Ray path does not return to the ground." << '\n' << '\n';
        }
    
        raypath.close();
        if(WriteCaustics) caustics.close();
        
        cout << "Continue plotting other ray paths? (y/n): ";
        cin >> keepgoing;
    }    
    GeoAc_DeleteSolutionArray(solution, length);
}

void GeoAc3D_RunEigSearch(char* inputs[], int count){
    double Source_Loc [3]   = {0.0, 0.0, 0.0};
    double Receiver_Loc [2] = {-250.0, 0.0};
    double theta_min = 0.5, theta_max = 45.0;
    int bnc_min = 0, bnc_max = 0;
    int iterations=25;
    double azimuth_err_lim=2.0;
    verbose_output=false;
    double freq = 0.1;
    char* ProfileFormat = "zTuvdp";
    char input_check;
    z_grnd=0.0;
    tweak_abs = 0.3;
    
    for(int i = 3; i < count; i++) if (strncmp(inputs[i], "profile_format=",15) == 0){ ProfileFormat = inputs[i]+15;}
    Spline_Single_G2S(inputs[2], ProfileFormat);
    
    for(int i = 3; i < count; i++){
        if (strncmp(inputs[i], "theta_min=",10) == 0){              theta_min = atof(inputs[i]+10);}
        else if (strncmp(inputs[i], "theta_max=",10) == 0){         theta_max = atof(inputs[i]+10);}
        else if (strncmp(inputs[i], "bnc_min=",8) == 0){            bnc_min = atoi(inputs[i]+8);}
        else if (strncmp(inputs[i], "bnc_max=",8) == 0){            bnc_max = atoi(inputs[i]+8);}
        else if (strncmp(inputs[i], "bounces=",8) == 0){
            bnc_min = atoi(inputs[i]+8);
            bnc_max = atoi(inputs[i]+8);
        }
        else if (strncmp(inputs[i], "z_src=",6) == 0){              Source_Loc[2] = atof(inputs[i]+6);}
        else if (strncmp(inputs[i], "x_rcvr=",7) == 0){             Receiver_Loc[0] = atof(inputs[i]+7);}
        else if (strncmp(inputs[i], "y_rcvr=",7) == 0){             Receiver_Loc[1] = atof(inputs[i]+7);}
        
        else if (strncmp(inputs[i], "Verbose=",8) == 0){            verbose_output = string2bool(inputs[i]+8);}
        else if (strncmp(inputs[i], "verbose=",8) == 0){            verbose_output = string2bool(inputs[i]+8);}
        else if (strncmp(inputs[i], "azimuth_err_lim=",16) == 0){   azimuth_err_lim=atof(inputs[i]+16);}
        else if (strncmp(inputs[i], "iterations=",11) == 0){        iterations=atof(inputs[i]+11);}

        else if (strncmp(inputs[i], "freq=",5) == 0){               freq = atof(inputs[i]+5);}
        else if (strncmp(inputs[i], "abs_coeff=",10) == 0){         tweak_abs = max(0.0, atof(inputs[i]+10));}
        else if (strncmp(inputs[i], "z_grnd=",7) == 0){             z_grnd = atof(inputs[i]+7);}
        else if (strncmp(inputs[i], "profile_format=",15) == 0){    ProfileFormat = inputs[i]+15;}
        else if (strncmp(inputs[i], "alt_max=",8) == 0){            GeoAc_vert_limit = atof(inputs[i]+8);}
        else if (strncmp(inputs[i], "rng_max=",8) == 0){            GeoAc_range_limit = atof(inputs[i]+8);}
        else{
            cout << "***WARNING*** Unrecognized parameter entry: " << inputs[i] << '\n';
            cout << "Continue? (y/n):"; cin >> input_check;
            if(input_check!='y' && input_check!='Y') return;
        }
    }
    Source_Loc[2]=max(z_grnd,Source_Loc[2]);
    
    // Extract the file name from the input and use it to distinguish the output
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
    
    results << "GeoAc3D - Eigenray Run Summary:" << '\n';
    results << '\t' << "Profile used: " << inputs[2] << '\n';
    results << '\t' << "Source Location (kilometers) : (" << Source_Loc[0] << ", " << Source_Loc[1] << ", " <<  Source_Loc[2] << ")." << '\n';
    results << '\t' << "Receiver Location (kilometers) : (" << Receiver_Loc[0] << ", " << Receiver_Loc[1] << ", " << z_grnd << ")." << '\n';
    results << '\t' << "Inclination range (degrees): " << theta_min << " - " << theta_max << "." << '\n';
    results << '\t' << "Ground reflection (bounce) limits: " << bnc_min << " - " << bnc_max << "." << '\n' << '\n';
    
    for(int n_bnc = bnc_min; n_bnc <= bnc_max; n_bnc++){
        cout << "Searching for " << n_bnc << " bounce eigenray(s) between " << theta_min << " and " << theta_max << "." << '\n';
        
        theta_start = theta_min;
        while(theta_start < theta_max){
            estimate_success = GeoAc_EstimateEigenray(Source_Loc, Receiver_Loc, theta_start, theta_max, theta_est, phi_est, theta_next, n_bnc, azimuth_err_lim);
            if(estimate_success)  GeoAc_3DEigenray_LM(Source_Loc, Receiver_Loc, theta_est, phi_est, freq, n_bnc, iterations, file_title);
            
            theta_start = theta_next;
        }
    }
    results.close();
    cout << '\t' << "Identified " << eigenray_count << " eigenray(s)." << '\n';
}

void GeoAc3D_RunEigDirect(char* inputs[], int count){
    double Source_Loc [3]   = {0.0, 0.0, 0.0};
    double Receiver_Loc [2] = {-250.0, 0.0};
    double theta_est = 0.5, phi_est=45.0;
    int bounces = 0, iterations=25;
    verbose_output=false;
    double freq = 0.1;
    char* ProfileFormat = "zTuvdp";
    char input_check;
    z_grnd=0.0;
    tweak_abs=0.3;
    
    for(int i = 3; i < count; i++) if (strncmp(inputs[i], "profile_format=",15) == 0){ ProfileFormat = inputs[i]+15;}
    Spline_Single_G2S(inputs[2], ProfileFormat);
    
    for(int i = 3; i < count; i++){
        if (strncmp(inputs[i], "theta_est=",10) == 0){              theta_est = atof(inputs[i]+10);}
        else if (strncmp(inputs[i], "bounces=",8) == 0){            bounces = atoi(inputs[i]+8);}
        else if (strncmp(inputs[i], "z_src=",6) == 0){              Source_Loc[2] = atof(inputs[i]+6);}
        else if (strncmp(inputs[i], "x_rcvr=",7) == 0){             Receiver_Loc[0] = atof(inputs[i]+7);}
        else if (strncmp(inputs[i], "y_rcvr=",7) == 0){             Receiver_Loc[1] = atof(inputs[i]+7);}
        else if (strncmp(inputs[i], "Verbose=",8) == 0){            verbose_output = string2bool(inputs[i]+8);}
        else if (strncmp(inputs[i], "verbose=",8) == 0){            verbose_output = string2bool(inputs[i]+8);}
        else if (strncmp(inputs[i], "iterations=",11) == 0){        iterations=atof(inputs[i]+11);}
        
        else if (strncmp(inputs[i], "freq=",5) == 0){               freq = atof(inputs[i]+5);}
        else if (strncmp(inputs[i], "abs_coeff=",10) == 0){         tweak_abs = max(0.0, atof(inputs[i]+10));}
        else if (strncmp(inputs[i], "z_grnd=",7) == 0){             z_grnd = atof(inputs[i]+7);}
        else if (strncmp(inputs[i], "profile_format=",15) == 0){    ProfileFormat = inputs[i]+15;}
        else if (strncmp(inputs[i], "alt_max=",8) == 0){            GeoAc_vert_limit = atof(inputs[i]+8);}
        else if (strncmp(inputs[i], "rng_max=",8) == 0){            GeoAc_range_limit = atof(inputs[i]+8);}
        else{
            cout << "***WARNING*** Unrecognized parameter entry: " << inputs[i] << '\n';
            cout << "Continue? (y/n):"; cin >> input_check;
            if(input_check!='y' && input_check!='Y') return;
        }
    }
    Source_Loc[2]=max(z_grnd,Source_Loc[2]);
    
    // Set the phi_start angle by the azimuth to the receiver and then change it if it's been input
    phi_est = 180.0/3.14159 * atan2(Receiver_Loc[1], Receiver_Loc[0]);
    for(int i = 3; i < count; i++){
        if (strncmp(inputs[i], "phi_est=",8) == 0){     phi_est = 90.0 - atof(inputs[i]+8);}
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
    if(argc < 3){
        GeoAc3D_Usage();
        return 0;
    }
    
    if (strncmp(argv[1], "-prop",5) == 0){
        GeoAc3D_RunProp(argv, argc);
    
    } else if (strncmp(argv[1], "-interactive",12) == 0){
        GeoAc3D_RunInteractive(argv, argc);
        
    } else if (strncmp(argv[1], "-eig_search",11) == 0){
        GeoAc3D_RunEigSearch(argv, argc);
        
    } else if (strncmp(argv[1], "-eig_direct",11) == 0){
        GeoAc3D_RunEigDirect(argv, argc);
        
    } else {
        cout << "Unrecognized option." << '\n';
    }
    
    return 0;
}
