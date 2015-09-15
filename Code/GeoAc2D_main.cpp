#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>

#include "GeoAc/GeoAc.Parameters.h"
#include "Atmo/Atmo_State.h"
#include "Atmo/G2S_Spline1D.h"
#include "GeoAc/GeoAc.EquationSets.h"
#include "GeoAc/GeoAc.Solver.h"
#include "GeoAc/GeoAc.Interface.h"

using namespace std;

void GeoAc2D_Usage(){
    cout << '\n';
    cout << '\t' << "#######################################################" << '\n';
    cout << '\t' << "####                    GeoAc2D                    ####" << '\n';
    cout << '\t' << "####          Two-Dimensional Ray Tracing          ####" << '\n';
    cout << '\t' << "#### Using the Effective Sound Speed Approximation ####" << '\n';
    cout << '\t' << "#######################################################" << '\n' << '\n';
    
    cout << "Usage: GeoAc2D [option] profile.met [parameters]" << '\n';
    cout << '\t' << '\t' << "Enter only 1 option." << '\n';
    cout << '\t' << '\t' << "Profile.met is expected to contain columns describing {z[km]  T[K]  u(zonal winds)[m/s]  v(meridional winds)[m/s]  density[g/cm^3]  p[mbar]} " << '\n';
    cout << '\t' << '\t' << '\t' << "Profile format can be modified, see manual document for details." << '\n';
    cout << '\t' << '\t' << "Parameter calls are expected using the format: parameter_name=value." << '\n' << '\n';
    
    cout << "Options and parameters are:" << '\n';
    cout << '\t' << "-prop (generate ray paths for propagations at fixed azimuth using multiple inclination angles)" << '\n';
    cout << '\t' << '\t' << "Parameter"  << '\t' << "Units/Options" << '\t' << "Default Value" << '\n';
    cout << '\t' << '\t' << "---------------------------------------------" << '\n';
    cout << '\t' << '\t' << "theta_min"         << '\t' << "degrees"            << '\t'  << '\t' << "0.5" << '\n';
    cout << '\t' << '\t' << "theta_max"         << '\t' << "degrees"            << '\t'  << '\t' << "45.0" << '\n';
    cout << '\t' << '\t' << "theta_step"        << '\t' << "degrees"            << '\t'  << '\t' << "0.5"  << '\n';
    cout << '\t' << '\t' << "azimuth"           << '\t' << '\t' << "degrees"    << '\t'  << '\t' << "-90.0"  << '\n';
    cout << '\t' << '\t' << "bounces"           << '\t' << '\t' << "integer"    << '\t'  << '\t' << "2" << '\n';
    cout << '\t' << '\t' << "z_src"             << '\t' << '\t' << "km"         << '\t'  << '\t' << "0.0" << '\n' << '\n';

    cout << '\t' << "-interactive (set a fixed source elevation and select individual ray paths to generate)" << '\n';
    cout << '\t' << '\t' << "Parameter"  << '\t' << "Units/Options" << '\t' << "Default Value" << '\n';
    cout << '\t' << '\t' << "---------------------------------------------" << '\n';
    cout << '\t' << '\t' << "z_src"             << '\t' << '\t' << "km"         << '\t'  << '\t' << "0.0" << '\n' << '\n';

    cout << '\t' << "Additional Parameters"  << '\t' << "Units/Options" << '\t' << "Default Value" << '\n';
    cout << '\t' << "---------------------------------------------" << '\n';
    cout << '\t' << "freq"              << '\t' << '\t' << '\t' << "Hz"         << '\t' << '\t' << "0.1" << '\n';
    cout << '\t' << "abs_coeff"         << '\t' << '\t' << "scalar"             << '\t' << '\t' << "0.3" << '\n';
    cout << '\t' << "z_grnd"            << '\t' << '\t' << '\t' << "km"         << '\t' << '\t' << "0.0" << '\n';
    cout << '\t' << "profile_format"    << '\t' << '\t' << "See Manual"         << '\t' << "zTuvdp" << '\n';
    cout << '\t' << "WriteCaustics"     << '\t' << '\t' << "True/False"         << '\t' << "False" << '\n';
    cout << '\t' << "CalcAmp"           << '\t' << '\t' << '\t' << "True/False" << '\t' << "True" << '\n';
    cout << '\t' << "alt_max"           << '\t' << '\t' << '\t' << "km"         << '\t' << '\t' << "Interpolation Max" << '\n';
    cout << '\t' << "rng_max"           << '\t' << '\t' << '\t' << "km"         << '\t' << '\t' << "10,000" << '\n' << '\n';
    
    cout << "Examples:" << '\n';
    cout << '\t' << "./GeoAc2D -prop ToyAtmo.met theta_step=1.0 bounces=0" << '\n';
    cout << '\t' << "./GeoAc2D -interactive ToyAtmo.met z_grnd=1.5 freq=0.2" << '\n' << '\n';
}

void GeoAc2D_RunProp(char* inputs[], int count){
    double theta_min = 0.5, theta_max=45.0, theta_step=0.5;
    double azimuth=-90.0;
    int bounces=2;
    double z_src = 0.0;
    bool CalcAmp=true, WriteCaustics=false;
    double freq=0.1;
    char* ProfileFormat = "zTuvdp";
    char input_check;
    z_grnd = 0.0;
    tweak_abs = 0.3;
    
    for(int i = 3; i < count; i++) if (strncmp(inputs[i], "profile_format=",15) == 0){ ProfileFormat = inputs[i]+15;}
    Spline_Single_G2S(inputs[2], ProfileFormat);    // Load profile into spline
    
    for(int i = 3; i < count; i++){
        if (strncmp(inputs[i], "theta_min=",10) == 0){              theta_min = atof(inputs[i]+10);}
        else if (strncmp(inputs[i], "theta_max=",10) == 0){         theta_max = atof(inputs[i]+10);}
        else if (strncmp(inputs[i], "theta_step=",11) == 0){        theta_step = atof(inputs[i]+11);}
        else if (strncmp(inputs[i], "azimuth=",8) == 0){            azimuth = atof(inputs[i]+8);}
        else if (strncmp(inputs[i], "bounces=",8) == 0){            bounces = atoi(inputs[i]+8);}
        else if (strncmp(inputs[i], "z_src=",6) == 0){              z_src = atof(inputs[i]+6);}
        
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
    z_src = max(z_src, z_grnd);
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
	double D, D_prev, travel_time_sum, attenuation, z_max;
	int k, length = GeoAc_ray_limit * int(1.0/(GeoAc_ds_min*10));
	bool BreakCheck;
	char output_buffer [60];
	
    GeoAc_WriteProfile("atmo.dat", 90.0 - azimuth);
    
	ofstream results;
    ofstream raypath;
    ofstream caustics;
    
    sprintf(output_buffer, "%s_results.dat", file_title);
    results.open(output_buffer);
    results << "# theta [deg]";
    results << '\t' << "phi [deg]";
    results << '\t' << "n_b";
    results << '\t' << "r_0 [km]";
    results << '\t' << "Travel Time [s]";
    results << '\t' << "Turning Height [km]";
    results << '\t' << "Inclination [deg]";
    results << '\t' << "Geo. Atten. [dB]";
    results << '\t' << "Atmo. Atten. [dB]";
    results << '\n';
    
    sprintf(output_buffer, "%s_raypaths.dat", file_title);
    raypath.open(output_buffer);
    
    raypath << "# r [km]";
    raypath << '\t' << "z [km]";
    raypath << '\t' << "Geo. Atten. [dB]";
    raypath << '\t' << "Atmo. Atten. [dB]";
    raypath << '\t' << "Travel Time [s]";
    raypath << '\n';
    
    if(WriteCaustics){
        for (int bnc = 0; bnc <= bounces; bnc++){
            sprintf(output_buffer, "%s_caustics-path%i.dat", file_title, bnc);
            caustics.open(output_buffer);
        
            caustics << "# r [km]";
            caustics << '\t' << "z [km]";
            caustics << '\t' << "Travel Time [s]";
            caustics << '\n';
        
            caustics.close();
        }
    }
    
	double** solution;
	GeoAc_BuildSolutionArray(solution,length);
    
	for(double theta = theta_min; theta <= theta_max; theta+=theta_step){
		cout << "Plotting ray path w/ theta = " << theta << ", phi = " << azimuth << '\n';
		GeoAc_theta = theta*Pi/180.0;
        GeoAc_phi = Pi/2.0 - azimuth*Pi/180.0;
        
		GeoAc_SetInitialConditions(solution, 0.0, z_src);
        travel_time_sum = 0.0;
        attenuation = 0.0;
        z_max = 0.0;
		
		for(int bnc_cnt = 0; bnc_cnt <= bounces; bnc_cnt++){
            
			k = GeoAc_Propagate_RK4(solution, BreakCheck);
            if(WriteCaustics) {
                sprintf(output_buffer, "%s_caustics-path%i.dat", file_title, bnc_cnt);
                caustics.open(output_buffer,fstream::app);
            }
            
            if(WriteCaustics) D_prev = GeoAc_Jacobian(solution,1);
            
            for(int m = 1; m < k ; m++){
                GeoAc_TravelTimeSegment(travel_time_sum, solution, m-1,m);
                GeoAc_SB_AttenSegment(attenuation, solution, m-1, m, freq);
                if(WriteCaustics) D = GeoAc_Jacobian(solution,m);
                
   				if(m % 25 == 0){
                    raypath << solution[m][0];
                    raypath << '\t' << max(solution[m][1],0.0);
                    if(CalcAmp){    raypath << '\t' << 20.0*log10(GeoAc_Amplitude(solution,m));}
                    else{           raypath << '\t' << 0.0;}
                    raypath << '\t' << -attenuation;
                    raypath << '\t' << travel_time_sum;
                    raypath << '\n';
                }
                if(WriteCaustics && D*D_prev < 0.0){
                    caustics << solution[m][0];
                    caustics << '\t' << solution[m][1];
                    caustics << '\t' << travel_time_sum << '\n';
                }
                if(WriteCaustics) D_prev = D;
            }
            
            if(WriteCaustics) caustics.close();
    		if(BreakCheck) break;
            for(int m = 0; m < k ; m++){ z_max = max (z_max, solution[m][1]);}
            
            results << theta;
			results << '\t' << azimuth;
            results << '\t' << bnc_cnt;
			results << '\t' << solution[k][0];
			results << '\t' << travel_time_sum;
            results << '\t' << z_max;
            results << '\t' << -theta;
            if(CalcAmp){    results << '\t' << 20.0*log10(GeoAc_Amplitude(solution,k));}
            else{           results << '\t' << 0.0;}
            results << '\t' << -attenuation;
			results << '\n';
            
			GeoAc_SetReflectionConditions(solution,k);
		}
		GeoAc_ClearSolutionArray(solution,k);
        raypath << '\n';
	}
	
	raypath.close();
	results.close();
	GeoAc_DeleteSolutionArray(solution, length);
}

void GeoAc2D_RunInteractive(char* inputs[], int count){
    double z_src=0.0, freq=0.1, D, D_prev;
    bool CalcAmp=true, WriteCaustics=false;
    char* ProfileFormat = "zTuvdp";
    char input_check;
    z_grnd = 0.0;
    tweak_abs = 0.3;
    
    for(int i = 3; i < count; i++) if (strncmp(inputs[i], "profile_format=",15) == 0){ ProfileFormat = inputs[i]+15;}
    Spline_Single_G2S(inputs[2], ProfileFormat);
    
    for(int i = 3; i < count; i++){
        if (strncmp(inputs[i], "z_src=",6) == 0){                   z_src = atof(inputs[i]+6);}

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
    if(z_src < z_grnd) z_src=z_grnd;
    if(WriteCaustics) CalcAmp=true;
    GeoAc_ConfigureCalcAmp(CalcAmp);
    
    // Define variables used for analysis
    double theta, azimuth;
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
        
        raypath << "# r [km]";
        raypath << '\t' << "z [km]";
        raypath << '\t' << "Geo. Atten. [dB]";
        raypath << '\t' << "Atmo. Atten. [dB]";
        raypath << '\t' << "Travel Time [s]";
        raypath << '\n';
        
        if(WriteCaustics){
            caustics.open("caustics.dat");
                
            caustics << "# r [km]";
            caustics << '\t' << "z [km]";
            caustics << '\t' << "Travel Time [s]";
            caustics << '\n';
        }
                
        cout << '\t' << "Enter inclination angle [degrees]: ";  cin >> theta;
        cout << '\t' << "Enter azimuth angle [degrees]: ";      cin >> azimuth;
        cout << '\t' << "Enter number of bounces: ";            cin >> bounces;
        cout << '\n';

        cout << '\t' << "Plotting ray path w/ theta = " << theta << ", phi = " << azimuth << '\n';
        GeoAc_theta = theta*Pi/180.0;
        GeoAc_phi = Pi/2.0 - azimuth*Pi/180.0;
        
        GeoAc_SetInitialConditions(solution, 0.0, z_src);
        travel_time_sum = 0.0;
        attenuation = 0.0;
        z_max = 0.0;
		
        for(int bnc_cnt = 0; bnc_cnt <= bounces; bnc_cnt++){
            k = GeoAc_Propagate_RK4(solution, BreakCheck);
            if(WriteCaustics) D_prev = GeoAc_Jacobian(solution,1);
            
            for(int m = 1; m < k ; m++){     // write profiles to data files and vector arrays
                GeoAc_TravelTimeSegment(travel_time_sum, solution, m-1,m);
                GeoAc_SB_AttenSegment(attenuation, solution, m-1, m, freq);
                z_max = max(z_max, solution[m][2]);
                if(WriteCaustics) D = GeoAc_Jacobian(solution,m);

                if(m % 25 == 0){
                    raypath << solution[m][0];
                    raypath << '\t' << max(solution[m][1],0.0);
                    if(CalcAmp){    raypath << '\t' << 20.0*log10(GeoAc_Amplitude(solution,m));}
                    else {          raypath << '\t' << 0.0;}
                    raypath << '\t' << -attenuation;
                    raypath << '\t' << travel_time_sum;
                    raypath << '\n';
                }
                if(WriteCaustics && D * D_prev < 0.0){
                    caustics << solution[m][0];
                    caustics << '\t' << solution[m][1];
                    caustics << '\t' << travel_time_sum << '\n';
                }
                if(WriteCaustics) D_prev = D;
            }
            if(BreakCheck) break;
            GeoAc_SetReflectionConditions(solution,k);
        }
    
        if(!BreakCheck){
            cout << '\t' << '\t' << "Arrival Range = " << solution[k][0] << '\n';
            if(CalcAmp){    cout << '\t' << '\t' << "Geometric Attenuation = " << 20.0*log10(GeoAc_Amplitude(solution,k)) << " dB." << '\n';}
            cout << '\t' << '\t' << "Atmospheric Attenuation = " << -attenuation << " dB." << '\n';
            cout << '\t' << '\t' << "Turning Height = " << z_max << "km." << '\n';
            cout << '\t' << '\t' << "Travel Time = " << travel_time_sum << "sec." << '\n';
            cout << '\t' << '\t' << "Arrival Celerity = " << solution[k][0]/travel_time_sum << "km/sec." << '\n' << '\n';
        } else { cout << '\t' << '\t' << "Ray path does not return to the ground." << '\n' << '\n';}
        
        raypath.close();
        if(WriteCaustics) caustics.close();

        cout << "Continue plotting other ray paths? (y/n): ";
        cin >> keepgoing;
    }
    
    GeoAc_DeleteSolutionArray(solution, length);
}


int main(int argc, char* argv[]){
    int skip = 3;
    
    if(argc < 3){
        GeoAc2D_Usage();
        return 0;
    }
    
    if (strncmp(argv[1], "-prop",5) == 0){
        GeoAc2D_RunProp(argv, argc);
    } else if (strncmp(argv[1], "-interactive",5) == 0){
        GeoAc2D_RunInteractive(argv, argc);
    } else {
        cout << "Unrecognized option." << '\n';
    }
    
    return 0;
}
