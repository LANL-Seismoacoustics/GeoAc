#ifndef _GEOAC_INTERFACE_CPP_
#define _GEOAC_INTERFACE_CPP_

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>

#include "GeoAc.Parameters.h"
#include "Atmo_State.h"
#include "GeoAc.EquationSets.h"
#include "GeoAc.Solver.h"

using namespace std;

//--------------------------------------//
//----Functions To Set Up The Solver----//
//--------------------------------------//
void GeoAc_SetEqCnt(){
	switch(GeoAc_dim){
		case(2):
			if(GeoAc_CalcAmp && GeoAc_AtmoStrat){               GeoAc_EqCnt = 6;}	// 2D propagation solving for rays and amplitude in stratified atmosphere
            else if(!(GeoAc_CalcAmp) && GeoAc_AtmoStrat){       GeoAc_EqCnt = 3;}	// 2D propagation solving for rays only in stratified atmosphere
            else if(GeoAc_CalcAmp && !(GeoAc_AtmoStrat)){       GeoAc_EqCnt = 8;}   // 2D propagation solving for rays and amplitude in horizontally varying atmosphere
            else if(!(GeoAc_CalcAmp) && !(GeoAc_AtmoStrat)){    GeoAc_EqCnt = 4;} 	// 2D propagation solving for rays only in horizontally varying atmosphere
			break;
            
		case(3):
			if(GeoAc_CalcAmp && GeoAc_AtmoStrat){               GeoAc_EqCnt = 12;}	// 3D propagation solving for rays and amplitude in stratified atmosphere
            else if(!(GeoAc_CalcAmp) && GeoAc_AtmoStrat){		GeoAc_EqCnt = 4;}   // 3D propagation solving for rays only in stratified atmosphere
            else if(GeoAc_CalcAmp && !(GeoAc_AtmoStrat)){		GeoAc_EqCnt = 18;}  // 3D propagation solving for rays and amplitude in horizontally varying atmosphere
			else if(!(GeoAc_CalcAmp) && !(GeoAc_AtmoStrat)){    GeoAc_EqCnt = 6;}   // 3D propagation solving for rays only in horizontally varying atmosphere
			break;

		default:
            cout << "***WARNING*** GeoAc Not Properly Configured" << '\n';
			break;
	}
}

void GeoAc_ConfigureCalcAmp(bool calc_amp){
	GeoAc_CalcAmp = calc_amp;
    GeoAc_SetSystem();    
	GeoAc_SetEqCnt();
}


//----------------------------------------------------------------//
//----Functions To Build, Clear, and Delete the Solution Array----//
//----------------------------------------------------------------//
void GeoAc_BuildSolutionArray(double ** & solution, int length){
	solution = new double* [length];
	for(int n = 0; n < length; n++){
		solution[n] = new double [GeoAc_EqCnt];
	}
}

void GeoAc_ClearSolutionArray(double ** & solution, int index){
	for(int n1 = 0; n1 < index; n1++){
		for(int n2 = 0; n2 < GeoAc_EqCnt; n2++){
			solution[n1][n2] = 0.0;
		}
	}
}

void GeoAc_DeleteSolutionArray(double ** & solution, int length){
	for(int n = 0; n < length; n++){
		delete [] solution[n];
	}
	delete [] solution;
}

//-------------------------------------------------//
//---------Functions To Output The Profile---------//
//-------------------------------------------------//
void GeoAc_WriteProfile(string file_name, double azimuth){
	ofstream file_out;	file_out.open(file_name.c_str() );
    
	if(!file_out.is_open()){
		cout << "Error opening file, check file name." << '\n';
	} else {
        for(int m = 0; m < 1400; m++){
            double z0 = m/10.0;
            file_out << z0 << '\t';
            file_out << pow(c(0, 0, z0), 2) / (R*gam) << '\t';
            file_out << u(0,0, z0) * 1000.0 << '\t';
            file_out << v(0,0, z0) * 1000.0 << '\t';
            file_out << rho(0,0, z0) << '\t';
            file_out << rho(0,0, z0) * pow(c(0,0,z0),2)/gam << '\t';
            file_out << c(0,0, z0) << '\t';
            file_out << c(0,0, z0) + cos(azimuth*Pi/180.0)*u(0,0,z0) + sin(azimuth*Pi/180.0)*v(0,0,z0) << '\t';
            file_out << '\n';
        }	
        file_out.close();
    }
}

void GeoAc_WriteProfile(string file_name, double x0, double y0, double azimuth){
	ofstream file_out;	file_out.open(file_name.c_str() );
    
	if(!file_out.is_open()){
		cout << "Error opening file, check file name." << '\n';
	} else {
        for(int m = 0; m < 1400; m++){
            double z0 = m/10.0;
            file_out << z0 << '\t';
            file_out << pow(c(x0, y0, z0), 2) / (R*gam) << '\t';
            file_out << u(x0, y0, z0) * 1000.0 << '\t';
            file_out << v(x0, y0, z0) * 1000.0 << '\t';
            file_out << rho(x0, y0, z0) << '\t';
            file_out << rho(x0, y0, z0) * pow(c(x0, y0, z0),2)/gam << '\t';
            file_out << c(x0, y0, z0) << '\t';
            file_out << c(x0, y0, z0) + cos(azimuth * Pi/180.0) * u(x0, y0, z0) + sin(azimuth * Pi/180.0) * v(x0, y0, z0) << '\t';
            file_out << '\n';
        }
        file_out.close();
    }
}


bool string2bool (std::string v){
    return !v.empty () &&
    (strcasecmp (v.c_str (), "true") == 0 || atoi (v.c_str ()) != 0);
}



#endif /* _GEOAC_INTERFACE_CPP_ */
