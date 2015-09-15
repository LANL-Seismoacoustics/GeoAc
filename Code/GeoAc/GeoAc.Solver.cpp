#ifndef GEOAC_SOLVER_CPP_
#define GEOAC_SOLVER_CPP_

#include <math.h>

#include "GeoAc.Parameters.h"
#include "Atmo_State.h"
#include "GeoAc.EquationSets.h"

using namespace std;

int GeoAc_Propagate_RK4(double ** & solution, bool & check){
	int k = 0;                                                      // Integer to track ending index of solution
	int step_limit = GeoAc_ray_limit * int(1.0/(GeoAc_ds_min*10));	// Limiting number of steps, defined by limiting ray length divided by step size (length/step)
	double s = 0, ds = GeoAc_ds_min;                                // Current ray length and (variable) ray length step

	double *temp0 = new double [GeoAc_EqCnt];		double *temp1 = new double [GeoAc_EqCnt];
	double *temp2 = new double [GeoAc_EqCnt];		double *temp3 = new double [GeoAc_EqCnt];
	double *temp4 = new double [GeoAc_EqCnt];

	double *partial1 = new double [GeoAc_EqCnt];	double *partial2 = new double [GeoAc_EqCnt];
	double *partial3 = new double [GeoAc_EqCnt];

	check = false;
	for(k = 0; k < (step_limit - 1); k++){
		for (int i = 0; i < GeoAc_EqCnt; i++){
			temp0[i] = solution[k][i];
		}
		GeoAc_UpdateSources(s, temp0);
		ds = GeoAc_Set_ds(temp0);

		s += ds;
		for (int i = 0; i < GeoAc_EqCnt; i++){
			temp1[i] = ds*GeoAc_EvalSrcEq(s,temp0,i);
			partial1[i] = solution[k][i] + temp1[i]/2.0;
		}		

		GeoAc_UpdateSources(s + ds/2, partial1);
		for (int i = 0; i < GeoAc_EqCnt; i++){
			temp2[i] = ds*GeoAc_EvalSrcEq(s + ds/2.0, partial1, i);
			partial2[i] = solution[k][i] + temp2[i]/2.0;

		}

		GeoAc_UpdateSources(s + ds/2, partial2);	
		for (int i = 0; i < GeoAc_EqCnt; i++){
			temp3[i] = ds*GeoAc_EvalSrcEq(s + ds/2.0, partial2,i);
			partial3[i] = solution[k][i] + temp3[i];
		}

		GeoAc_UpdateSources(s + ds, partial3);
		for (int i = 0; i < GeoAc_EqCnt; i++){
			temp4[i] = ds*GeoAc_EvalSrcEq(s+ds, partial3, i);
			solution[k+1][i] = solution[k][i] + temp1[i]/6.0 + temp2[i]/3.0 + temp3[i]/3.0 + temp4[i]/6.0;
		}

		if(GeoAc_BreakCheck(solution,k+1)){
			check = true;
			break;
		}
		if(GeoAc_GroundCheck(solution,k+1)){
			check = false;
			break;
		}
	}
	
	delete [] temp0;	delete [] temp1;	delete [] temp2;
	delete [] temp3;	delete [] temp4;
	delete [] partial1;	delete [] partial2;	delete [] partial3;
	
	return k+1;
}

#endif /* GEOAC_SOLVER_CPP_ */
