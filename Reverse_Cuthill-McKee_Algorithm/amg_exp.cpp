#include <stdlib.h>
#include <iostream>
#include "DataLoader.h"
#include "ReverseCuthillMckee.h"
extern "C" {
#include "fasp.h"
#include "fasp_functs.h"
}
using namespace std;

int main() {
	/*Loading parameters*/
	DataLoader loader = DataLoader();

	string path = "D:\\ProgramFiles\\NaViiX_fasp_AMG\\kvlcc2_AMG\\amg_par\\";
	string matrix_name = "Amg.m";
	string param_name = "fasp_amg_param.txt";
	string vector_name = "rhs.m";
	
	string out_mat_name = "reordered_matrix";
	string out_perm_name = "permutatio_matrix";
	string out_vector_name = "reoredered_vector";

	cout << "Loading file: " << path << endl;
	int* row_arr_ptr = NULL;
	int* col_arr_ptr = NULL;
	double* val_arr_ptr = NULL;

	dCSRmat orig_matrix;
	dvector orig_vector;
	dvector orig_x;

	dCSRmat ro_matrix;
	dvector ro_vector;
	dvector ro_x;

	AMG_param param;

	dCSRmat perm_matrix;


	/**
	* =====Executing AMG======
	*/


	//Loading parameters
	loader.LoadDiagonalMatrix(path + matrix_name, &orig_matrix, false);
	loader.LoadVector(path + vector_name, &orig_vector, false);
	loader.LoadFaspAMGParameters("MANUAL", &param);


	//Performing re-ordering
	ReverseCuthillMckee rcm;
	rcm.RCM(&orig_matrix, &ro_matrix, &perm_matrix);
	ro_vector.row = orig_vector.row;
	ro_vector.val = (double*)malloc(sizeof(double) * ro_vector.row);
	for (int i = 0; i < ro_vector.row; i++) {
		ro_vector.val[i] = orig_vector.val[i];
	}
	fasp_blas_dcsr_mxv(&perm_matrix,  orig_vector.val, ro_vector.val);

	double* orig_x_arr = (double*)malloc(sizeof(double) * orig_vector.row);
	for (int i = 0; i < orig_vector.row; i++) {
		orig_x_arr[i] = 0;
	}
	orig_x.row = orig_vector.row;
	orig_x.val = orig_x_arr;

	double* ro_x_arr = (double*)malloc(sizeof(double) * orig_vector.row);
	for (int i = 0; i < ro_vector.row; i++) {
		ro_x_arr[i] = 0;
	}
	ro_x.row = ro_vector.row;
	ro_x.val = ro_x_arr;
	
	//loader.WriteMatrix(path + out_mat_name, &ro_matrix, false);
	//loader.WriteMatrix(path + out_perm_name, &perm_matrix, false);
	//loader.WriteVector(path + out_vector_name, &ro_vector);

	//call amg func.
	fasp_solver_amg(&orig_matrix, &orig_vector, &orig_x, &param);

	//fasp_solver_amg(&ro_matrix, &ro_vector, &ro_x, &param);

	return 0;
}