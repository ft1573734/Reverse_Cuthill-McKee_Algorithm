#include <stdlib.h>
#include <iostream>
#include "DataLoader.h"
#include "ReverseCuthillMckee.h"
extern "C" {
#include "fasp.h"
#include "fasp_functs.h"
}

using namespace std;
int main1()
{
	DataLoader loader = DataLoader();

	string path = "D:\\ProgramFiles\\NaViiX_fasp_AMG\\kvlcc2_AMG\\amg_par\\";
	string matrix_name = "Amg.m";
	//string out_path = "D:\\ProgramFiles\\NaViiX_fasp_AMG\\kvlcc2_AMG\\reordered\\ro_0_level.dat";
	//string path = "D:\\ProgramFiles\\NaViiX_fasp_AMG\\kvlcc2_AMG\\1\\test.dat";
	cout << "Loading file: " << path << endl;
	int* row_arr_ptr = NULL;
	int* col_arr_ptr = NULL;
	double* val_arr_ptr = NULL;


	dCSRmat original_mat;
	dCSRmat result_mat;
	loader.LoadDiagonalMatrix(path+matrix_name, &original_mat, false);
	cout << "Loading complete" << endl;

//	for (int i = 0; i < original_mat.row; i++) {
//		cout << original_mat.IA[i] << endl;
//	}
	ReverseCuthillMckee rcm;
	dCSRmat perm_matrix;
	rcm.RCM(&original_mat, &result_mat, &perm_matrix);

	//Estimating how many arrays needed for the DIA format.
	unordered_set<int> distinctive_cols;
	for(int i=0; i< result_mat.row; i++){
		cout << "Processing row "<<i<<" ;" << endl;
		int start_col_ind_of_row = result_mat.IA[i];
		//traversing this row
		for (int j = 0; j < result_mat.JA[start_col_ind_of_row + 1] - result_mat.JA[start_col_ind_of_row]; j++) {
			int tmp_diag_index = result_mat.JA[j + start_col_ind_of_row] - i;
			distinctive_cols.insert(tmp_diag_index) ;
		}
	}


	//Print Info
	int* orig_mat_bandwidths = (int*)malloc(sizeof(int) * original_mat.row);
	int* result_mat_bandwidths = (int*)malloc(sizeof(int) * original_mat.row);
	int orig_mat_bandwidth = rcm.ComputeBandwidth(&original_mat, orig_mat_bandwidths);
	int result_mat_bandwidth = rcm.ComputeBandwidth(&result_mat, result_mat_bandwidths);
	int orig_mat_env = rcm.ComputeEnvelope(&original_mat, orig_mat_bandwidths);
	int result_mat_env = rcm.ComputeEnvelope(&result_mat, result_mat_bandwidths);
	int sum_orig_bandwidth = 0;
	int sum_result_bandwidth = 0;
	for (int i = 0; i < original_mat.row; i++) {
		sum_orig_bandwidth += orig_mat_bandwidths[i];
		sum_result_bandwidth += result_mat_bandwidths[i];
	}
	cout << "Original Matrix statistics:" << endl;
	cout << "Bandwidth is " << orig_mat_bandwidth << ". " << "Average bandwidth is " << sum_orig_bandwidth / original_mat.row << ". " << "Envelope is " << orig_mat_env <<". " << endl;
	cout << "========================================" << endl;
	cout << "Converted Matrix statistics:" << endl;
	cout << "Bandwidth is " << result_mat_bandwidth << ". " << "Average bandwidth is " << sum_result_bandwidth / original_mat.row << ". " << "Envelope is " << result_mat_env <<". " << endl;


	dCOOmat coo_matrix;
	fasp_format_dcsr_dcoo(&result_mat, &coo_matrix);
	coo_matrix.nnz = result_mat.nnz;

	string header = "Re-ordered matrix";

	string out_path = "D:\\ProgramFiles\\NaViiX_fasp_AMG\\kvlcc2_AMG\\reordered\\ro_0_level.dat";
	loader.WriteMatrix(out_path, &coo_matrix, true);


	system("pause");

	return 0;
}

