#include <stdlib.h>
#include <iostream>
#include "DataLoader.h"
#include "ReverseCuthillMckee.h"
extern "C" {
#include "fasp.h"
}

using namespace std;
int main()
{
	DataLoader loader = DataLoader();

	string path = "D:\\ProgramFiles\\NaViiX_fasp_AMG\\kvlcc2_AMG\\1\\0_level.dat";
	//string out_path = "D:\\ProgramFiles\\NaViiX_fasp_AMG\\kvlcc2_AMG\\reordered\\ro_0_level.dat";
	//string path = "D:\\ProgramFiles\\NaViiX_fasp_AMG\\kvlcc2_AMG\\1\\test.dat";
	cout << "Loading file: " << path << endl;
	int* row_arr_ptr = NULL;
	int* col_arr_ptr = NULL;
	double* val_arr_ptr = NULL;



	dCSRmat original_mat;
	dCOOmat result_mat;
	loader.LoadDiagonalMatrix(path, &original_mat);
	cout << "Loading complete" << endl;

//	for (int i = 0; i < original_mat.row; i++) {
//		cout << original_mat.IA[i] << endl;
//	}
	ReverseCuthillMckee rcm;
	rcm.RCM(&original_mat, &result_mat);
	//loader.WriteDiagonalMatrix(out_path, &result_mat);
	system("pause");
}

