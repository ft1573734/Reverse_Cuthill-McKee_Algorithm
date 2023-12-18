#include <stdlib.h>
#include <iostream>
#include "DataLoader.h"
#include "include\fasp.h"
#include "ReverseCuthillMckee.h"

using namespace std;
int main()
{
	DataLoader loader = DataLoader();

	string path = "D:\\ProgramFiles\\NaViiX_fasp_AMG\\kvlcc2_AMG\\1\\0_level.dat";

	int* row_arr_ptr = NULL;
	int* col_arr_ptr = NULL;
	double* val_arr_ptr = NULL;



	dCSRmat original_mat;
	dCSRmat result_mat;
	loader.LoadDiagonalMatrix(path, &original_mat);

//	for (int i = 0; i < original_mat.row; i++) {
//		cout << original_mat.IA[i] << endl;
//	}
	ReverseCuthillMckee rcm;
	rcm.RCM(&original_mat, &result_mat);

	system("pause");
}

