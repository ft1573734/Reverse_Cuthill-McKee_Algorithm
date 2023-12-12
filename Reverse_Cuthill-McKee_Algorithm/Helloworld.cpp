#include <stdlib.h>
#include <iostream>
#include "DataLoader.h"

using namespace std;
int main()
{
	DataLoader loader = DataLoader();

	string path = "D:\\ProgramFiles\\NaViiX_fasp_AMG\\kvlcc2_AMG\\1\\0_level.dat";


	int* row_arr_ptr = NULL;
	int* col_arr_ptr = NULL;
	double* val_arr_ptr = NULL;

	loader.LoadDiagonalMatrix(path, row_arr_ptr, col_arr_ptr, val_arr_ptr);
	for (int i = 0; i < 10; i++) {
		cout << row_arr_ptr[i] << endl;
	}
	system("pause");
}

