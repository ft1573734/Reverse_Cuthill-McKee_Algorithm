#pragma once
#ifndef DATA_LOADER_H
#define DATA_LOADER_H

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <list>
using namespace std;

class DataLoader {
public:
	DataLoader();

	void LoadDiagonalMatrix(string path, int* row_arr, int* col_arr, double* vals);
};

#endif
