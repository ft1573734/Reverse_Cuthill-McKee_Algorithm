#pragma once
#ifndef DATA_LOADER_H
#define DATA_LOADER_H

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <list>
#include "fasp.h"
using namespace std;

class DataLoader {
public:
	DataLoader();

	void LoadDiagonalMatrix(string path, dCSRmat* csr_matrix, bool has_header);

	void WriteMatrix(string path, dCOOmat* coo_matrix, bool need_header);

	void WriteMatrix(string path, dCSRmat* csr_matrix, bool need_header);

	void WriteVector(string path, dvector* v);

	void LoadFaspAMGParameters(string path, AMG_param* param);

	void LoadVector(string path, dvector* vector, bool has_header);
};

#endif
