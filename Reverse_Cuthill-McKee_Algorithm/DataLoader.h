#pragma once
#ifndef DATA_LOADER_H
#define DATA_LOADER_H

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <list>
#include "./base/include/fasp.h"
using namespace std;

class DataLoader {
public:
	DataLoader();

	void LoadDiagonalMatrix(string path, dCSRmat* csr_matrix);

};

#endif
