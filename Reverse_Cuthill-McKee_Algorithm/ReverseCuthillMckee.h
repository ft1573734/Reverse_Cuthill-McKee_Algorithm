#pragma once
#include <unordered_set>
#include "include/fasp.h"
using namespace std;

class ReverseCuthillMckee {
public:
	ReverseCuthillMckee();

	void RCM(dCSRmat* input, dCSRmat* output);

	int Peripheral_Node_Finder(dCSRmat* input);

	vector<unordered_set<int>> GenerateLevelStructure(dCSRmat* matrix, int initial_node);

private:
	int find_challenger(unordered_set<int> vector_set, dCSRmat* matrix, int schema);
};