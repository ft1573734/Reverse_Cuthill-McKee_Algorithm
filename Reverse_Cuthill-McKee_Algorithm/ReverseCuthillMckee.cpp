#include "ReverseCuthillMckee.h"
#include <iostream>
#include <vector>
#include <unordered_set>
#include <random>
#include <map>
#include <fstream>
extern "C" {
#include "fasp.h"
#include "fasp_functs.h"
}
using namespace std;
ReverseCuthillMckee::ReverseCuthillMckee(){

}
/**
 * RCM(.,.) implements the RCM algorithm.
 * The algorithm takes two inputs, a diagonal input matrix in the format of CSR, and a diagonal output matrix in the format of CSR.
 * The algorithm re-orders the input matrix in a more compact manner, so that the bandwidth of the diagnoal matrix is reduced.
 * The time complexity of the algorithm is high, further optimization is needed.
 * 
 * Input: 
 *		*input: a diagonal matrix in the format of CSR;
 *		*output: a diagonal matrix in the format of CSR.
 * Output: 
 *		void
 * 
 * Implemented by Xavier Wang on 12/15/2023.
 */
void ReverseCuthillMckee::RCM(dCSRmat* input, dCOOmat* output) {
	//Initialization:
	int row_count = input->row;
	int col_count = input->col;
	int nnz_count = input->nnz;
	int* IA = input->IA;
	int* JA = input->JA;
	double* vals = input->val;

	int chosen_node = -1;
	int challenger_node = -1;

	if (row_count != col_count) {
		cerr << "Error in RCM(.): the input matrix is not diagonal" << endl;
		return;
	}
	
	//Find the pseudo-peripheral node.
	int p_node = Peripheral_Node_Finder(input);
	
	int new_index = row_count - 1;

	vector<unordered_set<int>> level_structure = GenerateLevelStructure(input, p_node);

	//Generate the new index schema based on the level-strucure, in REVERSE (Reverse Cuthill Mckee).
	map<int, int> mp;
	for (unordered_set<int> i : level_structure) {
		for (int j : i) {
			mp.insert(pair<int, int>(j, new_index));
			new_index--;
		}
	}

	//Re-order the matrix using the new indices. A COO-formatted matrix is more suitable for the job.
	dCOOmat coo_matrix;
	fasp_format_dcsr_dcoo(input, &coo_matrix);
	
	//The fasp_format_dcsr_dcoo(.,.) might have a problem, I have to manually set the statistical values.
	coo_matrix.row = row_count;
	coo_matrix.col = col_count;
	coo_matrix.nnz = nnz_count;

	int* const reordered_rowind = new int[nnz_count];
	int* const reordered_colind = new int[nnz_count];


	for (int i = 0; i < coo_matrix.nnz; i++) {
		reordered_rowind[i] = mp.at(coo_matrix.rowind[i]);
		reordered_colind[i] = mp.at(coo_matrix.colind[i]);
	}

	coo_matrix.rowind = reordered_rowind;
	coo_matrix.colind = reordered_colind;
	
	string header = "Re-ordered matrix";

	string out_path = "D:\\ProgramFiles\\NaViiX_fasp_AMG\\kvlcc2_AMG\\reordered\\ro_0_level.dat";
	ofstream myfile(out_path);
	if (myfile.is_open())
	{
		string delimiter = "  ";
		myfile << header << endl;
		for (int i = 0; i < coo_matrix.nnz; i++) {
			myfile << coo_matrix.rowind[i] << delimiter << coo_matrix.colind[i] << delimiter << delimiter << coo_matrix.val[i] << delimiter << endl;
		}

		myfile.close();
	}
	else cout << "Unable to open file";

}



/**
 * Peripheral_Node_Finder(.) is a PUBLIC function that implements a pseudo-peripheral node finder algorithm, described in paper "An Implementation of a Pseudoperipheral Node Finder, 1979".
 * The algorithm takes one input, a diagonal input matrix in the format of CSR, representing the connectivity of a graph.
 * The algorithm returns the index of the Pseudo-peripheral node.
 * The time complexity of the algorithm is high, further optimization is needed.
 *
 * Input:
 *		*input: a diagonal matrix in the format of CSR;
 * Output:
 *		Index of the pseudo-peripheral node.
 *
 * Implemented by Xavier Wang on 12/15/2023.
 */

int ReverseCuthillMckee::Peripheral_Node_Finder(dCSRmat* input) {
	cout << "finding peripheral node" << endl;
	//Initialization:
	int row_count = input->row;
	int col_count = input->col;
	int nnz_count = input->nnz;
	int* IA = input->IA;
	int* JA = input->JA;
	double* vals = input->val;

	int chosen_node = -1;
	int challenger_node = -1;

	if (row_count != col_count) {
		cerr << "Error in Peripheral_Node_Finder (.), the input matrix is not diagonal" << endl;
		return 0;
	}
	//STEP 1: Since the chosen node is vacant, randomly set a node as the chosen node first;
	//Notice: there are other solutions for finding the starting node, a random approach is the most efficient, yet the randomness may cause problem when debugging due to randomness.
	//Other approaches include using the node with the minimum degree (which is more expensive), etc.

	if (chosen_node == -1) {
		random_device rd;     // Only used once to initialise (seed) engine
		mt19937 rng(rd());    // Random-number engine used (Mersenne-Twister in this case)
		uniform_int_distribution<int> uni(0, input->row-1); // Guaranteed unbiased

		chosen_node = uni(rng);
	}

	chosen_node = 0;
	
	//STEP 2: Generate the level structure of the chosen_node, which basically means scanning the graph layer by layer.
	vector<unordered_set<int>> level_structure = GenerateLevelStructure(input, chosen_node);
	//So far the level structure is generated, and the corresponding eccentricity is set.
	int chosen_node_eccentricity = level_structure.size();

	//STEP 3: Find the challenger in the level_structure
	challenger_node = find_challenger(level_structure.back(), input, 1);

	//STEP 4: Generate a level structure based on the new chosen node.
	//calculate the eccentricy of the challenger
	level_structure.clear();
	level_structure = GenerateLevelStructure(input, challenger_node);
	int challenger_node_eccentricity = level_structure.size();

	//If the eccentricity of the challenge is larger than the chosen_node, the challenger becomes the new chosen node.
	//The cycle continues until the chosen_node is not beaten.
	while (challenger_node_eccentricity > chosen_node_eccentricity) {
		chosen_node = challenger_node;
		chosen_node_eccentricity = challenger_node_eccentricity;
		
		challenger_node = find_challenger(level_structure.back(), input, 1);
		level_structure.clear();
		level_structure = GenerateLevelStructure(input, challenger_node);
		challenger_node_eccentricity = level_structure.size();
	}
	return chosen_node;
}

/**
 * GenerateLevelStructure(.,.) is a PUBLIC function used for generating a level-structure of a graph (adjacency matrix) based on a given node.
 * The algorithm sets the given node as layer 0, the adjacent nodes of layer0 as layer 1, the adjacent nodes of layer 1 as layer 2, etc. All the layers are exclusive.
 * The algorithm is frequently used in Peripheral_Node_Finder(.).
 * Input:
 *		*matrix: a pointer to the adjacency matrix of the graph;
 *		initial_node: the index of the given node. 
 * Output:
 *		(vector<unordered_set<int>>) A level-structure.
 *
 * Implemented by Xavier Wang on 12/18/2023.
 */
vector<unordered_set<int>> ReverseCuthillMckee::GenerateLevelStructure(dCSRmat* matrix, int initial_node) {
	cout << "Generating Level Structure" << endl;
	// The distance between nodes in the i-th layer and initial node is i.

	//Initialize the remaining nodes using an unordered set:
	int chosen_node_index = initial_node;
	unordered_set<int> remainder_nodes;
	for (int i = 0; i < matrix->row; i++) {
		remainder_nodes.insert(i); //The initial remainder_nodes is the whole data set.
	}

	unordered_set<int> tmp_level;
	//Start the sequence by setting the chosen node as level 0
	tmp_level.insert(chosen_node_index);
	remainder_nodes.erase(chosen_node_index);

	vector<unordered_set<int>> level_structure;
	level_structure.push_back(tmp_level);

	while (remainder_nodes.empty() == false) {
		unordered_set<int> previous_level = tmp_level;
		tmp_level.clear();
		//record the previous level
		//find the adjacent nodes of the previous layer (excluding visited nodes)
		for (int i : previous_level) {
			//get the i-th row of the input matrix. The NNZs within the row identifies the neighbours of the i-th node
			for (int j = matrix->IA[i]; j < matrix->IA[i + 1]; j++) {
				int tmp_node = matrix->JA[j];
				if (remainder_nodes.find(tmp_node) != remainder_nodes.end()) {
					//If the node exists in the remainder_nodes:
					if (tmp_node == i) {
						//Theoretically the program should NEVER reach here, since the i-th node is included in the previous_level, and should be subsequently removed from the remainder_nodes. This part of code is included only for debugging.
						cerr << "something is wrong" << endl;
					}
					tmp_level.insert(tmp_node);
					remainder_nodes.erase(tmp_node);
				}
				else {
					continue;
				}
			}
		}
		level_structure.push_back(tmp_level);
		cout << "tmp_level size is " << tmp_level.size() << ", including nodes:" << endl;
		cout << endl;
		cout << "Remainder Node Count: " << remainder_nodes.size() << endl;

	}
	cout << "The size of the level structure is " << level_structure.size() << endl;
	return level_structure;
}


/**
 * find_challenger(.) is a PRIVATE function used for finding the 'challenger node' within the last layer of the level-structure. There are three possible schemas:
 *	1. select the node with the smallest degree;
 *	2. TODO: select a node from the connected components within the last layer. 
 *	3. select a random point. (not recommended)
 * Input:
 *		vector_set: the nodes within the last layer.
 *		*matrix: a pointer to the matrix corresponding to the level-structure.
 *		schema: a user-defined parameter (1, 2 or 3) indicating which schema to use. So far only schema 1 is implemented
 * Output:
 *		(int) Index of the challenger point.
 *
 * Implemented by Xavier Wang on 12/18/2023.
 */
int ReverseCuthillMckee::find_challenger(unordered_set<int> vector_set, dCSRmat* matrix, int schema) {

	if (schema == 1) {
	//Schema 1: The challenger is the node in the last level of the level-structure, with the highest degree;
		int min_degree = INT_MAX;
		int tmp_node_index;
		int tmp_degree;
		int challenger_node;
		for (int i : vector_set) {
			tmp_node_index = i;
			tmp_degree = matrix->IA[i + 1] - matrix->IA[i];
			if (tmp_degree < min_degree) {
				min_degree = tmp_degree;
				challenger_node = tmp_node_index;
			}
		}
		return challenger_node;
	}
	return -1;
	
}