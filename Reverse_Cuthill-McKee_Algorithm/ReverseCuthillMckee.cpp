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
 * RCM(.,.,.) implements the RCM algorithm.
 * The algorithm takes three inputs, a diagonal input matrix in the format of CSR, a diagonal output matrix in the format of CSR and a permutation matrix responsible for the conversion.
 * The algorithm re-orders the input matrix in a more compact manner, so that the bandwidth of the diagnoal matrix is reduced.
 * The time complexity of the algorithm is high, further optimization is needed.
 *
 * Input:
 *		*input: a diagonal matrix in the format of CSR;
 *		*output: a diagonal matrix in the format of CSR;
 *		*perm_matrix: the permutation matrix used for converting the input to the output i.e. O = PIP^T.
 * Output:
 *		void
 *
 * Implemented & updated by Xavier Wang on 12/28/2023.
 */
void ReverseCuthillMckee::RCM(dCSRmat* input, dCSRmat* output, dCSRmat* perm_matrix) {
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
	
	/*
	* Generating re-ordering schema...
	*/

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


	/*
	* Constructing the permutation matrix...
	*/
	int* perm_matrix_row = (int*)malloc(sizeof(int) * row_count);
	int* perm_matrix_col = (int*)malloc(sizeof(int) * row_count);
	double* perm_matrix_vals = (double*)malloc(sizeof(double) * row_count);

	//initialize the permutation matrix
	for (int i = 0; i < row_count; i++) {
		perm_matrix_vals[i] = 1.0;
		perm_matrix_row[i] = i;
		perm_matrix_col[i] = i;
	}
	//reordering the rows of the perm_matrix
	for (int i = 0; i < row_count; i++) {
		int new_ind = mp.at(perm_matrix_row[i]);
		perm_matrix_row[i] = new_ind;
	}

	dCOOmat perm_matrix_coo;
	perm_matrix_coo.row = row_count;
	perm_matrix_coo.col = row_count; //There is no problem here, since perm matrix is a square matrix with 1 nnz per row/col.
	perm_matrix_coo.nnz = row_count; //Therefore, row, col and nnz counts are all row_count (or col_count, we are processing symmetic matrix here).
	perm_matrix_coo.rowind = perm_matrix_row;
	perm_matrix_coo.colind = perm_matrix_col;
	perm_matrix_coo.val = perm_matrix_vals;

	//Convert perm_matrix to CSR format, since most of the computations in fasp only supports CSR.
	fasp_format_dcoo_dcsr(&perm_matrix_coo, perm_matrix);


	//Compute PAP':

	//Construct the transpose of P i.e. P'
	dCSRmat perm_matrix_T;
	fasp_dcsr_trans(perm_matrix, &perm_matrix_T);
	//output = R*A*P
	fasp_blas_dcsr_rap(perm_matrix, input, &perm_matrix_T ,output);

	/*****
	* Deprecated: Conversion in old-fashioned way.
	*****/
	/*
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

	fasp_format_dcoo_dcsr(&coo_matrix, output);
	*/
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
		//cout << "tmp_level size is " << tmp_level.size() << ", including nodes:" << endl;
		//cout << endl;
		//cout << "Remainder Node Count: " << remainder_nodes.size() << endl;

	}
	cout << "The size of the level structure is " << level_structure.size() << endl;
	return level_structure;
}


/**
 * find_challenger(.) is a PRIVATE function used for finding the 'challenger node' within the last layer of the level-		structure. There are three possible schemas:
 *	1. select the node with the smallest degree;
 *	2. TODO: select nodes of the smallest degree from each of the connected components within the last layer.
 *	3. select a random point. (not recommended)
 * So far only schema 1 is implemented.
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
	//Schema 1: The challenger is the node in the last level of the level-structure, with the smallest degree;
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

/**
* Compute the bandwidth of a symmetric matrix, col by col / row by row.
* We denote the bandwidth of each col as b(i)
* Theoretically, the bandwidth of a symmetric matrix is max(b(i)).
* However, sometimes we want to study the bandwidth of each col.
* Therefore, this funtion returns an array corresponding to the bandwidth of each col.
* Input:
*		*matrix: a pointer to a symmetric matrix.
* Output:
*		int* a pointer to the 1st element of the bandwidth array.
* Implemented by Xavier Wang on 12/26/2023.
*/
int ReverseCuthillMckee::ComputeBandwidth(dCSRmat* matrix, int* bandwidths) {
 	int row_count = matrix->row;
	int max_bandwidth = 0;
	int* IA = matrix->IA;
	int* JA = matrix->JA;
	int sum_width = 0;
	for (int i = 0; i < row_count; i++) {
		int min_col_ind = INT_MAX;
		int max_col_ind = -1;
		for (int j = 0; j < (IA[i + 1] - IA[i]); j++) {
			if (JA[IA[i] + j] < min_col_ind) {
				min_col_ind = JA[IA[i] + j];
			}
			if (JA[IA[i] + j] > max_col_ind) {
				max_col_ind = JA[IA[i] + j];
			}
		}
		sum_width += max_col_ind - min_col_ind;
		bandwidths[i] = i - min_col_ind;
		if (bandwidths[i] > max_bandwidth) {
			max_bandwidth = bandwidths[i];
		}
	}
	cout << "Average width per row is " << sum_width / row_count << "." << endl;
	return max_bandwidth;
}

/**
* Compute the envelope of a symmetric matrix.
* The official definition of the envelope of a symmetric matrix is as follows:
* Env(A) = {(i,j)| 0 < j-i <= b(i), 1<=i,j<=n}[1]
* [1]The Reverse Cuthill-McKee Algorithm in Distributed-Memory, 2017.
*
* An intuitive understanding is that the envelope is the area surrounded by the edge-NNZs of the upper-triangle of the matrix.
*
* Input:
*		*matrix: a pointer to a symmetric matrix.
*		int*: a pointer to the 1st element of the bandwidth array.
* Output:
*		int: the envelope of the symmetric matrix.
* Implemented by Xavier Wang on 12/26/2023.
*/
int ReverseCuthillMckee::ComputeEnvelope(dCSRmat* matrix, int* bandwidths) {
	int envelope = 0;
	int row_count = matrix->row;
	for (int i = 0; i < row_count; i++) {
		envelope += bandwidths[i];
	}
	return envelope;
}