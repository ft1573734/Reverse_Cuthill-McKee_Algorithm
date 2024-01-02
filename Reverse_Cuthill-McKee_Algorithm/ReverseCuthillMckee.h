#pragma once
#include <unordered_set>
#include "fasp.h"
using namespace std;

class ReverseCuthillMckee {
public:
	ReverseCuthillMckee();

	/**
	 * RCM(.,.,.) implements the RCM algorithm.
	 * The algorithm takes three inputs, a diagonal input matrix in the format of CSR, a diagonal output matrix in the format of CSR and a permutation matrix responsible for the conversion.
	 * The algorithm re-orders the input matrix in a more compact manner, so that the bandwidth of the diagnoal matrix is reduced.
	 * The time complexity of the algorithm is high, further optimization is needed.
	 *
	 * Input:
	 *		*input: a diagonal matrix in the format of CSR;
	 *		*output: a diagonal matrix in the format of CSR.
	 *		*perm_matrix: the permutation matrix used for converting the input to the output i.e. O = PIP^T
	 * Output:
	 *		void
	 *
	 * Implemented & updated by Xavier Wang on 12/15/2023.
	 */
	void RCM(dCSRmat* input, dCSRmat* output, dCSRmat* perm_matrix);

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
	int Peripheral_Node_Finder(dCSRmat* input);

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
	vector<unordered_set<int>> GenerateLevelStructure(dCSRmat* matrix, int initial_node);

	/**
	* Compute the bandwidth of a symmetric matrix, col by col / row by row.
	* We denote the bandwidth of each col as b(i)
	* Theoretically, the bandwidth of a symmetric matrix is max(b(i)).
	* However, sometimes we want to study the bandwidth of each col.
	* Therefore, this funtion returns an array corresponding to the bandwidth of each col.
	* Input:
	*		*matrix: a pointer to a symmetric matrix.
	*		int*: a pointer to the 1st element of the bandwidth array.
	* Output:
	*		int: the maximum bandwidth of the cols/rows.
	* Implemented by Xavier Wang on 12/26/2023.
	*/
	int ComputeBandwidth(dCSRmat* matrix, int* bandwidths);

	/**
	* Compute the envelope of a symmetric matrix.
	* The formal definition of the envelope of a symmetric matrix is as follows:
	* Env(A) = {(i,j)| 0 < j-i <= b(i), 1<=i,j<=n}[1]
	* [1]The Reverse Cuthill-McKee Algorithm in Distributed-Memory, 2017.
	* 
	* An intuitive understanding is that the envelope is the area surrounded by the edge-NNZs of the upper-triangle of the matrix.
	* 
	* Input:
	*		*matrix: a pointer to a symmetric matrix.
	* Output:
	*		int: the envelope of the symmetric matrix.
	*/		
	int ComputeEnvelope(dCSRmat* matrix, int* bandwidths);

private:

	/**
	 * find_challenger(.) is a PRIVATE function used for finding the 'challenger node' within the last layer of the level-structure. There are three possible schemas:
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
	int find_challenger(unordered_set<int> vector_set, dCSRmat* matrix, int schema);
};
