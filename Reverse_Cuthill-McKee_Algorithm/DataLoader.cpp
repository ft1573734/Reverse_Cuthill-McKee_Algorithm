#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <list>
#include "DataLoader.h"
#include "./base/include/fasp.h"

using namespace std;
DataLoader::DataLoader() {

}

void DataLoader::LoadDiagonalMatrix(string path, dCSRmat* coo_matrix) {

    list<int> row_index_list;
    list<int> col_list;
    list<double> val_list;

    int num_row = 0;
    int num_col = 0;

    int counter = 0;
    int tmp_row_id = -1;

    string line;
    ifstream file(path);
    if (file.is_open()) {
        bool isheader = true;
        while (getline(file, line)) {
            if (isheader) {
                isheader = false;
                continue;
            }

            string delimiter = "  ";

            size_t pos = 0;
            string token;
            int token_count = 0;

            int row_id;
            int col_id;
            double val;

            

            while ((pos = line.find(delimiter)) != string::npos) {
                token = line.substr(0, pos);
                token_count++;
                switch (token_count) {
                case 1:
                    row_id = stoi(token);
                    break;
                case 2:
                    col_id = stoi(token);
                    break;
                case 3:
                    break;
                case 4:
                    val = stod(token);
                    token_count = 0;
                    break;
                }
                line.erase(0, pos + delimiter.length());
            }

            if (num_row < row_id) {
                num_row = row_id;
            }
            if (num_col < col_id) {
                num_col = col_id;
            }




            if (row_id != tmp_row_id) {
                row_index_list.push_back(counter);
                tmp_row_id = row_id;
            }
            else {
                counter++;
            }


            col_list.push_back(col_id);
            val_list.push_back(val);


        }
        file.close();

        int row_count = row_index_list.size();
        int col_count = col_list.size();
        int val_count = val_list.size();


        int* row_arr = new int[row_count];
        int* col_arr = new int[col_count];
        double* vals_arr = new double[val_count];

        copy(row_index_list.begin(), row_index_list.end(), row_arr);
        copy(col_list.begin(), col_list.end(), col_arr);
        copy(val_list.begin(), val_list.end(), vals_arr);

        


        coo_matrix->row = num_row;
        coo_matrix->col = num_col;
        coo_matrix->nnz = val_count;
        coo_matrix->IA = row_arr;
        coo_matrix->JA = col_arr;
        coo_matrix->val = vals_arr;



    }
    else {
        cerr << "Read file error" << endl;
    }


}
