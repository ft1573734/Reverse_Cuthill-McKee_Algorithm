#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <list>
#include "DataLoader.h"
#include "fasp.h"

using namespace std;
DataLoader::DataLoader() {

}

void DataLoader::LoadDiagonalMatrix(string path, dCSRmat* coo_matrix) {

    list<int> row_index_list;
    list<int> col_list;
    list<double> val_list;


    int counter = 0;
    int tmp_row_id = -1;

    int row_count = 0;
    int col_count = 0;

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
                    if (row_id > row_count) {
                        row_count = row_id;
                    }
                    break;
                case 2:
                    col_id = stoi(token);
                    if (col_id > col_count) {
                        col_count = col_id;
                    }
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


            if (row_id != tmp_row_id) {
                row_index_list.push_back(counter);
                tmp_row_id = row_id;
            }
            counter++;


            col_list.push_back(col_id);
            val_list.push_back(val);


        }
        file.close();
        
        //Since the index starts with 0, a "++" on the row/col count is necessary.
        row_count++;
        col_count++;


        int val_count = val_list.size();


        int* row_arr = new int[row_index_list.size()];
        int* col_arr = new int[col_list.size()];
        double* vals_arr = new double[val_list.size()];

        copy(row_index_list.begin(), row_index_list.end(), row_arr);
        copy(col_list.begin(), col_list.end(), col_arr);
        copy(val_list.begin(), val_list.end(), vals_arr);

        


        coo_matrix->row = row_count;
        coo_matrix->col = col_count;
        coo_matrix->nnz = val_count;
        coo_matrix->IA = row_arr;
        coo_matrix->JA = col_arr;
        coo_matrix->val = vals_arr;



    }
    else {
        cerr << "Read file error" << endl;
    }

}

void WriteDiagonalMatrix(string path, dCOOmat* coo_matrix){
    string header = "Re-ordered matrix";
    ofstream myfile(path);
    if (myfile.is_open())
    {
        string delimiter = "  ";
        for (int i = 0; i < coo_matrix->row; i++) {
            myfile << coo_matrix->rowind[i]<<delimiter<<coo_matrix->colind[i]<<delimiter<<delimiter<<coo_matrix->val[i]<<delimiter << endl;
        }

        myfile.close();
    }
    else cout << "Unable to open file";

}
