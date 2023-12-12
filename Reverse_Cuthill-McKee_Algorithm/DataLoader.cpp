#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <list>
#include "DataLoader.h"

using namespace std;
DataLoader::DataLoader() {

}

void DataLoader::LoadDiagonalMatrix(string path, int *row_arr_ptr, int *col_arr_ptr, double *vals_arr_ptr) {

    list<int> row_list;
    list<int> col_list;
    list<double> val_list;

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


            row_list.push_back(row_id);
            col_list.push_back(col_id);
            val_list.push_back(val);


        }
        file.close();

        int *const row_arr = new int[row_list.size()];
        int *const col_arr = new int[col_list.size()];
        double* const vals_arr = new double[val_list.size()];

        copy(row_list.begin(), row_list.end(), row_arr);
        copy(col_list.begin(), col_list.end(), col_arr);
        copy(val_list.begin(), val_list.end(), vals_arr);



        row_arr_ptr = &row_arr[0];
        col_arr_ptr = &col_arr[0];
        vals_arr_ptr = &vals_arr[0];


    }
    else {
        cerr << "Read file error" << endl;
    }
}
    

