#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <list>
#include "DataLoader.h"
extern "C" {
#include "fasp.h"
#include "fasp_functs.h"
}

using namespace std;
DataLoader::DataLoader() {

}

void DataLoader::LoadDiagonalMatrix(string path, dCSRmat* coo_matrix, bool has_header) {

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
        bool isheader = has_header;
        while (getline(file, line)) {
            if (isheader) {
                isheader = false;
                continue;
            }

            string delimiter = " ";

            size_t pos = 0;
            string token;
            int token_count = 0;

            int row_id;
            int col_id;
            double val;

            
            bool end_of_line = false;
            while (!end_of_line) {
                pos = line.find(delimiter);
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
                    end_of_line = true;
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
        //The size of the row_index_list is m+1, therefore another push_back and the end of the sequence is needed
        row_index_list.push_back(counter);
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

void DataLoader::WriteMatrix(string path, dCOOmat* coo_matrix, bool need_header){
    string header = "Re-ordered matrix";
    string delimiter = " ";
    ofstream myfile(path);
    if (myfile.is_open())
    {
        if (need_header) {
            myfile << header << endl;
        }
        for (int i = 0; i < coo_matrix->nnz; i++) {
            myfile << coo_matrix->rowind[i] << delimiter << coo_matrix->colind[i] << delimiter << delimiter << coo_matrix->val[i] << delimiter << endl;
        }

        myfile.close();
    }
    else cout << "Unable to open file";
}

void DataLoader::WriteMatrix(string path, dCSRmat* csr_matrix, bool need_header) {
    string header = "Re-ordered matrix";
    dCOOmat coo_matrix;
    fasp_format_dcsr_dcoo(csr_matrix, &coo_matrix);
    //The fasp_format_dcsr_dcoo(.,.) might have a problem, I have to manually set the statistical values.


    coo_matrix.nnz = csr_matrix->nnz;

    WriteMatrix(path, &coo_matrix, need_header);
}


void DataLoader::LoadFaspAMGParameters(string path, AMG_param* param){
    //Insert manually
    param->AMG_type = CLASSIC_AMG;
    param->print_level = PRINT_MIN;
    param->maxit = 100;
    param->tol = 1e-6;
    param->max_levels = 20;
    param->coarse_dof = 500;

    param->cycle_type = V_CYCLE;
    param->quality_bound = 10.0;
    param->smoother = SMOOTHER_GS;
    param->smooth_order = CF_ORDER;
    param->presmooth_iter = 1;
    param->postsmooth_iter = 1;
    param->relaxation = 1.0;
    param->polynomial_degree = 3;

    param->coarse_solver = 0;
    param->coarse_scaling = OFF;
    param->amli_degree = 2;
    param->amli_coef = NULL;
    param->nl_amli_krylov_type = SOLVER_GCG;
    param->coarsening_type = COARSE_RS;
    param->aggregation_type = PAIRWISE;
    param->aggregation_norm_type = -1; //not specified
    param->interpolation_type = 1;  
    param->strong_threshold = 0.3;
    param->max_row_sum = 0.9;

    param->truncation_threshold = 0.2;
    param->aggressive_level = 0;
    param->aggressive_path = 1;
    param->pair_number = 2;
    param->strong_coupled = 0.08;
    param->max_aggregation = 20;
    param->tentative_smooth = 0.67;
    param->smooth_filter = ON;
    param->smooth_restriction = ON;

    param->ILU_levels = 0;
    param->ILU_type = ILUk;
    param->ILU_lfil = 0;
    param->ILU_droptol = 0.001;
    param->ILU_relax = 0;
    param->ILU_permtol = 0.0;

    param->SWZ_levels = 0;
    param->SWZ_mmsize = 200;
    param->SWZ_maxlvl = 1;
    param->SWZ_type = 1;
    param->SWZ_blksolver = SOLVER_DEFAULT;
    param->theta = -1; //not specified
}

void DataLoader::LoadVector(string path, dvector* vector, bool has_header){
    list<double> val_list;
    string line;
    ifstream file(path);
    if (file.is_open()) {
        bool isheader = has_header;
        while (getline(file, line)) {
            if (isheader) {
                isheader = false;
                continue;
            }
            val_list.push_back(stod(line));
        }
        file.close();
        int row = val_list.size();
        double* vals_arr = (double*)malloc(sizeof(double) * row);
        copy(val_list.begin(), val_list.end(), vals_arr);

        vector->row = row;
        vector->val = vals_arr;
    }
    else {
        cerr << "Read file error" << endl;
    }
}

void DataLoader::WriteVector(string path, dvector* v) {
    ofstream myfile(path);
    if (myfile.is_open())
    {
        for (int i = 0; i < v->row; i++) {
            myfile << v->val[i] << endl;
        }

        myfile.close();
    }
    else cout << "Unable to open file";
}
