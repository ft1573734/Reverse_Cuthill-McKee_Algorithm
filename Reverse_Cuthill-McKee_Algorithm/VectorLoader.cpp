#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <list>

using namespace std;

class VectorLoader {
	void LoadDiagonalMatrix(string path, int* row_count, int* col_count, int* nnz_count, int* x_coor, int* y_coor, double* vals) {
        string tmp_line;
        ifstream myfile(path);
        if (myfile.is_open())
        {
            list<int> x_coordinates;
            list<int> y_coordinates;
            list<double> values;
            while (getline(myfile, tmp_line))
            {
                string delimiter = "  ";
                int pos = 0;
                string token;
                int token_counter = 0;
                while ((pos = tmp_line.find(delimiter)) <= tmp_line.length()) {
                    token = tmp_line.substr(0, pos);
                    token_counter++;
                    if (token_counter == 0) {
                        x_coordinates.push_back(stoi(token));
                    }
                    else if (token_counter == 1) {
                        y_coordinates.push_back(stoi(token));
                    }
                    else {
                        values.push_back(stod(token));
                    }
                    tmp_line.erase(0, pos + delimiter.length());
                }
                cout << tmp_line << endl;
            }
            myfile.close();
        }

        else cerr << "Unable to open file" << endl;


	}


};

