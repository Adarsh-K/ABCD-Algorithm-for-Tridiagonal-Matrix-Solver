#include <fstream>
#include <iostream>

using namespace std;

double* get_data(string file_name)
{
    std::ifstream file(file_name);
    int num_row, num_col, num_lines;

    // Ignore comments headers
    while (file.peek() == '%') file.ignore(2048, '\n');

    // Read number of rows and columns
    file >> num_row>> num_col >> num_lines;

    // Create 2D array and fill with zeros
    /* num_lines = 277774; */
    /* num_col = 57735; */
    /* num_row = 57735; */

    double* matrix;
    matrix = new double[num_row * num_col];
    // cout<<"hello";
    std::fill(matrix, matrix + num_row *num_col, 0.);
    // cout<<"hello";

    // fill the matrix with data
    for (int l = 0; l < num_lines; l++)
    {
        double data;
        int row, col;
        file >> row >> col >> data;
        matrix[(row -1) + (col -1) * num_col] = data; // changed to column major
    }

    file.close();

    return matrix;
}
