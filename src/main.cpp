#include "partition.h"
#include "get_matrices.h"
#include "read_matrix.h"
#include "utils.h"
#include "iostream"
#include "vector"
#include "armadillo"
#include "cstdlib"
#include "math.h"
#include "stdlib.h"

using namespace std;
using namespace arma;





int main() {
    double* matrix = get_data("../matrix_files/gre_512.mtx");

    // mat x(10, 1, fill::randu);
    // x = 0.0001465*x;
 // 8 x 8 sample matrix
    // [1,    3,  100, 0,   0,   0,    0,   0]
    // [2,   -2, -200, 0,   0,   0,    0,   0]
    // --------------------------------------
    // [0,   10,   -1, 4,  -9,   0,    0,   0]
    // [0, -0.3,    1, 2, 0.6,   0,    0,   0]
    // --------------------------------------
    // [0,    0,    0, 2,   9,  17,    8,   0]
    // [0,    0,    0, 2, -10, 0.1,    4,   0]
    // ---------------------------------------
    // [0,    0,    0, 0,   0, 0.5, -100, 2.3]
    // [0,    0,    0, 0,   0, 0.5,   -1, -34]
    // int num_row = 8;
    // double *matrix = new double[num_row * num_row];
    // std::fill(matrix, matrix + num_row * num_row, 0.);
    // matrix[0] = 1.;
    // matrix[1] = 2.;
    // matrix[8] = 3.;
    // matrix[9] = -2.;
    // matrix[10] = 10.;
    // matrix[11] = -0.3;
    // matrix[16] = 100;
    // matrix[17] = -200;
    // matrix[18] = -1.;
    // matrix[19] = 1.;
    // matrix[26] = 4.;
    // matrix[27] = 2.;
    // matrix[28] = 2.;
    // matrix[29] = 2.;
    // matrix[34] = -9.;
    // matrix[35] = 0.6;
    // matrix[36] = 9.;
    // matrix[37] = -10.;
    // matrix[44] = 17.;
    // matrix[45] = 0.1;
    // matrix[46] = 0.5;
    // matrix[47] = 0.5;
    // matrix[52] = 8.;
    // matrix[53] = 4.;
    // matrix[54] = -100.;
    // matrix[55] = -1.;
    // matrix[62] = 2.3;
    // matrix[63] = -34;
    // mat A(matrix, num_row, num_row, true, false);
    // A.print();
    // vector<int> partition_index{0, 2, 4, 6, num_row}; // first_row in .submat()
    // cout<<"Hekko";

    // hardcoding the size for now (will take it as input later along with file name)
    int num_row = 512;

    mat A(matrix, num_row, num_row, true, false);

    cout <<"Num rows in A = "<<A.n_rows<<endl;
    // get indices of the partition
    vector<int> partition_index = get_partition(num_row);
    int num_part = partition_index.size() - 1;
    int part_rows = round(num_row / num_part);

    // get matrices
    mat A_bar = get_augmented_A(A, partition_index);
    cout<<"Num rows in Augmented = "<<A_bar.n_rows<<endl;
    cout<<"Num cols in Augmented = "<<A_bar.n_cols<<endl;

    // A_bar.submat(0, 0, num_row - 1, num_row - 1).print();
    // cout<<endl;
    // A_bar.submat(0, num_row, num_row - 1, 2 * num_row -part_rows - 1).print();

////////////////////////////////////////////////////////////
    // mat P = get_p(A_bar, A_bar.n_cols, partition_index);
//     // cout<<"--------\n";
//     // P.print();
//     // cout<<P.n_rows<<endl;
//     // cout<<P.n_cols<<endl;
    // mat Y = get_y(A_bar.n_cols-num_row, num_row);
// // cout<<"--------\n";
//     // Y.print();
    // mat S = get_s(Y, P);

//     // cout<<S.n_rows<<endl;
//     // cout<<S.n_cols<<endl;

    // make test matrix
    
    mat x(A.n_rows, 1, fill::randu);
    mat b = A*x;
    mat newx = get_abcd(A, b);

//     // cout<<"--------\n";

//     // get the solution from abcd

    // mat w(A_bar.n_cols, 1, fill::zeros);
//     // mat tempa;
//     // mat tempb;
//     // for(int i=1; i < partition_index.size();i++){
//     //     tempa = A_bar.submat(partition_index[i-1], 0, partition_index[i]-1, A_bar.n_cols-1);
//     //     tempb = b.submat(partition_index[i-1], 0, partition_index[i]-1, 0);
//     //     w += pinv(tempa)*tempb;
//     // }
//      mat A_bpinv = pinv(A_bar); 
//      mat w = A_bpinv*b; 
//     // cout<<w.n_rows<<endl;
//     // cout<<w.n_cols<<endl;
//     mat f = -Y*w;

//     // cout<<f.n_rows<<endl;
//     // cout<<f.n_cols<<endl;

//     // cout<<"--------\n";

    // mat z = solve(S, b);
//     // cout<<z.n_rows<<endl;
//     // cout<<z.n_cols<<endl;
    // mat z_bar = Y.t()*z;
//     // cout<<z_bar.n_rows<<endl;
//     // cout<<z_bar.n_cols<<endl;


//     // new P
    // mat Abartrans = pinv(A_bar).t();
    // mat P_new = get_p(Abartrans, Abartrans.n_cols, partition_index);
    // mat I = eye(P_new.n_rows, P_new.n_rows);
    // mat u = (I-P_new)*z_bar;

//     // cout<<u.n_rows<<endl;
//     // cout<<u.n_cols<<endl;

    // mat xy = w+u;
//     // cout<<xy.n_rows<<endl;
//     // cout<<xy.n_cols<<endl;

    // mat X = xy.submat(0, 0, num_row-1, 0);

    // mat ch1 = A*X;


//     // cout<<"--------\n";
//     // all 0's => true

    double error = 0;
    for(int i=0;i<num_row;i++){
            // error += pow(pow(newx[i],2 )+pow(x[i], 2), 0.5);
        error += abs(newx[i]-x[i]);
            // cout<<newx[i]<<" "<<x[i]<<"\n";
    }
//     /* umat fin = (X == x); */
//     /* fin.print(); */
    cout<<"Error: "<<error<<" ";
    free(matrix);
    return 0;
}
