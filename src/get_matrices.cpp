#include "get_matrices.h"
#include "utils.h"
#include "iostream"
#include "partition.h"
// typedef std::vector<std::vector<double>> stdvecvec;
// typedef std::vector<double> stdvec;

using namespace std;
using namespace arma;

// Get Augmented matrix A_bar
mat get_augmented_A(mat A, vector<int> partition_index){

    // Get interconnections for all the Blocks
    // Correctly working, verified by Siddharth
    vector<pair<int, int>> interconnections = get_interconnections(A, partition_index);

    // for (int i = 0; i < interconnections.size(); i++) {
    //     cout << "Start = " << interconnections[i].first
    //         << ", End = " << interconnections[i].second << endl;
    // }

    int num_row = A.n_rows;
    int num_part = partition_index.size() - 1;
    int part_rows = round(num_row / num_part);
    cout<<"Num rows in partition = "<<part_rows<<endl;

    vector<mat> blocks = get_blocks(A, partition_index);
    mat I = -1. * arma::eye(part_rows, part_rows);
    mat C;

    // Construct rows of C from top to bottom, Block by Block
    for (int i = 0; i < interconnections.size(); i++) { // Last block done after the loop
        int start = interconnections[i].first; // Overlapping columns
        int end = interconnections[i].second;
        // blocks[i] is the upper matrix
        mat C_i = blocks[i].cols(start, end) * blocks[i + 1].cols(start, end).t();
        if(i > 0) {
            C_i.insert_cols(0, I); // -I & C_ij are always next to each other
            if(i > 1){
                C_i.insert_cols(0, part_rows * (i - 1)); // Zeros before -I
            }
        }
        if(i < interconnections.size() - 1) { // not for last interconnection
            C_i.insert_cols(part_rows * (i + 1), num_row - part_rows * (i + 2)); // Zeros after C_ij
        }
        if(C.size() == 0){ // i.e. C isn't initialised
            C = C_i;
        }
        else{
            C.insert_rows(part_rows * i, C_i);
        }
    }

    // C.print();
    // cout<<"Before last block\n";
    // WRONG!! For the final Block
    // I.insert_cols(0, part_rows * (num_part - 1));
    // C.insert_rows(part_rows * (num_part - 1), I);
    
    // Final -I block
    int done_cols = part_rows * (interconnections.size() - 1);
    mat last_aug_block(part_rows, done_cols, arma::fill::zeros);
    last_aug_block.insert_cols(done_cols, I);
    C.insert_rows(part_rows * interconnections.size(), last_aug_block); 
    // C.print();
    // cout<<endl;
    A.insert_cols(num_row, C); // A_bar

    return A;
}

// get the product of the matrix with its pseudo inverse
arma::mat get_pinv_product(arma::mat arma_matrix, int m, int n)
{
    /* arma::mat arma_matrix( mat, m, n, true, false ); */
    arma::mat pseudo_inverse = pinv(arma_matrix);
    return pseudo_inverse*arma_matrix;
}

arma::mat get_abcd(arma::mat A, arma::mat b){

    int num_row = A.n_rows;
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
    mat P = get_p(A_bar, A_bar.n_cols, partition_index);
    mat Y = get_y(A_bar.n_cols-num_row, num_row);
    mat S = get_s(Y, P);

//     // cout<<"--------\n";
//     // P.print();
//     // cout<<P.n_rows<<endl;
//     // cout<<P.n_cols<<endl;
// // cout<<"--------\n";
//     // Y.print();
    

//     // cout<<S.n_rows<<endl;
//     // cout<<S.n_cols<<endl;

//     // make test matrix
    // mat x(A.n_rows, 1, fill::randu);
    // mat b = A*x;

//     // cout<<"--------\n";

//     // get the solution from abcd

    mat w(A_bar.n_cols, 1, fill::zeros);
    mat v(num_row, 1, fill::randu);
    mat tempa;
    mat tempb;
    mat temp = pinv(A)*b;
    temp = temp+0.0001*v;
    for(int i=1; i < partition_index.size();i++){
        tempa = A_bar.submat(partition_index[i-1], 0, partition_index[i]-1, A_bar.n_cols-1);
        tempb = b.submat(partition_index[i-1], 0, partition_index[i]-1, 0);
        w += pinv(tempa)*tempb;
    }
//      mat A_bpinv = pinv(A_bar); 
//      mat w = A_bpinv*b; 
    cout<<"----------------------------------\n";
    cout<<w.n_rows<<endl;
    cout<<w.n_cols<<endl;
    mat f = -Y*w;
    cout<<f.n_rows<<endl;
    cout<<f.n_cols<<endl;

    cout<<"----------------------------------\n";

    mat z = solve(S, f);
    cout<<z.n_rows<<endl;
    cout<<z.n_cols<<endl;
    w=temp;
    mat z_bar = Y.t()*z;
    cout<<z_bar.n_rows<<endl;
    cout<<z_bar.n_cols<<endl;


    // new P
    mat Abartrans = pinv(A_bar).t();
    mat P_new = get_p(Abartrans, Abartrans.n_cols, partition_index);
    mat I = eye(P_new.n_rows, P_new.n_rows);
    mat u = (I-P_new)*z_bar;

    cout<<u.n_rows<<endl;
    cout<<u.n_cols<<endl;

    mat xy = join_cols(w, u);
    cout<<xy.n_rows<<endl;
    cout<<xy.n_cols<<endl;

    mat X = xy.submat(0, 0, num_row-1, 0);

    return X;
//     // cout<<"--------\n";
//     // all 0's => true

    // int error = 0;
    // for(int i=0;i<num_row;i++){
    //     if(ch1[i]!=b[i]){
    //         error += abs(X[i]-x[i]);
    //         cout<<ch1[i]<<" "<<b[i]<<"\n";
    //     }
    // }
//     /* umat fin = (X == x); */
//     /* fin.print(); */
    // cout<<"Error: "<<error;
}


// get the P matrix
arma::mat get_p(arma::mat A, int col, vector<int> partitions)
{
    arma::mat P(col, col, arma::fill::zeros);
    int rows;
    arma::mat temp;
    arma::mat test;
    for(int i=1; i < partitions.size();i++){
        rows = partitions[i]-partitions[i-1]; // 123
        /* double* temp_mat; */
        /* double* cmaj_temp; */
        /* temp_mat = new double[rows*col]; */
        /* cmaj_temp = new double[rows*col]; */

        temp = A.submat(partitions[i-1], 0, partitions[i]-1, col-1); // \bar{A_i}

        // creating the partioned arrays
        /* for(int j=col*partitions[i-1]; j<col*partitions[i];j++) */
        /*     temp_mat[j-col*partitions[i-1]] = A[j]; */

        //converting to column major
        /* int r, c; */
        /* for(int t=0;t<rows*col;t++){ */
        /*     r = int(t/col); */
        /*     c = int(t%col); */
        /*     cmaj_temp[c*rows + r] = temp_mat[t]; */
        /* } */

        // P = get_pinv(cmaj_temp, rows, col);
        test = temp.t()*((temp*temp.t()).i()); // Acc. to the equation

        // P = P + get_pinv_product(temp, rows, col);
        P = P+(test*temp);
        /* P.print(); */
        /* free(temp_mat); */
        /* free(cmaj_temp); */
    }
    return P;
}


// get Y matrix
// Y = [0 I]
arma::mat get_y(int columns_c, int columns_a)
{
    arma::mat O(columns_c, columns_a, fill::zeros);
    arma::mat I;
    I = eye(columns_c, columns_c);

    arma::mat Y = join_rows(O,I);

    return Y;
}


// S = Y(I-P)Y.T
arma::mat get_s(arma::mat Y, arma::mat P)
{
    arma::mat S;
    int rows = P.n_rows; // sqaure matrix

    arma::mat I;
    I = eye(rows, rows);

    S = Y*(I-P)*Y.t();

    return S;
}


// Testing code below

/* int main() */
/* { */
/*     arma::mat B; */

/*     double* mat; */
/*     mat = new double[4*4]; */

/*     int n_rows = 4; */
/*     int n_cols = 4; */

/*     for(int i=0;i<4;i++){ */
/*         for(int j=0;j<4;j++) */
/*             mat[j+4*i] = i*2+j*3; */
/*     } */

/*     double* cmat = new double[4*4]; */

/*     int r, c; */
/*     for(int i=0;i<16;i++){ */
/*         r = int(i/4); */
/*         c = int(i%4); */
/*         cmat[c*4 + r] = mat[i]; */
/*     } */
/*     vector<int> ps{0, 2,4}; */


/*     B = get_p(mat, 4, ps); */


/*     for(int i=0;i<4;i++){ */
/*         for(int j=0;j<4;j++) */
/*             cout<<mat[j+4*i]<<" "; */
/*         cout<<"\n"; */
/*     } */

    /* stdvecvec V(A.n_rows); */
    /* for (size_t i = 0; i < A.n_rows; ++i) { */
    /*      V[i] = arma::conv_to< stdvec >::from(A.row(i)); */
    /*  }; */

    /* for(int i=0;i<5;i++){ */
    /*     for(int i=0;i<5;i++) */
    /*         cout<<V[1][1]<<"\n"; */
    /* } */

    /* return 0; */
/* } */
