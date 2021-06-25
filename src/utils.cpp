#include "iostream"
#include "vector"
#include "armadillo"

using namespace std;
using namespace arma;

vector<mat> get_blocks(mat A, vector<int> partition_index){
    
    int num_row = A.n_rows;
    int num_part = partition_index.size() - 1;
    
    // Correctly getting Block matrices 
    vector<mat> blocks;
    for (int i = 0; i < num_part; i++) { 
        blocks.push_back(A.submat(partition_index[i], 0, partition_index[i + 1] - 1, num_row - 1));
    }

    // cout<<"Blocks\n\n";
    // for (int i = 0; i < num_row; i++) {
    //     blocks[i].print();
    //     cout<<"-------------------------------------------------------------"<<endl;
    // }

    return blocks;
}

// Get overlapping columns (interconnections) for all blocks of matrix A
// Working correctly
vector<pair<int, int>> get_interconnections(mat A, vector<int> partition_index) {

    // cout<<"Sample Matrix\n\n";
    // A.print(); // sample constructed correctly
    // cout<<"-------------------------------------------------------------"<<endl;

    vector<pair<int, int>> interconnections;
    int num_row = A.n_rows;
    int num_part = partition_index.size() - 1;
    
    vector<mat> blocks = get_blocks(A, partition_index); 

    // Get last non-0 column from upper Block & 1st non-0 column from lower Block
    // Working correctly
    for (int i = 0; i < num_part - 1; i++) { // For each overlap
        int start = 0, end = num_row - 1; // Starts according to the lower block
        for (int j = num_row - 1; j > -1; j--) { // blocks[i] is upper matrix
            if(!blocks[i].col(j).is_zero()) { 
                end = j;
                break;
            } 
        }
        for (int j = 0; j < num_row - 1; j++) { // blocks[i + 1] is lower matrix
            if(!blocks[i + 1].col(j).is_zero()){
                start = j;
                break;
            }
        }

        interconnections.push_back(make_pair(start, end));
    }

    return interconnections;
}