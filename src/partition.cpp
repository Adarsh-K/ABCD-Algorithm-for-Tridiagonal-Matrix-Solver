#include "iostream"
#include "vector"

using namespace std;


vector<int> get_partition(int size)
{

    int partitions = 0;

    // gre_1107
    if(size == 1107){
        partitions = 9;
    }
    // gre_512
    else if(size == 512){
        partitions = 8;
    }
    // gre_343
    else if(size == 343){
        partitions = 7;
    }
    // gre_216
    else if(size == 216){
        partitions = 6;
    }
    // bayer01
    else if(size == 57735){
        partitions = 15;
    }
    // 1hr36
    else{
        partitions = 16;
    }

    vector<int> index_arr;
    for(int i=0;i<partitions+1;i++)
        index_arr.push_back(i*size/partitions);

    return index_arr;
}
