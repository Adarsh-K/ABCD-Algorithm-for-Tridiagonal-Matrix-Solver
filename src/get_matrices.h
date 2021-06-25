#pragma once

#include "armadillo"
#include "vector"

arma::mat get_augmented_A(arma::mat A, std::vector<int> partition_index);

arma::mat get_p(arma::mat A, int col, std::vector<int> partitions);

arma::mat get_y(int columns_c, int columns_a);

arma::mat get_s(arma::mat Y, arma::mat P);
arma::mat get_abcd(arma::mat A, arma::mat b);