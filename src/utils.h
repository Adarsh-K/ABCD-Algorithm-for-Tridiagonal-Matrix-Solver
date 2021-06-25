#pragma once

#include "armadillo"

std::vector<std::pair<int, int>> get_interconnections(arma::mat A, std::vector<int> partition_index);

std::vector<arma::mat> get_blocks(arma::mat A, std::vector<int> partition_index);