//
// Created by Tirelli Andrea on 09/03/2020.
//
#include <vector>
#include <unordered_map>
#include <string>
#include <map>
#ifndef NANOPORE_DTW_HPP
#define NANOPORE_DTW_HPP

std::vector<std::vector<float>> compute_KL_matrix(const std::vector<float> &,
                                                        const std::vector<float> &,
                                                        const std::vector<float> &,
                                                        const std::vector<float> &);

float KL_gauss(const float &,const float &,const float &,const float &);

float KL_divergence_between_vectors(const std::vector<float> &v1, const std::vector<float> &v2,
                                    size_t start1, size_t end1, size_t start2, size_t end2);

void choose_barcode(std::vector<size_t>&, std::vector<std::vector<int>>&,
                    std::vector<std::vector<float>>&, std::vector<int>&,
                    size_t*, float*);

int compute_number_common_entries(std::vector<int> &, std::vector<int> &);

std::vector<int> compute_start_indices(std::vector<std::vector<int>>&kmer_barcode_matrix);
#endif //NANOPORE_DTW_HPP
