//
// Created by Tirelli Andrea on 09/03/2020.
//
#include <math.h>
#include <vector>
#include "signal_processing.hpp"
#include "AllToOneDtw.hpp"
#include <limits>
#include "dtw.hpp"

std::vector<std::vector<float>> compute_KL_matrix(const std::vector<float> &signal_means,
                                                  const std::vector<float> &signal_stds,
                                                  const std::vector<float> &all_kmer_means,
                                                  const std::vector<float> &all_kmer_stds) {
    std::vector<std::vector<float>> KL_matrix;
    size_t dim1 {signal_means.size()}, dim2 {all_kmer_means.size()};
    for (size_t i{0}; i<dim1;i++) {
        std::vector<float> row;
        for (size_t j{0}; j<dim2;j++) {
            row.push_back(KL_gauss(signal_means[i], signal_stds[i], all_kmer_means[j], all_kmer_stds[j]));
        }
        KL_matrix.push_back(row);
    }
    return KL_matrix;
}

float KL_gauss(const float &mean1, const float &std1, const float &mean2, const float &std2) {
    if (std1 == 0 or std2 == 0) {
        return std::numeric_limits<float>::infinity();
    } else {
        float first_piece {log(std2 / std1)};
        float num = pow(std1, 2) + pow(mean1 - mean2, 2);
        float den = 2 * pow(std2, 2);
        return first_piece + num / den - 0.5;
    }
}

float KL_divergence_between_vectors(const std::vector<float> &v1, const std::vector<float> &v2,
                                    size_t start1, size_t end1, size_t start2, size_t end2) {
    return KL_gauss(mean_(v1, start1, end1), std_(v1, start1, end1),
                    mean_(v2, start2, end2), std_(v2, start2, end2));
}

void choose_barcode(std::vector<size_t>& signal, std::vector<std::vector<int>>& kmer_barcode_matrix,
                    std::vector<std::vector<float>>& KL_matrix, std::vector<int>& start_indices,
                    size_t *index_barcode, float*distance) {

    AllToOneDtw d = AllToOneDtw(kmer_barcode_matrix, signal, KL_matrix);

    d.compute_all_to_one(start_indices, index_barcode, distance);
}

int compute_number_common_entries(std::vector<int> &v1, std::vector<int> &v2) {
    size_t i {0};
    size_t min_len = fmin(v1.size(), v2.size());
    while(v1[i] == v2[i] and i< min_len) {i+=1;}
    return i;
}

std::vector<int> compute_start_indices(std::vector<std::vector<int>>&kmer_barcode_matrix) {

    std::vector<int> start_indices;
    start_indices.push_back(0);
    for (size_t i{1};i<kmer_barcode_matrix.size(); i++) {
        start_indices.push_back(compute_number_common_entries(kmer_barcode_matrix.at(i),
                                                              kmer_barcode_matrix.at(i-1)));
    }
    return start_indices;
}
