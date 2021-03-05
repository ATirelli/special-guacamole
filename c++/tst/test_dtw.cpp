//
// Created by Tirelli Andrea on 16/05/2020.
//

#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "../src/dtw.hpp"
#include "../src/AllToOneDtw.hpp"
#include <math.h>

using ::testing::_;
using ::testing::Invoke;

TEST(test_compute_KL_matrix, zero_KL_matrix) {
    std::vector<float> means_signal (4, 1);
    std::vector<float> means_kmers (4, 1);
    std::vector<float> stds_signal (4, 3);
    std::vector<float> stds_kmers (4, 3);

    std::vector<std::vector<float>> expected_matrix (4, std::vector<float>{0, 0, 0, 0});
    auto KL_matrix = compute_KL_matrix(means_signal, stds_signal, means_kmers, stds_kmers);
    ASSERT_EQ(expected_matrix, KL_matrix);

}

TEST(test_compute_KL_matrix, non_trivial_matrix) {
    std::vector<float> means_signal (4, 1);
    std::vector<float> means_kmers (3, 1);
    std::vector<float> stds_signal (4, 3);
    std::vector<float> stds_kmers (3, 3);

    size_t expected_n_rows {4}, expected_n_columns {3};
    auto KL_matrix = compute_KL_matrix(means_signal, stds_signal, means_kmers, stds_kmers);
    ASSERT_EQ(expected_n_rows, KL_matrix.size());
    ASSERT_EQ(expected_n_columns, KL_matrix[0].size());

}

TEST(choose_barcode, trivial_solution) {

    //body of the test
    std::vector<size_t> signal {0,1,2,3};
    std::vector<std::vector<int>> kmer_barcode_matrix (3, std::vector<int> {1, 2, 3, 3});
    std::vector<std::vector<float>> KL_matrix (4, std::vector<float> {1, 1, 1, 1});
    std::vector<int> start_indices {0, 4, 4};
    size_t  index_barcode {0};
    float distance {0};
    choose_barcode(signal, kmer_barcode_matrix, KL_matrix, start_indices, &index_barcode, &distance);
    EXPECT_FLOAT_EQ(4, distance);
    EXPECT_EQ(0, index_barcode);

}

TEST(choose_barcode, all_different) {

    //body of the test
    std::vector<size_t> signal {0,1,2,3};
    size_t n_barcodes {3};
    std::vector<std::vector<int>> kmer_barcode_matrix;
    for (int i {0}; i<n_barcodes;i++) {
        kmer_barcode_matrix.push_back(std::vector<int> (3, n_barcodes-i));
    }
    std::vector<std::vector<float>> KL_matrix (4, std::vector<float> {1, 2, 3, 4});
    std::vector<int> start_indices {0, 0, 0};
    size_t  index_barcode {0};
    float distance {0};
    choose_barcode(signal, kmer_barcode_matrix, KL_matrix, start_indices, &index_barcode, &distance);
    EXPECT_FLOAT_EQ(6, distance);
    EXPECT_EQ(2, index_barcode);

}

TEST(test_kl_gauss, zero_stds) {
    float mean1 {1}, mean2 {3}, std1{0}, std2{0};
    auto kl = KL_gauss(mean1, std1, mean2, std2);
    ASSERT_EQ(kl, std::numeric_limits<float>::infinity());
}

TEST(test_kl_gauss, non_zero_stds) {
    float mean1 {1}, mean2 {3}, std1{1}, std2{1};
    auto KL = KL_gauss(mean1, std1, mean2, std2);
    float expected_KL {2};
    ASSERT_EQ(expected_KL, KL);
}

TEST(test_kl_gauss, equal_dist) {
    float mean1 {1}, mean2 {1}, std1{1}, std2{1};
    auto KL = KL_gauss(mean1, std1, mean2, std2);
    float expected_KL {0};
    ASSERT_EQ(expected_KL, KL);
}

TEST(test_kl_gauss, different_stds) {
    float mean1 {1}, mean2 {3}, std1{1}, std2{2};
    auto KL = KL_gauss(mean1, std1, mean2, std2);
    float expected_KL = log(2)+0.625-0.5;
    ASSERT_EQ(expected_KL, KL);
}

TEST(test_compute_n_common_entries, equal_vectors) {

    std::vector<int> v1 {1, 2, 3, 4, 5};
    std::vector<int> v2 = v1;
    auto n = compute_number_common_entries(v1, v2);
    ASSERT_EQ(n, 5);
}

TEST(test_compute_n_common_entries, different_from_start) {

    std::vector<int> v1 {1, 2, 3, 4, 5};
    std::vector<int> v2 {5, 4, 3, 2, 1};
    auto n = compute_number_common_entries(v1, v2);
    ASSERT_EQ(n, 0);
}

TEST(test_compute_n_common_entries, different_len_vectors) {

    std::vector<int> v1 {1, 2, 3, 4, 5};
    std::vector<int> v2 {1, 2, 3, 4, 5, 5, 4, 3, 2, 1};
    auto n = compute_number_common_entries(v1, v2);
    ASSERT_EQ(n, 5);
}

TEST(test_compute_start_indices, equal_vectors) {
    std::vector<int> v1 {1, 2, 3, 4, 5};
    std::vector<std::vector<int>> matrix {v1, v1};
    std::vector<int> expected_indices {0, 5};
    auto start_indices = compute_start_indices(matrix);

    ASSERT_EQ(expected_indices, start_indices);
}

TEST(test_compute_start_indices, easy_test) {
    std::vector<int> v1 {1, 2, 3, 4, 5};
    std::vector<int> v2 {1, 2, 4, 4, 5};
    std::vector<std::vector<int>> matrix {v1, v2};
    std::vector<int> expected_indices {0, 2};
    auto start_indices = compute_start_indices(matrix);

    ASSERT_EQ(expected_indices, start_indices);
}












