//
// Created by Tirelli Andrea on 02/05/2020.
//

#include "gtest/gtest.h"
#include "../src/signal_processing.hpp"
#include <string>


TEST(test_std_mean, prefefined_extremes) {
    std::vector<float> v {1, 1, 1, 1, 1};
    auto std = std_(v);
    auto mean = mean_(v);
    ASSERT_FLOAT_EQ(0, std);
    ASSERT_FLOAT_EQ(1, mean);
}

TEST(test_std_mean, non_constant_vector) {
    std::vector<float> v {1, 0, 1, 0, 1, 0};
    auto std = std_(v);
    auto mean = mean_(v);
    ASSERT_FLOAT_EQ(0.5, std);
    ASSERT_FLOAT_EQ(0.5, mean);
}

TEST(test_std_mean, non_standard_extremes) {
    std::vector<float> v {1, 0, 1, 0, 1, 0};
    auto std = std_(v, 1, 3);
    auto mean = mean_(v, 1, 3);
    ASSERT_FLOAT_EQ(0.4714045207910317, std);
    ASSERT_FLOAT_EQ(0.3333333333333333, mean);
}

TEST(long_const_test, no_th) {

    std::vector<float> v {1, 2, 2, 2, 5};
    size_t start{0}, end{0};
    find_longest_constant_subsequence(v, &start, &end);

    EXPECT_EQ (start,1);
    EXPECT_EQ (end, 3);

}

TEST(long_const_test, with_th) {

    std::vector<float> v{1, 5, 1, 1.2, 1.1, 0.9, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6};
    size_t start{0}, end{0};
    find_longest_constant_subsequence(v, &start, &end, 0.3);

    EXPECT_EQ (start, 2);
    EXPECT_EQ (end, 5);
}

TEST(test_find_first_increasing, test1) {
    std::vector<float> v {5, 4, 3, 2, 1, 3, 4, 7};
    int i;
    i = find_first_increasing(v);
    ASSERT_EQ(i, 4);
}

TEST(test_cut_first_part, strictly_increasing) {
    std::vector<float> v {5, 4,10, 20, 40, 39, 38, 25,  3, 2, 1, 3, 4, 7, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    std::vector<float> expected {1, 3, 4, 7, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    cut_first_part(v);
    ASSERT_EQ(v, expected);
}

TEST(test_cut_first_part, non_decreasing) {
    std::vector<float> v {5, 4,10, 20, 40, 39, 38, 25,  3, 2, 1, 1,  3, 4, 7, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    std::vector<float> expected {1,3, 4, 7, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    cut_first_part(v);
    ASSERT_EQ(v, expected);
}

TEST(test_segment_signal,easy_segmentation) {
    std::vector<float> v {9, 9, 1, 1};
    std::vector<size_t> expected_bkps {0, 2, 3};
    auto bkps = segment_signal(v, 2, 2);
    ASSERT_EQ(bkps, expected_bkps);
}

TEST(test_segment_signal,longer_segmentation) {
    std::vector<float> v {9, 9, 9, 3,3,3, 1, 1};
    std::vector<size_t> expected_bkps {0, 3, 6, 7};
    auto bkps = segment_signal(v, 3, 1);
    ASSERT_EQ(bkps, expected_bkps);
}

TEST(test_segment_signal,difficult_segmentation) {
    std::vector<float> v {9, 9.1, 8.9, 8.8, 9, 9, -3,-2.5,-3, 1, 1.6, 1.2, 0.9, 1};
    std::vector<size_t> expected_bkps {0, 6, 9, 13};
    auto bkps = segment_signal(v, 3, 1);
    ASSERT_EQ(bkps, expected_bkps);
}

TEST(test_quantize_signal,easy_quantization) {
    std::vector<float> v {2, 2, 1, 1};
    std::vector<size_t> bkps {0, 2, 3};
    auto quantized = quantize_signal(v, bkps);
    ASSERT_EQ(v, quantized);
}

TEST(test_quantize_signal,non_trivial_quantization) {
    std::vector<float> v {3, 2,3, 2,  1, 1.5, 9, 10, 9};
    std::vector<size_t> expected_bkps={0, 4, 6, 8};
    std::vector<float> expected_quantized {2.5, 2.5,2.5, 2.5,  1.25, 1.25, 9.3333333, 9.3333333, 9.3333333};
    auto quantized = quantize_signal(v, expected_bkps);
    ASSERT_EQ(expected_quantized, quantized);
}

TEST(test_compute_direction, forward_direction) {
    std::vector<float> quantized {1, 2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 5,4,5,4,5,4,5,4,5,4,5,4,5,4,5,4,5,4,5,4,5,4,5,4,5,4,5};
    size_t start, end, direction;
    compute_direction(quantized, &start, &end, &direction);
    size_t expected_start = 4;
    size_t expected_end = 10;
    size_t expected_direction = 0;

    ASSERT_EQ(expected_start, start);
    ASSERT_EQ(expected_end, end);
    ASSERT_EQ(expected_direction, direction);
}

TEST(test_compute_direction, backward_direction) {
    std::vector<float> quantized {1, 2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 5,4,5,4,5,4,5,4,5,4,5,4,5,4,5,4,5,4,5,4,5,4,5,4,5,4,5};
    std::reverse(quantized.begin(), quantized.end());
    size_t start, end, direction;
    compute_direction(quantized, &start, &end, &direction);
    size_t expected_start = 27;
    size_t expected_end = 33;
    size_t expected_direction = 1;

    ASSERT_EQ(expected_start, start);
    ASSERT_EQ(expected_end, end);
    ASSERT_EQ(expected_direction, direction);
}

TEST(test_compute_direction, compute_direction_with_th) {
    std::vector<float> quantized {1, 5, 1, 1.2, 1.1, 0.9, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6};
    size_t start, end, direction;
    float th=0.3;
    compute_direction(quantized, &start, &end, &direction, th);
    size_t expected_start = 2;
    size_t expected_end = 5;
    size_t expected_direction = 0;

    ASSERT_EQ(expected_start, start);
    ASSERT_EQ(expected_end, end);
    ASSERT_EQ(expected_direction, direction);
}

TEST(test_cut_subsignal, forward_direction) {
    std::vector<float> signal {1,2,4,5,4,5,4,5,4,5,4,5,4};
    size_t strand {0}, start {3}, end {6};
    cut_subsignal(signal, strand, start, end);
    std::vector<float> expected_cut_signal {1, 2, 4};

    ASSERT_EQ(expected_cut_signal, signal);
}

TEST(test_cut_subsignal, backward_direction) {
    std::vector<float> signal {1,2,4,5,4,5,4,5,4,5,4,5,4};
    size_t strand {1}, start {3}, end {6};
    cut_subsignal(signal, strand, start, end);
    std::vector<float> expected_cut_signal {5,4,5,4,5,4};

    ASSERT_EQ(expected_cut_signal, signal);
}

TEST(test_relevant_bkps, forward_direction_simple) {
    std::vector<size_t> bkps {0, 3, 7, 9, 12, 19};
    size_t strand {0}, start{8}, end {17};
    auto relevant_bkps = select_bkps(bkps, strand, start, end);
    std::vector<size_t> expected {0, 3, 7};
    ASSERT_EQ(expected, relevant_bkps);
}

TEST(test_relevant_bkps, forward_direction_non_trivial) {
    std::vector<size_t> bkps {0, 3, 4, 5, 9, 12, 19};
    size_t strand {0}, start{8}, end {17};
    auto relevant_bkps = select_bkps(bkps, strand, start, end);
    std::vector<size_t> expected {0, 3, 4, 5, 7};
    ASSERT_EQ(expected, relevant_bkps);
}

TEST(test_relevant_bkps, backward_direction_easy) {
    std::vector<size_t> bkps {0, 3, 4, 5, 9, 12, 19};
    size_t strand {1}, start{8}, end {15};
    auto relevant_bkps = select_bkps(bkps, strand, start, end);
    std::vector<size_t> expected {0, 3};
    ASSERT_EQ(expected, relevant_bkps);
}

TEST(test_relevant_bkps, backward_direction_non_trivial) {
    std::vector<size_t> bkps {0, 3, 4, 5, 9, 12, 19, 22, 27};
    size_t strand {1}, start{8}, end {12};
    auto relevant_bkps = select_bkps(bkps, strand, start, end);
    std::vector<size_t> expected {0, 6, 9, 14};
    ASSERT_EQ(expected, relevant_bkps);
}

TEST(test_get_kmers_info, trivial) {
    std::vector<float> signal {1, 1, 1, 1, 3, 3, 3, 3, 3, 5, 5};
    std::vector<size_t> bkps {0, 4, 9, 10};
    std::vector<float> expexted_means {1, 3, 5};
    std::vector<float> expexted_stds (3, 0);

    std::vector<float> means, stds;
    get_kmers_info(signal, bkps, &means, &stds);

    ASSERT_EQ(expexted_means, means);
    ASSERT_EQ(expexted_stds, stds);
}

TEST(test_get_kmers_info, with_bkps_search) {
    std::vector<float> signal {1, 1, 1, 1, 3, 3, 3, 3, 3, 5, 5};
    std::vector<size_t> bkps =segment_signal(signal, 3, 1);
    std::vector<float> expexted_means {1, 3, 5};
    std::vector<float> expexted_stds (3, 0);

    std::vector<float> means, stds;
    get_kmers_info(signal, bkps, &means, &stds);

    ASSERT_EQ(expexted_means, means);
    ASSERT_EQ(expexted_stds, stds);
}