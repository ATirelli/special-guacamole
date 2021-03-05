//
// Created by Andrea Tirelli on 02/05/2020.
//

#include "gtest/gtest.h"
#include "../src/utils.hpp"
#include <string>

TEST(rev_comp_test, easy_tests) {
    //arrange
    //act
    //assert
    EXPECT_EQ (reverse_complement("AAA"), "TTT");
    EXPECT_EQ (reverse_complement("CCCCCC"), "GGGGGG");

}
TEST(rej_cond_test, test1) {
    //arrange
    //act
    //assert
    std::string s {"Hello"};
    size_t n_kmers, strand, len_read,start, end, upper_bound,lower_bound;
    n_kmers = 50;
    strand = 0;
    len_read = 500;
    start = 0;
    end = 50;
    upper_bound = 200;
    lower_bound = 40;
    bool rej_cond = compute_rejection_condition(n_kmers, strand, len_read,start,
                           end, upper_bound,lower_bound);
    EXPECT_EQ(rej_cond, false);



}