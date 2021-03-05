//
// Created by Tirelli Andrea on 03/03/2020.
//

#ifndef NANOPORE_SIGNAL_PROCESSING_HPP
#define NANOPORE_SIGNAL_PROCESSING_HPP

#include <vector>
#include <iostream>


void cut_first_part(std::vector<float> &);

size_t find_first_increasing(const std::vector<float> &vec);

std::vector<size_t> segment_signal(std::vector<float> &, size_t, size_t);

std::vector<float> quantize_signal(std::vector<float> &,const std::vector<size_t> &);//, std::vector<float> *);

float std_(const std::vector<float> &, size_t=0, size_t=0);

float mean_(const std::vector<float> &, size_t=0, size_t=0);

void compute_direction(const std::vector<float>&, size_t*, size_t*, size_t*, float=0);

void cut_subsignal(std::vector<float>&, size_t &, size_t &, size_t &);

std::vector<size_t> select_bkps(std::vector<size_t>&, size_t&, size_t&, size_t&);

void get_kmers_info(const std::vector<float>&, const std::vector<size_t>&,
                    std::vector<float>*, std::vector<float>*);

void find_longest_constant_subsequence(const std::vector<float> &, size_t*, size_t*, float=0);

#endif //NANOPORE_SIGNAL_PROCESSING_HPP
