//
// Created by Tirelli Andrea on 24/03/2020.
//
#include <vector>
#include <unordered_map>
#include <string>
#ifndef NANOPORE_ALLTOONEDTW_HPP
#define NANOPORE_ALLTOONEDTW_HPP

class AllToOneDtw {

private:
    std::vector<std::vector<float>> cost_matrix;
    std::vector<std::vector<int>> kmer_barcode_matrix;
    std::vector<size_t> seq;
    size_t len_subseq;
    size_t len_seq;
    std::vector<std::vector<float>> KL_matrix;
    size_t start {0}, end{0}, reps{0};

public:
    AllToOneDtw(const std::vector<std::vector<int>>&, const std::vector<size_t>&,const std::vector<std::vector<float>>&);

    void compute_accumulated_cost_matrix(int, int &);

    std::vector<std::pair<size_t, size_t>> return_subsequence_path();

    std::vector<std::vector<float>> get_cost_matrix();
    float get_distance();

    void compute_all_to_one(std::vector<int> &, size_t *, float*);

};


#endif //NANOPORE_ALLTOONEDTW_HPP
