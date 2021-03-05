//
// Created by Tirelli Andrea on 24/03/2020.
//
#include <limits>
#include <vector>
#include <algorithm>
#include <iostream>
#include "AllToOneDtw.hpp"

AllToOneDtw::AllToOneDtw(const std::vector<std::vector<int>>& kmer_barcode_matrix, const std::vector<size_t>& seq,
                         const std::vector<std::vector<float>>& KL_matrix)  {
    this->kmer_barcode_matrix = kmer_barcode_matrix;
    this->seq = seq;
    this->len_subseq = kmer_barcode_matrix.at(0).size();
    this->len_seq = seq.size();
    this->KL_matrix = KL_matrix;
    for (size_t i{0}; i<len_subseq+1; i++) {
        cost_matrix.emplace_back(len_seq+1, std::numeric_limits<float>::infinity());
    }
    cost_matrix.at(0) = std::vector<float>(len_seq + 1, 0);
}

void AllToOneDtw::compute_accumulated_cost_matrix(int index, int &start_index) {

    for (int i{start_index}; i < len_subseq; i++) {
        for (size_t j{0}; j < len_seq; j++) {
            auto val = std::min(cost_matrix.at(i).at(j + 1),
                                std::min(cost_matrix.at(i + 1).at(j),cost_matrix.at(i).at(j)));
            cost_matrix.at(i + 1).at(j + 1) = KL_matrix.at(j).at(kmer_barcode_matrix.at(index)[i]) + val;
        }
    }

}

std::vector<std::vector<float>> AllToOneDtw::get_cost_matrix() {
    return cost_matrix;
}
float AllToOneDtw::get_distance() {
    return  *std::min_element(cost_matrix.at(len_subseq).begin(), cost_matrix.at(len_subseq).end());
}
std::vector<std::pair<size_t, size_t>> AllToOneDtw::return_subsequence_path() {
    size_t dim1{cost_matrix.size()}, dim2{cost_matrix.at(0).size()};
    size_t rep1{0}, rep2{0}, end_vertical_path{0}, start_vertical_path{0}, current_horizontal_position{dim1-1};
    std::vector<float> last_row = cost_matrix.at(dim1-1);
    std::vector<std::pair<size_t, size_t>> path;

    end_vertical_path = std::distance(last_row.begin(),
                                      std::min_element(last_row.begin()+dim2-2, last_row.end()))+dim2-2;

    start_vertical_path = end_vertical_path;

    while (current_horizontal_position!=1 and start_vertical_path!=1) {
        std::vector<float> v {cost_matrix.at(current_horizontal_position-1).at(start_vertical_path-1),
                              cost_matrix.at(current_horizontal_position-1).at(start_vertical_path),
                              cost_matrix.at(current_horizontal_position).at(start_vertical_path-1)};

        auto min_pos {std::distance(v.begin(), std::min_element(v.begin(), v.end()))};
        if (min_pos==0) {
            current_horizontal_position--;
            start_vertical_path--;
            path.emplace_back(std::make_pair(current_horizontal_position-1, start_vertical_path-1));

        }
        else if (min_pos==1) {
            current_horizontal_position--;
            rep2++;
            path.emplace_back(std::make_pair(current_horizontal_position-1, start_vertical_path));
        }
        else {
            start_vertical_path--;
            rep1++;
            path.emplace_back(std::make_pair(current_horizontal_position, start_vertical_path-1));
        }

    }
    start = start_vertical_path;
    end = end_vertical_path;
    reps = std::max(rep1, rep2);

    std::reverse(path.begin(), path.end());
    return path;

}

void AllToOneDtw::compute_all_to_one(std::vector<int> &start_indices, size_t *index_barcode, float*distance) {
    size_t size1 {this->kmer_barcode_matrix.size()};
    std::vector<float> dtw (size1, std::numeric_limits<float>::infinity());

    //first barcode
    compute_accumulated_cost_matrix(0, start_indices.at(0));
    dtw.at(0)=get_distance();
    auto min_dist = get_distance();
    int min_index = 0;
    int k=1;
    //all the others
    while (k<size1) {
        for (size_t j{0}; j < len_seq; j++) {
            auto val = std::min(cost_matrix.at(start_indices[k]).at(j + 1),
                                std::min(cost_matrix.at(start_indices[k] + 1).at(j),cost_matrix.at(start_indices[k]).at(j)));
            cost_matrix.at(start_indices[k] + 1).at(j + 1) = KL_matrix.at(j).at(kmer_barcode_matrix.at(k)[start_indices[k]]) + val;
        }
        if (*std::min_element(cost_matrix.at(start_indices[k] + 1).begin(), cost_matrix.at(start_indices[k] + 1).end())>=min_dist) {
            auto l {k};
            while (start_indices[k]>=start_indices[l]) {k++;}
        }
        else {
            for (int i{start_indices[k] + 1}; i < len_subseq; i++) {
                for (size_t j{0}; j < len_seq; j++) {
                    auto val = std::min(cost_matrix.at(i).at(j + 1),
                                        std::min(cost_matrix.at(i + 1).at(j), cost_matrix.at(i).at(j)));
                    cost_matrix.at(i + 1).at(j + 1) = KL_matrix.at(j).at(kmer_barcode_matrix.at(k)[i]) + val;
                }
            }
            auto temp_min_dist = get_distance();
            if (temp_min_dist<min_dist) {min_dist = temp_min_dist; min_index = k;}
            k++;
        }
        *index_barcode = min_index;
        *distance = min_dist;

    }
}

