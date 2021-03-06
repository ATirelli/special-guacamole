//
// Created by Tirelli Andrea on 06/03/2020.
//
#include <algorithm> 
#include <iostream>
#include <vector>
#include <cmath>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include "dtw.hpp"
using namespace boost::accumulators;

float std_(const std::vector<float> &vec, size_t start, size_t end) {
    if (end == 0 and start==0) {end=vec.size()-1;}
    if (start==end) {
        return 0;
    }
    else {
        accumulator_set<double, stats<tag::variance>> acc;
        for (size_t i{start}; i <=end; i++) { acc(vec.at(i)); }
        return sqrt(variance(acc));
    }
}

float mean_(const std::vector<float> &vec, size_t start, size_t end) {
    if (end == 0 and start==0) {end=vec.size()-1;}
    if (start == end) {
        return vec[start];
    }
    else {
    accumulator_set<double, stats<tag::mean>> acc;
    for (size_t i {start}; i<=end;i++ ) {acc(vec.at(i));}
    return mean(acc);
    }
}

size_t find_first_increasing(const std::vector<float> &vec) {
    bool found = false;
    size_t i {0}, k {0}, len_vec {vec.size()};
    while (!found and i<len_vec-1) {
        if (vec.at(i+1) > vec.at(i)) {
            k = i;
            found = true;

        }
        i++;
    }
    return k;
}

void cut_first_part(std::vector<float> &signal) {
    size_t one_third = signal.size()/3;
    int max_height = std::distance(signal.begin(), std::max_element(signal.begin(), signal.begin()+one_third));
    signal.erase(signal.begin(), signal.begin()+max_height);
    size_t start {find_first_increasing(signal)};
    signal.erase(signal.begin(), signal.begin()+start);

}

std::vector<size_t> segment_signal(std::vector<float> &signal, size_t n_bkps, size_t min_len) {
    size_t len_signal {signal.size()};
    size_t multiple {len_signal%min_len};
    float last_value {signal.at(len_signal-1)};
    if (multiple!=0) {
        for (auto i {0}; i<min_len-multiple; i++) { signal.push_back(last_value);}
        len_signal += min_len-multiple;
    }

    size_t half_len {len_signal/min_len};
    std::vector<std::tuple<size_t, size_t>> segments (half_len, std::make_pair(0, 0));
    size_t j {0};
    size_t i {0};
    while (j+min_len-1<len_signal) {
        segments.at(i) = std::make_pair(j, j+min_len-1);
        i++;
        j+=min_len;
    }

    std::vector<float> distances;

    for (size_t i {0};i<half_len-1;i++) {
        distances.push_back(std_(signal,
                           std::get<0>(segments.at(i)),
                            std::get<1>(segments.at(i + 1))));

    }

    size_t n_segs{half_len};
    while (n_segs>n_bkps) {
        int min_index = std::distance(distances.begin(), std::min_element(distances.begin(), distances.end()));
        std::get<1>(segments.at(min_index)) = std::get<1>(segments.at(min_index+1));
        segments.erase(segments.begin()+min_index+1);
        if (min_index!=0) {
            distances.at(min_index - 1) = std_(signal,
                                               std::get<0>(segments.at(min_index - 1)),
                                               std::get<1>(segments.at(min_index)));
        }
        if (min_index!=distances.size()-1) {
            distances.at(min_index) = std_(signal,
                                           std::get<0>(segments.at(min_index)),
                                           std::get<1>(segments.at(min_index + 1)));
            distances.erase(distances.begin()+min_index+1);
        }
        else {
            distances.erase(distances.begin()+min_index);
        }

        n_segs--;
    }
    std::vector<size_t> bkps;
    for (auto val: segments) {bkps.push_back(std::get<0>(val));}
    if (bkps[bkps.size()-1]!=len_signal-1) {bkps.push_back(len_signal-1);}
    return bkps;
}

void find_longest_constant_subsequence(const std::vector<float> &segments, size_t* startptr, size_t* endptr, float th) {
    size_t i{0}, j{0}, start{0}, end{0}, longest{0};
    size_t len_seg {segments.size()};

    while (i<len_seg-1) {
        j = i+1;
        while (std::fabs(segments.at(i)-segments.at(j))<=th and j<len_seg-1) {j++;}

        if (j-i >= longest) {
            start = i;
            end = j;
            longest = j-1;

        }
        i = j;
    }
    *startptr = start;
    *endptr = end-1;
}

std::vector<float> quantize_signal(std::vector<float> &signal,const std::vector<size_t> &bkps) {
    std::vector<float> quantized;

    size_t dwell {0}, start{0}, end{0};
    float mean_sig {0};
    for (size_t i{0};i<bkps.size()-1; i++) {
        start = bkps.at(i);
        end = bkps.at((i+1))+(i==bkps.size()-2?0:-1);
        dwell = end-start+1;
        mean_sig = mean_(signal, start, end);
        for (size_t j{0}; j< dwell; j++) {quantized.push_back(mean_sig);}
    }
    return quantized;
}

void compute_direction(const std::vector<float>& signal, size_t* start, size_t* end, size_t* direction, float th) {
    find_longest_constant_subsequence(signal, start, end, th);
    if (*start < signal.size()-*end) {
        *direction = 0;
    }
    else {*direction = 1;}
}

void cut_subsignal(std::vector<float> &signal, size_t &strand, size_t &start, size_t &end) {
    if (strand == 0) {
        signal.erase(signal.begin()+start,signal.end());
    }
    else {
        signal.erase(signal.begin(),signal.begin()+end+1);
    }
}

std::vector<size_t> select_bkps(std::vector<size_t>& bkps, size_t& strand,size_t& start,size_t& end) {
    std::vector<size_t> relevant_bkps;
    size_t i {0};
    if (strand == 0) {

        while (bkps.at(i) < start) {relevant_bkps.push_back(bkps.at(i)); i++;}
        if (relevant_bkps[relevant_bkps.back()]!=start-1) {
            relevant_bkps.push_back(start - 1);
        }
    }
    else {
        i = bkps.size()-1;
        while (bkps.at(i) > end) { relevant_bkps.push_back(bkps.at(i)-end-1); i--;}
        if (relevant_bkps[relevant_bkps.size()-1]!=0) {
            relevant_bkps.push_back(0);
        }
        std::reverse(relevant_bkps.begin(), relevant_bkps.end());
    }
    return relevant_bkps;
}

void get_kmers_info(const std::vector<float>& signal,const  std::vector<size_t>& bkps,
        std::vector<float>* means, std::vector<float>* stds) {
    for (size_t i{0};i<bkps.size()-1; i++) {
        (*means).push_back(mean_(signal, bkps.at(i), bkps.at(i+1)-(i+1==bkps.size()-1?0:1)));
        (*stds).push_back(std_(signal, bkps.at(i), bkps.at(i+1)-(i+1==bkps.size()-1?0:1)));
    }
}
