//
// Created by Tirelli Andrea on 05/03/2020.
//

#include <string.h>
#include <algorithm>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include "utils.hpp"
#include "SmithWaterman.hpp"

std::vector<std::vector<std::string>> read_tsv(std::string fname, char sep) {
    std::vector<std::vector<std::string> > items;
    std::ifstream ifs(fname);
    if (ifs.fail()) {
        std::cerr << "error" << std::endl;
        exit(-1);
    }
    std::string line;
    while (getline(ifs, line)) {
        std::stringstream ss(line);
        std::vector<std::string> item;
        std::string tmp;
        while(getline(ss, tmp, sep)) {
            item.push_back(tmp);
        }
        items.push_back(item);
    }
    return items;
}

std::unordered_map<std::string, std::unordered_map<std::string, float>> get_hash_kmers(std::vector<std::vector<std::string>> &kmers_vec) {
    std::unordered_map<std::string, std::unordered_map<std::string, float>>hash_kmers{};
    for (const auto& row: kmers_vec) {
        std::unordered_map<std::string, float> stats;
        std::string kmer {row.at(0)};
        stats.insert(std::make_pair("level_mean", stof(row.at(1))));
        stats.insert(std::make_pair("level_std", stof(row.at(2))));
        stats.insert(std::make_pair("sd_mean", stof(row.at(3))));
        stats.insert(std::make_pair("sd_std", stof(row.at(4))));
        hash_kmers.insert(make_pair(kmer, stats));
    }
    return hash_kmers;
}

std::unordered_map<std::string, size_t> get_pos_hash_vector(std::vector<std::string>&vec) {
    std::unordered_map<std::string, size_t> hash {};
    for (size_t i{0}; i<vec.size(); i++) {
        hash.insert(std::make_pair(vec.at(i), i));
    }
    return hash;
}

std::unordered_map<size_t, std::map<std::string, std::vector<int>>>* get_barcode_data(std::string &fname,
        std::string &adapter, std::vector<std::string>&kmer_list, std::vector<std::string>*bc_list) {
    std::vector<std::string> barcodes=read_barcode_file(fname);
    std::unordered_map<std::string, size_t> kmer_pos {get_pos_hash_vector(kmer_list)};
    std::unordered_map<std::string, size_t> barcode_pos {get_pos_hash_vector(barcodes)};
    std::map<std::string, std::vector<int>> fwd_bc_kmers {}, bkw_bc_kmers{};
    std::unordered_map<size_t, std::map<std::string, std::vector<int>>> * out_hash_ptr;
    out_hash_ptr = new std::unordered_map<size_t, std::map<std::string, std::vector<int>>>;
    for (auto &bc: barcodes) {
        (*bc_list).push_back(bc);
        std::string bc_adapter = adapter+bc;
        fwd_bc_kmers.insert(std::make_pair(bc, get_sequence_kmers(bc_adapter, 6, kmer_pos)));
        std::vector<int >seq_kmer {get_sequence_kmers(reverse_complement(bc_adapter), 6, kmer_pos)};
        std::reverse(seq_kmer.begin(), seq_kmer.end());
        bkw_bc_kmers.insert(std::make_pair(bc, seq_kmer));
    }
    out_hash_ptr->insert(std::make_pair(0, fwd_bc_kmers));
    out_hash_ptr->insert(std::make_pair(1, bkw_bc_kmers));
    return out_hash_ptr;
}

std::vector<int> get_sequence_kmers(const std::string& seq, size_t kmer_len, std::unordered_map<std::string, size_t> &kmer_pos) {
    size_t len{seq.length()};
    std::vector<int> pos_seq {};
    for (size_t i {0}; i<=len-kmer_len; i++) {
        pos_seq.push_back(kmer_pos[seq.substr(i, kmer_len)]);
    }
    return pos_seq;
}

std::vector<std::string> read_barcode_file(std::string fname){
    std::ifstream infile(fname);
    std::vector<std::string> barcodes;
    std::string barcode;
    while(infile>>barcode) {
        barcodes.push_back(barcode);
    }
    return barcodes;
}

std::vector<std::vector<int>> get_barcode_matrix(const std::map<std::string, std::vector<int>>&bc_map) {
    std::vector<std::vector<int>> barcode_matrix;
    auto it {bc_map.begin()};
    while (it!=bc_map.end()) {
        barcode_matrix.push_back(it->second);
        it++;
    }
    return barcode_matrix;
}

std::string reverse_complement(std::string seq) {
    std::unordered_map<char, char> conversion_table {{'A', 'T'},
                                                     {'T', 'A'},
                                                     {'C', 'G'},
                                                     {'G', 'C'}};

    std::string reverse_complemented {};
    std::reverse(seq.begin(), seq.end());
    for (auto c: seq) {
        reverse_complemented+=conversion_table[c];
    }
    return reverse_complemented;
}

bool compute_rejection_condition(size_t &n_kmers, size_t &strand, size_t &len_read,
                                 size_t &start, size_t &end, size_t &upper_bound, size_t &lower_bound) {
    bool rej_cond1 {start>len_read/3 and strand == 0};
    bool rej_cond2 {end<2*len_read/3 and strand == 1};
    return  n_kmers>= upper_bound or n_kmers <=lower_bound or rej_cond1 or rej_cond2;
}

void double_strand_smith_waterman(std::string&query, std::string& ref,
        size_t* strand, uint16_t* sw_score) {
    int32_t maskLen = strlen(query.c_str())/2;
    maskLen = maskLen < 15 ? 15 : maskLen;
    // Declares a default Aligner
    StripedSmithWaterman::Aligner aligner {StripedSmithWaterman::Aligner(1, 1, 1, 1)};
    // Declares a default filter
    StripedSmithWaterman::Filter filter;
    // Declares an alignment that stores the result
    StripedSmithWaterman::Alignment alignment;
    StripedSmithWaterman::Alignment alignment_rev_comp;
    // Aligns the query to the ref
    aligner.Align(query.c_str(), ref.c_str(), ref.size(), filter, &alignment, maskLen);
    aligner.Align(query.c_str(), reverse_complement(ref).c_str(), ref.size(), filter, &alignment_rev_comp, maskLen);
    if (alignment.sw_score >= alignment_rev_comp.sw_score) {
        *strand=0;
        *sw_score=alignment.sw_score;
    }
    else {
        *strand=1;
        *sw_score=alignment_rev_comp.sw_score;
    }


}
