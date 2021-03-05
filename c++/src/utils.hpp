//
// Created by Tirelli Andrea on 05/03/2020.
//

#ifndef NANOPORE_UTILS_HPP
#define NANOPORE_UTILS_HPP
#include <string>
#include <vector>
#include <unordered_map>

std::vector<std::vector<std::string>> read_tsv(std::string, char);

std::vector<std::string> read_barcode_file(std::string);

std::unordered_map<std::string, std::unordered_map<std::string, float>> get_hash_kmers(std::vector<std::vector<std::string>> &);

std::unordered_map<std::string, size_t> get_pos_hash_vector(std::vector<std::string>&);

std::unordered_map<size_t, std::map<std::string, std::vector<int>>>* get_barcode_data(std::string&, std::string&, std::vector<std::string>&, std::vector<std::string>*);

std::vector<int> get_sequence_kmers(const std::string&, size_t, std::unordered_map<std::string, size_t> &);

std::vector<std::vector<int>> get_barcode_matrix(const std::map<std::string, std::vector<int>>&);

std::string reverse_complement(std::string);

void double_strand_smith_waterman(std::string&, std::string&,size_t*, uint16_t*);

bool compute_rejection_condition(size_t &, size_t &, size_t &,size_t &, size_t &,size_t &, size_t &);


#endif //NANOPORE_UTILS_HPP
