//
// Created by Tirelli Andrea on 05/04/2020.
//
#include <iostream>
#include <vector>
#include <numeric>
#include <unordered_map>
#include <tuple>
#include <regex>
#include <fstream>
#include "signal_processing.hpp"
#include "dtw.hpp"
#include "utils.hpp"
#include "AllToOneDtw.hpp"
#include "nanopolish_fast5_io.h"
#include "boost/algorithm/string.hpp"


#ifndef NANOPORE_DEMULTIPLEX_HPP
#define NANOPORE_DEMULTIPLEX_HPP


void demultiplex_file(const std::string &, const std::string &, std::string &, std::string &, std::string &);

void process_read(fast5_file &, std::string& ,
                  std::unordered_map<size_t, std::map<std::string, std::vector<int>>>&,
std::vector<std::string> &,
std::vector<float> &,std::vector<float>& ,
std::string* , float* );

void write_output(std::string &, std::vector<std::string> &,std::vector<std::string> &, std::vector<float>&);

void assign_barcode(std::vector<float>&, size_t& ,std::map<std::string, std::vector<int>>&,
                    std::vector<std::string> &, std::vector<float> &,std::vector<float>& ,
                    size_t &len_read,std::string*, float*);
#endif //NANOPORE_DEMULTIPLEX_HPP
