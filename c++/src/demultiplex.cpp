//
// Created by Tirelli Andrea on 05/04/2020.
//

#include "demultiplex.hpp"
#include <regex>
#include <boost/filesystem.hpp>

void demultiplex_file(const std::string &fname,const std::string &kmer_file,  std::string &barcode_file,
                      std::string & umi_dir, std::string & adapter, std::string & out_dir) {
    fast5_file f_p = fast5_open(fname);
    std::vector<std::string> read_ids = fast5_get_multi_read_groups(f_p);
    size_t n_reads {read_ids.size()};

    auto kmers = read_tsv(kmer_file, '\t');
    std::vector<std::string> kmer_list {};
    std::vector<float> kmer_means, kmer_stds;
    kmer_list.reserve(kmers.size());
    for (const auto &kmer:kmers) {
        kmer_list.push_back(kmer.at(0));
        kmer_means.push_back(stof(kmer.at(1)));
        kmer_stds.push_back(stof(kmer.at(3)));
    }

    std::vector<std::string> bc_list;
    auto barcodes = get_barcode_data(barcode_file, adapter, kmer_list, &bc_list);

    std::vector<std::string> assigned_barcodes (read_ids.size(), "");
    std::vector<std::string> assigned_umis (read_ids.size(), "");
    std::vector<float> computed_distances (read_ids.size(), 0);

    for (size_t i{0}; i<n_reads; i++) {
        std::string read_id {read_ids.at(i)};
        read_id.erase(0, 5);
        std::cout<< "Processing read number " << i << std::endl;
        process_read(f_p, read_id, *barcodes, bc_list, umi_dir, adapter, kmer_list, kmer_means, kmer_stds,
                     &(assigned_barcodes.at(i)),&(assigned_umis.at(i)), &(computed_distances.at(i)));
    }
    boost::filesystem::path p(fname);
    std::string basename = p.filename().string();;
    std::string fname_output=out_dir+std::regex_replace(basename, std::regex(".fast5"), ".txt");
    write_output(fname_output,read_ids, assigned_barcodes, assigned_umis,computed_distances);

}

void process_read(fast5_file &f_p, std::string& read_id,
                  std::unordered_map<size_t, std::map<std::string, std::vector<int>>>&barcodes,
                  std::vector<std::string> &bc_list,
                  std::string &umi_dir,
                  std::string &adapter,
                  std::vector<std::string>& kmer_list,
                  std::vector<float> &kmer_means,std::vector<float>& kmer_stds,
                  std::string* assigned_barcode,std::string* assigned_umi, float* computed_dist) {

    // load data to retrieve signal and read
    fast5_raw_scaling scaling{fast5_get_channel_params(f_p, read_id)};
    raw_table sig{fast5_get_raw_samples(f_p, read_id, scaling)};
    std::string seq{fast5_get_sequence(f_p, read_id)};

    // load signal and read
    std::vector<float> signal(sig.raw, sig.raw + sig.n);
    std::string read{fast5_get_sequence(f_p, read_id)};
    size_t len_read (read.size());

    // find poly-A tail and strand with Smith-Waterman algorithm
    uint16_t sw_score {0};
    size_t sw_strand {0};
    std::string query (15, 'T');
    double_strand_smith_waterman(query, read, &sw_strand, &sw_score);

    if (sw_score>=12) {
        //compute barcode
        std::map<std::string, std::vector<int>> bc_dir{barcodes[sw_strand]};
        assign_barcode(signal, sw_strand, bc_dir, bc_list, umi_dir,adapter, kmer_list, kmer_means, kmer_stds,
                       len_read, assigned_barcode,assigned_umi, computed_dist);
    }
    else {
        // no poly-A tail was found so we discard the read
        *assigned_barcode = "no_polyA_tail";
        *computed_dist = 0;
    }
}

void assign_barcode(std::vector<float>&signal, size_t& sw_strand,std::map<std::string, std::vector<int>>&barcodes,
                    std::vector<std::string> &bc_list, std::string &umi_dir, std::string &adapter,
                    std::vector<std::string>& kmer_list, std::vector<float> &kmer_means,std::vector<float>& kmer_stds,
                    size_t &len_read,std::string*assigned_barcode, std::string*assigned_umi, float*computed_distance) {
    // signal processing
    cut_first_part(signal);

    // segment signal and find breakpoints
    std::vector<size_t> bkps (segment_signal(signal, len_read, 2));

    //quantize signal
    std::vector<float> quantized_signal (quantize_signal(signal,bkps));

    float differences[quantized_signal.size()];
    std::adjacent_difference(quantized_signal.begin(), quantized_signal.end(), differences);
    for (auto &val:differences) {val = val < 0 ? -val : val;}
    std::vector<float> diff;
    for (auto &val :differences) {
        if (val!=0) { diff.push_back(val);}
    }
    float std {std_(diff, 1)};

    //compute direction
    size_t start {0}, end{0}, strand{0};
    compute_direction(quantized_signal, &start,&end,&strand, std/4);

    if (sw_strand!=strand or (strand==0 and start<=bkps.front()) or (strand==1 and end >=bkps.back())) {
        *assigned_barcode = "";
        *computed_distance = 0;
    }
    else {
        //cut subsignal
        cut_subsignal(signal, strand, start, end);
        auto relevant_bkps = select_bkps(bkps, strand,start,end);

        //get segments means and stds
        std::vector<float> means(relevant_bkps.size()-1), stds(relevant_bkps.size()-1);
        get_kmers_info(signal,relevant_bkps,&means,&stds);
        size_t n_kmers {means.size()}, upper_bound{200}, lower_bound{50};
        if (strand == 1) {reverse(means.begin(), means.end());reverse(stds.begin(), stds.end());}

        bool rej_cond = compute_rejection_condition(n_kmers,strand, len_read,
                                                    start, end, upper_bound, lower_bound);
        if (!rej_cond) {
            *assigned_barcode = "";
            *computed_distance = 0;
        }
        else {
            auto KL_matrix = compute_KL_matrix(means,stds,kmer_means,kmer_stds);
            std::vector<std::vector<int>> bc_matrix = get_barcode_matrix(barcodes);
            std::vector<size_t> seq_;
            for (size_t i{0}; i<means.size(); i++) {seq_.push_back(i);}
            auto start_indices = compute_start_indices(bc_matrix);
            size_t barcode_index {0};
            choose_barcode(seq_, bc_matrix, KL_matrix, start_indices, &barcode_index, computed_distance);
            *assigned_barcode = bc_list[barcode_index];

            //open file with corresponding UMI list
            std::vector<std::string> umi_list;
            std::string umi_file{umi_dir+*assigned_barcode+".txt"};
            std::string ad_barcode = adapter+*assigned_barcode;
            auto umis = get_barcode_data(umi_file, ad_barcode, kmer_list, &umi_list);
            auto umi_dir{(*umis)[sw_strand]};
	    
            //repeat above with adapter+barcode+umi
            std::vector<std::vector<int>> umi_matrix = get_barcode_matrix(umi_dir);
            auto start_indices_umi = compute_start_indices(umi_matrix);
            size_t umi_index {0};
            choose_barcode(seq_, umi_matrix, KL_matrix, start_indices_umi, &umi_index, computed_distance);
            *assigned_umi = umi_list[umi_index];
            delete [] umis;


        }
    }
}

void write_output(std::string &fname, std::vector<std::string> &read_ids,
                  std::vector<std::string> &barcodes, std::vector<std::string> &umis, std::vector<float>&distances) {
    std::ofstream outfile (fname);
    size_t n_reads {read_ids.size()};
    for (size_t i{0}; i<n_reads; i++) {
        outfile << read_ids[i] << "\t" << barcodes[i] << "\t"<<  umis[i] << "\t" <<distances[i] << std::endl;
    }



}
