//
// Created by Tirelli Andrea on 02/03/2020.
//

#include "demultiplex.hpp"
#include <boost/program_options.hpp>

#include <string>
using namespace boost::program_options;

int main(int argc, const char *argv[])
{
    try {
        std::string fast5, kmer, barcodes, adapter, outdir, umi_dir;
        options_description desc{"Options"};
        desc.add_options()
                ("help,h", "Help screen")
                ("fast5", value<std::string>(&fast5)->default_value(""), "fast5 file to be demultiplexed")
                ("kmer", value<std::string>(&kmer)->default_value(""), "kmers model file")
                ("barcode", value<std::string>(&barcodes)->default_value(""), "barcode list")
                ("umi_dir", value<std::string>(&umi_dir)->default_value(""), "UMIs directory")
                ("adapter", value<std::string>(&adapter)->default_value(""), "adapter before barcode")
                ("out", value<std::string>(&outdir)->default_value(""), "output directory");

        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);
        notify(vm);
        demultiplex_file(fast5, kmer, barcodes, umi_dir, adapter, outdir);

    }
    catch (const error &ex)
    {
        std::cerr << ex.what() << '\n';
    }


}