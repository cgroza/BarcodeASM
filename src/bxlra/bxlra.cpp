#include "BxBamWalker.h"
#include "LocalAssemblyWindow.h"
#include "RegionFileReader.h"
#include "SeqLib/BamRecord.h"
#include <algorithm>
#include <iostream>
#include <iterator>
#include <unistd.h>

namespace opt {
std::string bam_path;
std::string bx_bam_path;
std::string regions_path;
bool weird_reads_only = true;
} // namespace opt

int main(int argc, char **argv) {
  opterr = 0;
  int c;
  while ((c = getopt(argc, argv, "ab:B:r:")) != -1)
    switch (c) {
    case 'b':
      opt::bam_path = optarg;
      break;
    case 'B':
      opt::bx_bam_path = optarg;
      break;
    case 'r':
      opt::regions_path = optarg;
      break;
    case 'a':
      opt::weird_reads_only = false;
      break;
    default:
      abort();
    }

  SeqLib::BamReader bam_reader = SeqLib::BamReader();
  bam_reader.Open(opt::bam_path);

  BxBamWalker bx_bam_walker = BxBamWalker(opt::bx_bam_path, "0000", opt::weird_reads_only);

  RegionFileReader region_reader(opt::regions_path, bam_reader.Header());

  for (auto &region : region_reader.getRegions()) {
    LocalAssemblyWindow local_win(region, bam_reader, bx_bam_walker);
    // std::cerr << local_win.retrieveGenomewideReads() << std::endl;
    local_win.assembleReads();
    std::cerr << "Reads: " << local_win.getReads().size() << std::endl;
    std::cerr << "Contigs: " << local_win.getContigs().size() << std::endl;
    for(auto &contig : local_win.getContigs())
        {
            std::cerr << contig.Seq.length() << std::endl;
        }
  }
}
