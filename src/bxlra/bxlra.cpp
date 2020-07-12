#include "BxBamWalker.h"
#include "RegionFileReader.h"
#include "LocalAssemblyWindow.h"
#include "SeqLib/BamRecord.h"
#include <algorithm>
#include <iostream>
#include <iterator>

int main(int argc, char **argv) {
  SeqLib::BamReader bam_reader = SeqLib::BamReader();
  bam_reader.Open(
      "/Users/cgroza/sv_assemble/NA12878_WGS_v2_phased_possorted_10p2.bam");
  BxBamWalker bx_bam_walker =
      BxBamWalker("/Users/cgroza/sv_assemble/"
                  "NA12878_WGS_v2_phased_possorted_10p2_bxindex.bam",
                  "0000", false);

  RegionFileReader region_reader(
      "/Users/cgroza/git/bx_local_read_assembler/src/bxlra/regions.bed",
      bam_reader.Header());

  for (auto &region : region_reader.getRegions()) {
      LocalAssemblyWindow local_win(region, bam_reader, bx_bam_walker);
      // std::cerr << local_win.retrieveGenomewideReads() << std::endl;
      local_win.assembleReads();
      std::cerr << "Reads: " << local_win.getReads().size() << std::endl;
      std::cerr << "Contigs: " << local_win.getContigs().size() << std::endl;
  }
}
