#include "BxBamWalker.h"
#include "CTPL/ctpl_stl.h"
#include "LocalAlignment.h"
#include "LocalAssemblyWindow.h"
#include "RegionFileReader.h"
#include "SeqLib/BamRecord.h"
#include "SeqLib/RefGenome.h"
#include "SeqLib/UnalignedSequence.h"
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <mutex>
#include <stdexcept>
#include <unistd.h>
#include <vector>
#include <future>

#ifdef DEBUG_MINIMAP
std::string ref("GTCAAAGTAAAGAAAACAAAACGAATCAAAATCCCTGTGCGTTTTCCAGTGAGCTGTAGCCAAGCCGAGTGAATTTTCTGTCTCTTGAATTCTAGGTGCTGAGAACTGATTGGCATCATCATCTTTTGAGAAATTAAACTTTAGTCATA"
                "AAATTTCTATTTGAACTAAAATTATATGACTCAGATATCTTTTGATTTAAACCTAATGACAGCAAATTGCTTCATAGTGGGAATACAGTAATTCTGGATAGGCATGAGGTTTGTATTTTAAAATGAAACGCAGTTTTTGAATTGCAATG"
                "AACAAGAGAATGTAATTAGTCCCAATTGCCAAATCAGAAATTCTGCCAGGAATTTCAAGAAAGGAAAAAAAGAATGGAATCATGTTCAACATCAATGCCAGGGAAAGAATGCTTCCAGAGGCTCGAGGACTCTCCATGGGAAATCAGTT"
                "CCCTCCAGTTCCTTAAATCTCCACAGCCGTACCTTTCATCAACACTGAGCATTTCCTCGATTTTAAAATTCCTGGTTTGGTATTCTGTTTTCCTTTATTAAAACTCCTGTCATTGCTCTAAAATCCCCACCTCCGTTTATAAAGTCTGA"
                "ATGTATGAAGGTCGTGAGGATTTTGTTTTCCAATGGTATCCGATTTTTCCTTTAACGTAACCAGTGAGGAGGAATAGAATAAGGGCTGGACTGGGGCGGCTGGCATCCGAAGTGCTTGCTTGTGAGCTTTTCTTTATCGCTTGTTTGTTTC"
                "CTGGATGTCTGTCTAAGGAGATGCCAGAGTCTGTGTCCTGATAACTTACTGGGTCTGTAAAGTGACTGCTGATGAGAAGCTTCCCCACCACCCCACCATGATTACACTGAAATTCAGAATTCTACCAGAAGACAATTTTCTCAGTTCATCA"
                "CCCTTCCAAGGACATTACATAGGTCTCACTTTCTCAGAGTGAGCAACAGCAAAGCAACATAGAAAGAAGGCACATCAGGCTTCCATCCCTTCCAATTTATATTACTTTGTTTCTCAACTGGCTTTTAAGCTATTGCAAGTATGTATTAATA"
                "TCGTCTCTCTTTGCAAAACTGCTTAGGTCTGAAGACTTGAAAACAAAAACATAGGAATATATATATTTCTTCCAAAGCAAAAGATACATATGTTTTCTTCTCAAGCAATAAAATATTGGGGATAGGTGTCAAAGCATCACATGATTCTAGT"
                "GTCTTATTAAAATTGGTATTTTTCTAGATCCCAGGAGTTCTATGAATAGACAGTAGCATTTCTCTTCATTCTCAATGCCAGCTATGTCTACCTGGAAGAAGAGATATACACTTAAAGCTGGTGGCAAAGGCTACTGAATTGGTTTATAGG"
                "AGAATTCCACTTTGATTATTTTGTTTAGACTTTTTTTGAAAATATTTTATAATCTGTTAAGATTCAACAACATCAGGGAAAAAAAGATCTCAAAAGCGCTATGCTTTACTCAATTTCAGATAATTATGTTAAAAATTGTTTAGAACTACA"
                "ATCAGAGGAAGGTTCTAACTAATTTTTAACAAAATAAAGTCCAACAAGTAGAATTCCATAGAACTGGCTCATTGGCATTTAGCATATGTTAAGAGACAGAATTGTTTTATTTGAAATTGTCATATTTCACCATCCTCAACATTAGGTGAA"
                "TTTTTTTTAAAGTTTCAGTGATCCCAAATGAATTTCTGATTTTCAGTAGAAAATTGCTAAACCTTCATGTTGTCTAGA");

std::string seq("GTCAAAGTAAAGAAAACAAAACGAATCAAAATCCCTGTGCGTTTTCCAGTGAGCTGTAGCCAAGCCGAGTGAATTTTCTGTCTCTTGAATTCTAGGTGCTGAGAACTGATTGGCATCATCATCTTTTGAGAAATTAAACTTTAGTCATA"
                "AAATTTCTATTTGAACTAAAATTATATGACTCAGATATCTTTTGATTTAAACCTAATGACAGCAAATTGCTTCATAGTGGGAATACAGTAATTCTGGATAGGCATGAGGTTTGTATTTTAAAATGAAACGCAGTTTTTGAATTGCAATG"
                "AACAAGAGAATGTAATTAGTCCCAATTGCCAAATCAGAAATTCTGCCAGGAATTTCAAGAAAGGAAAAAAAGAATGGAATCATGTTCAACATCAATGCCAGGGAAAGAATGCTTCCAGAGGCTCGAGGACTCTCCATGGGAAATCAGTT"
                "CCCTCCAGTTCCTTAAATCTCCACAGCCGTACCTTTCATCAACACTGAGCATTTCCTCGATTTTAAAATTCCTGGTTTGGTATTCTGTTTTCCTTTATTAAAACTCCTGTCATTGCTCTAAAATCCCCACCTCCGTTTATAAAGTCTGA"
                "ATGTATGAAGGTCGTGAGGATTTTGTTTTCCAATGGTATCCGATTTTTCCTTTAACGTAACCAGTGAGGAGGAATAGAATAAGGGCTGGACTGGGGCGGCTGGCATCCGAAGTGCTTGCTTGTGAGCTTTTCTTTATCGCTTGTTTGTTTC"
                //insertion
                "TCAAAGTAAAGAAAACAAAACGAATCAAAATCCCTGTGCGTTTTCCAGTGAGCTGTAGCCAAGCCGAGTGAATTTTCTGTCTCTTGAATTCTAGGTGCTGAGAACTGATTGGCATCATCATCTTTTGAGAAATTAAACTTTAGTCATAAAA"
                "TTTCTATTTGAACTAAAATTATATGACTCAGATATCTTTTGATTTAAACCTAATGACAGCAAATTGCTTCATAGTGGGAATACAGTAATTCTGGATAGGCATGAGGTTTGTATTTTAAAATGAAACGCAGTTTTTGAATTGCAATGAAC"
                //reference
                "CTGGATGTCTGTCTAAGGAGATGCCAGAGTCTGTGTCCTGATAACTTACTGGGTCTGTAAAGTGACTGCTGATGAGAAGCTTCCCCACCACCCCACCATGATTACACTGAAATTCAGAATTCTACCAGAAGACAATTTTCTCAGTTCATCA"
                "CCCTTCCAAGGACATTACATAGGTCTCACTTTCTCAGAGTGAGCAACAGCAAAGCAACATAGAAAGAAGGCACATCAGGCTTCCATCCCTTCCAATTTATATTACTTTGTTTCTCAACTGGCTTTTAAGCTATTGCAAGTATGTATTAATA"
                "TCGTCTCTCTTTGCAAAACTGCTTAGGTCTGAAGACTTGAAAACAAAAACATAGGAATATATATATTTCTTCCAAAGCAAAAGATACATATGTTTTCTTCTCAAGCAATAAAATATTGGGGATAGGTGTCAAAGCATCACATGATTCTAGT"
                "GTCTTATTAAAATTGGTATTTTTCTAGATCCCAGGAGTTCTATGAATAGACAGTAGCATTTCTCTTCATTCTCAATGCCAGCTATGTCTACCTGGAAGAAGAGATATACACTTAAAGCTGGTGGCAAAGGCTACTGAATTGGTTTATAGG"
                "AGAATTCCACTTTGATTATTTTGTTTAGACTTTTTTTGAAAATATTTTATAATCTGTTAAGATTCAACAACATCAGGGAAAAAAAGATCTCAAAAGCGCTATGCTTTACTCAATTTCAGATAATTATGTTAAAAATTGTTTAGAACTACA"
                "ATCAGAGGAAGGTTCTAACTAATTTTTAACAAAATAAAGTCCAACAAGTAGAATTCCATAGAACTGGCTCATTGGCATTTAGCATATGTTAAGAGACAGAATTGTTTTATTTGAAATTGTCATATTTCACCATCCTCAACATTAGGTGAA"
                "TTTTTTTTAAAGTTTCAGTGATCCCAAATGAATTTCTGATTTTCAGTAGAAAATTGCTAAACCTTCATGTTGTCTAGA");
#endif

namespace opt {
std::string bam_path;
std::string bx_bam_path;
std::string regions_path;
std::string reference_path;
size_t min_overlap = 90;
size_t num_threads = 1;
bool weird_reads_only = true;
} // namespace opt

int main(int argc, char **argv) {

  #ifdef DEBUG_MINIMAP
    std::cerr << "DEBUG_MINIMAP" << std::endl;
    LocalAlignment test_aln(ref, "test");
    SeqLib::UnalignedSequence insertion_seq("insertion", seq);
    SeqLib::UnalignedSequenceVector test_seqs;
    test_seqs.push_back(insertion_seq);
    test_aln.align(test_seqs);
    test_aln.writeAlignments(std::cerr);
    return 0;
  #endif

  opterr = 0;
  int c;
  while ((c = getopt(argc, argv, "at:b:B:r:g:o:")) != -1)
    switch (c) {
    case 't':
        try {
            opt::num_threads = std::stoi(optarg);
        }
        catch(std::invalid_argument){
            std::cerr << "Number of threads -t must be integer!" << std::endl;
            return -1;
        }
        break;
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
    case 'g':
      opt::reference_path = optarg;
      break;
    case 'o':
      opt::min_overlap = std::stoi(optarg);
      break;
    default:
      abort();
    }
  AssemblyParams params;
  params.min_overlap = opt::min_overlap;

  // Storage for thread pooled resources
  // These not be guarded by mutex, since they assigned to individual thread IDs
  std::vector<SeqLib::BamReader*> bam_readers(opt::num_threads);
  std::vector<BxBamWalker*> bx_bam_walkers(opt::num_threads);
  std::vector<SeqLib::RefGenome*> ref_genomes(opt::num_threads);

  // initialize pooled bam, bx_bam, and genome readers
  for(size_t i = 0; i < opt::num_threads; i++) {
    // one reference genome reader for each thread
    SeqLib::RefGenome *ref_genome = new SeqLib::RefGenome();
    bool load_ref = ref_genome -> LoadIndex(opt::reference_path);
    std::cerr << "Loaded " << opt::reference_path << ": " << load_ref << std::endl;
    ref_genomes[i] = ref_genome;

    // one BamReader for each thread
    SeqLib::BamReader *bam_reader = new SeqLib::BamReader();
    bam_reader -> Open(opt::bam_path);
    bam_readers[i] = bam_reader;

    // and one BX_BamReader for each thread
    BxBamWalker *bx_bam_walker = new BxBamWalker(opt::bx_bam_path, "0000", opt::weird_reads_only);
    bx_bam_walkers[i] = bx_bam_walker;
  }

  // Thread pool to run all the regions
  ctpl::thread_pool thread_pool(opt::num_threads);

  // Regions to be locally assembled
  RegionFileReader region_reader(opt::regions_path, bam_readers[0]->Header());

  // Stores output from each thread, must be protected by mutex
  std::vector<std::pair<LocalAssemblyWindow*, LocalAlignment*>> output;
  std::mutex output_mutex;

  // file to write contig sequences in
  std::ofstream fasta("contigs.fa");
  std::mutex fasta_mutex;

  // file to write alignments
  std::ofstream alns("alignments.tsv");
  std::mutex alns_mutex;
  // output alignments
  alns << LocalAlignment::getAlignmentHeader() << std::endl;

  for (auto region : region_reader.getRegions()) {
      auto future = thread_pool.push([region, &output, &output_mutex,
                                      &fasta, &fasta_mutex,
                                      &alns, &alns_mutex, &params,
                                      &ref_genomes, &bam_readers, &bx_bam_walkers](int id) {

      std::cerr << "ID " << id << std::endl;

      LocalAssemblyWindow *local_win = new LocalAssemblyWindow(region, *bam_readers[id], *bx_bam_walkers[id], params);

      local_win -> assembleReads();
      std::cerr << "Reads: " << local_win -> getReads().size() << std::endl;
      std::cerr << "Contigs: " << local_win -> getContigs().size() << std::endl;

      // MUTEX: only one thread must write to the fasta file at a time
      fasta_mutex.lock();
      local_win -> writeContigs(fasta);
      fasta_mutex.unlock();

      LocalAlignment *local_alignment = new LocalAlignment(region.ChrName(bam_readers[id] -> Header()),
                                                           region.pos1, region.pos2, *ref_genomes[id]);
      local_alignment -> align(local_win -> getContigs());

      // MUTEX: only one thread must write to alignments
      alns_mutex.lock();
      local_alignment -> writeAlignments(alns);
      alns_mutex.unlock();

      // MUTEX: only one thread must push to std::vector
      output_mutex.lock();
      output.push_back(std::pair<LocalAssemblyWindow*, LocalAlignment*>(local_win, local_alignment));
      output_mutex.unlock();
    });
  }

  thread_pool.stop(true);
  fasta.close();
  alns.close();

  // cleanup allocated memory
  for (auto &region_output : output) {
      delete region_output.first;
      delete region_output.second;
  }

  output.clear();
}
