#include "BxBamWalker.h"
#include "CTPL/ctpl_stl.h"
#include "ContigAlignment.h"
#include "LocalAlignment.h"
#include "LocalAssemblyWindow.h"
#include "RegionFileReader.h"
#include "SeqLib/BamRecord.h"
#include "SeqLib/RefGenome.h"
#include "SeqLib/UnalignedSequence.h"
#include <ContigAlignment.h>
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <future>
#include <iostream>
#include <iterator>
#include <mutex>
#include <stdexcept>
#include <unistd.h>
#include <vector>
#include "SeqLib/FastqReader.h"


namespace opt {
std::string bam_path;
std::string bx_bam_path;
std::string regions_path;
std::string reference_path;
size_t min_overlap = 90;
size_t num_threads = 1;
bool weird_reads_only = true;
bool aggressive_bubble_pop = false;
bool split_reads_by_phase = false;
bool write_gfa = false;
std::string detect_seqs_fa;
} // namespace opt

int main(int argc, char **argv) {
  opterr = 0;
  int c;
  while ((c = getopt(argc, argv, "GSPat:b:B:r:g:o:F:")) != -1)
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
    case 'S':
      opt::split_reads_by_phase = true;
      break;
    case 'P':
      opt::aggressive_bubble_pop = true;
      break;
    case 'a':
      opt::weird_reads_only = false;
      break;
    case 'g':
      opt::reference_path = optarg;
      break;
    case 'G':
      opt::write_gfa = true;
      break;
    case 'o':
      opt::min_overlap = std::stoi(optarg);
      break;
    case 'F' :
      opt::detect_seqs_fa = optarg;
      break;
    default:
      abort();
    }

  AssemblyParams params;
  params.min_overlap = opt::min_overlap;
  params.aggressive_bubble_pop = opt::aggressive_bubble_pop;
  params.split_reads_by_phase = opt::split_reads_by_phase;
  params.write_gfa = opt::write_gfa;

  std::cerr << "Params o: " << params.min_overlap << std::endl
            << "Params P: " << params.aggressive_bubble_pop << std::endl
            << "Params a: " << !opt::weird_reads_only << std::endl
            << "Params S: " << opt::split_reads_by_phase << std::endl
            << "Params G: " << opt::write_gfa << std::endl;

  // Storage for thread pooled resources
  // These not be guarded by mutex, since they assigned to individual thread IDs
  std::vector<SeqLib::BamReader*> bam_readers(opt::num_threads);
  std::vector<BxBamWalker*> bx_bam_walkers(opt::num_threads);
  std::vector<SeqLib::RefGenome*> ref_genomes(opt::num_threads);

  // load sequences to be detected in contigs
  SeqLib::UnalignedSequenceVector detect_seqs;
  if(opt::detect_seqs_fa.size() > 0) {
      SeqLib::FastqReader seq_fa(opt::detect_seqs_fa);
      SeqLib::UnalignedSequence s;
      while(seq_fa.GetNextSequence(s))
          detect_seqs.push_back(s);
  }

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

  // file to write contig sequences in
  std::ofstream fasta("contigs.fa");
  std::mutex fasta_mutex;

  // file to write TE hits in
  std::ofstream hits("hits.tsv");
  std::mutex hits_mutex;

  // file to write alignments
  std::ofstream alns("alignments.tsv");
  std::mutex alns_mutex;
  // output alignments
  alns << LocalAlignment::getAlignmentHeader() << std::endl;

  for (auto region : region_reader.getRegions()) {
    std::string chrom = region.ChrName(bam_readers[0]->Header());
    std::cerr << "Running " << chrom << " " << region.pos1 << " " << region.pos2 << std::endl;
    auto future = thread_pool.push([region, &fasta, &fasta_mutex,
                                    &alns, &alns_mutex,
                                    &hits, &hits_mutex,
                                    &params, &detect_seqs,
                                    &ref_genomes, &bam_readers,
                                    &bx_bam_walkers](int id) {

      std::cerr << "ID " << id << std::endl;
      LocalAssemblyWindow local_win(region, *bam_readers[id], *bx_bam_walkers[id], params);

      local_win.assembleReads();

      ContigAlignment read_aln(local_win.getContigs(), local_win.getPrefix());
      read_aln.alignReads(local_win.getReads());

      hits_mutex.lock();
      read_aln.detectSequences(detect_seqs, hits);
      hits_mutex.unlock();

      std::cerr << "Reads: " << local_win.getReads().size() << std::endl;
      local_win.clearReads();
      std::cerr << "Contigs: " << local_win.getContigs().size() << std::endl;

      // MUTEX: only one thread must write to the fasta file at a time
      fasta_mutex.lock();
      local_win.writeContigs(fasta);
      fasta_mutex.unlock();

      LocalAlignment local_alignment(region.ChrName(bam_readers[0]->Header()),
                                      region.pos1, region.pos2, *ref_genomes[id]);
      local_alignment.align(local_win.getContigs());

      // MUTEX: only one thread must write to alignments
      alns_mutex.lock();
      local_alignment.writeAlignments(alns);
      alns_mutex.unlock();
    });
  }

  thread_pool.stop(true);
  fasta.close();
  alns.close();
}
