#ifndef LOCAL_ALGIGNMENT_H
#define LOCAL_ALGIGNMENT_H

#include "SeqLib/BWAWrapper.h"
#include "SeqLib/BamRecord.h"
#include "SeqLib/GenomicRegion.h"
#include "SeqLib/RefGenome.h"
#include "SeqLib/UnalignedSequence.h"
#include "minimap2/minimap.h"
#include <cstring>
#include <ostream>
#include <stdlib.h>
#include <unordered_map>
#include <sstream>
#include "AlignmentCommon.h"

struct LocalAlignmentParams {
  int max_join_long = 20000;
  int max_join_short = 10000;
  int min_join_flank_sc = 10;
  float min_join_flank_ratio = 0.1f;
  int max_gap = 10000;
  int bw = 2000;
  float pri_ratio = 0.8f;
  int max_chain_skip = 25;
  int end_bonus = 10;
  float chain_gap_scale = 1.0f;
  int min_chain_score = 40;
  int e = 4;  int e2 = 0;
  int q = 10;  int q2 = 300;
  int max_chain_iter = 5000;
  float max_clip_ratio = 1.0f;
  int zdrop = 10000;
  int zdrop_inv = 1000;
  int max_gap_ref = -1;
  int a = 2;  int b = 4;
  // minimizer parameters
  int minimizer_w = 10;
  int minimizer_k = 15;
  int bucket_bits = 64;
  int is_hpc = 1;
};

class LocalAlignment {
public:
  LocalAlignment(std::string chr, size_t start, size_t end,
                 const SeqLib::RefGenome &genome);
  LocalAlignment(std::string target_sequence, std::string target_name);

  ~LocalAlignment();
  void align(const SeqLib::UnalignedSequenceVector &seqs);
  size_t writeAlignments(std::ostream &out);

  // default minimap2 parameters

  static std::string getAlignmentHeader() {
    return "TName TLength TStart TEnd QName QLength QStart QEnd Hit Strand "
           "CIGAR";
  }

private:
  void setupIndex(std::string target_sequence);

  mm_idx_t *m_minimap_index;
  mm_idxopt_t m_index_opt;
  mm_mapopt_t m_map_opt;

  char *m_local_sequence;

  std::string m_target_name;

  std::unordered_map<SeqLib::UnalignedSequence, MinimapAlignment,
                     UnalignedSequenceHash, UnalignedSequenceEqualsTo>
      m_alignments;

  LocalAlignmentParams m_params;
};

#endif
