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
  int max_join_short = 2000;
  int min_join_flank_sc = 1000;
  float min_join_flank_ratio = 0.5f;
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
    const int MINIMIZER_K = 15;
    const int MINIMIZER_W = 10;
    const int BUCKET_BITS = 64;
    const int IS_HPC = 0;

    static std::string getAlignmentHeader(){
        return "TName TLength TStart TEnd QName QLength QStart QEnd Hit CIGAR";
    }

  private:
    void setupIndex(std::string target_sequence);

    mm_idx_t* m_minimap_index;
    mm_idxopt_t m_index_opt;
    mm_mapopt_t m_map_opt;

    char *m_local_sequence;

    std::string m_target_name;

    std::unordered_map<SeqLib::UnalignedSequence,
                       MinimapAlignment,
                       UnalignedSequenceHash,
                       UnalignedSequenceEqualsTo> m_alignments;

    LocalAlignmentParams m_params;
};

#endif
