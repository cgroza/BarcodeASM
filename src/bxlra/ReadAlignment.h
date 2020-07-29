#ifndef READ_ALIGNMENT_H
#define READ_ALIGNMENT_H

#include "AlignmentCommon.h"
#include "BxBamWalker.h"
#include "SeqLib/UnalignedSequence.h"
#include "minimap2/minimap.h"
#include <vector>
#include <string>

class ReadAlignment {
    public:
    ReadAlignment(const SeqLib::UnalignedSequenceVector &contigs);
    ~ReadAlignment();

    void alignReads(const BamReadVector &reads);
    // default minimap2 parameters
    const int MINIMIZER_K = 15;
    const int MINIMIZER_W = 10;
    const int BUCKET_BITS = 64;
    const int IS_HPC = 0;

  private:
    char** m_sequences;         // contigs to be aligned
    char** m_names;             // names of contigs
    size_t m_num_seqs;

    mm_idx_t *m_minimap_index;
    mm_idxopt_t m_index_opt;
    mm_mapopt_t m_map_opt;

};

#endif
