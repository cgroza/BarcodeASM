#ifndef LOCAL_ALGIGNMENT_H
#define LOCAL_ALGIGNMENT_H

#include "SeqLib/BWAWrapper.h"
#include "SeqLib/BamRecord.h"
#include "SeqLib/GenomicRegion.h"
#include "SeqLib/RefGenome.h"
#include "SeqLib/UnalignedSequence.h"
#include "minimap.h"
#include <stdlib.h>
#include <cstring>

class LocalAlignment {
public:
    LocalAlignment(std::string chr, size_t start, size_t end,
                   const SeqLib::RefGenome &genome);

    ~LocalAlignment();
    void align(const SeqLib::UnalignedSequenceVector &seqs);
    size_t writeAlignments();

    // default minimap2 parameters
    const int MINIMIZER_K = 15;
    const int MINIMIZER_W = 10;
    const int BUCKET_BITS = 64;
    const int IS_HPC = 0;

private:
    mm_idx_t* m_minimap_index;
    mm_idxopt_t m_index_opt;
    mm_mapopt_t m_map_opt;

    char *m_local_sequence;

    SeqLib::BamRecordVector m_local_alignments;
    std::string m_chr;
    size_t m_start, m_end;
};

#endif
