#ifndef MINIGRAPH_ALIGNMENT_H
#define MINIGRAPH_ALIGNMENT_H
#include "AlignmentCommon.h"
#include "SeqLib/UnalignedSequence.h"
#include "minigraph/gfa.h"
#include "minigraph/mgpriv.h"
#include "minigraph/minigraph.h"
#include <string>
#include <unordered_map>

class MinigraphAlignment {
public:
    MinigraphAlignment(std::string gfa_file);
    void alignSequence(const SeqLib::UnalignedSequenceVector seqs);
    ~MinigraphAlignment();

private:
    mg_mapopt_t m_opt;
    mg_idxopt_t m_ipt;
    mg_ggopt_t m_gpt;

    gfa_t* m_gfa;
    mg_idx_t *m_gfa_index;

    std::unordered_map<SeqLib::UnalignedSequence, mg_gchains_t*,
        UnalignedSequenceHash, UnalignedSequenceEqualsTo> m_alignments;
};

#endif
