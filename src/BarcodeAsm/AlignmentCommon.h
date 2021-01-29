#ifndef ALIGNMENT_COMMON_H
#define ALIGNMENT_COMMON_H

#include "SeqLib/UnalignedSequence.h"
#include "minimap2/minimap.h"
#include <string>

struct MinimapAlignment {
  mm_reg1_t *reg;
  int num_hits;
};

struct UnalignedSequenceHash {
  std::size_t operator()(const SeqLib::UnalignedSequence &k) const {
    return std::hash<std::string>()(k.Seq);
  }
};

struct UnalignedSequenceEqualsTo {
  bool operator()(const SeqLib::UnalignedSequence &a,
                  const SeqLib::UnalignedSequence b) const {
    return std::equal_to<std::string>()(a.Seq, b.Seq) &&
           std::equal_to<std::string>()(a.Name, b.Name);
  }
};

struct UnitigHit {
    std::string unitig_name;
    size_t ql, qs, qe;
    size_t tl, ts, te;
    std::string cigar;
    char strand;
};

#endif
