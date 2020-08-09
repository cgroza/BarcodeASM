#ifndef READ_ALIGNMENT_H
#define READ_ALIGNMENT_H

#include "AlignmentCommon.h"
#include "BxBamWalker.h"
#include "SeqLib/UnalignedSequence.h"
#include "minimap2/minimap.h"
#include <ostream>
#include <vector>
#include <string>
#include <unordered_map>

typedef std::unordered_map<std::string, std::unordered_set<std::string>> MatePairContigMap;
typedef std::pair<std::string, std::string> Edge;
typedef std::vector<UnitigHit> UnitigHits;

struct EdgeHash {
  // Property: return the same value for <s1, s2> and <s2, s1>.
  // Desirable since edges are undirected in ContigMatePairGraph
  std::size_t operator()(const Edge &k) const {
    size_t h1 = std::hash<std::string>()(k.first);
    size_t h2 = std::hash<std::string>()(k.second);
    return std::hash<size_t>()(h1 + h2);
  }
};

class ContigMatePairGraph {
public:
    ContigMatePairGraph(std::unordered_set<std::string> &contigs,
                        MatePairContigMap &mate_contig_map);

    void writeGFA(std::ostream &out);
private:
    std::unordered_set<std::string> m_segments;
    std::unordered_set<Edge, EdgeHash> m_edges;
};

class ContigAlignment {
public:
    ContigAlignment(const SeqLib::UnalignedSequenceVector &contigs, const std::string &prefix);
    ~ContigAlignment();

    ContigMatePairGraph alignReads(const BamReadVector &reads);
    UnitigHits alignSequence(SeqLib::UnalignedSequence seq);
    void detectTEs(std::ostream &out);

    // default minimap2 parameters
    const int MINIMIZER_K = 15;
    const int MINIMIZER_W = 10;
    const int BUCKET_BITS = 64;
    const int IS_HPC = 0;

    // Reference ALU sequence, used for detecting contigs with potential Alu
    // inserts
    static const std::string ALU_REF;
private:
    std::string m_prefix;
    char** m_sequences;         // contigs to be aligned
    char** m_names;             // names of contigs
    size_t m_num_seqs;

    mm_idx_t *m_minimap_index;
    mm_idxopt_t m_index_opt;
    mm_mapopt_t m_map_opt;

};

#endif
