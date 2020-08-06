#ifndef MINIGRAPH_ALIGNMENT_H
#define MINIGRAPH_ALIGNMENT_H
#include "minigraph/minigraph.h"
#include <string>

class MinigraphAlignment {
    public:

    MinigraphAlignment(std::string gfa);
    ~MinigraphAlignment();

  private:
    mg_mapopt_t m_opt;
    mg_idxopt_t m_ipt;
    mg_ggopt_t m_gpt;

    gfa_t* m_gfa;
    mg_idx_t *m_gfa_index;
};

#endif
