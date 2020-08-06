#include "MinigraphAlignment.h"
#include "minigraph/gfa.h"
#include "minigraph/minigraph.h"

MinigraphAlignment::MinigraphAlignment(std::string gfa) {
  mg_opt_set(0, &m_ipt, &m_opt, &m_gpt);

  m_gfa = gfa_read(gfa.c_str());
  m_gfa_index = mg_index(m_gfa, &m_ipt, 1, &m_opt); // combine mg_index_core() and mg_opt_update()
}

MinigraphAlignment::~MinigraphAlignment() {
    gfa_destroy(m_gfa);
    mg_idx_destroy(m_gfa_index);
}
