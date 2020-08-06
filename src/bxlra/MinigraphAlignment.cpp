#include "MinigraphAlignment.h"

MinigraphAlignment::MinigraphAlignment(std::string gfa_file) {
  mg_opt_set(0, &m_ipt, &m_opt, &m_gpt);

  m_gfa = gfa_read(gfa_file.c_str());
  // 1 for only one thread
  m_gfa_index = mg_index(m_gfa, &m_ipt, 1, &m_opt); // combine mg_index_core() and mg_opt_update()
}

void MinigraphAlignment::alignSequence(const SeqLib::UnalignedSequenceVector seqs) {
    mg_tbuf_t *tbuf = mg_tbuf_init();

    for(auto &seq : seqs) {
       m_alignments[seq] = mg_map(m_gfa_index, seq.Seq.length(), seq.Seq.c_str(), tbuf, &m_opt, seq.Name.c_str());
    }
    mg_tbuf_destroy(tbuf);
}

MinigraphAlignment::~MinigraphAlignment() {
  // eventually free the chain to prevent memory leak
  for(auto& pair : m_alignments)
      mg_gchain_free(pair.second);

  // destroy gfa graph and its index
  gfa_destroy(m_gfa);
  mg_idx_destroy(m_gfa_index);
}
