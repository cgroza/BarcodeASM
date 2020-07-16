#include "LocalAlignment.h"

LocalAlignment::LocalAlignment(std::string chr, size_t start, size_t end,
                               const SeqLib::RefGenome &genome)
{
    std::string region = genome.QueryRegion(chr, start, end);
    setupIndex(region);
}

LocalAlignment::LocalAlignment(std::string target_sequence) {
  setupIndex(target_sequence);
}

void LocalAlignment::setupIndex(std::string target_sequence) {
  m_local_sequence = new char[target_sequence.size() + 1];
  memcpy(m_local_sequence, target_sequence.c_str(), target_sequence.size() + 1);

  mm_set_opt(0, &m_index_opt, &m_map_opt);
  m_map_opt.flag |= MM_F_CIGAR; // perform alignment

  m_minimap_index = mm_idx_str(MINIMIZER_W, MINIMIZER_W, IS_HPC, BUCKET_BITS, 1,
                               (const char **)&m_local_sequence, NULL);
  // update the mapping options
  mm_mapopt_update(&m_map_opt, m_minimap_index);
  mm_idx_stat(m_minimap_index);
  }

LocalAlignment::~LocalAlignment() {
    // free allocated memory
    mm_idx_destroy(m_minimap_index);
    free(m_local_sequence);
}

void LocalAlignment::align(const SeqLib::UnalignedSequenceVector &seqs) {
  for (auto &seq : seqs) {
    std::cerr << ">Query " << seq.Name.c_str() << " " << seq.Seq.c_str() << " "
              << seq.Seq.length() << std::endl;
    std::cerr << ">Target " << m_local_sequence << " " << m_minimap_index->seq->len << " "
              << m_minimap_index->b << " " << m_minimap_index->w << " "
              << m_minimap_index->k << " " << m_minimap_index->flag
              << std::endl;

    mm_reg1_t *reg;
    int num_hits;
    mm_tbuf_t *thread_buf = mm_tbuf_init();

    reg = mm_map(m_minimap_index, seq.Seq.length(), seq.Seq.c_str(), &num_hits,
                 thread_buf, &m_map_opt, seq.Name.c_str());

    std::cerr << "Number of hits: " << num_hits << std::endl;
    for (int j = 0; j < num_hits; ++j) { // traverse hits and print them out
        std::cerr << "HIT: " << j << std::endl;
      mm_reg1_t *r = &reg[j];
      assert(r->p); // with MM_F_CIGAR, this should not be NULL
      for (int i = 0; i < r->p->n_cigar; ++i)
        std::cerr << (r->p->cigar[i] >> 4) << ("MIDNSH"[r->p->cigar[i] & 0xf]);
      std::cerr << std::endl;
      std::cerr << "Q: " << (r->qs) << " " << (r->qe) << std::endl;
      std::cerr << "R: " << (r->rs) << " " << (r->re) << std::endl;
      free(r->p);
    }
    std::cerr << std::endl;
    free(reg);
    mm_tbuf_destroy(thread_buf);
  }
}

size_t LocalAlignment::writeAlignments() {
    return m_local_alignments.size();
}

// TODO: extend the local alignment window since contigs sometimes go beyond it.
// TODO: test insertion calling by aligning the reference sequence with an artificial insertion
