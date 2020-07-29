#include "ReadAlignment.h"

ReadAlignment::ReadAlignment(const SeqLib::UnalignedSequenceVector &contigs) {
    m_num_seqs = contigs.size();
    m_sequences = new char*[m_num_seqs];

    for(size_t i = 0; i < m_num_seqs; i++) {
        size_t seq_size = contigs.at(i).Seq.length();
        m_sequences[i] = new char[seq_size + 1];
        memcpy(m_sequences[i], contigs.at(i).Seq.c_str(), seq_size + 1);
    }

    mm_set_opt(0, &m_index_opt, &m_map_opt);
    m_map_opt.flag |= MM_F_CIGAR; // perform alignment

    m_minimap_index = mm_idx_str(MINIMIZER_W, MINIMIZER_W, IS_HPC, BUCKET_BITS,
                                 1, (const char**) m_sequences , NULL);
    // update the mapping options
    mm_mapopt_update(&m_map_opt, m_minimap_index);
    mm_idx_stat(m_minimap_index);
}

ReadAlignment::~ReadAlignment() {
  // free allocated memory
  mm_idx_destroy(m_minimap_index);
  for(size_t i = 0; i < m_num_seqs; i++)
      delete m_sequences[i];
  delete m_sequences;

}

void ReadAlignment::alignReads(const BamReadVector &reads) {
  for (auto &read : reads) {
    mm_tbuf_t *thread_buf = mm_tbuf_init();
    int num_hits;
    mm_reg1_t *reg =
        mm_map(m_minimap_index, read.Sequence().length(), read.Sequence().c_str(),
               &num_hits, thread_buf, &m_map_opt, read.Qname().c_str());
    for (int j = 0; j < num_hits; ++j) { // traverse hits and print them out
      mm_reg1_t *r = &reg[j];
      assert(r->p); // with MM_F_CIGAR, this should not be NULL
      std::cerr << read.Qname() << " " << read.Sequence().length() << " " <<
          r->qs << " " << r->qe << " " << "+-"[r->rev] << " " <<
          m_minimap_index->seq[r->rid].name << " " <<
          m_minimap_index->seq[r->rid].len << " " <<
          r->rs << " " << r->re  << " " << r->mlen << " " <<
          r->blen << " " << r->mapq << std::endl;
      // for (i = 0; i < r->p->n_cigar; ++i)
      //   printf("%d%c", r->p->cigar[i] >> 4, "MIDNSH"[r->p->cigar[i] & 0xf]);
      free(r->p);
      mm_tbuf_destroy(thread_buf);
    }
  }
}
