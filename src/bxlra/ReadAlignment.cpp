#include "ReadAlignment.h"
#include <algorithm>

ReadAlignment::ReadAlignment(const SeqLib::UnalignedSequenceVector &contigs) {
    m_num_seqs = contigs.size();
    m_sequences = new char*[m_num_seqs];
    m_names = new char*[m_num_seqs];

    for(size_t i = 0; i < m_num_seqs; i++) {
        size_t seq_size = contigs.at(i).Seq.length();
        size_t name_size = contigs.at(i).Name.length();

        m_sequences[i] = new char[seq_size + 1];
        m_names[i] = new char[name_size + 1];

        memcpy(m_sequences[i], contigs.at(i).Seq.c_str(), seq_size + 1);
        memcpy(m_names[i], contigs.at(i).Name.c_str(), name_size + 1);
    }

    mm_set_opt(0, &m_index_opt, &m_map_opt);
    m_map_opt.flag |= MM_F_CIGAR; // perform alignment

    m_minimap_index = mm_idx_str(MINIMIZER_W, MINIMIZER_W, IS_HPC, BUCKET_BITS,
                                 m_num_seqs, (const char**) m_sequences , (const char**) m_names);
    // update the mapping options
    mm_mapopt_update(&m_map_opt, m_minimap_index);
    mm_idx_stat(m_minimap_index);
}

ReadAlignment::~ReadAlignment() {
  // free allocated memory
  mm_idx_destroy(m_minimap_index);
  for(size_t i = 0; i < m_num_seqs; i++) {
      delete m_sequences[i];
      delete m_names[i];
  }
  delete m_sequences;
  delete m_names;
}

void ReadAlignment::alignReads(const BamReadVector &reads) {
  std::unordered_map<std::string, std::unordered_set<std::string>> read_contig_map;
  mm_tbuf_t *thread_buf = mm_tbuf_init();
  for (auto &read : reads) {
    read_contig_map.emplace(read.Qname(), std::unordered_set<std::string>());

    int num_hits;
    mm_reg1_t *reg = mm_map(m_minimap_index, read.Sequence().length(), read.Sequence().c_str(),
               &num_hits, thread_buf, &m_map_opt, read.Qname().c_str());

    if (num_hits > 0) { // traverse hits and print them out
      mm_reg1_t *r = &reg[0];
      assert(r->p); // with MM_F_CIGAR, this should not be NULL
      read_contig_map[read.Qname()].insert(m_minimap_index->seq[r->rid].name);
          // std::cerr << read.Qname() << " " << read.Sequence().length() << " "
          // <<
          //     r->qs << " " << r->qe << " " << "+-"[r->rev] << " " <<
          //     m_minimap_index->seq[r->rid].name << " " <<
          //     m_minimap_index->seq[r->rid].len << " " <<
          //     r->rs << " " << r->re  << " " << r->mlen << " " <<
          //     r->blen << " " << r->mapq << " " << j << std::endl;
          // for (i = 0; i < r->p->n_cigar; ++i)
          //   printf("%d%c", r->p->cigar[i] >> 4, "MIDNSH"[r->p->cigar[i] &
          //   0xf]);
      free(r->p);
    }
  }
  mm_tbuf_destroy(thread_buf);

  // clear empty entries or that are fully within one contig
  for(auto it = read_contig_map.begin(); it != read_contig_map.end(); it++) {
      if(it -> second.size() <= 1)
          read_contig_map.erase(it);
  }
  for(auto &r : read_contig_map) {
      std::cerr << r.first << " ";
      for (auto &c : r.second) std::cerr << c << " ";
      std::cerr << std::endl;
  }
}
