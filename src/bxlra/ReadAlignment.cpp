#include "ReadAlignment.h"
#include <ostream>
#include <utility>

ContigMatePairGraph::ContigMatePairGraph(std::unordered_set<std::string> &contigs,
                                         MatePairContigMap &mate_contig_map) :
    m_segments(contigs) {

    for(auto& b : mate_contig_map) {
        std::vector<std::string> keys;
        std::copy(b.second.begin(), b.second.end(), std::back_inserter(keys));
        m_edges.emplace(keys.at(0), keys.at(1));
    }
}

void ContigMatePairGraph::writeGFA(std::ostream &out){
    // Header
    out << "H\tVN:Z:1.0" << std::endl;

    // Segments
    for(auto &s : m_segments)
        out << "S\t" << s << "\t*" << std::endl;
    // Links
    for(auto &e : m_edges)
        out << "L\t" << e.first << "\t+\t" << e.second << "\t+\t*" << std::endl;
}

ReadAlignment::ReadAlignment(const SeqLib::UnalignedSequenceVector &contigs, const std::string &prefix)
    : m_prefix(prefix) {
  m_num_seqs = contigs.size();
  m_sequences = new char *[m_num_seqs];
  m_names = new char *[m_num_seqs];

  for (size_t i = 0; i < m_num_seqs; i++) {
    size_t seq_size = contigs.at(i).Seq.length();
    size_t name_size = contigs.at(i).Name.length();

    m_sequences[i] = new char[seq_size + 1];
    m_names[i] = new char[name_size + 1];

    memcpy(m_sequences[i], contigs.at(i).Seq.c_str(), seq_size + 1);
    memcpy(m_names[i], contigs.at(i).Name.c_str(), name_size + 1);
  }

  mm_set_opt(0, &m_index_opt, &m_map_opt);
  m_map_opt.flag |= MM_F_CIGAR; // perform alignment

  m_minimap_index =
      mm_idx_str(MINIMIZER_W, MINIMIZER_W, IS_HPC, BUCKET_BITS, m_num_seqs,
                 (const char **)m_sequences, (const char **)m_names);
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

ContigMatePairGraph ReadAlignment::alignReads(const BamReadVector &reads) {
  MatePairContigMap read_contig_map;
  mm_tbuf_t *thread_buf = mm_tbuf_init();
  for (auto &read : reads) {
    read_contig_map.emplace(read.Qname(), std::unordered_set<std::string>());

    int num_hits;
    mm_reg1_t *reg = mm_map(m_minimap_index, read.Sequence().length(), read.Sequence().c_str(),
               &num_hits, thread_buf, &m_map_opt, read.Qname().c_str());

    if (num_hits > 0) { // include first hit
      mm_reg1_t *r = &reg[0];
      assert(r->p); // with MM_F_CIGAR, this should not be NULL
      read_contig_map[read.Qname()].emplace(std::string(m_minimap_index->seq[r->rid].name));
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
    free(reg);
  }

  // clear empty entries or that are fully within one contig
  for(auto it = read_contig_map.begin(); it != read_contig_map.end(); ) {
      if(it -> second.size() != 2)
          it = read_contig_map.erase(it);
      else it++;
  }

  #ifdef DEBUG_READ_ALIGNMENT
  for(auto &r : read_contig_map) {
      std::cerr << r.first << " ";
      for (auto &c : r.second) std::cerr << c << " ";
      std::cerr << std::endl;
  }
  #endif

  mm_tbuf_destroy(thread_buf);

  std::unordered_set<std::string> contigs;
  for(size_t i = 0; i < m_num_seqs; i++)
      contigs.emplace(m_names[i]);

  ContigMatePairGraph contig_mate_pair_grah(contigs, read_contig_map);
  std::ofstream gfa_out(m_prefix + "_mates.gfa");
  contig_mate_pair_grah.writeGFA(gfa_out);
  gfa_out.close();
  return contig_mate_pair_grah;
}

UnitigHits ReadAlignment::alignSequence(SeqLib::UnalignedSequence seq) {
  UnitigHits unitig_hits;
  mm_tbuf_t *thread_buf = mm_tbuf_init();
  int num_hits;

  mm_reg1_t *reg = mm_map(m_minimap_index, seq.Seq.length(), seq.Seq.c_str(),
                          &num_hits, thread_buf, &m_map_opt, seq.Name.c_str());
  for(size_t i = 0; i < num_hits; i++) {
    mm_reg1_t *r = &reg[i];
    assert(r->p); // with MM_F_CIGAR, this should not be NULL

    unitig_hits.emplace_back(m_minimap_index->seq[r->rid].name);

    free(r->p);
  }
  free(reg);
  mm_tbuf_destroy(thread_buf);

  return unitig_hits;
}
