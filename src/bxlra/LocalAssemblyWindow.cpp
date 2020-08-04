#include "LocalAssemblyWindow.h"
#include "BxBamWalker.h"

LocalAssemblyWindow::LocalAssemblyWindow(SeqLib::GenomicRegion region,
                                         SeqLib::BamReader bam,
                                         BxBamWalker bx_bam,
                                         AssemblyParams params)
    : m_params(params), m_region(region), m_bam(bam), m_bx_bam(bx_bam) {

  std::stringstream prefix_ss;
  prefix_ss << region.ChrName(bam.Header()) << "_" << region.pos1 << "_" << region.pos2;
  m_prefix = prefix_ss.str();

}

size_t LocalAssemblyWindow::retrieveGenomewideReads() {
  collectLocalBarcodes();
  std::cerr << "Pre barcode collection: " << m_reads.size() << std::endl;

  // Barcode frequency in assembly window
  std::cerr << "Barcode frequency" <<  std::endl;
  size_t total = 0;
  for (const auto &b : m_barcode_count) {
    std::cerr << b.first << " " << b.second << std::endl;
    total += b.second;
  }
  std::cerr << "Total reads " << total << std::endl;
  std::cerr << "Total barcodes " << m_barcode_count.size() << std::endl;

  std::cerr << "Barcode phase sets" << std::endl;
  for (const auto &b : m_barcode_phase) {
    std::cerr << b.first << " " << b.second << std::endl;
    total += b.second;
  }
  std::cerr << "Total phase sets " << m_barcode_phase.size() << std::endl;

  std::cerr << "Barcode haplotype" << std::endl;
  for (const auto &b : m_barcode_hap) {
    std::cerr << b.first << " " << b.second << std::endl;
    total += b.second;
  }
  std::cerr << "Phased barcodes " << m_barcode_hap.size() << std::endl;

  // make sure to only import unique reads
  // tally already imported reads
  std::unordered_set<std::string> seqs;
  for(auto& r : m_reads)
      seqs.insert(r.Sequence());
  // add the genome wide reads to assembly
  std::vector<BxBarcode> barcodes;
  for (auto &b : m_barcode_count)
      barcodes.push_back(b.first);

  BamReadVector genomewide_reads = m_bx_bam.fetchReadsByBxBarcode(barcodes);
  for(auto &r : genomewide_reads) {
      // add this record if it is new
      if(seqs.count(r.Sequence()) == 0) {
          seqs.insert(r.Sequence());
          m_reads.push_back(r);
      }
  }
  std::cerr << "Post barcode collection: " << m_reads.size() << std::endl;
  return m_reads.size();
}

size_t LocalAssemblyWindow::assembleReads() {
  retrieveGenomewideReads();

  // Use the phased reads to do haploid assembly of the region
  if(m_params.split_reads_by_phase) {
      PhaseSplit split = separateReadsByPhase();

      std::cerr << "Phase 1 reads " << split.first.size() << std::endl;
      std::cerr << "Phase 2 reads " << split.second.size() << std::endl;

      size_t h1 = assemblePhase(split.first, "1");
      size_t h2 = assemblePhase(split.second, "2");
      return h1 + h2;
  }
  // Assemble the diploid assembly of the region
  else {
      return assemblePhase(m_reads, "0");
  }
}

size_t LocalAssemblyWindow::assemblePhase(BamReadVector &phased_reads, std::string phase) {
  SeqLib::FermiAssembler fermi;

  fermi.SetMinOverlap(m_params.min_overlap);
  fermi.AddReads(phased_reads);
  if (m_params.aggressive_bubble_pop)
    fermi.SetAggressiveTrim();
  // fermi.CorrectAndFilterReads();
  fermi.PerformAssembly();

  std::ofstream gfa_out(m_prefix + "_p" + phase + ".gfa");
  fermi.WriteGFA(gfa_out);
  gfa_out.close();

  size_t count = 0;
  for (auto contig : fermi.GetContigs()) {
    std::stringstream ss;
    ss << m_prefix << "_p" << phase << "_" << count;
    m_contigs.push_back(SeqLib::UnalignedSequence(ss.str(), contig));
    ++count;
  }
  return count;
}


void LocalAssemblyWindow::collectLocalBarcodes() {
  std::cerr << m_region.ToString(m_bam.Header()) << std::endl;
  m_bam.SetRegion(m_region);

  while (true) {
    // Retrieve all reads within this region and their barcode frequencies and
    // phase sets
    SeqLib::BamRecord bam_record;

    if (m_bam.GetNextRecord(bam_record)) {
      m_reads.push_back(bam_record);

      // collect barcode and its phase set
      std::string bx_tag;
      // barcode tag may not always be present
      if (bam_record.GetZTag("BX", bx_tag)) {
          if (m_barcode_count.find(bx_tag) == m_barcode_count.end()) {

              // barcode init
              m_barcode_count[bx_tag] = 1;
              fillPhasingData(bam_record, bx_tag);
          }
          else {
              m_barcode_count[bx_tag] = m_barcode_count[bx_tag] + 1;
              fillPhasingData(bam_record, bx_tag);
          }
      }
    } else
      break;
  }
}

void LocalAssemblyWindow::fillPhasingData(SeqLib::BamRecord &bam_record, std::string &bx_tag) {
  // Do nothing if already filled
  if(m_barcode_hap.count(bx_tag) == 1)
    return;
  // barcode phase set init
  int ps;
  if (bam_record.GetIntTag("PS", ps)) {
    m_barcode_phase[bx_tag] = ps;

    // barcode haplotype init
    int hp;
    if (bam_record.GetIntTag("HP", hp))
      m_barcode_hap[bx_tag] = hp;
  }
}

SeqLib::UnalignedSequenceVector LocalAssemblyWindow::getContigs() const {
  return m_contigs;
}

BamReadVector LocalAssemblyWindow::getReads() const { return m_reads; }

void LocalAssemblyWindow::sortContigs() {
    // sort contigs in decreasing sequence length order
    std::sort(m_contigs.begin(), m_contigs.end(),
              [](SeqLib::UnalignedSequence &a, SeqLib::UnalignedSequence &b) {
                  return a.Seq.length() > b.Seq.length();
              });
}

PhaseSplit LocalAssemblyWindow::separateReadsByPhase() {
    PhaseSplit phase_split;
    // run through the reads and split according to barcode/phase association
    for(auto &r : m_reads) {
        std::string bx_tag;
        if(r.GetZTag("BX", bx_tag)) {
            // check if we have a phasing for this barcode
            if(m_barcode_hap.count(bx_tag) == 0) {
                std::cerr << "Unphased barcode " << bx_tag << std::endl;
                // add read to both phases if read is unphased
                phase_split.first.push_back(r);
                phase_split.second.push_back(r);
                continue;
            }
            // inspect the haplotype tag for the phase
            switch(m_barcode_hap[bx_tag]) {
            case 1:
                phase_split.first.push_back(r);
                break;
            case 2:
                phase_split.second.push_back(r);
                break;
            }
        }
    }
    return phase_split;
}


void LocalAssemblyWindow::clearReads() {
    m_reads.clear();
}

void LocalAssemblyWindow::writeContigs(std::ostream &out) {
    for(auto &contig : m_contigs) {
        out << ">" << contig.Name << std::endl;
        out << contig.Seq << std::endl;
    }
}

std::string LocalAssemblyWindow::getPrefix() const {
    return m_prefix;
}
