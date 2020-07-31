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
  std::cerr << "Sum " << total << std::endl;

  std::cerr << "Barcode phase sets" << std::endl;
  for (const auto &b : m_barcode_phase) {
    std::cerr << b.first << " " << b.second << std::endl;
    total += b.second;
  }

  std::cerr << "Barcode haplotype" << std::endl;
  for (const auto &b : m_barcode_hap) {
    std::cerr << b.first << " " << b.second << std::endl;
    total += b.second;
  }

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
  SeqLib::FermiAssembler fermi;

  fermi.SetMinOverlap(m_params.min_overlap);
  fermi.AddReads(m_reads);
  if(m_params.aggressive_bubble_pop)
      fermi.SetAggressiveTrim();
  // fermi.CorrectAndFilterReads();
  fermi.PerformAssembly();

  std::ofstream gfa_out(m_prefix + ".gfa");
  fermi.WriteGFA(gfa_out);
  gfa_out.close();

  size_t count = 0;
  for(auto contig : fermi.GetContigs()) {
      std::stringstream ss;
      ss << m_prefix << "_" << count;
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

              // barcode phase set init
              int ps;
              if(bam_record.GetIntTag("PS", ps)) {
                  m_barcode_phase[bx_tag] = ps;

                  // barcode haplotype init
                  int hp;
                  if(bam_record.GetIntTag("HP", hp))
                      m_barcode_hap[bx_tag] = hp;
              }
          }
        else
          m_barcode_count[bx_tag] = m_barcode_count[bx_tag] + 1;
      }
    } else
      break;
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
