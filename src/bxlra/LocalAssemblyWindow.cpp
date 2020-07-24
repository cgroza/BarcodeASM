#include "LocalAssemblyWindow.h"
#include "SeqLib/FermiAssembler.h"
#include <algorithm>
#include <iterator>
#include <vector>

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
  BxBarcodeCounts barcodes = collectLocalBarcodes();
  std::cerr << "Pre barcode collection: " << m_reads.size() << std::endl;

  // Barcode frequency in assembly window
  // std::cerr << std::endl;
  size_t total = 0;
  for (const auto &b : barcodes) {
    std::cerr << b.first << " " << b.second << std::endl;
    total += b.second;
  }
  // std::cerr << "Sum " << total << std::endl;

  // add the genome wide reads to assembly
  BamReadVector genomewide_reads = m_bx_bam.fetchReadsByBxBarcode(barcodes);
  std::copy(genomewide_reads.begin(), genomewide_reads.end(), std::back_inserter(m_reads));
  // std::cerr << "Post barcode collection: " << m_reads.size() << std::endl;
  return m_reads.size();
}

size_t LocalAssemblyWindow::assembleReads() {
  retrieveGenomewideReads();
  SeqLib::FermiAssembler fermi;

  fermi.SetMinOverlap(m_params.min_overlap);
  fermi.AddReads(m_reads);
  fermi.CorrectAndFilterReads();
  fermi.PerformAssembly();

  size_t count = 0;
  for(auto contig : fermi.GetContigs()) {
      std::stringstream ss;
      ss << m_prefix << "_" << count;
      m_contigs.push_back(SeqLib::UnalignedSequence(ss.str(), contig));
      ++count;
  }
  return count;
}

BxBarcodeCounts LocalAssemblyWindow::collectLocalBarcodes() {
  std::cerr << m_region.ToString(m_bam.Header()) << std::endl;
  m_bam.SetRegion(m_region);

  BamReadVector read_vector;

  while (true) {
    // Retrieve all reads within this region.
    SeqLib::BamRecord bam_record;

    if (m_bam.GetNextRecord(bam_record)) {
      read_vector.push_back(bam_record);
    } else
      break;
  }
  // Add these records to the assembly
  m_reads = read_vector;

  std::cerr << "Local reads: " << read_vector.size() << std::endl;
  return BxBamWalker::collectBxBarcodes(read_vector);
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
