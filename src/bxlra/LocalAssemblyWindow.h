#ifndef LOCAL_ASSEMBLY_WINDOW_H
#define LOCAL_ASSEMBLY_WINDOW_H

#include "BxBamWalker.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/BamRecord.h"
#include "SeqLib/FermiAssembler.h"
#include "SeqLib/GenomicRegion.h"
#include "SeqLib/SeqLibUtils.h"
#include "SeqLib/UnalignedSequence.h"
#include <algorithm>
#include <fstream>
#include <iterator>
#include <ostream>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

struct AssemblyParams {
    size_t min_overlap = 90;
    size_t min_contig_length = 200;
    size_t aggressive_bubble_pop = false;
    bool split_reads_by_phase = false;
};

typedef std::unordered_map<BxBarcode, int> BxBarcodeCounts;
typedef std::unordered_map<BxBarcode, int> BxBarcodePS;
typedef std::unordered_map<BxBarcode, int> BxBarcodeHap;
typedef std::pair<BamReadVector, BamReadVector> PhaseSplit;

class LocalAssemblyWindow {
public:
    LocalAssemblyWindow(SeqLib::GenomicRegion region, SeqLib::BamReader bam, BxBamWalker bx_bam, AssemblyParams params);
    size_t retrieveGenomewideReads();
    size_t assembleReads();
    void collectLocalBarcodes();
    SeqLib::UnalignedSequenceVector getContigs() const;
    BamReadVector getReads() const;
    void clearReads();
    void writeContigs(std::ostream &out);
    std::string getPrefix() const;

  private:
    void sortContigs();
    size_t assemblePhase(BamReadVector &phased_reads, std::string phase);
    PhaseSplit separateReadsByPhase();
    void fillPhasingData(SeqLib::BamRecord &bam_record, std::string &bx_tag);

    AssemblyParams m_params;
    SeqLib::GenomicRegion m_region;
    SeqLib::BamReader m_bam;
    BxBamWalker m_bx_bam;
    BamReadVector m_reads;
    std::string m_prefix;
    SeqLib::UnalignedSequenceVector m_contigs;

    // keep track of barcode frequency and their phase set
    BxBarcodeCounts m_barcode_count;
    BxBarcodePS m_barcode_phase;
    BxBarcodeHap m_barcode_hap;
};

#endif
