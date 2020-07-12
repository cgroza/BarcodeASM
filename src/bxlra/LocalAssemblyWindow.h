#ifndef LOCAL_ASSEMBLY_WINDOW_H
#define LOCAL_ASSEMBLY_WINDOW_H

#include "ASQG.h"
#include "OverlapCommon.h"
#include "ReadTable.h"
#include "SGACommon.h"
#include "SeqLib/UnalignedSequence.h"
#include "overlap.h"
#include "SGUtil.h"
#include "SGVisitors.h"
#include "EncodedString.h"
#include "Util.h"
#include "SGSearch.h"
#include "CorrectionThresholds.h"

#include "SeqLib/BamReader.h"
#include "SeqLib/GenomicRegion.h"

#include "BxBamWalker.h"
#include <sstream>

struct AssemblyParams {
    double error_rate = 0;
    int seed_length = 0;
    int seed_stride = 0;
    bool irr_only = true;
    bool get_components = true;
    size_t min_overlap = 90;
    bool perform_trim = true;
    size_t trim_length_threshold = 100;
    size_t trim_rounds = 1;
    bool validate = true;
};


class LocalAssemblyWindow {

public:
    LocalAssemblyWindow(SeqLib::GenomicRegion region, SeqLib::BamReader& bam, BxBamWalker& bx_bam);
    size_t retrieveGenomewideReads();
    size_t assembleReads();
    BxBarcodeCounts collectLocalBarcodes();
    SeqLib::UnalignedSequenceVector getContigs() const;
    BamReadVector getReads() const;

  private:
    void createReadTable();
    void assembleFromGraph(std::stringstream &asqg_stream, bool exact);

    SeqLib::GenomicRegion m_region;
    SeqLib::BamReader m_bam;
    BxBamWalker m_bx_bam;
    BamReadVector m_reads;
    AssemblyParams m_params;
    ReadTable m_read_table;
    std::string m_prefix;
    SeqLib::UnalignedSequenceVector m_contigs;
};

#endif
