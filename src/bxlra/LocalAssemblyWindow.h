#ifndef LOCAL_ASSEMBLY_WINDOW_H
#define LOCAL_ASSEMBLY_WINDOW_H

#include "ASQG.h"
#include "BxBamWalker.h"
#include "CorrectionThresholds.h"
#include "EncodedString.h"
#include "OverlapAlgorithm.h"
#include "OverlapCommon.h"
#include "RLBWT.h"
#include "ReadTable.h"
#include "SGACommon.h"
#include "SGSearch.h"
#include "SGUtil.h"
#include "SGVisitors.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/GenomicRegion.h"
#include "SeqLib/SeqLibUtils.h"
#include "SeqLib/UnalignedSequence.h"
#include "SuffixArray.h"
#include "Util.h"
#include "overlap.h"
#include <algorithm>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <ostream>

struct AssemblyParams {
    double error_rate = 0;
    int seed_length = 0;
    int seed_stride = 0;
    bool irr_only = true;
    bool get_components = true;
    size_t min_overlap = 90;
    bool perform_trim = true;
    size_t trim_length_threshold = 200;
    size_t trim_rounds = 3;
    bool validate = true;
    size_t min_contig_length = 200;
    size_t walk_max_distance = 100000;
    size_t min_repeat_size = 20;
};


class LocalAssemblyWindow {

public:
    LocalAssemblyWindow(SeqLib::GenomicRegion region, SeqLib::BamReader bam, BxBamWalker bx_bam, AssemblyParams params);
    size_t retrieveGenomewideReads();
    size_t assembleReads();
    BxBarcodeCounts collectLocalBarcodes();
    SeqLib::UnalignedSequenceVector getContigs() const;
    BamReadVector getReads() const;
    void clearReads();
    void writeContigs(std::ostream &out);


  private:
    void createReadTable();
    void assembleFromGraph(std::stringstream &asqg_stream, bool exact);
    void sortContigs();
    void walkAssemblyGraph(StringGraph *str_graph);

    AssemblyParams m_params;
    SeqLib::GenomicRegion m_region;
    SeqLib::BamReader m_bam;
    BxBamWalker m_bx_bam;
    BamReadVector m_reads;
    ReadTable m_read_table;
    std::string m_prefix;
    SeqLib::UnalignedSequenceVector m_contigs;
};

#endif
