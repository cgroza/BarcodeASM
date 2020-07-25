#ifndef LOCAL_ASSEMBLY_WINDOW_H
#define LOCAL_ASSEMBLY_WINDOW_H

#include "BxBamWalker.h"
#include "SeqLib/BamReader.h"
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
#include <vector>

struct AssemblyParams {
    size_t min_overlap = 90;
    size_t min_contig_length = 200;
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
    void sortContigs();

    AssemblyParams m_params;
    SeqLib::GenomicRegion m_region;
    SeqLib::BamReader m_bam;
    BxBamWalker m_bx_bam;
    BamReadVector m_reads;
    std::string m_prefix;
    SeqLib::UnalignedSequenceVector m_contigs;
};

#endif
