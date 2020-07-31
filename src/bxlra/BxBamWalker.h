#ifndef BX_BAM_WALKER_H
#define BX_BAM_WALKER_H

#include "SeqLib/BFC.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/BamRecord.h"
#include "SeqLib/GenomicRegion.h"
#include <algorithm>
#include <exception>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <system_error>
#include <unordered_map>
#include <vector>

/* Let's distinguish regular strings from BxBarcodes in the source code */
typedef std::string BxBarcode;
typedef std::vector<SeqLib::BamRecord> BamReadVector;

class BxBamWalker : public SeqLib::BamReader {
    /* Reads a BAM file that was produced by the lariat aligner. This file must be
       prepared by flipping the chromosome and BX tag with bxtools convert. Then
       this file must be sorted and indexed by the BX tag using samtools. This
       allows fast retrieval of reads having a specific bx tag.
    */

    public:
    /* bx_bam_path: BAM file indexed by bx tag. */
    BxBamWalker();
    BxBamWalker(const std::string &bx_bam_path,
                     const std::string _prefix = "0000",
                     bool _weird_reads_only = true);

    /*  */
    BamReadVector fetchReadsByBxBarcode(const BxBarcode &bx_barcode);
    BamReadVector fetchReadsByBxBarcode(const std::vector<BxBarcode> &bx_barcodes);
    std::string prefix;

    static bool isBxReadWeird(SeqLib::BamRecord &r);
    static const int POOR_ALIGNMENT_MAX_MAPQ = 10;

    private:
    bool weird_reads_only;

};

#endif
