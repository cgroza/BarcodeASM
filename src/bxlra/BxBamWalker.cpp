#include "BxBamWalker.h"
#include <iostream>
#include <stdexcept>
#include <system_error>

BxBamWalker::BxBamWalker(const std::string &bx_bam_path,
                         const std::string _prefix,
                         bool _weird_reads_only)
    : prefix(_prefix), weird_reads_only(_weird_reads_only) {

    BamReader();
    Open(bx_bam_path);
}

BxBamWalker::BxBamWalker() : weird_reads_only(true) { BamReader(); }

BamReadVector
    BxBamWalker::fetchReadsByBxBarcode(const BxBarcode &bx_barcode) {
    BamReadVector read_vector;
    // We must convert the stirng barcode into an index ID. We retrieve this ID
    // from the BAM header.
    SeqLib::BamHeader header = Header();
    SeqLib::GenomicRegion bx_region(bx_barcode, "1", "2", header);

    // Check if this barcode exists in this BxBamWalker
    if (header.Name2ID(bx_barcode) < 0) {
        return read_vector;
    }

    SetRegion(bx_region);

    // BX tags are sorted and indexed in large contiguous blocks within the BAM.
    // Each BX tag is a region. Therefore, we stop when we exhaust the region.
    while (true) {
        // Careful. We cannot reuse BamRecords since we must avoid pushing shallow
        // copies into the BamReadVector.
        SeqLib::BamRecord bx_record;

        // Are we still within the same BX block?
        if (GetNextRecord(bx_record)) {
            // do we only want weird reads?
            if (weird_reads_only) {
                if(isBxReadWeird(bx_record))
                    read_vector.push_back(bx_record);
            } else
                read_vector.push_back(bx_record);
        } else
            break;
    }
    return read_vector;
}

BamReadVector BxBamWalker::fetchReadsByBxBarcode(
                                              const std::set<BxBarcode> &bx_barcodes) {
    BamReadVector all_reads;
    for (const BxBarcode &barcode : bx_barcodes) {
        std::string barcode_copy = barcode;
        // "-" gets translated to "_" during bxtools convert. We must correct this.
        std::replace(barcode_copy.begin(), barcode_copy.end(), '-', '_');
        BamReadVector reads = fetchReadsByBxBarcode(barcode_copy);
        // move the reads into all_reads
        std::move(reads.begin(), reads.end(), std::back_inserter(all_reads));
    }
    return all_reads;
}

std::set<BxBarcode>
    BxBamWalker::collectBxBarcodes(const BamReadVector &reads) {
    std::set<BxBarcode> barcodes;
    for (const SeqLib::BamRecord &r : reads) {
        std::string bx_tag;
        // tag may not always be present
        if (r.GetZTag("BX", bx_tag)) {
            barcodes.insert(bx_tag);
        }
    }
    return barcodes;
}

bool BxBamWalker::isBxReadWeird(SeqLib::BamRecord &r) {
    // look for unmapped reads, unpaired reads and poor alignments
    // NOTE: the ProperPair flag is part of the SAM specification and is set by
    // the aligner Includes information on insert size, mate read orientation
    // etc...
    bool isWeird = !r.PairMappedFlag() || !r.MappedFlag() ||
        (r.MapQuality() <= POOR_ALIGNMENT_MAX_MAPQ); // ||
    // !r.ProperPair();

    // report which flag was set
#ifdef DEBUG_BX_BAM_WALKER
    if (isWeird)
        std::cerr << !r.PairMappedFlag() << !r.MappedFlag() << !r.ProperPair()
                  << (r.MapQuality() < POOR_ALIGNMENT_MAX_MAPQ) << std::endl;
#endif
    return isWeird;
}
