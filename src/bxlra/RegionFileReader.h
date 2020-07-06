#ifndef REGION_FILE_READER_H
#define REGION_FILE_READER_H

#include "SeqLib/BamHeader.h"
#include "SeqLib/GenomicRegion.h"
#include <fstream>
#include <string>

class RegionFileReader {

public:
    RegionFileReader(const std::string &path, SeqLib::BamHeader header);

    SeqLib::GenomicRegionVector getRegions() const;

private:
    std::string m_path;
    SeqLib::GenomicRegionVector m_regions;
    SeqLib::BamHeader m_header;
};

#endif
