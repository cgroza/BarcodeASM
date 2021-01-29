#include "RegionFileReader.h"

RegionFileReader::RegionFileReader(const std::string &path, SeqLib::BamHeader header) : m_path(path), m_header(header) {
    // assume the file is a BED file with chrom, start, end fields
    std::ifstream file;
    file.open(path);
    std::string chrom;
    std::string start;
    std::string end;

    while(file >> chrom >> start >> end) {
        m_regions.push_back(SeqLib::GenomicRegion(chrom, start, end, header));
    }

}

SeqLib::GenomicRegionVector RegionFileReader::getRegions() const {
    return m_regions;
}
