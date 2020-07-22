#include "LocalAssemblyWindow.h"

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
  std::cerr << std::endl;
  size_t total = 0;
  for (const auto &b : barcodes) {
    std::cerr << b.first << " " << b.second << std::endl;
    total += b.second;
  }
  std::cerr << "Sum " << total << std::endl;

  // add the genome wide reads to assembly
  BamReadVector genomewide_reads = m_bx_bam.fetchReadsByBxBarcode(barcodes);
  std::copy(genomewide_reads.begin(), genomewide_reads.end(), std::back_inserter(m_reads));
  std::cerr << "Post barcode collection: " << m_reads.size() << std::endl;
  return m_reads.size();
}

size_t LocalAssemblyWindow::assembleReads() {
  retrieveGenomewideReads();
  createReadTable();
  // forward suffix index
  SuffixArray suffix_array_f(&m_read_table, 1, false);
  RLBWT BWT_f(&suffix_array_f, &m_read_table);

  // reverse suffix index
  m_read_table.reverseAll();
  SuffixArray suffix_array_r(&m_read_table, 1, false);
  RLBWT BWT_r(&suffix_array_r, &m_read_table);

  // reverse the read table back to its forward state
  m_read_table.reverseAll();

  suffix_array_f.writeIndex();
  suffix_array_r.writeIndex();

  // TODO: learn read overlap parameters. 

  OverlapAlgorithm overlapper(&BWT_f, &BWT_r, m_params.error_rate,
                              m_params.seed_length, m_params.seed_stride,
                              m_params.irr_only);

  bool exact = m_params.error_rate < 0.0001;
  overlapper.setExactModeOverlap(exact);
  overlapper.setExactModeIrreducible(exact);

  // write ASQG in memory, dump later.
  std::stringstream asqg_writer;

  // Write ASGQ header
  ASQG::HeaderRecord header_record;
  header_record.setOverlapTag(m_params.min_overlap);
  header_record.setErrorRateTag(m_params.error_rate);
  header_record.setContainmentTag(true);
  header_record.setTransitiveTag(!m_params.irr_only);
  header_record.write(asqg_writer);

  // reset read table iterator
  m_read_table.setZero();

  size_t workid = 0;
  SeqItem si;

  std::stringstream hits_stream;
  // find read overlaps for each read
  while (m_read_table.getRead(si)) {
    SeqRecord read;
    read.id = si.id;
    read.seq = si.seq;
    OverlapBlockList obl;

    OverlapResult rr = overlapper.overlapRead(read, m_params.min_overlap, &obl);
    overlapper.writeOverlapBlocks(hits_stream, workid, rr.isSubstring, &obl);

    // write a vertex for this read
    ASQG::VertexRecord record(read.id, read.seq.toString());
    record.setSubstringTag(rr.isSubstring);
    record.write(asqg_writer);

    ++workid;
  }

  // parse hits and write edges for each read vertex
  std::string line;
  bool b_is_self_compare = true;
  ReadInfoTable query_read_info_table(&m_read_table);

  while (std::getline(hits_stream, line)) {
    size_t read_idx;
    size_t total_entries;
    bool is_substring;
    OverlapVector ov;
    OverlapCommon::parseHitsString(line, &query_read_info_table,
                                   &query_read_info_table, &suffix_array_f,
                                   &suffix_array_r, b_is_self_compare, read_idx,
                                   total_entries, ov, is_substring);

    for (OverlapVector::iterator iter = ov.begin(); iter != ov.end(); ++iter) {
      ASQG::EdgeRecord edge_record(*iter);
      edge_record.write(asqg_writer);
    }
  }

  // dump raw overlap graph to file
  std::ofstream asqg_file_writer;
  asqg_file_writer.open(m_prefix + "_graph.asqg");
  asqg_file_writer << asqg_writer.str();
  asqg_file_writer.close();

  // pre-process the overlap graph, retrieve contigs
  assembleFromGraph(asqg_writer, exact);
  sortContigs();
  return m_contigs.size();
}

void LocalAssemblyWindow::assembleFromGraph(std::stringstream &asqg_stream,
                                            bool exact) {
  // configure string graph object
  StringGraph *str_graph =
      SGUtil::loadASQG(asqg_stream, m_params.min_overlap, true);
  str_graph->m_get_components = m_params.get_components;
  str_graph->setExactMode(exact);

  // The returned contigs are enumerations of walks in the ASGQ graph
  // The following preprocessing steps remove unwanted walks

  SGTransitiveReductionVisitor trans_visitor;
  SGGraphStatsVisitor stat_visitor;
  SGTrimVisitor trim_visitor(m_params.trim_length_threshold);
  SGContainRemoveVisitor contain_visitor;
  SGValidateStructureVisitor validate_visitor;

  // ?????
  while (str_graph->hasContainment())
    str_graph->visit(contain_visitor);

  // removes redundant paths from the graph
  str_graph->visit(trans_visitor);


  str_graph->simplify(); // merges vertices by removing transitive edges

  if (m_params.validate)
    str_graph->visit(validate_visitor);

  // Remove branches that do not merge into to form a bubble
  if (m_params.perform_trim) {
    for (size_t i = 0; i < m_params.trim_rounds; i++)
      str_graph->visit(trim_visitor);
  }

  // identify these vertices with a unique prefix
  str_graph->renameVertices(m_prefix + "_");

  // walk graph to retrieve contigs
  walkAssemblyGraph(str_graph);

  str_graph -> writeASQG(m_prefix + "_pruned_graph.asqg");
  delete str_graph;
}

void LocalAssemblyWindow::walkAssemblyGraph(StringGraph *str_graph) {
    // this code is adapted from SGA
    typedef std::vector<VertexPtrVec> ComponentVector;
    VertexPtrVec all_vertices = str_graph -> getAllVertices();
    ComponentVector components;
    SGSearchTree::connectedComponents(all_vertices, components);

    // Select the largest component
    int selected_idx = -1;
    size_t largestSize = 0;

    for(size_t i = 0; i < components.size(); ++i)
    {
        std::cerr << "Component " << i << ": " << components[i].size() << " vertices" << std::endl;
        if(components[i].size() > largestSize)
        {
            selected_idx = i;
            largestSize = components[i].size();
        }
    }

    assert(selected_idx != -1);
    VertexPtrVec selected_component = components[selected_idx];

    std::cerr << "component-walk: selected component of size " << selected_component.size() << std::endl;

    // Build a vector of the terminal vertices
    VertexPtrVec terminals;
    for(size_t i = 0; i < selected_component.size(); ++i)
    {
        Vertex* vertex = selected_component[i];
        size_t antisense_count = vertex->getEdges(ED_ANTISENSE).size();
        size_t sense_count = vertex->getEdges(ED_SENSE).size();

        if(antisense_count == 0 || sense_count == 0)
            terminals.push_back(vertex);
    }

    std::cout << "selected component has " << terminals.size() << " terminal vertices" << std::endl;

    // Find walks between all-pairs of terminal vertices
    SGWalkVector temp_walks;
    for(size_t i = 0; i < terminals.size(); ++i)
    {
        for(size_t j = i + 1; j < terminals.size(); j++)
        {
            Vertex* x = terminals[i];
            Vertex* y = terminals[j];
            SGSearch::findWalks(x, y, ED_SENSE, m_params.walk_max_distance, 1000000, false, temp_walks);
            SGSearch::findWalks(x, y, ED_ANTISENSE, m_params.walk_max_distance, 1000000, false, temp_walks);
        }
    }

    // Remove duplicate walks
    std::map<std::string, SGWalk> walk_map;
    for(size_t i = 0; i < temp_walks.size(); ++i)
    {
        std::string walk_string = temp_walks[i].getString(SGWT_START_TO_END);
        walk_map.insert(std::make_pair(walk_string, temp_walks[i]));
    }

    // Copy unique walks to the output
    size_t count = 0;
    for(std::map<std::string, SGWalk>::iterator map_iter = walk_map.begin(); map_iter != walk_map.end(); ++map_iter) {
        std::stringstream ss;
        ss << m_prefix << "_" << count;
        m_contigs.push_back(SeqLib::UnalignedSequence(ss.str(), map_iter -> first));
        ++count;
    }
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

void LocalAssemblyWindow::createReadTable() {
  size_t count = 0;
  for (auto &i : m_reads) {
    std::string seq = i.Sequence();
    // ignore uncalled bases or reads that are too short
    // IMPORTANT: cannot build suffix array/index when reads have uncalled bases
    if (seq.find("N") != std::string::npos || seq.length() < m_params.min_overlap)
        continue;
    std::string seq_id = std::to_string(++count);

    if (!i.MappedFlag() && !i.MateReverseFlag())
      SeqLib::rcomplement(seq);

    SeqItem seq_item;
    seq_item.seq = seq;
    seq_item.id = seq_id;
    m_read_table.addRead(seq_item);
  }

  // duplicate removal

  //construct an index and check each read for duplicates
  // against this index
  SuffixArray forward_suffix_array(&m_read_table, 1, false); 
  RLBWT forward_BWT(&forward_suffix_array, &m_read_table);

  // reverse
  m_read_table.reverseAll();
  SuffixArray reverse_suffix_array(&m_read_table, 1, false);
  RLBWT reverse_RBWT(&reverse_suffix_array, &m_read_table);
  m_read_table.reverseAll();

  OverlapAlgorithm dup_overlapper(&forward_BWT, &reverse_RBWT, 0, 0, 0, false);

  m_read_table.setZero();
  ReadTable undup_read_table;
  SeqItem sir;
  while (m_read_table.getRead(sir)) {
    OverlapBlockList overlap_block;
    SeqRecord read;
    read.id = sir.id;
    read.seq = sir.seq;
    // check duplicate
    OverlapResult rr = dup_overlapper.alignReadDuplicate(read, &overlap_block);

    if (!rr.isSubstring)
      undup_read_table.addRead(sir);
  }
  m_read_table = undup_read_table;
}


void LocalAssemblyWindow::sortContigs() {
    // sort contigs in decreasing sequence length order
    std::sort(m_contigs.begin(), m_contigs.end(),
              [](SeqLib::UnalignedSequence &a, SeqLib::UnalignedSequence &b) {
                  return a.Seq.length() > b.Seq.length();
              });
}

void LocalAssemblyWindow::writeContigs(std::ostream &out) {
    for(auto &contig : m_contigs) {
        out << ">" << contig.Name << std::endl;
        out << contig.Seq << std::endl;
    }
}
