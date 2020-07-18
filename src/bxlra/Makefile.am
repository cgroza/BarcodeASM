bin_PROGRAMS = bxlra

bxlra_CPPFLAGS = \
	-I$(top_srcdir)/sga/src/Util \
	-I$(top_srcdir)/sga/src/Bigraph \
	-I$(top_srcdir)/sga/src/SuffixTools \
	-I$(top_srcdir)/sga/src/StringGraph \
	-I$(top_srcdir)/sga/src/Algorithm \
	-I$(top_srcdir)/sga/src/SQG \
	-I$(top_srcdir)/sga/src/SGA \
	-I$(top_srcdir)/SeqLib \
	-I$(top_srcdir)/SeqLib/htslib -Wno-sign-compare

bxlra_LDADD = \
	$(top_builddir)/minimap2/libminimap2.a \
	$(top_builddir)/sga/src/SGA/libsga.a \
	$(top_builddir)/sga/src/StringGraph/libstringgraph.a \
	$(top_builddir)/sga/src/Algorithm/libalgorithm.a \
	$(top_builddir)/sga/src/SuffixTools/libsuffixtools.a \
	$(top_builddir)/sga/src/Bigraph/libbigraph.a \
	$(top_builddir)/sga/src/Util/libutil.a \
	$(top_builddir)/sga/src/SQG/libsqg.a \
	$(top_builddir)/SeqLib/src/libseqlib.a \
	$(top_builddir)/SeqLib/bwa/libbwa.a \
	$(top_builddir)/SeqLib/htslib/libhts.a \
	$(top_builddir)/SeqLib/fermi-lite/libfml.a \
	-llzma -lbz2 -lz

bxlra_SOURCES = bxlra.cpp BxBamWalker.cpp RegionFileReader.cpp LocalAssemblyWindow.cpp LocalAlignment.cpp

install:
	mkdir -p ../../bin && mv bxlra ../../bin
