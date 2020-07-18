bin_PROGRAMS = bxlra

bxlra_CPPFLAGS = \
	-I$(top_srcdir)/src/SGA/Util \
	-I$(top_srcdir)/src/SGA/Bigraph \
	-I$(top_srcdir)/src/SGA/SuffixTools \
	-I$(top_srcdir)/src/SGA/StringGraph \
	-I$(top_srcdir)/src/SGA/Algorithm \
	-I$(top_srcdir)/src/SGA/SQG \
	-I$(top_srcdir)/src/SGA/SGA \
	-I$(top_srcdir)/SeqLib \
	-I$(top_srcdir)/SeqLib/htslib -Wno-sign-compare

bxlra_LDADD = \
	$(top_builddir)/minimap2/libminimap2.a \
	$(top_builddir)/src/SGA/SGA/libsga.a \
	$(top_builddir)/src/SGA/StringGraph/libstringgraph.a \
	$(top_builddir)/src/SGA/Algorithm/libalgorithm.a \
	$(top_builddir)/src/SGA/SuffixTools/libsuffixtools.a \
	$(top_builddir)/src/SGA/Bigraph/libbigraph.a \
	$(top_builddir)/src/SGA/Util/libutil.a \
	$(top_builddir)/src/SGA/SQG/libsqg.a \
	$(top_builddir)/SeqLib/src/libseqlib.a \
	$(top_builddir)/SeqLib/bwa/libbwa.a \
	$(top_builddir)/SeqLib/htslib/libhts.a \
	$(top_builddir)/SeqLib/fermi-lite/libfml.a \
	-llzma -lbz2 -lz

bxlra_SOURCES = bxlra.cpp BxBamWalker.cpp RegionFileReader.cpp LocalAssemblyWindow.cpp LocalAlignment.cpp

install:
	mkdir -p ../../bin && mv bxlra ../../bin
