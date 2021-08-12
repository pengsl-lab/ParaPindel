include ../Makefile.local

CXX = mpicxx
SOURCES=paraPindel.cpp reader.cpp reporter.cpp searcher.cpp parameter.cpp \
    refreader.cpp control_state.cpp search_deletions_nt.cpp \
    search_inversions.cpp search_inversions_nt.cpp bam2depth.cpp\
    search_tandem_duplications.cpp search_tandem_duplications_nt.cpp \
    output_sorter.cpp farend_searcher.cpp search_variant.cpp \
    searchshortinsertions.cpp searchdeletions.cpp output_file_data.cpp \
    bddata.cpp shifted_vector.cpp read_buffer.cpp \
    line_reader.cpp ifstream_line_reader.cpp gz_line_reader.cpp\
    pindel_read_reader.cpp user_defined_settings.cpp \
    fn_parameters.cpp logstream.cpp search_MEI_util.cpp \
    search_MEI.cpp assembly.cpp genotyping.cpp

 
all:	paraPindel paraPindel-debug
-include $(SOURCES:.cpp=.d)
-include $(SOURCES:.cpp=-debug.d)
test:	paraPindel $(addprefix cppcheck-result-,$(SOURCES:.cpp=.xml))

pindel-debug: $(SOURCES:.cpp=-debug.o)
	$(CXX) $(LDFLAGS) $^ -fopenmp -g -pg -fprofile-arcs -ftest-coverage \
	  $(HTSLIB_LDFLAGS) -lhts -lm -lz -o $@

paraPindel:	$(SOURCES:.cpp=.o)
	$(CXX) $(LDFLAGS) $(HTSLIB_LDFLAGS) $^ -O3 -fopenmp -fPIE -pg -lhts -lm -lz -o $@
	$(CXX) -O3 -w -fopenmp -pg  SAM_2_PINDEL_cin_2011Sept26.cpp -o  sam2pindel
	$(CXX) -O3 -w -fopenmp -pg  pindel2vcf.cpp -o pindel2vcf
	$(CXX) -O3 -w -fopenmp -pg  pindel2vcf4tcga.cpp -o pindel2vcf4tcga
clean:
	$(RM) paraPindel pindel2vcf pindel2vcf4tcga pindel-debug $(OBJECTS) $(OBJECTS_DEBUG) $(CPPCHECKS) \
	    *.gcov *.gcda *.gcno gmon.out *.d *.o


%.o: %.cpp
	$(CXX) $(MPICXXFLAGS) $(HTSLIB_CPPFLAGS) -Wall -g -c -O3 -fopenmp   $< -o $@

%-debug.o: %.cpp
	$(CXX) $(MPICXXFLAGS) -Wall -g -c -fopenmp -g -pg  \
	  -fprofile-arcs -ftest-coverage $(HTSLIB_CPPFLAGS) $< -o $@


cppcheck-result-%.xml: %.cpp
	cppcheck --quiet --enable=all --xml $(HTSLIB_CPPFLAGS) $(realpath $<) 2> $@


%.d: %.cpp
	@set -e; rm -f $@; \
	    $(CXX) -MM $(CPPFLAGS) $(HTSLIB_CPPFLAGS) $< > $@.$$$$; \
	    sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	    $(RM) $@.$$$$

%-debug.d: %.cpp
	@set -e; rm -f $@; \
	    $(CXX) -MM $(CPPFLAGS) $(HTSLIB_CPPFLAGS) $< > $@.$$$$; \
	    sed 's,\($*\)\.o[ :]*,\1-debug.o $@ : ,g' < $@.$$$$ > $@; \
	    $(RM) $@.$$$$

.PHONY: all clean test
.IGNORE: clean
