# Main make file for the pindel project
# Include the local configuration
-include Makefile.local

default: paraPindel

all: paraPindel cppcheck functional-tests coverage-tests acceptance-tests \
	regression-tests
test: paraPindel cppcheck functional-tests

paraPindel: Makefile.local
	make -C src paraPindel

paraPindel-debug: Makefile.local
	make -C src paraPindel-debug

cppcheck: Makefile.local
	make -C src test

acceptance-tests: Makefile.local paraPindel
	make -C test acceptance-tests

coverage-tests: Makefile.local paraPindel-debug
	make -C test coverage-tests

functional-tests: Makefile.local paraPindel
	make -C test functional-tests

regression-tests: Makefile.local paraPindel
	make -C test regression-tests

clean:
	make -C src clean
	make -C test clean

Makefile.local:
	@echo '# Local configuration' > $@
	@echo '# Location of HTSlib' >> $@
	@if [ -z "$(HTSLIB)" ]; then \
	     echo "HTSLIB_CPPFLAGS=" >> $@; \
	 elif [ -d $(HTSLIB)/include ]; then \
	     echo "HTSLIB_CPPFLAGS=-I$(realpath $(HTSLIB)/include)" >> $@; \
	 else \
	     echo "HTSLIB_CPPFLAGS=-I$(realpath $(HTSLIB))" >> $@; \
	 fi
	@if [ -z "$(HTSLIB)" ]; then \
	     echo "HTSLIB_LDFLAGS=" >> $@; \
	 elif [ -d $(HTSLIB)/lib ]; then \
	     echo "HTSLIB_LDFLAGS=-L$(HTSLIB)/lib -Wl,-rpath $(HTSLIB)/lib" >> $@; \
	 else \
	     echo "HTSLIB_LDFLAGS=-L$(HTSLIB) -Wl,-rpath $(HTSLIB)" >> $@; \
	 fi
	@echo '' >> $@
	@echo '# Number of threads for functional tests, set to 2 or more, recommended to match number of cores' >> $@
	@(if [ -e /proc/cpuinfo ] ; then THREADS=`fgrep -c processor /proc/cpuinfo` ; echo "THREADS=${THREADS}" ; else echo 'THREADS=2' ; fi) >> $@
	@echo '' >> $@
	@echo '# Acceptance test tuning variables (seconds), set to realistic values for your system' >> $@
	@echo '# Numbers based on running in CI on Intel i7 2.8GHz, 8 cores, 24GB RAM' >> $@
	@echo 'COLOUSINGBD15_TIME=60' >> $@
	@echo 'COLOWOBD15_TIME=80' >> $@
	@echo 'SIM1CHRVS20305_TIME=60' >> $@
	@false

# Pseudo targets for configuration
.PHONY: default all clean test paraPindel paraPindel-debug cppcheck acceptance-tests \
	coverage-tests functional-tests regression-tests
.PRECIOUS: Makefile.local
.IGNORE: clean
