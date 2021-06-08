/*
 * This File is part of Pindel; a program to locate genomic variation.
 * https://trac.nbic.nl/pindel/
 *
 *   Copyright (C) 2011 Kai Ye
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

// System header files
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <getopt.h>
#include <omp.h>
#include <set>
#include <map>
#include <sstream>
#include <ctime>
#include <mpi.h>
#include <vector>
#include <memory>
#include <assert.h>
#include <cstdlib>

#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <unistd.h>

//#include <mpi.h>
// Pindel header files
#include "logstream.h"
#include "pindel.h"
#include "fn_parameters.h"
#include "bddata.h"
#include "reader.h"
#include "refreader.h"
#include "searcher.h"
#include "reporter.h"
#include "parameter.h"
#include "control_state.h"
#include "search_deletions_nt.h"
#include "search_inversions.h"
#include "search_inversions_nt.h"
#include "search_MEI.h"
#include "search_tandem_duplications.h"
#include "search_tandem_duplications_nt.h"
#include "read_buffer.h"
#include "line_reader.h"
#include "pindel_read_reader.h"
#include "farend_searcher.h"
#include "searcher.h"
#include "search_variant.h"
#include "searchshortinsertions.h"
#include "searchdeletions.h"
#include "user_defined_settings.h"
#include "logdef.h"
#include "assembly.h"
//#include "genotyping.h"
#include "ifstream_line_reader.h"
#include "gz_line_reader.h"

#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "htslib/kseq.h"
#include "htslib/khash.h"
#include "htslib/ksort.h"

#define NUM_CHR 3000
#define MAX_NUM_BAM_FILES 10
#define LEN_BAMFILENAME 200
#define LEN_BAMTAG 40

typedef struct WindowInfo {
	std::string chrName;
	int bed_index;
	int win_start;
	int win_end;
	unsigned int global_start;
	unsigned int global_end;
	int isused;
} WindowInfo;

std::vector<WindowInfo> windowSet;
typedef struct {
	int tasks[NUM_CHR];
	int size;
} OneNode;

#define handle_error(msg) \
    do { perror(msg); exit(EXIT_FAILURE); } while (0)

clock_t time_start, time_end;

/*v Kai Ye update 0.2.4h, Oct 31 2011, update for MOSAIK */
/*v EW update 0.2.4j, Pindel will now abort when insert size is set too small. */
/*v Kai/EW update 0.2.4k, pindel will now give the consensus inserted sequence instead of always the first one. */
/*v EW update 0.2.4l, -Werror removed for user convenience, also made some typecasts in reader.cpp explicit; also reset -l and -k options to true by default. */
/*v Kai update 0.2.4m, removed an error in invertion calling. */
/* EW update 0.2.4n, improved VCF creator: less memory, more speed, removed bug. */
/* EW/Kai/Matthijs: update 0.2.4o: does not report short inversions as deletions anymore, also displays the reference for short inversions correctly, slightly changed sensitivity, small memory leak fixed. */
/* Kai/EW: update 0.2.4p: sorts reads in output, number of unique calls should be correct now, pindel now gives error if using config file that does not exist or has other problems */
/* EW: update 0.2.4q: also works with -c all instead of -c ALL */
/* Kai: update 0.2.4q: use ChrName and ChrSeq for clarity; start to include assembly */
/* Kai: min_perfect_match_around_BP to 5 */
/* EW: update 0.2.4s: bugfix for -p option of Pindel0.2.4r */
/* EW: update 0.2.4t, updates now shown in RELEASE document in trunk directory */

const std::string Pindel_Version_str = "Pindel version 0.2.5b9, 20160729.";

const Chromosome g_dummyChromosome("", "");
Genome g_genome;
std::ofstream g_logFile;
std::vector<ChrNameAndSizeAndIndex> g_ChrNameAndSizeAndIndex;

int g_binIndex = -1; // global variable for the bin index, as I cannot easily pass an extra parameter to the diverse functions
unsigned int g_maxPos = 0; // to calculate which is the last position in the chromosome, and hence to calculate the number of bins
short g_MinClose = 8;
std::set<std::string> g_sampleNames;
std::map<std::string, unsigned> g_SampleName2Index;
std::map<std::string, unsigned> g_ChrName2Ploidy;

short Before, After;
std::vector<RefCoveragePerPosition> g_RefCoverageRegion;
BDData g_bdData;

const int alphs = 4;
const char alphabet[alphs] = { 'A', 'C', 'G', 'T' };
unsigned int BoxSize = 10000; // 10k is fine
const double Min_Filter_Ratio = 0.5;
unsigned int NumberOfSIsInstances = 0;
unsigned int NumberOfDeletionsInstances = 0;
unsigned int NumberOfDIInstances = 0;
unsigned int g_numberOfInvInstances = 0;
unsigned int NumberOfTDInstances = 0;
short g_reportLength = 1;
unsigned g_RegionStart, g_RegionEnd;
char Match2N[256];
char Match[256];
char Convert2RC[256];
char Convert2RC4N[256];
char Cap2LowArray[256];
bool FirstChr = true;
unsigned int DSizeArray[15];
int g_maxInsertSize = 0;
unsigned g_NumberOfGapAlignedReads = 0;
std::vector<Parameter*> parameters;
std::string f = "";
std::string CurrentChrMask = "";
UserDefinedSettings *userSettings;
std::string Spacer = "";
double loadChrs = 0, readingAndMatching = 0, searchingFarEnds = 0, searchingSVs = 0, total_time = 0;
// #########################################################

int NumRead2ReportCutOff_BP = 2;
const int g_MAX_RANGE_INDEX = 9; // 5 or 6 or 7 or maximum 8      //# // user
unsigned int WINDOW_SIZE = 10000000;
unsigned int CROSS_WINDOW_SIZE = 5000;
const unsigned int AROUND_REGION_BUFFER = 10000; // how much earlier reads should be selected if only a region of the chromosome needs be specified.
// #########################################################

const short MaxDI = 30;

int cabs(int a) {
	if (a >= 0)
		return a;
	else
		return (-1) * a;
}

// Note: in case one needs to handle differnt line delimuters (like crlf)
// from http://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf
void safeGetline(std::istream &is, std::string &t) {
	t.clear();
	getline(is, t);
}

void SPLIT_READ::setUnmatchedSeq(const std::string &unmatchedSeq) {
	//const double EPSILON = 0.00001; // to compensate for downrounding off errors (0.04 = 0.03999999 on some computers)

	UnmatchedSeq = unmatchedSeq;
	if (UnmatchedSeq.size() > 0) {
		/*for (unsigned int x=0; x<UnmatchedSeq.size();x++ ) {
		 std::cout << "seq[" << x << "]="<< int(UnmatchedSeq[ x ]) << "('" << UnmatchedSeq[ x ] << "')\n";
		 }*/
		unsigned int lastCharIndex = UnmatchedSeq.size() - 1;
		while (!isalnum(UnmatchedSeq[lastCharIndex])) {
			UnmatchedSeq.resize(lastCharIndex);
			lastCharIndex--;
		}
		/*for (unsigned int x=0; x<UnmatchedSeq.size();x++ ) {
		 std::cout << "resseq[" << x << "]="<< int(UnmatchedSeq[ x ]) << "('" << UnmatchedSeq[ x ] << "')\n";
		 }*/
	}

	ReadLength = UnmatchedSeq.size();
	ReadLengthMinus = ReadLength - 1;

	//UserDefinedSettings* userSettings = UserDefinedSettings::Instance();

	MAX_SNP_ERROR = g_maxMismatch[ReadLength];
	TOTAL_SNP_ERROR_CHECKED_Minus = MAX_SNP_ERROR + userSettings->ADDITIONAL_MISMATCH;
	TOTAL_SNP_ERROR_CHECKED = TOTAL_SNP_ERROR_CHECKED_Minus + 1;
}

const Chromosome* Genome::addChromosome(Chromosome *newChromosome) {
	for (unsigned int i = 0; i < m_chromosomes.size(); i++) {
		if (m_chromosomes[i]->getName() == newChromosome->getName()) {
			delete m_chromosomes[i];
			//std::cout << "delete one chromsome." << std::endl;
			m_chromosomes[i] = newChromosome;
			return newChromosome;
		}
	}
	m_chromosomes.push_back(newChromosome);
	return newChromosome;
}

const Chromosome* Genome::getChr(unsigned int index) const {
	if (index < m_chromosomes.size()) {
		return m_chromosomes[index];
	} else {
		return NULL;
	}
}

const Chromosome* Genome::getChr(const std::string &chromosomeName) const {
	for (unsigned int i = 0; i < m_chromosomes.size(); i++) {
		if (m_chromosomes[i]->getName() == chromosomeName) {
			return m_chromosomes[i];
		}
	}
	return NULL;
}

short Genome::getChrID(const std::string &chromosomeName) {
	for (unsigned int i = 0; i < m_chromosomes.size(); i++) {
		if (m_chromosomes[i]->getName() == chromosomeName) {
			return i;
		}
	}
	return -1;
}

void Genome::clear() {
	for (unsigned int i = 0; i < m_chromosomes.size(); i++) {
		delete m_chromosomes[i];
	}
	m_fullMode = false;
	m_currentChromosomeIndex = -1;
	m_referenceFile.close();
}

void Genome::reset() {
	m_referenceFile.clear();
	m_referenceFile.seekg(0, std::ios::beg);
	m_currentChromosomeIndex = -1;
}

void Genome::load(const std::string &referenceFileName) {
	clear();
	m_referenceFile.open(referenceFileName.c_str());
}

/*void Genome::loadAll(const std::string& referenceFileName)
 {
 load( referenceFileName );
 while (!m_referenceFile.eof() ) {
 loadChromosome();
 //std::cout << "m_chromosomes.size() " << m_chromosomes.size() << " " << m_chromosomes[m_chromosomes.size() - 1]->getName() << std::endl;
 }
 m_fullMode = true;
 }*/

void Genome::loadAll(const std::string &referenceFileName) {
	load(referenceFileName);
	while (!m_referenceFile.eof()) {
		loadChromosome();
	}
	m_fullMode = true;

	//
	/*char *addr;
	 int fd;
	 struct stat sb;
	 off_t offset, pa_offset;
	 size_t length;
	 int mmap_flag = 1;

	 fd = open(referenceFileName.c_str(), O_RDONLY);
	 if (fd == -1){
	 //handle_error("open");
	 mmap_flag = 0;
	 }

	 if (fstat(fd, &sb) == -1){    // To obtain reference file size
	 //handle_error("fstat");
	 mmap_flag = 0;
	 }
	 offset = 0;
	 pa_offset = offset & ~(sysconf(_SC_PAGE_SIZE) - 1); //offset for mmap() must be page aligned
	 if (offset >= sb.st_size) {
	 fprintf(stderr, "offset is past end of file\n");
	 mmap_flag = 0;
	 //exit(EXIT_FAILURE);
	 }
	 length = sb.st_size - offset;
	 addr = (char*) mmap(NULL, length + offset - pa_offset, PROT_READ, MAP_PRIVATE, fd, pa_offset);
	 if (addr == MAP_FAILED)
	 handle_error("mmap");

	 unsigned int i, j = 0, len_chrName = 0, len_chrSeq = 0;
	 for (i = 0; i < length + offset - pa_offset; i++) {
	 if (*(addr + i) == '>') {
	 std::string chromosomeName = "";
	 std::string chromosomeSeq = "";
	 j = i + 1;
	 while (*(addr + j) != '\n' && *(addr + j) != ' ') {
	 chromosomeName += *(addr + j++);
	 len_chrName++;
	 }

	 while (*(addr + j++) != '\n') {
	 len_chrName++;
	 }
	 while (*(addr + j) != '>') {
	 char ch = *(addr + j);
	 if (ch != '\n') {
	 ch = toupper(ch);
	 if (ch != 'A' && ch != 'C' && ch != 'G' && ch != 'T') {
	 ch = 'N';
	 }
	 chromosomeSeq += ch;
	 }
	 len_chrSeq++;
	 j++;
	 }
	 //std::cout << chromosomeName << std::endl;
	 chromosomeSeq = Spacer + chromosomeSeq + Spacer;
	 Chromosome *newChromosome = new Chromosome(chromosomeName, chromosomeSeq);
	 addChromosome(newChromosome);
	 }
	 i = i + len_chrName + len_chrSeq;
	 len_chrName = 0;
	 len_chrSeq = 0;
	 }
	 munmap(addr, length + offset - pa_offset);*/

}

const Chromosome* Genome::getNextChromosome() {
	if (m_fullMode) {
		// all chromosomes loaded
		if ((unsigned int) (m_currentChromosomeIndex + 1) < m_chromosomes.size()) {
			m_currentChromosomeIndex++;
			return m_chromosomes[m_currentChromosomeIndex];
		} else {
			return NULL;
		}
	} else {
		// need to load next chromosome
		if (m_currentChromosomeIndex >= 0) {
			// there is a previous chromosome
			const std::string oldName = m_chromosomes[m_currentChromosomeIndex]->getName();
			delete m_chromosomes[m_currentChromosomeIndex];
			m_chromosomes[m_currentChromosomeIndex] = new Chromosome(oldName, "");
		}
		m_currentChromosomeIndex++;
		return loadChromosome();
	}
}

const Chromosome* Genome::loadChromosome() {
	std::string chromosomeName = "";
	std::string chromosomeSeq = "";
	if (m_referenceFile.eof()) {
		return NULL;
	}
	char chrIndicatorChar; // '>' preceeds all chromosome names in a FASTA file
	m_referenceFile >> chrIndicatorChar;
	//std::cout << "-------------" << chrIndicatorChar << std::endl;
	if (chrIndicatorChar != '>') {
		std::cout << "Error: fasta line starts with " << chrIndicatorChar << " instead of '>'. Aborting.\n";
		exit(EXIT_FAILURE);
		return NULL;
	}
	m_referenceFile >> chromosomeName;
	//std::cout << "-----------chrName" << chromosomeName << std::endl;
	std::string restOfChromosomeHeadline;
	safeGetline(m_referenceFile, restOfChromosomeHeadline);
	//std::cout << "--------rest" << restOfChromosomeHeadline << std::endl;
	char ch;
	do {
		m_referenceFile >> ch;
		if (ch == '>') {
			// next chromosome
			m_referenceFile.putback(ch);
			break;
		} else {
			ch = toupper(ch);
			if (ch != 'A' && ch != 'C' && ch != 'G' && ch != 'T') {
				ch = 'N';
			}
			chromosomeSeq += ch;
		}
	} while (!m_referenceFile.eof());

	chromosomeSeq = Spacer + chromosomeSeq + Spacer;
	Chromosome *newChromosome = new Chromosome(chromosomeName, chromosomeSeq);
	//std::cout << "chromosomeName = " <<  chromosomeName << std::endl;
	//if(chromosomeName == "chr1")
//	std::cout << "seq = " << chromosomeSeq << std::endl;
	return addChromosome(newChromosome);
}

UniquePoint::UniquePoint(const Chromosome *chromosome_ptr, const short lengthStr, const unsigned int absLoc, const char direction, const char strand, const short mismatches) :
		chromosome_p(chromosome_ptr), LengthStr(lengthStr), AbsLoc(absLoc), Direction(direction), Strand(strand), Mismatches(mismatches) {
}

SearchWindow& SearchWindow::operator=(const SearchWindow &other) {
	if (this != &other) {
		this->SearchWindow::~SearchWindow(); // explicit non-virtual destructor
		new (this) SearchWindow(other); // placement new
	}
	return *this;
}

bool SearchWindow::encompasses(const std::string &chromosomeName, const unsigned int position) const {
	return ((m_chromosome->getName() == chromosomeName) && (position >= m_currentStart) && (position <= m_currentEnd));
}

SearchWindow::SearchWindow(const Chromosome *chromosome, const unsigned int regionStart, const unsigned int regionEnd) :
		m_chromosome(chromosome) {
	m_currentStart = regionStart;
	m_currentEnd = regionEnd;
}

SearchWindow::SearchWindow(const SearchWindow &other) :
		m_chromosome(other.m_chromosome), m_currentStart(other.m_currentStart), m_currentEnd(other.m_currentEnd) {
}

LoopingSearchWindow::LoopingSearchWindow(const SearchRegion *region, const Chromosome *chromosome, const unsigned int binSize, int Chr_start, int Chr_end) :
		SearchWindow(chromosome, 0, chromosome->getBiolSize()), m_BIN_SIZE(binSize), chr_start(Chr_start), chr_end(Chr_end) {

	bool noOverlapWithChromosome = false;
	if (region->isStartDefined()) {
		m_officialStart = region->getStart();
		if (m_officialStart > chromosome->getBiolSize()) {
			std::cout << "m_officialStart > chromosome->getBiolSize() :" << m_officialStart << ">" << chromosome->getBiolSize() << std::endl;
			noOverlapWithChromosome = true;
		}
		// if the user defines a region, you need to start with reads before that, but not before the start of the chromosome
		m_globalStart = (region->getStart() >= AROUND_REGION_BUFFER ? region->getStart() - AROUND_REGION_BUFFER : 0);
	} else {
		m_officialStart = m_globalStart = 0;
	}

	if (region->isEndDefined()) {
		m_officialEnd = region->getEnd();
		if (m_officialEnd < 0) {
			std::cout << "m_officialEnd < 0 :" << m_officialEnd << "<" << "0" << std::endl;
			noOverlapWithChromosome = true;
		}
		// if the user defines a region, you need to end with reads after that, but not after the end of the chromosome
		m_globalEnd = std::min(chromosome->getBiolSize(), region->getEnd() + AROUND_REGION_BUFFER);
	} else {
		m_officialEnd = m_globalEnd = chromosome->getBiolSize();
	}

	if (noOverlapWithChromosome) {
		std::cout << "Error: the region to scan (" << m_officialStart << ", " << m_officialEnd << ") does not overlap with the " << "chromosome (positions 0 to " << chromosome->getBiolSize() << std::endl;
		exit(EXIT_FAILURE);
	}

	m_currentStart = m_globalStart;
	m_displayedStart = m_officialStart;
	updateEndPositions();

}

/*LoopingSearchWindow::LoopingSearchWindow(const SearchRegion *region, const Chromosome *chromosome, const int binSize, const unsigned Bed_start, const unsigned Bed_end, int step) :
 SearchWindow(chromosome, 0, chromosome->getBiolSize()), m_BIN_SIZE(binSize) {
 m_officialStart = Bed_start;
 m_globalStart = (Bed_start >= AROUND_REGION_BUFFER ? Bed_start - AROUND_REGION_BUFFER : 0);

 m_officialEnd = Bed_end;
 m_globalEnd = std::min(chromosome->getBiolSize(), Bed_end + AROUND_REGION_BUFFER);

 m_currentStart = m_globalStart;
 m_displayedStart = m_officialStart;
 updateEndPositionsbyStep(step);

 }
 LoopingSearchWindow::LoopingSearchWindow(const SearchRegion *region, const Chromosome *chromosome, const int binSize, const unsigned Bed_start, const unsigned Bed_end) :
 SearchWindow(chromosome, 0, chromosome->getBiolSize()), m_BIN_SIZE(binSize) {
 m_officialStart = Bed_start;
 m_globalStart = (Bed_start >= AROUND_REGION_BUFFER ? Bed_start - AROUND_REGION_BUFFER : 0);

 m_officialEnd = Bed_end;
 m_globalEnd = std::min(chromosome->getBiolSize(), Bed_end + AROUND_REGION_BUFFER);

 m_currentStart = m_globalStart;
 m_displayedStart = m_officialStart;
 updateEndPositions();
 }*/
void LoopingSearchWindow::updateEndPositionsbyStep(int step) {
	m_currentEnd = m_currentStart + m_BIN_SIZE * step;
	if (m_currentEnd > m_globalEnd) {
		m_currentEnd = m_globalEnd;
	}
	m_displayedEnd = m_displayedStart + m_BIN_SIZE * step;
	if (m_displayedEnd > m_officialEnd) {
		m_displayedEnd = m_officialEnd;
	}
}
void LoopingSearchWindow::updateEndPositions() {
	m_currentEnd = m_currentStart + m_BIN_SIZE;
	if (m_currentEnd > m_globalEnd) {
		m_currentEnd = m_globalEnd;
	}
	m_displayedEnd = m_displayedStart + m_BIN_SIZE;
	if (m_displayedEnd > m_officialEnd) {
		m_displayedEnd = m_officialEnd;
	}
}

void LoopingSearchWindow::next() {
	m_currentStart = m_currentStart + m_BIN_SIZE;
	m_displayedStart = m_displayedStart + m_BIN_SIZE;
	updateEndPositions();
}

std::string LoopingSearchWindow::display() const {
	std::stringstream ss;
	if (m_displayedStart < m_displayedEnd) {
		ss << "\nLooking at chromosome " << m_chromosome->getName() << " bases " << m_displayedStart << " to " << m_displayedEnd << " of the bed region: chromosome " << m_chromosome->getName() << ":" << chr_start << "-" << chr_end << " \n";
	} else {
		ss << "Checking out reads near the borders of the specified regions for extra evidence.\n";
	}
	return ss.str();
}

bool LoopingSearchWindow::finished() const {
	// ugly hack for speed purposes when using Pindel-formatted input
	if (userSettings->pindelFilesAsInput() && m_currentStart >= g_maxPos) {
		return true;
	}
	return (m_currentStart > m_globalEnd);
}

unsigned int LoopingSearchWindow::get_chr_start() {
	return chr_start;
}
unsigned int LoopingSearchWindow::get_chr_end() {
	return chr_end;
}

unsigned int SPLIT_READ::getLastAbsLocCloseEnd() const {
	return UP_Close[UP_Close.size() - 1].AbsLoc;
}

bool SPLIT_READ::goodFarEndFound() const {
	return ((UP_Far.MaxLen() + UP_Close.MaxLen() >= UnmatchedSeq.size()));
}

bool SPLIT_READ::hasCloseEnd() const {
	return !UP_Close.empty();
}

unsigned int SortedUniquePoints::MaxLen() const {
	if (m_positions.size() == 0) {
		return 0;
	} else {
		int lastElementIndex = m_positions.size() - 1;
		return m_positions[lastElementIndex].LengthStr;
	}
}

unsigned int SortedUniquePoints::NumMismatch() const {
	if (m_positions.size() == 0) {
		return 1000;
	} else {
		int lastElementIndex = m_positions.size() - 1;
		return m_positions[lastElementIndex].Mismatches;
	}
}

unsigned int SPLIT_READ::MaxLenFarEnd() const {
	return UP_Far.MaxLen();
}

unsigned int SPLIT_READ::MaxLenCloseEnd() const {
	return UP_Close.MaxLen();
}

bool doOutputBreakdancerEvents() {
	return (userSettings->breakdancerOutputFilename != "" && parameters[findParameter("-b", parameters)]->isSet());
}

void outputBreakDancerEvent(const std::string &chromosomeName, const int leftPosition, const int rightPosition, const int svSize, const std::string &svType, const int svCounter) {
	std::ofstream breakDancerOutputFile(userSettings->breakdancerOutputFilename.c_str(), std::ios::app);
	breakDancerOutputFile << chromosomeName << "\t" << leftPosition << "\t" << rightPosition << "\t" << svSize << "\t" << svType << "\t" << svCounter << "\n";
}

void reportBreakDancerEvent(const std::string &chromosomeName, const int leftPosition, const int rightPosition, const int svSize, const std::string &svType, const int svCounter) {
	if (doOutputBreakdancerEvents() && g_bdData.isBreakDancerEvent(leftPosition, rightPosition)) {
		outputBreakDancerEvent(chromosomeName, leftPosition, rightPosition, svSize, svType, svCounter);
	}
}

struct Region {
	Region() {
		start = 0;
		end = 0;
	}
	unsigned start;
	unsigned end;
};

short WhetherRemoveDuplicates;

std::string TempLine_DB_Unique;

std::vector<Region>
Merge(const std::vector<Region> &AllRegions);

bool readTransgressesBinBoundaries(SPLIT_READ &read, const unsigned int &upperBinBorder) {
	return (read.BPRight > upperBinBorder - 2 * read.InsertSize);
}

/** 'readInSpecifiedRegion' if a region is specified, check if the read is in it. */
bool readInSpecifiedRegion(const SPLIT_READ &read, // in: the read
		const SearchRegion *region) {
	bool passesFilter = true;

	// if a start position has been defined, and the breakpoint is before it
	if (region->isStartDefined() && (read.BPLeft + 1 < (unsigned int) region->getStart())) {
		passesFilter = false;
	}

	// if an end of the region has been specified
	if (region->isEndDefined() && (read.BPLeft + 1 > (unsigned int) region->getEnd())) {
		passesFilter = false;

	}
	// no region specified, so all positions are okay
	return passesFilter;
}

void saveReadForNextCycle(SPLIT_READ &read, std::vector<SPLIT_READ> &futureReads) {
	futureReads.push_back(read);
	read.Used = true; // as it cannot be used for this round of analyses anymore
}
void expandCurrentWindow(LoopingSearchWindow currentWindow) {

}
std::ostream& operator<<(std::ostream &os, const UniquePoint &up) {
	os << "LengthStr: " << up.LengthStr << " * Absloc: " << up.AbsLoc << " * Direction: " << up.Direction << " * Strand: " << up.Strand;
	os << " * Mismatches: " << up.Mismatches << std::endl;
	return os;
}

std::ostream& operator<<(std::ostream &os, const SPLIT_READ &splitRead) {
	os << "\nName: " << splitRead.Name << std::endl;
	os << "Fragname: " << splitRead.FragName << std::endl;
	os << "FarFragname: " << splitRead.FarFragName << std::endl;

	os << "UnmatchedSeq: " << splitRead.getUnmatchedSeq() << std::endl;
	os << "MatchedD: " << splitRead.MatchedD << " * MatchedRelPos: " << splitRead.MatchedRelPos << " * MS: " << splitRead.MS << " * ";
	os << "InsertSize: " << splitRead.InsertSize << std::endl;
	os << "Tag: " << splitRead.Tag << std::endl;
	os << "ReadLength: " << splitRead.getReadLength() << std::endl;
	os << "ReadLengthMinus: " << splitRead.getReadLengthMinus() << std::endl;
	os << "MAX_SNP_ERROR:" << splitRead.getMAX_SNP_ERROR() << std::endl;
	os << "TOTAL_SNP_ERROR_CHECKED:" << splitRead.getTOTAL_SNP_ERROR_CHECKED() << std::endl;
	os << "getTOTAL_SNP_ERROR_CHECKED_Minus():" << splitRead.getTOTAL_SNP_ERROR_CHECKED_Minus() << std::endl;
	os << "BP:" << splitRead.BP << std::endl;
	os << "Left:" << splitRead.Left << std::endl;
	os << "Right:" << splitRead.Right << std::endl;
	os << "BPLeft:" << splitRead.BPLeft << std::endl;
	os << "BPRight:" << splitRead.BPRight << std::endl;
	os << "IndelSize:" << splitRead.IndelSize << std::endl;
	os << "UniqueRead:" << splitRead.UniqueRead << std::endl;
	os << "NT_str:" << splitRead.NT_str << std::endl;
	os << "NT_size:" << splitRead.NT_size << std::endl;
	os << "UP_Close: ";
	for (unsigned int i = 0; i < splitRead.UP_Close.size(); i++) {
		os << "[" << i << "]=" << splitRead.UP_Close[i] << " ";
	}
	os << std::endl;
	os << "UP_Far: ";
	for (unsigned int i = 0; i < splitRead.UP_Far.size(); i++) {
		os << "[" << i << "]=" << splitRead.UP_Far[i] << " ";
	}
	os << std::endl;
	//os << std::endl;

	return os;
}

bool fileExists(const std::string &filename) {
	std::ifstream test(filename.c_str());
	bool exists = (bool) test;
	test.close();
	return exists;
}

/** 'convertToXBaiFileName' converts the name of the BAM file "bamFileName", for example x.bam, into x.bai
 (or at least returns "x.bai" as a result" **/
//convertToXBaiFilename
std::string convertToXBaiFilename(const std::string &bamFileName) {
	std::string outputName = bamFileName;
	int nameLength = outputName.length();
	outputName[nameLength - 1] = 'i'; // "x.bam" -> "x.bai"
	return outputName;
}

void readBamConfigFile(std::string &bamConfigFilename, ControlState &currentState) {
	int sampleCounter = 0;
	std::ifstream BamConfigFile(bamConfigFilename.c_str());
	if (BamConfigFile) {
		while (BamConfigFile.good()) {
			bam_info tempBamInfo;
			BamConfigFile >> tempBamInfo.bamFileName >> tempBamInfo.insertSize;
			if (!BamConfigFile.good()) {
				break;
			}
			tempBamInfo.tag = "";
			BamConfigFile >> tempBamInfo.tag;
			if (tempBamInfo.tag == "") {
				*logStream << "Missing tag in line '" << tempBamInfo.bamFileName << "\t" << tempBamInfo.insertSize << "' in configuration file " << bamConfigFilename << "\n";
				exit(EXIT_FAILURE);
			}
			g_sampleNames.insert(tempBamInfo.tag);
			if (!fileExists(tempBamInfo.bamFileName)) {
				*logStream << "I cannot find the file '" << tempBamInfo.bamFileName << "'. referred to in configuration file '" << bamConfigFilename << "'. Please change the BAM configuration file.\n\n";
				exit(EXIT_FAILURE);
			}
			if (!fileExists(tempBamInfo.bamFileName + ".bai") && !fileExists(convertToXBaiFilename(tempBamInfo.bamFileName))) {
				*logStream << "I cannot find the bam index-file '" << tempBamInfo.bamFileName << ".bai' that should accompany the file " << tempBamInfo.bamFileName << " mentioned in the configuration file " << bamConfigFilename << ". Please run samtools index on " << tempBamInfo.bamFileName << ".\n\n";
				exit(EXIT_FAILURE);
			}

			//copy kai and throw crap into useless variable
			std::string restOfLine;
			safeGetline(BamConfigFile, restOfLine);
			currentState.bams_to_parse.push_back(tempBamInfo);
			sampleCounter++;
		} // while
		if (sampleCounter == 0) {
			*logStream << "Could not find any samples in the sample file '" << bamConfigFilename << "'. Please run Pindel again with a config-file of the specified type (format 'A.bam	<insert-size>	sample_label)\n\n";
			exit( EXIT_FAILURE);
		}
	} else {
		// no config-file defined
		*logStream << "BAM configuration file '" << bamConfigFilename << "' does not exist. Please run Pindel again with an existing config-file (format 'A.bam	insert-size	sample_label')\n\n";
		exit( EXIT_FAILURE);
	}
}

void readPindelConfigFile(std::string &pindelConfigFilename, std::vector<std::string> &pindelfilesToParse) {
	int sampleCounter = 0;
	std::ifstream pindelConfigFile(pindelConfigFilename.c_str());
	if (pindelConfigFile) {
		while (pindelConfigFile.good()) {
			std::string pindelFilename;
			pindelConfigFile >> pindelFilename;
			if (!pindelConfigFile.good()) {
				break;
			}

			if (!fileExists(pindelFilename)) {
				*logStream << "I cannot find the file '" << pindelFilename << "'. referred to in configuration file '" << pindelConfigFilename << "'. Please change the Pindel-configurationfile.\n\n";
				exit(EXIT_FAILURE);
			}

			std::string restOfLine;
			safeGetline(pindelConfigFile, restOfLine);
			pindelfilesToParse.push_back(pindelFilename);
			sampleCounter++;
		} // while
		if (sampleCounter == 0) {
			*logStream << "Could not find any samples in the sample file '" << pindelConfigFilename << "'. Please run Pindel again with a config-file of the specified type (format 'coloA.pindel_input<newline>colo_tumor.pindel_input<newline>...')\n\n";
			exit(EXIT_FAILURE);
		}
	} else {
		// no config-file defined
		*logStream << "Pindel configuration file '" << pindelConfigFilename << "' does not exist. Please run Pindel again with an existing config-file (format 'coloA.pindel_input<newline>colo_tumor.pindel_input<newline>...')\n\n";
		exit(EXIT_FAILURE);
	}
}

LineReader* getLineReaderByFilename(const char *filename) {
	const std::string strFilename = filename;
	const size_t len = strFilename.length();
	const std::string suffix = strFilename.substr(len - 3, 3);

	if (len > 3 && suffix == ".gz") {
		return new GZLineReader(filename);
	} else {
		return new IfstreamLineReader(filename);
	}
}

void TestFileForOutput(const char *filename) {
	std::ofstream outputfileTest(filename);
	if (!outputfileTest) {
		LOG_ERROR(*logStream << "Sorry, cannot write to the file: " << filename << std::endl);
		exit(EXIT_FAILURE);
	}
	outputfileTest.close();
}

void TestFileForOutput(const std::string &filename) {
	TestFileForOutput(filename.c_str());
}

void CheckWhetherFasta(const std::string &filename) {
	std::ifstream FastaFile(filename.c_str());
	char FirstCharOfFasta;
	FastaFile >> FirstCharOfFasta;
	if (FirstCharOfFasta != '>') {
		*logStream << "The reference genome must be in fasta format!" << std::endl;
		exit(EXIT_FAILURE);
	}
	FastaFile.close();
}

/** 'probOfReadWithTheseErrors' gives the probability that a read of length 'length' will have number of errors 'numberOfErrors' if the sequencing error rate is 'seqErrorRate'. */
double probOfReadWithTheseErrors(const unsigned int length, const unsigned int numberOfErrors, const double seqErrorRate) {
	double chanceOfCorrect = 1.0 - seqErrorRate;
	unsigned int numberOfCorrectBases = length - numberOfErrors;
	double matchedPart = pow(chanceOfCorrect, numberOfCorrectBases);

	double mismatchedPart = 1.0;
	for (unsigned int i = 0; i < numberOfErrors; i++) {
		mismatchedPart *= (((length - i) * seqErrorRate) / (numberOfErrors - i));
	}
	return matchedPart * mismatchedPart;
}

std::vector<unsigned int> g_maxMismatch;

/**  'createProbTable' fills the g_probMismatch[ length ] table with the maximum amount of differences between read and reference are acceptable.
 For example: if sequencing error rate is 1%, and SNP rate is 0.1%, total error rate is 1.1%. Throwing away all 100-base reads with 1 error would throw away over 1.1% of reads,
 which is acceptable if the sensitivity is set to 95% (0.95), but not if it is set to 1%. Uses binomial formula. */
void createProbTable(const double seqErrorRate, const double sensitivity) {
	const unsigned int MAX_LENGTH = 500;
	g_maxMismatch.assign(MAX_LENGTH, 0);
	for (unsigned int length = 0; length < MAX_LENGTH; length++) {
		double totalErrorProb = 0.0;
		for (unsigned int numberOfErrors = 0; numberOfErrors <= length; numberOfErrors++) {
			totalErrorProb += probOfReadWithTheseErrors(length, numberOfErrors, seqErrorRate);
			if (totalErrorProb > sensitivity) {
				g_maxMismatch[length] = numberOfErrors + 1;
				//g_maxMisMatch.push_back( numberOfErrors );
				//std::cout << length << " bases has max errors \t" << g_maxMismatch[length] << "\n";
				break;// break out of this length, up to the next
			}
		}
	}
	g_maxMismatch[0] = 0;
	g_maxMismatch[1] = 0;
	g_maxMismatch[2] = 0;
	g_maxMismatch[3] = 0;
}

void hello() {
	*logStream << "                                                                               " << std::endl;
	*logStream << "       ___                   ___         _____         ___                     " << std::endl;
	*logStream << "      /  /\\     ___         /__/\\       /  /::\\       /  /\\                    " << std::endl;
	*logStream << "     /  /::\\   /  /\\        \\  \\:\\     /  /:/\\:\\     /  /:/_                   " << std::endl;
	*logStream << "    /  /:/\\:\\ /  /:/         \\  \\:\\   /  /:/  \\:\\   /  /:/ /\\   ___     ___    " << std::endl;
	*logStream << "   /  /:/~/://__/::\\     _____\\__\\:\\ /__/:/ \\__\\:| /  /:/ /:/_ /__/\\   /  /\\   " << std::endl;
	*logStream << "  /__/:/ /:/ \\__\\/\\:\\__ /__/::::::::\\\\  \\:\\ /  /://__/:/ /:/ /\\\\  \\:\\ /  /:/   " << std::endl;
	*logStream << "  \\  \\:\\/:/     \\  \\:\\/\\\\  \\:\\~~\\~~\\/ \\  \\:\\  /:/ \\  \\:\\/:/ /:/ \\  \\:\\  /:/    " << std::endl;
	*logStream << "   \\  \\::/       \\__\\::/ \\  \\:\\  ~~~   \\  \\:\\/:/   \\  \\::/ /:/   \\  \\:\\/:/     " << std::endl;
	*logStream << "    \\  \\:\\       /__/:/   \\  \\:\\        \\  \\::/     \\  \\:\\/:/     \\  \\::/      " << std::endl;
	*logStream << "     \\  \\:\\      \\__\\/     \\  \\:\\        \\__\\/       \\  \\::/       \\__\\/       " << std::endl;
	*logStream << "      \\__\\/                 \\__\\/                     \\__\\/                    " << std::endl;
	*logStream << "                                                                               " << std::endl;
	*logStream << std::endl << std::endl << std::endl;

	*logStream << "           #******************************************************#" << std::endl;
	*logStream << "           #                     PINDEL 0.3                       #" << std::endl;
	*logStream << "           #                                                      #" << std::endl;
	*logStream << "           #         This is a parallel version of Pindel         #" << std::endl;
	*logStream << "           #                                                      #" << std::endl;
	*logStream << "           #                Author: Yaning Yang                   #" << std::endl;
	*logStream << "           #                  (C) 2021-, CSEE                     #" << std::endl;
	*logStream << "           #                  HUNAN UNIVERSITY                    #" << std::endl;
	*logStream << "           #******************************************************#" << std::endl;
	*logStream << std::endl << std::endl;
}

void init(int argc, char *argv[], ControlState &currentState, int mpirank) {
	logStream = &std::cout;
	if (mpirank == 0) {
		hello();
		std::cout << "Initializing parameters..." << std::endl;
	}
	defineParameters(parameters);
	readParameters(argc, argv, parameters);
	if (userSettings->logFilename != "") {
		g_logFile.open(userSettings->logFilename.c_str());
		logStream = &g_logFile;
	}
	if (argc <= 1) { // the user has not given any parameters
		printHelp(parameters);
		exit(EXIT_FAILURE);
	}
	if (!checkParameters(parameters)) {
		exit(EXIT_FAILURE);
	}
	for (unsigned i = 0; i < g_SpacerBeforeAfter; i++) {
		Spacer += "N";
	}
	createProbTable(0.001 + userSettings->Seq_Error_Rate, userSettings->sensitivity);
	std::string fastaFilename(userSettings->referenceFilename.c_str());
	if (mpirank == 0)
		std::cout << "Loading reference genome ..." << std::endl;
	time_start = clock();
	g_genome.loadAll(fastaFilename);
	time_end = clock();
	double endtime = (double) (time_end - time_start) / CLOCKS_PER_SEC;
	if (mpirank == 0) {
		std::cout << "Loading reference genome done." << std::endl;
		std::cout << "Time of loading reference genome:" << endtime << " (s)" << std::endl;
	}

	bool BreakDancerDefined = parameters[findParameter("-b", parameters)]->isSet();
	if (BreakDancerDefined) {
		g_bdData.loadBDFile(userSettings->breakdancerFilename);
	}
	if (userSettings->FLOAT_WINDOW_SIZE > 5000.0) {
		LOG_ERROR(*logStream << "Window size of " << userSettings->FLOAT_WINDOW_SIZE << " million bases is too large" << std::endl);
		exit(EXIT_FAILURE);
	} else if (userSettings->FLOAT_WINDOW_SIZE > 100.0) {
		LOG_ERROR(*logStream << "Window size of " << userSettings->FLOAT_WINDOW_SIZE << " million bases is rather large; this may produce bad::allocs or segmentation faults. If that happens, either try to reduce the window size or deactivate the searching for breakpoints and long insertions by adding the command-line options \"-l false -k false\"." << std::endl);
	}
	WINDOW_SIZE = (unsigned int) (1000000 * userSettings->FLOAT_WINDOW_SIZE);
	if (userSettings->FLOAT_CROSS_WINDOW_SIZE > 100.0) {
		LOG_ERROR(*logStream << "Cross window size of " << userSettings->FLOAT_CROSS_WINDOW_SIZE << " kb is too large" << std::endl);
		exit(EXIT_FAILURE);
	} else if (userSettings->FLOAT_CROSS_WINDOW_SIZE > 50.0) {
		LOG_ERROR(*logStream << "Cross window size of " << userSettings->FLOAT_CROSS_WINDOW_SIZE << " kb is rather large; "
				"This may take longer to process." << std::endl);
	}
	CROSS_WINDOW_SIZE = (unsigned int) (1000 * userSettings->FLOAT_CROSS_WINDOW_SIZE);
	//userSettings->Analyze_LI = true;
	if (userSettings->singlePindelFileAsInput()) {
		//      std::cout << "userSettings->singlePindelFileAsInput : " << userSettings->singlePindelFileAsInput() << std::endl;
		currentState.lineReader = getLineReaderByFilename(userSettings->pindelFilename.c_str());
		currentState.inf_Pindel_Reads = new PindelReadReader(*currentState.lineReader);
	}

	if (userSettings->pindelConfigFileAsInput()) {
		readPindelConfigFile(userSettings->pindelConfigFilename, currentState.pindelfilesToParse);
	}
	if (userSettings->bamFilesAsInput()) {
		readBamConfigFile(userSettings->bamConfigFilename, currentState);
	}
	bool AssemblyInputDefined = parameters[findParameter("-z", parameters)]->isSet();
	if (AssemblyInputDefined) {
		std::ifstream inf_AssemblyInput;
		inf_AssemblyInput.open(userSettings->inf_AssemblyInputFilename.c_str());
	}

	bool GenotypingInputDefined = parameters[findParameter("-g", parameters)]->isSet();
	if (GenotypingInputDefined) {
		std::ifstream inf_GenotypingInput;
		inf_GenotypingInput.open(userSettings->inf_GenotypingInputFilename.c_str());
	}
	omp_set_num_threads(userSettings->numThreads);
	g_MinClose = userSettings->minClose;
	//std::cout << "minClose = " << g_MinClose << std::endl;
//std::cout << "14" << std::endl;
	if (userSettings->MaxRangeIndex > g_MAX_RANGE_INDEX) {
		LOG_ERROR(*logStream << "Maximal range index (-x) exceeds the allowed value (" << g_MAX_RANGE_INDEX << ") - resetting to " << g_MAX_RANGE_INDEX << ".\n");
		userSettings->MaxRangeIndex = g_MAX_RANGE_INDEX;
	}
//std::cout << "15" << std::endl;
	if (userSettings->ADDITIONAL_MISMATCH < 1) {
		LOG_ERROR(*logStream << "Number of additional mismatches (-a) is less than the allowed value (1) - resetting to 1.\n");
		userSettings->ADDITIONAL_MISMATCH = 1;
	}

	TestFileForOutput(userSettings->getSIOutputFilename(mpirank));
	TestFileForOutput(userSettings->getDOutputFilename(mpirank));
	TestFileForOutput(userSettings->getTDOutputFilename(mpirank));
	TestFileForOutput(userSettings->getINVOutputFilename(mpirank));
	TestFileForOutput(userSettings->getLIOutputFilename(mpirank));
	TestFileForOutput(userSettings->getBPOutputFilename(mpirank));
	TestFileForOutput(userSettings->getCloseEndOutputFilename(mpirank));
	CheckWhetherFasta(userSettings->referenceFilename);

//std::cout << "16" << std::endl;

	if (userSettings->breakdancerOutputFilename != "") {
		TestFileForOutput(userSettings->breakdancerOutputFilename);
	}
//std::cout << "17" << std::endl;

	Match[(short) 'A'] = 'A';
	Match[(short) 'C'] = 'C';
	Match[(short) 'G'] = 'G';
	Match[(short) 'T'] = 'T';
	Match[(short) 'N'] = 'X';
	Match[(short) '$'] = '$';
	Match2N[(short) 'A'] = 'N';
	Match2N[(short) 'C'] = 'N';
	Match2N[(short) 'G'] = 'N';
	Match2N[(short) 'T'] = 'N';
	Match2N[(short) 'N'] = 'X';
	Match2N[(short) '$'] = '$';
	Convert2RC[(short) 'A'] = 'T';
	Convert2RC[(short) 'C'] = 'G';
	Convert2RC[(short) 'G'] = 'C';
	Convert2RC[(short) 'T'] = 'A';
	Convert2RC[(short) 'N'] = 'X';
	Convert2RC[(short) '$'] = '$';
	Convert2RC4N[(short) 'A'] = 'T';
	Convert2RC4N[(short) 'C'] = 'G';
	Convert2RC4N[(short) 'G'] = 'C';
	Convert2RC4N[(short) 'T'] = 'A';
	Convert2RC4N[(short) 'N'] = 'N';
	Cap2LowArray[(short) 'A'] = 'a';
	Cap2LowArray[(short) 'C'] = 'c';
	Cap2LowArray[(short) 'G'] = 'g';
	Cap2LowArray[(short) 'T'] = 't';
	Cap2LowArray[(short) 'N'] = 'n';
	Cap2LowArray[(short) '$'] = 'n';

//std::cout << "18" << std::endl;

	/*   for (unsigned int i = 0; i < g_SpacerBeforeAfter; i++) {
	 Spacer += "N";
	 }*/
//std::cout << "19" << std::endl;
	//Distance = 300;
	DSizeArray[0] = 0;
	DSizeArray[1] = 128;
	for (int dIndex = 2; dIndex < 15; dIndex++) {
		DSizeArray[dIndex] = DSizeArray[dIndex - 1] * 4;
	}
	if (mpirank == 0)
		std::cout << "Initializing parameters done." << std::endl;
	//if (!userSettings->getRegion()->isTargetChromosomeDefined() && AssemblyInputDefined == false && GenotypingInputDefined == false) {
	//    *logStream << "Looping over all chromosomes." << std::endl;
	//}
//std::cout << "20" << std::endl;
}

void SearchFarEnd(const std::string &chromosome, SPLIT_READ &read, const Chromosome &currentChromosome) {

	const int START_SEARCH_SPAN = 64;
	//std::cout << "getting searchCluster" << std::endl;
	const std::vector<SearchWindow> &searchCluster = g_bdData.getCorrespondingSearchWindowCluster(read);
	//std::cout << "searchCluster size: " << searchCluster.size() << std::endl;
	if (searchCluster.size() != 0) {
		//std::cout << "Breakdancer input is not empty " << searchCluster.size() << std::endl;
		//for (unsigned index = 0; index < searchCluster.size(); index++)
		//	searchCluster[index].display();

		// chromosome: currentChromosome->getSeq(),  read: currentState.Reads_SR[i]
		//std::cout << "SearchFarEndAtPos(chromosome, read, searchCluster)" << std::endl;

		SearchFarEndAtPos(chromosome, read, searchCluster); // SearchFarEndAtPos
		//std::cout << "SearchFarEndAtPos(chromosome, read, searchCluster)  done....." << std::endl;
		//std::cout << "finished" << std::endl;

		if (read.goodFarEndFound()) {
			//read.Investigated = true;
			return;
		}
		//else SearchFarEndAtPosPerfect( chromosome, read, searchCluster);
	}
	//std::cout << "SearchFarEnd	2" << std::endl;
	//UserDefinedSettings* userSettings = UserDefinedSettings::Instance();

	// if breakdancer does not find the event, or not find an event we trust, we turn to regular pattern matching
	unsigned int searchSpan = START_SEARCH_SPAN;
	unsigned int centerOfSearch = read.getLastAbsLocCloseEnd();
	//std::cout << "SearchFarEnd	3" << std::endl;
	// userSettings->MaxRangeIndex + 1 = 3
	//std::cout << "read.getUnmatchedSeq().size() = " << read.getUnmatchedSeq().size() << std::endl;
	//std::cout << "centerOfSearch = " << centerOfSearch << std::endl;
	//std::cout <<"userSettings->MaxRangeIndex + 1 = " << userSettings->MaxRangeIndex + 1 << std::endl;

	for (int rangeIndex = 1; rangeIndex <= userSettings->MaxRangeIndex + 1; rangeIndex++) {
		//std::cout << "rangeIndex " << rangeIndex << "\t";
		// note: searching the central range again and again may seem overkill,
		//but since Pindel at the moment wants an unique best point, you can't skip the middle part
		// may be stuff for future changes/refactorings though
		std::vector<SearchWindow> aroundCESearchCluster;
		unsigned int Start, End;
		if (centerOfSearch > searchSpan + g_SpacerBeforeAfter) {
			Start = centerOfSearch - searchSpan;
		} else {
			Start = g_SpacerBeforeAfter;
		}
		if (centerOfSearch + searchSpan + g_SpacerBeforeAfter < chromosome.size()) {
			End = centerOfSearch + searchSpan;
		} else {
			End = chromosome.size() - g_SpacerBeforeAfter;
		}
		//std::cout << "Start is (abs) " << Start << " (rel): " << Start - g_SpacerBeforeAfter << "\n";
		//std::cout << "End is (abs) " << End << " (rel): " << End - g_SpacerBeforeAfter << " Size chrom= " << chromosome.size() << "\n";
		//std::cout <<"Start is " << Start << "\t" << "FirstStart " << centerOfSearch+searchSpan-10000000 << "<COS"
				//<< centerOfSearch-10000000 << " span " << searchSpan<< std::endl;
		SearchWindow regularWindow(&currentChromosome, Start, End);

		aroundCESearchCluster.clear();
		aroundCESearchCluster.push_back(regularWindow);
		//std::cout << rangeIndex << "\tSearchFarEndAtPos" << std::endl;
		//std::cout << "Before searching, read.UP_Far.MaxLen()  = " << read.UP_Far.MaxLen() << ", " << "read.UP_Close.MaxLen() = " << read.UP_Close.MaxLen() << std::endl;
		SearchFarEndAtPos(chromosome, read, aroundCESearchCluster); // SearchFarEndAtPosPerfect
		//std::cout << "end\tSearchFarEndAtPos" << std::endl;
		//std::cout << "After searching, read.UP_Far.MaxLen()  = " << read.UP_Far.MaxLen() << ", " << "read.UP_Close.MaxLen() = " << read.UP_Close.MaxLen() << std::endl;
		//std::cout << "UnmatchedSeq.size() = " << read.getUnmatchedSeq().size() << std::endl;
		if (read.goodFarEndFound()) {
			//read.Investigated = true;
			//std::cout << "read.goodFarEndFound()........" << std::endl;
			//std::cout << std::endl<< std::endl<< std::endl << std::endl;
			return;
		}

		//else SearchFarEndAtPosPerfect( chromosome, read, searchCluster);
		searchSpan *= 4;
	}
	//std::cout << std::endl<< std::endl<< std::endl << std::endl;
}

void ReportCloseMappedReads(const std::vector<SPLIT_READ> &reads) {
	std::ofstream CloseEndMappedOutput(userSettings->getCloseEndOutputFilename().c_str(), std::ios::app);
	int TotalNumReads = reads.size();
	for (int Index = 0; Index < TotalNumReads; Index++) {
		const SPLIT_READ &currentRead = reads[Index];
		CloseEndMappedOutput << currentRead.Name << "\n" << currentRead.getUnmatchedSeq() << "\n" << currentRead.MatchedD << "\t" << currentRead.FragName << "\t" << currentRead.MatchedRelPos << "\t" << currentRead.MS << "\t" << currentRead.InsertSize << "\t" << currentRead.Tag << "\n";
	}
}

void ReportCloseAndFarEndCounts(const std::vector<SPLIT_READ> &reads) {
	unsigned Count_Far = 0;
	unsigned Count_Used = 0;
	unsigned Count_Unused = 0;

	for (int Index = reads.size() - 1; Index >= 0; Index--) {
		const SPLIT_READ &currentRead = reads[Index];
		if (!currentRead.UP_Far.empty()) {
			Count_Far++;
		}
		if (currentRead.Used) {
			Count_Used++;
		}
	}
	Count_Unused = reads.size() - Count_Far;
	*logStream << "Total: " << reads.size() << ";\tClose_end_found " << reads.size() << ";\tFar_end_found " << Count_Far << ";\tUsed\t" << Count_Used << ".\n\n";
	*logStream << "For LI and BP: " << Count_Unused << "\n\n";
}

void SearchFarEnds(const std::string &chromosomeSeq, std::vector<SPLIT_READ> &reads, const Chromosome &currentChromosome) {
	//std::cout << "report per 1k reads, not sequential due to openmp" << std::endl;
	//std::cout << "vector<SPLIT_READ> &reads size = " << reads.size() << std::endl;
#pragma omp parallel default(shared)
	{
#pragma omp for
		for (int readIndex = 0; readIndex < (int) reads.size(); readIndex++) {
			if (reads[readIndex].MapperSplit == false) {
				SearchFarEnd(chromosomeSeq, reads[readIndex], currentChromosome);

			}
		}
	}

	*logStream << "Far end searching completed for this window." << std::endl;
}

void SearchSVs(ControlState &currentState, const int NumBoxes, const SearchWindow &currentWindow, int mpirank) {
	//std::cout << "1" << std::endl;
	//UserDefinedSettings* userSettings = UserDefinedSettings::Instance();

	SearchDeletions searchD;
	searchD.SearchByMpirank(g_bdData, currentState, NumBoxes, currentWindow, mpirank);
//std::cout << "2" << std::endl;
	searchIndels(currentState, NumBoxes, currentWindow, mpirank);
//std::cout << "3" << std::endl;
	if (userSettings->Analyze_TD) {
		searchTandemDuplications(currentState, NumBoxes, currentWindow, mpirank);
		searchTandemDuplicationsNT(currentState, NumBoxes, currentWindow, mpirank);
	}
//std::cout << "4" << std::endl;
	if (userSettings->Analyze_INV) {
		searchInversions(currentState, NumBoxes, currentWindow, mpirank);
		searchInversionsNT(currentState, NumBoxes, currentWindow, mpirank);
	}
//std::cout << "5" << std::endl;

	SearchShortInsertions searchSI;
	searchSI.SearchByMpirank(g_bdData, currentState, NumBoxes, currentWindow, mpirank);
//std::cout << "6" << std::endl;
	ReportCloseAndFarEndCounts(currentState.Reads_SR);
//std::cout << "7" << std::endl;
	if (userSettings->Analyze_LI) {
		SortOutputLI(currentState, currentWindow.getChromosome()->getSeq(), currentState.Reads_SR, currentWindow, userSettings->getLIOutputFilename(mpirank));
	}
//std::cout << "8" << std::endl;
	//if (userSettings->Analyze_BP) {
	//SortOutputRest(currentState, currentWindow.getChromosome()->getSeq(), currentState.Reads_SR, currentWindow, userSettings->getBPOutputFilename());
	//}
//std::cout << "9" << std::endl;
}

class TimerItem {

public:
	TimerItem(const std::string &id);
	void stop();
	void restart();
	const std::string getReport() const;
	const std::string& getId() const {
		return m_id;
	}

private:
	std::string m_id;
	time_t m_countSoFar;
	time_t m_lastStart;
};

TimerItem::TimerItem(const std::string &id) {
	m_id = id;
	m_countSoFar = 0;
	m_lastStart = time(NULL);
}

void TimerItem::stop() {
	m_countSoFar += (time(NULL) - m_lastStart);
}

void TimerItem::restart() {
	m_lastStart = time(NULL);
}

const std::string TimerItem::getReport() const {
	std::stringstream ss;
	ss << m_id << " " << m_countSoFar << " seconds.\n";
	return ss.str();
}

class Timer {

public:
	void switchTo(const std::string &itemName);
	void reportAll(std::ostream &os);
	Timer() :
			m_currentItemIndex(-1) {
	}

private:
	std::vector<TimerItem> m_timerItems;
	int m_currentItemIndex;
};

void Timer::switchTo(const std::string &itemName) {
	if (m_currentItemIndex != -1) {
		m_timerItems[m_currentItemIndex].stop();
	}
	for (unsigned int itemIter = 0; itemIter < m_timerItems.size(); itemIter++) {
		if (m_timerItems[itemIter].getId() == itemName) {
			m_timerItems[itemIter].restart();
			m_currentItemIndex = (int) itemIter;
			return;
		}
	}
	// no existing element found? new element needs to be constructed
	TimerItem newItem(itemName);
	m_timerItems.push_back(newItem);
	m_currentItemIndex = m_timerItems.size() - 1;
}

void Timer::reportAll(std::ostream &os) {
	if (m_currentItemIndex != -1) {
		m_timerItems[m_currentItemIndex].stop();
	}
	for (std::vector<TimerItem>::iterator itemIter = m_timerItems.begin(); itemIter != m_timerItems.end(); itemIter++) {
		os << itemIter->getReport();
	}
}

void UpdateFarFragName(std::vector<SPLIT_READ> &input) {
	for (unsigned index = 0; index < input.size(); index++) {
		if (input[index].UP_Far.size()) {
			input[index].FarFragName = input[index].UP_Far[0].chromosome_p->getName();
			input[index].MatchedFarD = input[index].UP_Far[0].Strand;
		}
	}
}

short UpdateRefReadCoverage(ControlState &currentState, const SearchWindow &currentWindow) {
	//std::set<std::string> g_sampleNames;
	//std::map<std::string,unsigned> g_SampleName2Index;
	std::cout << "There are " << g_sampleNames.size() << " samples." << std::endl;
	g_SampleName2Index.clear();
	std::set<std::string>::iterator it;
	unsigned index = 0;
	for (it = g_sampleNames.begin(); it != g_sampleNames.end(); it++) {
		g_SampleName2Index.insert(std::pair<std::string, unsigned>(*it, index));
		index++;
	}
	std::cout << "SampleName2Index done" << std::endl;
	// std::vector <RefCoveragePerPosition> g_RefCoverageRegion;
	unsigned Start = currentWindow.getStart();
	unsigned End = currentWindow.getEnd();
	unsigned Length = End - Start + 1;
	unsigned NumberOfSamples = g_sampleNames.size();
	std::cout << "declaring g_RefCoverageRegion for " << g_sampleNames.size() << " samples and " << Length << " positions." << std::endl;
	g_RefCoverageRegion.clear();
	RefCoveragePerPosition OnePos;
	for (unsigned SampleIndex = 0; SampleIndex < NumberOfSamples; SampleIndex++) {
		OnePos.RefCoveragePerSample.push_back(0);
	}
	for (unsigned index = 0; index < Length; index++) {
		g_RefCoverageRegion.push_back(OnePos);
	}
	std::cout << "declare g_RefCoverageRegion done" << std::endl;

#pragma omp parallel default(shared)
	{
#pragma omp for
		for (int readIndex = 0; readIndex < (int) currentState.RefSupportingReads.size(); readIndex++) {
			if (currentState.RefSupportingReads[readIndex].Pos < Start || currentState.RefSupportingReads[readIndex].Pos + currentState.RefSupportingReads[readIndex].ReadLength > End) {
				continue;
			}
			unsigned SampleID = g_SampleName2Index.find(currentState.RefSupportingReads[readIndex].Tag)->second;
			unsigned PosStart = currentState.RefSupportingReads[readIndex].Pos - Start;
#pragma omp critical
			{
				for (unsigned posIndex = 1; posIndex < (unsigned) currentState.RefSupportingReads[readIndex].ReadLength - 1; posIndex++) {
					g_RefCoverageRegion[PosStart + posIndex].RefCoveragePerSample[SampleID]++;
				}
			}
		}   //RefSupportingReads
	}
	//std::cout << "update g_RefCoverageRegion done" << std::endl;
	/*
	 for (unsigned PosIndex = 0; PosIndex < Length; PosIndex++) {
	 std::cout << Start + PosIndex;
	 it=g_sampleNames.begin();
	 for (unsigned SampleIndex = 0; SampleIndex < NumberOfSamples; SampleIndex++) {
	 std::cout << "\t" << *it << ":" << g_RefCoverageRegion[PosIndex].RefCoveragePerSample[SampleIndex];
	 it++;
	 }
	 std::cout << "\n";
	 }*/
	return 0;
}

short init_g_ChrNameAndSizeAndIndex(std::string RefIndexFileName) {
	ChrNameAndSizeAndIndex OneChr;
	std::string TempStr;
	std::ifstream RefIndexInput(RefIndexFileName.c_str());
	if (!RefIndexInput) {
		return 1;
	}
	short ChrCount = 0;
	while (RefIndexInput >> OneChr.ChrName >> OneChr.ChrSize) {
		getline(RefIndexInput, TempStr);
		OneChr.ChrIndex = ChrCount;
		ChrCount++;
		g_ChrNameAndSizeAndIndex.push_back(OneChr);
	}
	return 0;
}

unsigned getChrIndex(std::string &ChrName) {
	for (unsigned index = 0; index < g_ChrNameAndSizeAndIndex.size(); index++) {
		if (g_ChrNameAndSizeAndIndex[index].ChrName == ChrName) {
			return index;
		}
	}
	return g_ChrNameAndSizeAndIndex.size();
}

unsigned getChrSize(std::string &ChrName) {
	for (unsigned index = 0; index < g_ChrNameAndSizeAndIndex.size(); index++) {
		if (g_ChrNameAndSizeAndIndex[index].ChrName == ChrName) {
			return g_ChrNameAndSizeAndIndex[index].ChrSize;
		}
	}
	return 0;
}

short CheckChrName(std::string ChrName) {
	for (unsigned index = 0; index < g_ChrNameAndSizeAndIndex.size(); index++) {
		if (g_ChrNameAndSizeAndIndex[index].ChrName == ChrName) {
			return 1;
		}
	}
	return 0;
}

void CleanUpBedRecord(std::vector<BED> &include, std::vector<BED> &exclude) {
	if (exclude.size() == 0) {
		return;
	}
	//std::cout << "before include.size() " << include.size() << std::endl;
	for (unsigned include_index = 0; include_index < include.size(); include_index++) {
		//std::cout << "include " << include[include_index].ChrName << "\t" << include[include_index].Start << "\t" << include[include_index].End << std::endl;
		for (unsigned exclude_index = 0; exclude_index < exclude.size(); exclude_index++) {
			//std::cout << "exclude " << include[include_index].ChrName << "\t" << include[include_index].Start << "\t" << include[include_index].End << std::endl;
			if (include[include_index].Start == include[include_index].End) {
				break;
			}
			if (include[include_index].ChrName != exclude[exclude_index].ChrName) {
				continue;
			} else { // same chromosome
				if (include[include_index].Start > exclude[exclude_index].End || exclude[exclude_index].Start > include[include_index].End) {
					continue;
				}
				// "exclude" contains current "include", set start = end
				if (exclude[exclude_index].Start <= include[include_index].Start && include[include_index].End <= exclude[exclude_index].End) { // "exclude" contains current "include", set start = end
					include[include_index].End = include[include_index].Start;
				}
				// "include" contains current "exclude", break into two regions
				else if (include[include_index].Start < exclude[exclude_index].Start && exclude[exclude_index].End < include[include_index].End) {

					BED NewOne;
					NewOne.ChrName = include[include_index].ChrName;
					NewOne.Start = exclude[exclude_index].End;
					NewOne.End = include[include_index].End;
					include.push_back(NewOne);
					include[include_index].End = exclude[exclude_index].Start;
					//std::cout << include[include_index].ChrName << " " << include[include_index].Start << " " << include[include_index].End << "\t" << NewOne.ChrName << " " << NewOne.Start << " " << NewOne.End << std::endl;
				}
				// intersect I
				else if (exclude[exclude_index].Start <= include[include_index].Start && include[include_index].Start < exclude[exclude_index].End && exclude[exclude_index].End < include[include_index].End) {
					include[include_index].Start = exclude[exclude_index].End;
				}
				// intersect II
				else if (include[include_index].Start < exclude[exclude_index].Start && exclude[exclude_index].Start < include[include_index].End && include[include_index].End < exclude[exclude_index].End) {
					include[include_index].End = exclude[exclude_index].Start;
				}
			}
		}
	}
	//std::cout << "after include.size() " << include.size() << std::endl;
	//for (unsigned index = 0; index < include.size(); index++) {
	//	std::cout << include[index].ChrName << " " << include[index].Start << " " << include[index].End << std::endl;
	//}
	std::vector<BED> result;
	//std::cout << "a list of regions after excluding regions: " << std::endl;
	for (unsigned include_index = 0; include_index < include.size(); include_index++) {
		if (include[include_index].Start != include[include_index].End) {
			result.push_back(include[include_index]);
			//std::cout << include[include_index].ChrName << "\t" << include[include_index].Start << "\t" << include[include_index].End << std::endl;
		}
	}
	//std::cout << "merging the list if there is overlap: " << std::endl;
	for (unsigned first = 0; first < result.size() - 1; first++) {
		for (unsigned second = first + 1; second < result.size(); second++) {
			if (result[first].ChrName != result[second].ChrName) {
				continue;
			} else { // same chromosome

				// non overlap
				if (result[first].Start > result[second].End || result[second].Start > result[first].End) {
					continue;
				}
				// "second" contains current "first", set start = end
				if (result[second].Start <= result[first].Start && result[first].End <= result[second].End) { // "second" contains current "first", set start = end
					result[first].End = result[first].Start;
					break;
				}
				// "first" contains current "second", discard second
				else if (result[first].Start <= result[second].Start && result[second].End <= result[first].End) {
					result[second].Start = result[second].End;
					break;
				}
				// intersect I
				else if (result[second].Start <= result[first].Start && result[first].Start <= result[second].End && result[second].End <= result[first].End) {
					result[first].Start = result[second].Start;
					result[second].Start = result[second].End;
				}
				// intersect II
				else if (result[first].Start <= result[second].Start && result[second].Start <= result[first].End && result[first].End <= result[second].End) {
					result[first].End = result[second].End;
					result[second].Start = result[second].End;
				}
			}
		}
	}
	std::vector<BED> final;
	for (unsigned index = 0; index < result.size(); index++)
		if (result[index].Start != result[index].End) {
			final.push_back(result[index]);
		}
	// sorting
	//std::cout << "final regions" << std::endl;
	bool exchange;
	for (unsigned first = 0; first < final.size() - 1; first++) {
		for (unsigned second = first + 1; second < final.size(); second++) {
			exchange = false;
			unsigned first_ChrIndex = getChrIndex(final[first].ChrName);
			unsigned second_ChrIndex = getChrIndex(final[second].ChrName);
			if (first_ChrIndex < second_ChrIndex) {
				continue;
			} else if (first_ChrIndex > second_ChrIndex) {
				exchange = true;
			} else { // same chr
				if (final[first].Start > final[second].Start) {
					exchange = true;
				}
			}
			if (exchange) {
				BED temp = final[first];
				final[first] = final[second];
				final[second] = temp;
			}
		}
		unsigned CurrentSize = getChrSize(final[first].ChrName);
		if (CurrentSize < final[first].End) {
			final[first].End = CurrentSize;
		}
		//std::cout << final[first].ChrName << "\t" << final[first].Start << "\t" << final[first].End << std::endl;
	}
	exclude.clear();
	std::cout << "\nfinal list of regions:" << std::endl;
	for (unsigned first = 0; first < final.size(); first++) {
		std::cout << "\t" << final[first].ChrName << "\t" << final[first].Start << "\t" << final[first].End << std::endl;
	}
	std::cout << std::endl;
	include = final;
}

struct InterChrCall {
	char AnchorD;
	char FirstD;
	std::string FirstChrName;
	unsigned FirstPos;
	char SecondD;
	std::string SecondChrName;
	unsigned SecondPos;
	unsigned NumSupport;
	std::string InsertedSequence;
};

void MergeInterChr(ControlState &currentState, UserDefinedSettings *usersettings) {
	unsigned cutoff = 2;
	std::ifstream INT_input(usersettings->getINTOutputFilename().c_str());
	InterChrCall one;
	std::vector<InterChrCall> All;
	std::string tempstr;
	while (INT_input >> tempstr >> one.AnchorD >> one.FirstChrName >> one.FirstPos >> one.FirstD >> one.SecondChrName >> one.SecondPos >> one.SecondD >> one.InsertedSequence >> tempstr >> one.NumSupport) {
		//if (one.FirstPos != 0 && one.SecondPos != 0)
		All.push_back(one);
		//std::cout << "getting " << one.FirstChrName << "\t" << one.FirstPos << "\t" << one.SecondChrName << "\t" << one.SecondPos << "\t" << one.NumSupport << std::endl;
	}
	std::ofstream INToutputfile((usersettings->getINTOutputFilename() + "_final").c_str());
	if (All.size() == 0) {
		return;
	} else if (All.size() < 2) {
		if (All[0].NumSupport >= cutoff * 2) {
			INToutputfile << All[0].FirstChrName << "\t" << All[0].FirstPos << "\t" << All[0].SecondChrName << "\t" << All[0].SecondPos << "\t" << All[0].InsertedSequence << "\t" << All[0].NumSupport << "\t" << All[0].AnchorD << "\t" << All[0].FirstChrName << "\t" << All[0].FirstPos << "\t" << All[0].FirstD << "\t" << All[0].SecondChrName << "\t" << All[0].SecondPos << "\t" << All[0].SecondD << "\t" << All[0].InsertedSequence << "\t" << All[0].NumSupport << std::endl;
		}
	}
	bool reported;
	for (unsigned index_a = 0; index_a < All.size(); index_a++) {
		reported = false;
		for (unsigned index_b = index_a; index_b < All.size(); index_b++) {
			if (index_a == index_b) {
				continue;
			}
			if (All[index_a].FirstChrName == All[index_b].FirstChrName && All[index_a].SecondChrName == All[index_b].SecondChrName) {
				if (cabs(All[index_a].FirstPos - All[index_b].FirstPos) < 10 && cabs(All[index_a].SecondPos - All[index_b].SecondPos) < 10 && All[index_a].NumSupport + All[index_b].NumSupport >= cutoff) {

					INToutputfile << "chr\t" << All[index_a].FirstChrName << "\tpos\t" << unsigned((All[index_a].FirstPos + All[index_b].FirstPos) / 2) << "\tchr\t" << All[index_a].SecondChrName << "\tpos\t" << unsigned((All[index_a].SecondPos + All[index_b].SecondPos) / 2) << "\tseq\t" << All[index_a].InsertedSequence << "\tsupport\t" << All[index_a].NumSupport + All[index_b].NumSupport << "\tINFOR\t" << All[index_a].AnchorD << "\t" << All[index_a].FirstChrName << "\t" << All[index_a].FirstPos
							<< "\t" << All[index_a].FirstD << "\t" << All[index_a].SecondChrName << "\t" << All[index_a].SecondPos << "\t" << All[index_a].SecondD << "\t" << All[index_a].InsertedSequence << "\t" << All[index_a].NumSupport << "\t" << All[index_b].AnchorD << "\t" << All[index_b].FirstChrName << "\t" << All[index_b].FirstPos << "\t" << All[index_b].FirstD << "\t" << All[index_b].SecondChrName << "\t" << All[index_b].SecondPos << "\t" << All[index_b].SecondD << "\t"
							<< All[index_b].InsertedSequence << "\t" << All[index_b].NumSupport << std::endl;
					reported = true;
					break;
				}
			}
		}
		if (reported == false) {
			if (All[index_a].NumSupport >= cutoff * 2) {
				INToutputfile << "chr\t" << All[index_a].FirstChrName << "\tpos\t" << All[index_a].FirstPos << "\tchr\t" << All[index_a].SecondChrName << "\tpos\t" << All[index_a].SecondPos << "\tseq\t" << All[index_a].InsertedSequence << "\tsupport\t" << All[index_a].NumSupport << "\tINFOR\t" << All[index_a].AnchorD << "\t" << All[index_a].FirstChrName << "\t" << All[index_a].FirstPos << "\t" << All[index_a].FirstD << "\t" << All[index_a].SecondChrName << "\t" << All[index_a].SecondPos << "\t"
						<< All[index_a].SecondD << "\t" << All[index_a].InsertedSequence << "\t" << All[index_a].NumSupport << std::endl;
			}
		}
	}
}

void readUserSetting(ControlState &currentState, int mpirank) {
	if (init_g_ChrNameAndSizeAndIndex(userSettings->getRefFilename() + ".fai") == 1) {
		if (mpirank == 0)
			std::cout << "Please use samtools to index your reference file.\n .fai is missing.\n" << std::endl;
		exit(1);
	}
	if (userSettings->loopOverAllChromosomes()) { // WGS
		if (parameters[findParameter("-j", parameters)]->isSet()) { // if a bed file is provided to processing
			if (userSettings->inf_InclusiveBedFileName != "") {
				BED OneBedRecord;
				std::ifstream inf_IncludeBed(userSettings->inf_InclusiveBedFileName.c_str());
				std::string restofbedline;
				while (inf_IncludeBed >> OneBedRecord.ChrName >> OneBedRecord.Start >> OneBedRecord.End) {
					getline(inf_IncludeBed, restofbedline);
					if (OneBedRecord.Start > OneBedRecord.End) {
						unsigned tempint = OneBedRecord.Start;
						OneBedRecord.Start = OneBedRecord.End;
						OneBedRecord.End = tempint;
					}
					unsigned ChrSize = g_genome.getChr(OneBedRecord.ChrName)->getBiolSize();
					if (OneBedRecord.End > ChrSize) {
						OneBedRecord.End = ChrSize;
					}
					currentState.IncludeBed.push_back(OneBedRecord);
					//std::cout << currentState.IncludeBed.size() << "\t" << OneBedRecord.ChrName << "\t" << OneBedRecord.Start << "\t" << OneBedRecord.End << std::endl;
				}
			}
		} else {
			// g_ChrNameAndSizeAndIndex.size() == 2779
			for (unsigned index = 0; index < g_ChrNameAndSizeAndIndex.size(); index++) {
				BED OneBedRecord;
				OneBedRecord.ChrName = g_ChrNameAndSizeAndIndex[index].ChrName;
				OneBedRecord.Start = 1;
				OneBedRecord.End = g_ChrNameAndSizeAndIndex[index].ChrSize;
				currentState.IncludeBed.push_back(OneBedRecord);
				//std::cout << currentState.IncludeBed.size() << "\t" << OneBedRecord.ChrName << "\t" << OneBedRecord.Start << "\t" << OneBedRecord.End << std::endl;
			}
		}
	} else { // one region
		if (parameters[findParameter("-j", parameters)]->isSet()) { // if a bed file is provided to processing
			if (userSettings->inf_InclusiveBedFileName != "") {

				std::string ChrName = userSettings->getRegion()->getTargetChromosomeName();
				unsigned Start = userSettings->getRegion()->getStart();
				unsigned End = userSettings->getRegion()->getEnd();
				unsigned ChrSize = g_genome.getChr(ChrName)->getBiolSize();
				if (End > ChrSize) {
					End = ChrSize;
				}
				BED OneBedRecord;
				std::ifstream inf_IncludeBed(userSettings->inf_InclusiveBedFileName.c_str());
				std::string restofbedline;
				while (inf_IncludeBed >> OneBedRecord.ChrName >> OneBedRecord.Start >> OneBedRecord.End) {
					getline(inf_IncludeBed, restofbedline);
					if (OneBedRecord.Start > OneBedRecord.End) {
						unsigned tempint = OneBedRecord.Start;
						OneBedRecord.Start = OneBedRecord.End;
						OneBedRecord.End = tempint;
					}
					if (OneBedRecord.ChrName != ChrName) {
						continue;   // not in the same chromosome
					}
					if (OneBedRecord.Start > End) {
						continue;   // no overlap
					}
					if (Start > OneBedRecord.End) {
						continue;   // no overlap
					}
					if (OneBedRecord.Start < Start) {
						OneBedRecord.Start = Start;
					}
					if (OneBedRecord.End > End) {
						OneBedRecord.End = End;
					}
					currentState.IncludeBed.push_back(OneBedRecord);
					//std::cout << currentState.IncludeBed.size() << "\t" << OneBedRecord.ChrName << "\t" << OneBedRecord.Start << "\t" << OneBedRecord.End << std::endl;
				}
			}
		} else {
			BED OneBedRecord;
			OneBedRecord.ChrName = userSettings->getRegion()->getTargetChromosomeName();
			OneBedRecord.Start = userSettings->getRegion()->getStart();
			OneBedRecord.End = userSettings->getRegion()->getEnd();

			unsigned ChrSize = g_genome.getChr(OneBedRecord.ChrName)->getBiolSize();
			if (OneBedRecord.End > ChrSize) {
				OneBedRecord.End = ChrSize;
			}
			currentState.IncludeBed.push_back(OneBedRecord);
		}
	}

	//std::cout << "currentState.IncludeBed.size() " << currentState.IncludeBed.size() << std::endl;
	//for (unsigned index = 0; index < currentState.IncludeBed.size(); index++) {
	//	std::cout << currentState.IncludeBed.size() << "\t" << currentState.IncludeBed[index].ChrName << "\t" << currentState.IncludeBed[index].Start << "\t" << currentState.IncludeBed[index].End << std::endl;
	//}

	//std::cout << "inf_ExclusiveBedFileName" << userSettings->inf_ExclusiveBedFileName << std::endl;
	if (parameters[findParameter("-J", parameters)]->isSet()) { // if a bed file is provided to exclude regions for processing
		if (userSettings->inf_ExclusiveBedFileName != "") {
			BED OneBedRecord;
			std::ifstream inf_ExcludeBed(userSettings->inf_ExclusiveBedFileName.c_str());
			std::string restofbedline;
			while (inf_ExcludeBed >> OneBedRecord.ChrName >> OneBedRecord.Start >> OneBedRecord.End) {
				getline(inf_ExcludeBed, restofbedline);
				if (OneBedRecord.Start > OneBedRecord.End) {
					unsigned tempint = OneBedRecord.Start;
					OneBedRecord.Start = OneBedRecord.End;
					OneBedRecord.End = tempint;
				}
				currentState.ExcludeBed.push_back(OneBedRecord);
				//std::cout << currentState.ExcludeBed.size() << "\t" << OneBedRecord.ChrName << "\t" << OneBedRecord.Start << "\t" << OneBedRecord.End << std::endl;
			}
		}
	}

	CleanUpBedRecord(currentState.IncludeBed, currentState.ExcludeBed);

	//return 0;

	if (!userSettings->loopOverAllChromosomes()) {
		if (CheckChrName(userSettings->getRegion()->getTargetChromosomeName()) == 0) {
			if (mpirank == 0)
				std::cout << "Please check chromosome name in the reference file, BAM files and command line. \n Make sure that they are consistent.\n" << std::endl;
			exit(1);
		}
	}
	/* Start of shortcut to genotyping */ //
	std::ifstream inf_AssemblyInput;
	inf_AssemblyInput.open(userSettings->inf_AssemblyInputFilename.c_str());
	bool GenotypingInputDefined = parameters[findParameter("-g", parameters)]->isSet();

	//std::ifstream FastaFile(userSettings->getRefFilename().c_str());
	if (GenotypingInputDefined) {
		//doGenotyping(currentState, userSettings );
		exit(EXIT_SUCCESS);
	}

	bool AssemblyInputDefined = parameters[findParameter("-z", parameters)]->isSet();
	if (AssemblyInputDefined) {
		//doAssembly(currentState, userSettings );
		exit(EXIT_SUCCESS);
	}

	// If -q parameter given, search for mobile element insertions and quit.
	if (parameters[findParameter("-q", parameters)]->isSet()) {
		exit(searchMEImain(currentState, g_genome, userSettings));
	}

	if (parameters[findParameter("-Y", parameters)]->isSet()) {
		std::string ChrName, TempStr;
		unsigned Ploidy;
		std::ifstream PloidyFile(userSettings->PloidyFileName.c_str());
		std::map<std::string, unsigned>::iterator it;
		while (PloidyFile >> ChrName >> Ploidy) {
			std::cout << ChrName << "\t" << Ploidy << std::endl;
			getline(PloidyFile, TempStr);
			it = g_ChrName2Ploidy.find(ChrName);
			if (it != g_ChrName2Ploidy.end()) {
				std::cout << "This chromosome " << ChrName << " already seen. Please check " << userSettings->PloidyFileName << std::endl;
				exit(1);
			} else {
				g_ChrName2Ploidy.insert(std::pair<std::string, unsigned>(ChrName, Ploidy));
			}
		}
	}

	std::ofstream RPoutputfile(userSettings->getRPOutputFilename().c_str());
}
const Chromosome* getBedInfo(std::string &Bed_ChrName, int &Bed_start, int &Bed_end, int index, ControlState &currentState) {
	Bed_ChrName = currentState.IncludeBed[index].ChrName;
	Bed_start = currentState.IncludeBed[index].Start;
	Bed_end = currentState.IncludeBed[index].End;
	const Chromosome *currentChromosome = g_genome.getChr(Bed_ChrName);
	return currentChromosome;
}
void clearChrInfo(ControlState &currentState) {
	currentState.InterChromosome_SR.clear();
	currentState.InputReads_SR.clear();
	currentState.Reads_SR.clear();
	currentState.FutureReads_SR.clear();
	currentState.OneEndMappedReads.clear();
	currentState.Reads_RP.clear();
	currentState.Reads_RP_Discovery.clear();
	currentState.RefSupportingReads.clear();
	currentState.Reads_RP_Discovery_InterChr.clear();
}

void mpiVecWinSend(std::vector<WindowInfo> &winSet, int &vWin_size, int dest, int tag, MPI_Comm comm) {
	MPI_Send(&vWin_size, 1, MPI_INT, dest, tag, comm);
	for (int v = 0; v < vWin_size; v++) {
		WindowInfo tempWin = winSet[v];
		int len_chrName = tempWin.chrName.length();
		MPI_Send(&len_chrName, 1, MPI_INT, dest, tag, comm);
		if (len_chrName != 0) {
			MPI_Send(tempWin.chrName.data(), len_chrName, MPI_CHAR, dest, tag, comm);
		}
		MPI_Send(&tempWin.bed_index, 1, MPI_INT, dest, tag, comm);
		MPI_Send(&tempWin.win_start, 1, MPI_INT, dest, tag, comm);
		MPI_Send(&tempWin.win_end, 1, MPI_INT, dest, tag, comm);
		MPI_Send(&tempWin.global_start, 1, MPI_INT, dest, tag, comm);
		MPI_Send(&tempWin.global_end, 1, MPI_INT, dest, tag, comm);
		MPI_Send(&tempWin.isused, 1, MPI_INT, dest, tag, comm);

	}
}

void mpiVecWinRecv(std::vector<WindowInfo> &winSet, int &vWin_size, int src, int tag, MPI_Comm comm, MPI_Status status) {
	MPI_Recv(&vWin_size, 1, MPI_INT, src, tag, comm, &status);
	for (int v = 0; v < vWin_size; v++) {
		WindowInfo tempWin;
		int len_chrName;
		int tempBed_index;
		int tempStart;
		int tempEnd;
		int tempGlobalStart;
		int tempGlobalEnd;
		int tempIsUsed;
		MPI_Recv(&len_chrName, 1, MPI_INT, src, tag, comm, &status);

		std::string tempChrName;
		if (len_chrName != 0) {
			std::vector<char> tmp(len_chrName);
			MPI_Recv(tmp.data(), len_chrName, MPI_CHAR, src, tag, comm, &status);
			tempChrName.assign(tmp.begin(), tmp.end());
		} else {
			tempChrName.clear();
		}
		tempWin.chrName = tempChrName;
		MPI_Recv(&tempBed_index, 1, MPI_INT, src, tag, comm, &status);
		tempWin.bed_index = tempBed_index;
		MPI_Recv(&tempStart, 1, MPI_INT, src, tag, comm, &status);
		tempWin.win_start = tempStart;
		MPI_Recv(&tempEnd, 1, MPI_INT, src, tag, comm, &status);
		tempWin.win_end = tempEnd;
		MPI_Recv(&tempGlobalStart, 1, MPI_INT, src, tag, comm, &status);
		tempWin.global_start = tempGlobalStart;
		MPI_Recv(&tempGlobalEnd, 1, MPI_INT, src, tag, comm, &status);
		tempWin.global_end = tempGlobalEnd;
		MPI_Recv(&tempIsUsed, 1, MPI_INT, src, tag, comm, &status);
		tempWin.isused = tempIsUsed;
		winSet.push_back(tempWin);
	}
}

void searchOneWindow(ControlState &currentState, int win, int mpirank) {
	time_t reading_start, reading_end, searchingFar_start, searchingFar_end, searchingSV_start, searchingSV_end;
	WindowInfo currWindow = windowSet.at(win);
	std::string CurrentChrName = currWindow.chrName;
	int current_window_bedIndex = currWindow.bed_index;
	int current_window_start = currWindow.win_start;
	int current_window_end = currWindow.win_end;
	int current_global_start = currWindow.global_start;
	int current_global_end = currWindow.global_end;

	int Bed_start, Bed_end;
	const Chromosome *currentChromosome = getBedInfo(CurrentChrName, Bed_start, Bed_end, current_window_bedIndex, currentState);
	if (currentChromosome == NULL) {
		std::cout << "There is no " << CurrentChrName << " in the reference file." << std::endl;
		exit(1);
	}
	*logStream << "Process: " << mpirank << ", processing region: " << CurrentChrName << "\t" << current_window_start << "\t" << current_window_end << std::endl;
	g_maxPos = 0;
	CurrentChrMask.resize(currentChromosome->getCompSize());
	for (unsigned int i = 0; i < currentChromosome->getCompSize(); i++) {
		CurrentChrMask[i] = 'N';
	}
	BoxSize = currentChromosome->getCompSize() / 30000;
	if (BoxSize == 0) {
		BoxSize = 1;
	}
	unsigned NumBoxes = (unsigned) (currentChromosome->getCompSize() * 2 / BoxSize) + 1; // box size
	userSettings->getRegion()->SetRegion(CurrentChrName, current_window_start, current_window_end);
	LoopingSearchWindow currentWindow(userSettings->getRegion(), currentChromosome, WINDOW_SIZE, current_global_start, current_global_end);

	g_RegionStart = currentWindow.getStart();
	g_RegionEnd = currentWindow.getEnd();
	SearchWindow currentWindow_cs = currentWindow.makePindelCoordinateCopy();

	time(&reading_start);
	if (userSettings->bamFilesAsInput() && userSettings->SearchDiscordantReadPair) {
		get_RP_Cross_Reads_Discovery(currentState, currentWindow, current_global_start, current_global_end);
		g_bdData.UpdateBD(currentState);
		std::cout << "external BD events: " << g_bdData.GetBDSize_external() << " Added BreakDancer-like events: " << (g_bdData.GetBDSize_total() - g_bdData.GetBDSize_external()) / 2 << std::endl;

	}
	g_bdData.loadRegion(currentWindow_cs);

	get_SR_Cross_Reads(currentState, currentWindow, current_global_start, current_global_end);

	std::cout << "There are " << currentState.RefSupportingReads.size() << " reads supporting the reference allele." << std::endl;
	UpdateRefReadCoverage(currentState, currentWindow);
	if (currentState.Reads_SR.size()) {
		*logStream << "There are " << currentState.Reads_SR.size() << " split-reads for this chromosome region.\n" << std::endl; // what region?
		std::cout << "There are " << g_NumberOfGapAlignedReads << " split-reads mapped by aligner." << std::endl;
		g_NumberOfGapAlignedReads = 0;
		if (userSettings->reportCloseMappedReads()) {
			*logStream << "report closeMappedReads" << std::endl;
			ReportCloseMappedReads(currentState.Reads_SR);
		}
		//Time_Load_E = time(NULL);
		time(&reading_end);
		readingAndMatching += difftime(reading_end, reading_start);

		if (!userSettings->reportOnlyCloseMappedReads) {
			//timer.switchTo("Searching far ends");

			time(&searchingFar_start);
			*logStream << "search far ends" << std::endl;
			SearchFarEnds(currentChromosome->getSeq(), currentState.Reads_SR, *currentChromosome);
			*logStream << "update FarFragName" << std::endl;

			UpdateFarFragName(currentState.Reads_SR);
			*logStream << "update FarFragName done" << std::endl;

			if (userSettings->reportInterchromosomalEvents) {
				*logStream << "save interchromsome SR" << std::endl;
				for (unsigned index = 0; index < currentState.Reads_SR.size(); index++) {
					if (currentState.Reads_SR[index].UP_Far.size()) {
						//std::cout << currentState.Reads_SR[index].FragName << " " << currentState.Reads_SR[index].FarFragName << std::endl;
						if (currentState.Reads_SR[index].FragName != currentState.Reads_SR[index].FarFragName) {
							currentState.InterChromosome_SR.push_back(currentState.Reads_SR[index]);
						}
					}
				}
			}
			time(&searchingFar_end);
			searchingFarEnds += difftime(searchingFar_end, searchingFar_start);
			//timer.switchTo("Searching and reporting variations");

			time(&searchingSV_start);
			*logStream << "Searching and reporting variations" << std::endl;
			SearchSVs(currentState, NumBoxes, currentWindow, mpirank);
			time(&searchingSV_end);
			searchingSVs += difftime(searchingSV_end, searchingSV_start);
		}
		*logStream << "There are " << currentState.FutureReads_SR.size() << " split-reads saved for the next cycle.\n" << std::endl;
	} else {
		*logStream << "no reads " << std::endl;
	}
	//std::cout << "InterChromosome_SR.size(): " << currentState.InterChromosome_SR.size() << std::endl;
	if (userSettings->reportInterchromosomalEvents && currentState.InterChromosome_SR.size()) {
		SortAndReportInterChromosomalEvents(currentState, g_genome, userSettings);
	}
	{
		currentState.InterChromosome_SR.clear();
		currentState.InputReads_SR.clear();
		currentState.Reads_SR.clear();
		currentState.FutureReads_SR.clear();
		currentState.OneEndMappedReads.clear();
		currentState.Reads_RP.clear();
		currentState.Reads_RP_Discovery.clear();
		currentState.RefSupportingReads.clear();
		currentState.Reads_RP_Discovery_InterChr.clear();
	}

	if (parameters[findParameter("-q", parameters)]->isSet()) {
		int DDresult = searchMEImain(currentState, g_genome, userSettings);
		if (DDresult != 0) {
			exit(DDresult);
		}
	}
	MergeInterChr(currentState, userSettings);
}

int main(int argc, char *argv[]) {

	time_t total_start, total_end;
	int mpirc, mpisize, mpirank;
	mpirc = MPI_Init(&argc, (char***) &argv);
	if (mpirc != MPI_SUCCESS) {
		fprintf(stderr, "MPI cannot be initialized. Terminating.\n");
		MPI_Abort(MPI_COMM_WORLD, mpirc);
	}
	mpirc = MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
	mpirc = MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
	MPI_Status mpiStatus, mpiStatus4task, mpiStatus4finish;
	int winTag = 1, vWin_size = 0;
	userSettings = UserDefinedSettings::Instance();

	time(&total_start);
	ControlState currentState;
	init(argc, argv, currentState, mpirank);
	readUserSetting(currentState, mpirank);

	if (mpirank == 0) {
		*logStream << "The number of chromosome is " << currentState.IncludeBed.size() << std::endl;
		for (size_t index = 0; index < currentState.IncludeBed.size(); index++) {
			std::string currChrName = currentState.IncludeBed[index].ChrName;
			unsigned int win = currentState.IncludeBed[index].Start;
			while (win < currentState.IncludeBed[index].End) {
				WindowInfo tempWin;
				tempWin.chrName = currChrName;
				tempWin.global_start = currentState.IncludeBed[index].Start;
				tempWin.global_end = currentState.IncludeBed[index].End;
				tempWin.bed_index = index;

				tempWin.win_start = win > CROSS_WINDOW_SIZE ? (win - CROSS_WINDOW_SIZE) : win;
				//tempWin.win_start = win - 2 * insertSize < 0 ? win : win-2 * insertSize;
				win += WINDOW_SIZE;
				tempWin.win_end = win < tempWin.global_end ? win : tempWin.global_end;
				tempWin.isused = 0;
				windowSet.push_back(tempWin);
			}
		}
		//std::cout << "WindowSet.size = " << windowSet.size() <<std::endl;;

		*logStream << "The number of processes is " << mpisize << std::endl;
		*logStream << "WINDOW SIZE = " << WINDOW_SIZE << ", the number of window is " << windowSet.size() << std::endl;
		vWin_size = windowSet.size();
		for (int dest = 1; dest < mpisize; dest++) {
			mpiVecWinSend(windowSet, vWin_size, dest, winTag, MPI_COMM_WORLD);
		}
	} else {
		mpiVecWinRecv(windowSet, vWin_size, 0, winTag, MPI_COMM_WORLD, mpiStatus);
	}

	bool firstTurn = true;
	int tag4newtask = 444, tag4finishedtask = 555, tag4finish = 6, isFinshed = 0;

	unsigned int nextWindowId = 0, finishedWin = 0, receiveCount = 0;
	int max_num_processes = windowSet.size() + 1;
	if (mpisize > max_num_processes && mpirank == 0) {
		std::cout << "In order to save computing resources, please use less than or equal to " << windowSet.size() + 1 << " processes" << std::endl;
		exit(1);
	} else if (mpisize == 1) {
		for (unsigned int winId = 0; winId < windowSet.size(); winId++) {
			searchOneWindow(currentState, winId, mpirank);
		}
	} else {
		while (true) {
			if (mpirank == 0) {
				for (int dest = 1; dest < mpisize; dest++) {
					MPI_Send(&isFinshed, 1, MPI_INT, dest, tag4finish, MPI_COMM_WORLD);
				}
				if (isFinshed == 1) {
					std::cout << "All windows is finished." << std::endl;
					break;
				}
				while (isFinshed == 0) {
					if (firstTurn == true) {
						for (int i = 1; i < mpisize; i++) {
							int winId = i - 1;
							MPI_Send(&winId, 1, MPI_INT, i, tag4newtask, MPI_COMM_WORLD);
							std::cout << "The window " << windowSet.at(winId).chrName << ":" << windowSet.at(winId).win_start << "-" << windowSet.at(winId).win_end << " is send to process " << i << std::endl;
							windowSet.at(winId).isused = 1;
							firstTurn = false;
						}
					} else {
						MPI_Recv(&finishedWin, 1, MPI_INT, MPI_ANY_SOURCE, tag4finishedtask, MPI_COMM_WORLD, &mpiStatus4task);
						receiveCount++;
						//std::cout << "receiveCount = " << receiveCount << std::endl;
						int finishedRank = mpiStatus4task.MPI_SOURCE;
						if (finishedWin < windowSet.size()) {
							windowSet.at(finishedWin).isused = 1;
							std::cout << "The window " << windowSet.at(finishedWin).chrName << ":" << windowSet.at(finishedWin).win_start << "-" << windowSet.at(finishedWin).win_end << " is finished by process " << finishedRank << std::endl;
							std::cout << std::endl;
						}
						for (nextWindowId = 0; nextWindowId < windowSet.size(); ++nextWindowId) {
							if (windowSet.at(nextWindowId).isused == 0) {
								break;
							}
						}
						if (receiveCount == windowSet.size()) {
							isFinshed = 1;
							break;
						}
						if (nextWindowId < windowSet.size()) {
							if (windowSet.at(nextWindowId).isused == 0) {
								MPI_Send(&nextWindowId, 1, MPI_INT, finishedRank, tag4newtask, MPI_COMM_WORLD);
								//MPI_Isend(&nextWindowId, 1, MPI_INT, finishedRank, tag4newtask, MPI_COMM_WORLD, &mpiRequest);
								windowSet.at(nextWindowId).isused = 1;
								std::cout << "The next window " << windowSet.at(nextWindowId).chrName << ":" << windowSet.at(nextWindowId).win_start << "-" << windowSet.at(nextWindowId).win_end << " is sent to " << finishedRank << std::endl;
							}
						}
					}
				}

			} else {
				while (isFinshed == 0) {
					if (isFinshed == 0) {
						MPI_Recv(&nextWindowId, 1, MPI_INT, 0, tag4newtask, MPI_COMM_WORLD, &mpiStatus4task);
						std::cout << "Process " << mpirank << " received next window " << windowSet.at(nextWindowId).chrName << ":" << windowSet.at(nextWindowId).win_start << "-" << windowSet.at(nextWindowId).win_end << std::endl;
						searchOneWindow(currentState, nextWindowId, mpirank);
						finishedWin = nextWindowId;
						MPI_Send(&finishedWin, 1, MPI_INT, 0, tag4finishedtask, MPI_COMM_WORLD);
						windowSet.at(finishedWin).isused = 1;
					}
				}
				MPI_Recv(&isFinshed, 1, MPI_INT, 0, tag4finish, MPI_COMM_WORLD, &mpiStatus4finish);
				//MPI_Irecv(&slaveFinished, 1, MPI_INT, 0, tag4finish, MPI_COMM_WORLD, &mpiRequest);
				if (isFinshed == 1) {
					std::cout << "slave " << mpirank << " is finished." << std::endl;
					MPI_Finalize();
					//break;
				}
			}

		}
	}

	time(&total_end);
	total_time = difftime(total_end, total_start);

//		std::cout << "Loading chromosomes: " << loadChrs << " (s)" << std::endl;
//		std::cout << "Reading in reads + matching close ends: " << readingAndMatching << " (s)" << std::endl;
//		std::cout << "Searching far ends: " << searchingFarEnds << " (s)" << std::endl;
//		std::cout << "Searching and reporting variations: " << searchingSVs << " (s)" << std::endl;
	std::cout << "The total running time: " << total_time << " (s)" << std::endl;

	//MPI_Finalize();
	exit(EXIT_SUCCESS);

} //main

std::vector<std::string> ReverseComplement(const std::vector<std::string> &InputPatterns) {
	std::vector<std::string> OutputPatterns; // = InputPatterns;
	unsigned int NumPattern = InputPatterns.size();
	OutputPatterns.resize(NumPattern);
	for (unsigned int i = 0; i < NumPattern; i++) {
		OutputPatterns[i] = ReverseComplement(InputPatterns[i]);
	}
	return OutputPatterns;
}

std::string Reverse(const std::string &InputPattern) {
	std::string OutputPattern = InputPattern;
	unsigned int LenPattern = InputPattern.size();
	for (unsigned int j = 0; j < LenPattern; j++) {
		OutputPattern[j] = InputPattern[j];
	}
	return OutputPattern;
}

std::string ReverseComplement(const std::string &InputPattern) {
	std::string OutputPattern = InputPattern;

	unsigned int LenPattern = InputPattern.size();

	for (unsigned int j = 0; j < LenPattern; j++) {
		OutputPattern[j] = Convert2RC4N[(unsigned int) InputPattern[LenPattern - j - 1]];
	}

	return OutputPattern;
}

std::string Cap2Low(const std::string &input) {
	std::string output = input;
	for (unsigned int i = 0; i < output.size(); i++) {
		output[i] = Cap2LowArray[(short) input[i]];
	}
	return output;
}

bool ReportEvent(const std::vector<SPLIT_READ> &Deletions, const unsigned int &S, const unsigned int &E) {
//return true;
	short ReadLength = Deletions[S].getReadLength() - Deletions[S].NT_size;
	short Min_Length = (short) ((ReadLength * Min_Filter_Ratio) + 0.5) - 1;
	short Max_Length = (short) (ReadLength * (1 - Min_Filter_Ratio) - 0.5) - 1;
	bool LeftMin = false;
	bool LeftMax = false;
	bool RightMin = false;
	bool RightMax = false;
	for (unsigned i = S; i <= E; i++) {
		ReadLength = Deletions[i].getReadLength() - Deletions[i].NT_size;
		Min_Length = (short) ((ReadLength * Min_Filter_Ratio) + 0.5) - 1;
		Max_Length = (short) (ReadLength * (1 - Min_Filter_Ratio) - 0.5) - 1;
		if (Deletions[i].BP <= Min_Length) {
			LeftMin = true;
		}
		if (Deletions[i].getReadLength() - Deletions[i].BP - Deletions[i].NT_size <= Min_Length) {
			RightMin = true;
		}
		if (Deletions[i].BP >= Max_Length) {
			LeftMax = true;
		}
		if (Deletions[i].getReadLength() - Deletions[i].BP - Deletions[i].NT_size >= Max_Length) {
			RightMax = true;
		}
	}

	if (LeftMin && LeftMax && RightMin && RightMax) {
		return true;
	} else {
		return false;
	}
}

void GetRealStart4Deletion(const std::string &TheInput, unsigned int &RealStart, unsigned int &RealEnd) {

	if (TheInput.size() < RealStart || TheInput.size() < RealEnd) {
		return;
	}
	unsigned int PosIndex = RealStart + g_SpacerBeforeAfter;
	unsigned int Start = PosIndex + 1;
	unsigned int End = RealEnd + g_SpacerBeforeAfter - 1;
	while (TheInput[PosIndex] == TheInput[End] && TheInput[PosIndex] != 'N') {
		--PosIndex;
		--End;
	}
	RealStart = PosIndex - g_SpacerBeforeAfter;
	PosIndex = RealEnd + g_SpacerBeforeAfter;
	while (TheInput[PosIndex] == TheInput[Start] && TheInput[PosIndex] != 'N') {
		++PosIndex;
		++Start;
	}
	RealEnd = PosIndex - g_SpacerBeforeAfter;
}

// rotates a string back: KAI -> IKA
void rotateBack(std::string &str) {
	char lastChar = str[str.size() - 1];
	str = lastChar + str.substr(0, str.size() - 1);
}

void rotateForward(std::string &str) {
	char firstChar = str[0];
	str = str.substr(1) + firstChar;
}

// the true principle here would be that you can move the insertion backwards,
void GetRealStart4Insertion(const std::string &chromosomeSeq, std::string &InsertedStr, unsigned int &RealStart, unsigned int &RealEnd) {

	if (chromosomeSeq.size() < RealStart || chromosomeSeq.size() < RealEnd) {
		return;
	}
	unsigned int lastPosAfterInsertion_comp = RealEnd + g_SpacerBeforeAfter;
//std::cout << "AInsertedStr: " << InsertedStr << ", start= " << RealStart  << chromosomeSeq[ g_SpacerBeforeAfter + RealStart] << ", end= " << RealEnd << chromosomeSeq[ g_SpacerBeforeAfter + RealEnd] << std::endl;
//for (int x=-5; x<5; x++ ) { std::cout << chromosomeSeq[ g_SpacerBeforeAfter + RealStart + x] ; }
//std::cout << "\n";
	while (chromosomeSeq[lastPosAfterInsertion_comp] == InsertedStr[0] && chromosomeSeq[lastPosAfterInsertion_comp] != 'N') {
		rotateForward(InsertedStr);
		lastPosAfterInsertion_comp++;
	}
	RealEnd = lastPosAfterInsertion_comp - g_SpacerBeforeAfter;
//std::cout << "BInsertedStr: " << InsertedStr << ", start= " << RealStart << ", end= " << RealEnd << std::endl;
//for (int x=-5; x<5; x++ ) { std::cout << chromosomeSeq[ g_SpacerBeforeAfter + RealStart + x] ; }
//std::cout << "\n";
	unsigned int lastPosBeforeInsertion_comp = lastPosAfterInsertion_comp - 1;
	while (chromosomeSeq[lastPosBeforeInsertion_comp] == InsertedStr[InsertedStr.size() - 1] && chromosomeSeq[lastPosBeforeInsertion_comp] != 'N') {
		rotateBack(InsertedStr);
		lastPosBeforeInsertion_comp--;
	}
	RealStart = lastPosBeforeInsertion_comp - g_SpacerBeforeAfter;
//std::cout << "CInsertedStr: " << InsertedStr << ", start= " << RealStart << ", end= " << RealEnd << std::endl;
//for (int x=-5; x<5; x++ ) { std::cout << chromosomeSeq[ g_SpacerBeforeAfter + RealStart + x] ; }
//std::cout << "\n";

}
//std::cout << "RS: " << RealStart << " RE: " << RealEnd << " IS: " <<  InsertedStr << "\n";
//std::cout << "----\n";

/*unsigned int IndelSize = InsertedStr.size();
 unsigned int PosIndex = RealStart + g_SpacerBeforeAfter;
 unsigned int original_RealStart = RealStart;

 for (int i = IndelSize - 1; i >= 0; i--) {
 if (TheInput[PosIndex] == InsertedStr[i]) {
 PosIndex--;
 }
 else {
 break;
 }
 }
 if (PosIndex == RealStart + g_SpacerBeforeAfter - IndelSize) {
 while (TheInput[PosIndex] == TheInput[PosIndex + IndelSize]) {
 PosIndex--;
 }
 }
 RealStart = PosIndex - g_SpacerBeforeAfter;
 PosIndex = RealEnd + g_SpacerBeforeAfter;
 for (unsigned int i = 0; i < IndelSize; i++) {
 if (TheInput[PosIndex] == InsertedStr[i]) {
 PosIndex++;
 }
 else {
 break;
 }
 }
 if (PosIndex == RealEnd + g_SpacerBeforeAfter + IndelSize) {
 while (TheInput[PosIndex] == TheInput[PosIndex - IndelSize]) {
 PosIndex++;
 }
 }
 RealEnd = PosIndex - g_SpacerBeforeAfter;
 unsigned DIFF = RealStart - original_RealStart;
 InsertedStr = InsertedStr.substr(0, IndelSize - DIFF) + InsertedStr.substr(IndelSize, DIFF);*/
//}
/*void GetRealStart4Insertion(const std::string & TheInput,
 std::string & InsertedStr, unsigned int &RealStart,
 unsigned int &RealEnd)
 {
 unsigned int IndelSize = InsertedStr.size();
 unsigned int PosIndex = RealStart + g_SpacerBeforeAfter;
 unsigned int original_RealStart = RealStart;

 for (int i = IndelSize - 1; i >= 0; i--) {
 if (TheInput[PosIndex] == InsertedStr[i]) {
 PosIndex--;
 }
 else {
 break;
 }
 }
 if (PosIndex == RealStart + g_SpacerBeforeAfter - IndelSize) {
 while (TheInput[PosIndex] == TheInput[PosIndex + IndelSize]) {
 PosIndex--;
 }
 }
 RealStart = PosIndex - g_SpacerBeforeAfter;
 PosIndex = RealEnd + g_SpacerBeforeAfter;
 for (unsigned int i = 0; i < IndelSize; i++) {
 if (TheInput[PosIndex] == InsertedStr[i]) {
 PosIndex++;
 }
 else {
 break;
 }
 }
 if (PosIndex == RealEnd + g_SpacerBeforeAfter + IndelSize) {
 while (TheInput[PosIndex] == TheInput[PosIndex - IndelSize]) {
 PosIndex++;
 }
 }
 RealEnd = PosIndex - g_SpacerBeforeAfter;
 unsigned DIFF = RealStart - original_RealStart;
 InsertedStr = InsertedStr.substr(0, IndelSize - DIFF) + InsertedStr.substr(IndelSize, DIFF);
 }*/

std::vector<Region> Merge(const std::vector<Region> &AllRegions) {
	return AllRegions;
}

void GetCloseEndInner(const std::string &CurrentChrSeq, SPLIT_READ &Temp_One_Read, int RangeIndex) {
	std::string CurrentReadSeq;
//std::vector<unsigned int> PD[Temp_One_Read.getTOTAL_SNP_ERROR_CHECKED()];
	std::vector<PosVector> PD;       //typedef std::vector <unsigned> PosVector;
	PosVector emptyPosVector;
	PD.assign(Temp_One_Read.getTOTAL_SNP_ERROR_CHECKED(), emptyPosVector);
	if (Temp_One_Read.InsertSize > g_maxInsertSize) {
		g_maxInsertSize = Temp_One_Read.InsertSize;
	}
	short total_snp_error = Temp_One_Read.getTOTAL_SNP_ERROR_CHECKED();   // 6
	for (int CheckIndex = 0; CheckIndex < total_snp_error; CheckIndex++) {
		PD[CheckIndex].reserve(4 * Temp_One_Read.InsertSize);
	}
	SortedUniquePoints UP;
	int Start, End;
	short BP_Start; // = MinClose;
	short BP_End; // = ReadLength - MinClose;

	Temp_One_Read.UP_Close.clear();
	BP_Start = g_MinClose;
	BP_End = Temp_One_Read.getReadLengthMinus();

//std::cout << "CurrentChrSeq = " << CurrentChrSeq<< std::endl;
	if (Temp_One_Read.MatchedD == Plus) {
		CurrentReadSeq = ReverseComplement(Temp_One_Read.getUnmatchedSeq());

		Start = Temp_One_Read.MatchedRelPos + g_SpacerBeforeAfter - RangeIndex * Temp_One_Read.InsertSize; //	RangeIndex = 0, 1
		End = Start + (2 * RangeIndex + 1) * Temp_One_Read.InsertSize;
//std::cout << "1+" << Start << " " << End << std::endl;
//Start = Temp_One_Read.MatchedRelPos + g_SpacerBeforeAfter - 2 * RangeIndex * Temp_One_Read.getReadLength(); /////////////
//End = Start + 2 * (RangeIndex + 1) * Temp_One_Read.getReadLength();
//std::cout << "CurrentReadSeq: " << CurrentReadSeq <<std::endl;
		char LeftChar;
		LeftChar = CurrentReadSeq[0];
		if (LeftChar != 'N') {
			{
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChrSeq[pos] == LeftChar) {
						PD[0].push_back(pos);
					}
				}
			}
		}

		CheckLeft_Close(Temp_One_Read, CurrentChrSeq, CurrentReadSeq, PD, BP_Start, BP_End, 1, UP); // LengthStr
		if (UP.empty()) {
		} else {
			Temp_One_Read.Used = false;
			Temp_One_Read.UP_Close.swap(UP);
			UP.clear();
		}
	} else if (Temp_One_Read.MatchedD == Minus) {

		CurrentReadSeq = Temp_One_Read.getUnmatchedSeq();
		End = Temp_One_Read.MatchedRelPos + g_SpacerBeforeAfter + (RangeIndex) * Temp_One_Read.InsertSize; /////////
		Start = End - (2 * RangeIndex + 1) * Temp_One_Read.InsertSize;
//std::cout << "1-" << Start << " " << End << std::endl;
//End = Temp_One_Read.MatchedRelPos + g_SpacerBeforeAfter + 2 * (RangeIndex) * Temp_One_Read.getReadLength(); /////////
//Start = End - 2 * (RangeIndex + 1) * Temp_One_Read.getReadLength();
//std::cout << "CurrentReadSeq: " << CurrentReadSeq <<std::endl;
		char RightChar;
		RightChar = CurrentReadSeq[Temp_One_Read.getReadLengthMinus()];
//std::cout << "Starting to fit the close end with character" << RightChar << "\n";
		if (RightChar != 'N') {
			for (int pos = Start; pos < End; pos++) {
				if (CurrentChrSeq[pos] == RightChar) {
					PD[0].push_back(pos);
				}
			}
		}
//std::cout << "1\t" << PD[0].size() << "\t" << PD[1].size() << std::endl;
//        LOG_DEBUG(*logStream << "1\t" << PD[0].size() << "\t" << PD[1].size() << std::endl);
		CheckRight_Close(Temp_One_Read, CurrentChrSeq, CurrentReadSeq, PD, BP_Start, BP_End, 1, UP);
//        LOG_DEBUG(*logStream << UP.size() << std::endl);
		if (UP.empty()) {
		} else {
			Temp_One_Read.Used = false;
			Temp_One_Read.UP_Close.swap(UP);
			UP.clear();
		}
	}
	/*
	 if (Temp_One_Read.UP_Close.size() == 0) {
	 if (Temp_One_Read.MatchedD == Plus) {
	 CurrentReadSeq = ReverseComplement(Temp_One_Read.getUnmatchedSeq());
	 Start = Temp_One_Read.MatchedRelPos + g_SpacerBeforeAfter - (RangeIndex - 1) * Temp_One_Read.InsertSize; /////////////
	 End = Start + (RangeIndex + 1) * Temp_One_Read.InsertSize;
	 //std::cout << "2+" << Start << " " << End << std::endl;
	 //Start = Temp_One_Read.MatchedRelPos + g_SpacerBeforeAfter - 2 * (RangeIndex - 1) * Temp_One_Read.getReadLength(); /////////////
	 //End = Start + 2 * (RangeIndex + 1) * Temp_One_Read.getReadLength();
	 if (Start + 10000 < End) {
	 std::cout << "warning: in GetCloseEndInner Start + 10000 < End, slow" << std::endl;
	 }
	 char LeftChar;
	 LeftChar = CurrentReadSeq[0];
	 if (LeftChar != 'N') {
	 {
	 for (int pos = Start; pos < End; pos++) {
	 if (CurrentChrSeq[pos] == LeftChar) {
	 PD[0].push_back(pos);
	 }
	 }
	 }
	 }

	 CheckLeft_Close(Temp_One_Read, CurrentChrSeq, CurrentReadSeq, PD, BP_Start, BP_End, 1, UP); // LengthStr
	 if (UP.empty()) {}
	 else {
	 Temp_One_Read.Used = false;
	 Temp_One_Read.UP_Close.swap(UP);
	 UP.clear();
	 }
	 }
	 else if (Temp_One_Read.MatchedD == Minus) {

	 CurrentReadSeq = Temp_One_Read.getUnmatchedSeq();
	 End = Temp_One_Read.MatchedRelPos + g_SpacerBeforeAfter + (RangeIndex - 1) * Temp_One_Read.InsertSize; /////////
	 Start = End - (RangeIndex + 1) * Temp_One_Read.InsertSize;
	 //std::cout << "2-" << Start << " " << End << std::endl;
	 //End = Temp_One_Read.MatchedRelPos + g_SpacerBeforeAfter + 2 * (RangeIndex - 1) * Temp_One_Read.getReadLength(); /////////
	 //Start = End - 2 * (RangeIndex + 1) * Temp_One_Read.getReadLength();
	 if (Start + 10000 < End) {
	 std::cout << "warning: in GetCloseEndInner Start + 10000 < End, slow" << std::endl;
	 }
	 char RightChar;
	 RightChar = CurrentReadSeq[Temp_One_Read.getReadLengthMinus()];
	 //std::cout << "Starting to fit the close end with character" << RightChar << "\n";
	 if (RightChar != 'N') {
	 for (int pos = Start; pos < End; pos++) {
	 if (CurrentChrSeq[pos] == RightChar) {
	 PD[0].push_back(pos);
	 }
	 }
	 }
	 //std::cout << "1\t" << PD[0].size() << "\t" << PD[1].size() << std::endl;
	 //        LOG_DEBUG(*logStream << "1\t" << PD[0].size() << "\t" << PD[1].size() << std::endl);
	 CheckRight_Close(Temp_One_Read, CurrentChrSeq, CurrentReadSeq, PD, BP_Start, BP_End, 1, UP);
	 //        LOG_DEBUG(*logStream << UP.size() << std::endl);
	 if (UP.empty()) {}
	 else {
	 Temp_One_Read.Used = false;
	 Temp_One_Read.UP_Close.swap(UP);
	 UP.clear();
	 }
	 }
	 }
	 */

	return;
}

void GetCloseEndInnerPerfectMatch(const std::string &CurrentChrSeq, SPLIT_READ &Temp_One_Read, unsigned RangeIndex) {
	std::string CurrentReadSeq;
//std::vector<unsigned int> PD[Temp_One_Read.getTOTAL_SNP_ERROR_CHECKED()];
	std::vector<PosVector> PD;
	PosVector emptyPosVector;
	PD.assign(Temp_One_Read.getTOTAL_SNP_ERROR_CHECKED(), emptyPosVector);
	if (Temp_One_Read.InsertSize > g_maxInsertSize) {
		g_maxInsertSize = Temp_One_Read.InsertSize;
	}
	for (int CheckIndex = 0; CheckIndex < Temp_One_Read.getTOTAL_SNP_ERROR_CHECKED(); CheckIndex++) {
		PD[CheckIndex].reserve(3 * Temp_One_Read.InsertSize);
	}
	SortedUniquePoints UP;
	int Start, End;
	short BP_Start; // = MinClose;
	short BP_End; // = ReadLength - MinClose;

	Temp_One_Read.UP_Close.clear();
	BP_Start = g_MinClose;
	BP_End = Temp_One_Read.getReadLengthMinus();
	if (Temp_One_Read.MatchedD == Plus) {
		CurrentReadSeq = ReverseComplement(Temp_One_Read.getUnmatchedSeq());
		Start = Temp_One_Read.MatchedRelPos + g_SpacerBeforeAfter - RangeIndex * Temp_One_Read.InsertSize; /////////////
		End = Start + (RangeIndex + 1) * Temp_One_Read.InsertSize;
		char LeftChar;
		LeftChar = CurrentReadSeq[0];
		if (LeftChar != 'N') {
			{
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChrSeq[pos] == LeftChar) {
						PD[0].push_back(pos);
					}
				}
			}
		}

		if (PD[0].size()) {
			CheckLeft_Close_Perfect(Temp_One_Read, CurrentChrSeq, CurrentReadSeq, PD, BP_Start, BP_End, 1, UP); // LengthStr
		}
		if (UP.empty()) {
		} else {
			Temp_One_Read.Used = false;
			Temp_One_Read.UP_Close.swap(UP);
			UP.clear();
		}
	} else if (Temp_One_Read.MatchedD == Minus) {

		CurrentReadSeq = Temp_One_Read.getUnmatchedSeq();
		End = Temp_One_Read.MatchedRelPos + g_SpacerBeforeAfter + (RangeIndex) * Temp_One_Read.InsertSize; /////////
		Start = End - (RangeIndex + 1) * Temp_One_Read.InsertSize;
		char RightChar;
		RightChar = CurrentReadSeq[Temp_One_Read.getReadLengthMinus()];
//std::cout << "Starting to fit the close end with character" << RightChar << "\n";
		if (RightChar != 'N') {
			for (int pos = Start; pos < End; pos++) {
				if (CurrentChrSeq[pos] == RightChar) {
					PD[0].push_back(pos);
				}
			}
		}
//std::cout << "1\t" << PD[0].size() << "\t" << PD[1].size() << std::endl;
//        LOG_DEBUG(*logStream << "1\t" << PD[0].size() << "\t" << PD[1].size() << std::endl);
		if (PD[0].size()) {
			CheckRight_Close_Perfect(Temp_One_Read, CurrentChrSeq, CurrentReadSeq, PD, BP_Start, BP_End, 1, UP);
		}
//        LOG_DEBUG(*logStream << UP.size() << std::endl);
		if (UP.empty()) {
		} else {
			Temp_One_Read.Used = false;
			Temp_One_Read.UP_Close.swap(UP);
			UP.clear();
		}
	}

	if (Temp_One_Read.UP_Close.size()) {
		if (Temp_One_Read.MatchedD == Plus) {
			CurrentReadSeq = ReverseComplement(Temp_One_Read.getUnmatchedSeq());
			Start = Temp_One_Read.MatchedRelPos + g_SpacerBeforeAfter - (RangeIndex - 1) * Temp_One_Read.InsertSize; /////////////
			End = Start + (RangeIndex + 1) * Temp_One_Read.InsertSize;
			char LeftChar;
			LeftChar = CurrentReadSeq[0];
			if (LeftChar != 'N') {
				{
					for (int pos = Start; pos < End; pos++) {
						if (CurrentChrSeq[pos] == LeftChar) {
							PD[0].push_back(pos);
						}
					}
				}
			}

			if (PD[0].size()) {
				CheckLeft_Close_Perfect(Temp_One_Read, CurrentChrSeq, CurrentReadSeq, PD, BP_Start, BP_End, 1, UP); // LengthStr
			}
			if (UP.empty()) {
			} else {
				Temp_One_Read.Used = false;
				Temp_One_Read.UP_Close.swap(UP);
				UP.clear();
			}
		} else if (Temp_One_Read.MatchedD == Minus) {

			CurrentReadSeq = Temp_One_Read.getUnmatchedSeq();
			End = Temp_One_Read.MatchedRelPos + g_SpacerBeforeAfter + (RangeIndex - 1) * Temp_One_Read.InsertSize; /////////
			Start = End - (RangeIndex + 1) * Temp_One_Read.InsertSize;
			char RightChar;
			RightChar = CurrentReadSeq[Temp_One_Read.getReadLengthMinus()];
			//std::cout << "Starting to fit the close end with character" << RightChar << "\n";
			if (RightChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChrSeq[pos] == RightChar) {
						PD[0].push_back(pos);
					}
				}
			}
			//std::cout << "1\t" << PD[0].size() << "\t" << PD[1].size() << std::endl;
			//        LOG_DEBUG(*logStream << "1\t" << PD[0].size() << "\t" << PD[1].size() << std::endl);
			if (PD[0].size()) {
				CheckRight_Close_Perfect(Temp_One_Read, CurrentChrSeq, CurrentReadSeq, PD, BP_Start, BP_End, 1, UP);
			}
			//        LOG_DEBUG(*logStream << UP.size() << std::endl);
			if (UP.empty()) {
			} else {
				Temp_One_Read.Used = false;
				Temp_One_Read.UP_Close.swap(UP);
				UP.clear();
			}
		}
	}
	return;
}

void GetCloseEnd(const std::string &CurrentChrSeq, SPLIT_READ &Temp_One_Read) {
	const int MaxRange = 2; // 3 insert size away
//if (Temp_One_Read.Name == "@M01144:44:000000000-A6N99:1:1104:14364:9012/2") {
//    std::cout << "found @M01144:44:000000000-A6N99:1:1104:14364:9012/2" << std::endl;
//}
	for (int RangeIndex = 0; RangeIndex < MaxRange; RangeIndex++) {
		GetCloseEndInner(CurrentChrSeq, Temp_One_Read, RangeIndex);
//std::cout << "\nfirst: " << Temp_One_Read.Name << " " << Temp_One_Read.UP_Close.size() << std::endl;
		if (Temp_One_Read.UP_Close.size() == 0) { // no good close ends found
			//if (Temp_One_Read.Name == "@M01144:44:000000000-A6N99:1:1104:14364:9012/2") {
			//    std::cout << "found step 1 @M01144:44:000000000-A6N99:1:1104:14364:9012/2" << std::endl;
			//}
			//std::cout << Temp_One_Read.Name << " " << Temp_One_Read.UP_Close.size() << std::endl;
			Temp_One_Read.setUnmatchedSeq(ReverseComplement(Temp_One_Read.getUnmatchedSeq()));
			GetCloseEndInner(CurrentChrSeq, Temp_One_Read, RangeIndex);
			//std::cout << "second: " << Temp_One_Read.Name << " " << Temp_One_Read.UP_Close.size() << std::endl;
		}
//if (Temp_One_Read.UP_Close.size()==0) {
		/*if (Temp_One_Read.Name == "@M01144:44:000000000-A6N99:1:1104:14364:9012/2") {
		 std::cout << "found step 2 @M01144:44:000000000-A6N99:1:1104:14364:9012/2" << std::endl;
		 }
		 */
//GetCloseEndInnerPerfectMatch( CurrentChrSeq, Temp_One_Read, RangeIndex );
//std::cout << "third: " << Temp_One_Read.Name << " " << Temp_One_Read.UP_Close.size() << std::endl;
//}
//if (Temp_One_Read.UP_Close.size()==0) { // no good close ends found
		/*if (Temp_One_Read.Name == "@M01144:44:000000000-A6N99:1:1104:14364:9012/2") {
		 std::cout << "found step 3 @M01144:44:000000000-A6N99:1:1104:14364:9012/2" << std::endl;
		 }
		 */
//Temp_One_Read.setUnmatchedSeq( ReverseComplement( Temp_One_Read.getUnmatchedSeq() ) );
//GetCloseEndInnerPerfectMatch( CurrentChrSeq, Temp_One_Read, RangeIndex );
//std::cout << "fourth: " << Temp_One_Read.Name << " " << Temp_One_Read.UP_Close.size() << "\n" <<  std::endl;
//}
		/*
		 if (Temp_One_Read.Name == "@M01144:44:000000000-A6N99:1:1104:14364:9012/2") {
		 std::cout << "found step 4 @M01144:44:000000000-A6N99:1:1104:14364:9012/2" << std::endl;
		 }
		 */
		if (Temp_One_Read.hasCloseEnd()) {
			break;
		}
	}

	/*
	 GetCloseEndInnerPerfectMatch( CurrentChrSeq, Temp_One_Read );

	 SortedUniquePoints First_UP = Temp_One_Read.UP_Close; // backup

	 Temp_One_Read.setUnmatchedSeq( ReverseComplement( Temp_One_Read.getUnmatchedSeq() ) );

	 GetCloseEndInnerPerfectMatch( CurrentChrSeq, Temp_One_Read );



	 if (First_UP.size() + Temp_One_Read.UP_Close.size()) {
	 std::cout << First_UP.size() << " " << Temp_One_Read.UP_Close.size() << std::endl;
	 if (First_UP.size() > Temp_One_Read.UP_Close.size())
	 Temp_One_Read.UP_Close = First_UP;
	 }
	 else {
	 GetCloseEndInner( CurrentChrSeq, Temp_One_Read );
	 First_UP = Temp_One_Read.UP_Close;
	 Temp_One_Read.setUnmatchedSeq( ReverseComplement( Temp_One_Read.getUnmatchedSeq() ) );
	 GetCloseEndInner( CurrentChrSeq, Temp_One_Read );
	 if (First_UP.size() > Temp_One_Read.UP_Close.size())
	 Temp_One_Read.UP_Close = First_UP;
	 }
	 */
}

unsigned CountElements(const FarEndSearchPerRegion *OneRegionSearchResult_input, int levels) {
	unsigned Sum = 0;

	for (int LevelIndex = 0; LevelIndex < levels; LevelIndex++) {
		Sum += OneRegionSearchResult_input->PD_Plus[LevelIndex].size() + OneRegionSearchResult_input->PD_Minus[LevelIndex].size();
	}
	return Sum;
}

void ExtendMatchPerfect(SPLIT_READ &read, const std::string &readSeq, const std::vector<FarEndSearchPerRegion*> &WholeGenomeSearchResult_input, const short minimumLengthToReportMatch, const short BP_End, const short CurrentLength, SortedUniquePoints &UP) {
//UserDefinedSettings *userSettings = UserDefinedSettings::Instance();
	const char CurrentChar = readSeq[CurrentLength];
	const char CurrentCharRC = Convert2RC4N[(short) CurrentChar];
	bool AllEmpty = true;
	std::vector<FarEndSearchPerRegion*> WholeGenomeSearchResult_output;
	for (unsigned IndexOfRegion = 0; IndexOfRegion < WholeGenomeSearchResult_input.size(); IndexOfRegion++) {
		const FarEndSearchPerRegion *CurrentRegion_input = WholeGenomeSearchResult_input[IndexOfRegion];
		unsigned int Max_size = 0;
		for (int CheckedIndex = 0; CheckedIndex <= userSettings->ADDITIONAL_MISMATCH; CheckedIndex++) {
			if (Max_size < CurrentRegion_input->PD_Plus[CheckedIndex].size()) {
				Max_size = CurrentRegion_input->PD_Plus[CheckedIndex].size();
			}
			if (Max_size < CurrentRegion_input->PD_Minus[CheckedIndex].size()) {
				Max_size = CurrentRegion_input->PD_Minus[CheckedIndex].size();
			}
		}
		const std::string &chromosomeSeq = CurrentRegion_input->CurrentChromosome->getSeq();
		FarEndSearchPerRegion *CurrentRegion_output = new FarEndSearchPerRegion(CurrentRegion_input->CurrentChromosome, read.getTOTAL_SNP_ERROR_CHECKED(), Max_size);
		for (int i = 0; i <= read.getTOTAL_SNP_ERROR_CHECKED_Minus(); i++) {
			CategorizePositions(CurrentChar, chromosomeSeq, CurrentRegion_input->PD_Plus, CurrentRegion_output->PD_Plus, i, 1, userSettings->ADDITIONAL_MISMATCH);
			CategorizePositions(CurrentCharRC, chromosomeSeq, CurrentRegion_input->PD_Minus, CurrentRegion_output->PD_Minus, i, -1, userSettings->ADDITIONAL_MISMATCH);
		}
		if (CountElements(CurrentRegion_output, 1)) {
			AllEmpty = false;
			WholeGenomeSearchResult_output.push_back(CurrentRegion_output);
		} else {
			delete CurrentRegion_output;
		}
	}
	/*std::cout << "Matching " << CurrentLength << " length and char " << CurrentChar << ", " << CurrentCharRC << std::endl;
	 for (unsigned IndexOfRegion = 0; IndexOfRegion < WholeGenomeSearchResult_input.size(); IndexOfRegion++) {
	 std::cout << "Region index: " << IndexOfRegion << "\n";
	 for (int CheckedIndex = 0; CheckedIndex < read.getTOTAL_SNP_ERROR_CHECKED(); CheckedIndex++) {
	 std::cout << "R["<< CheckedIndex << "]=" << WholeGenomeSearchResult_input[IndexOfRegion].PD_Plus[CheckedIndex].size() << "/" << WholeGenomeSearchResult_input[IndexOfRegion].PD_Minus[CheckedIndex].size() << std::endl;
	 }
	 }
	 std::cout << "UP-size:" << UP.size() << "\n";*/
// this loop looks familiar; candidate for factoring out mini-function?
	if (AllEmpty == false) {
		const short CurrentLengthOutput = CurrentLength + 1;
		CheckBothPerfect(read, readSeq, WholeGenomeSearchResult_output, minimumLengthToReportMatch, BP_End, CurrentLengthOutput, UP);
	} else {
	} // else-if Sum
	for (unsigned int i = 0; i < WholeGenomeSearchResult_output.size(); i++) {
		delete WholeGenomeSearchResult_output[i];
	}
}

void ExtendMatch(SPLIT_READ &read, const std::string &readSeq, const std::vector<FarEndSearchPerRegion*> &WholeGenomeSearchResult_input, const short minimumLengthToReportMatch, const short BP_End, const short CurrentLength, SortedUniquePoints &UP) {
	const char CurrentChar = readSeq[CurrentLength];
	const char CurrentCharRC = Convert2RC4N[(short) CurrentChar];
	bool AllEmpty = true;
	std::vector<FarEndSearchPerRegion*> WholeGenomeSearchResult_output;
	for (unsigned IndexOfRegion = 0; IndexOfRegion < WholeGenomeSearchResult_input.size(); IndexOfRegion++) {
		const FarEndSearchPerRegion *CurrentRegion_input = WholeGenomeSearchResult_input[IndexOfRegion];
		unsigned int Max_size = 0;
		for (int CheckedIndex = 0; CheckedIndex < read.getTOTAL_SNP_ERROR_CHECKED(); CheckedIndex++) {
			if (Max_size < CurrentRegion_input->PD_Plus[CheckedIndex].size()) {
				Max_size = CurrentRegion_input->PD_Plus[CheckedIndex].size();
			}
			if (Max_size < CurrentRegion_input->PD_Minus[CheckedIndex].size()) {
				Max_size = CurrentRegion_input->PD_Minus[CheckedIndex].size();
			}
		}
		const std::string &chromosomeSeq = CurrentRegion_input->CurrentChromosome->getSeq();
		FarEndSearchPerRegion *CurrentRegion_output = new FarEndSearchPerRegion(CurrentRegion_input->CurrentChromosome, read.getTOTAL_SNP_ERROR_CHECKED(), Max_size);
		for (int i = 0; i <= read.getTOTAL_SNP_ERROR_CHECKED_Minus(); i++) {
			CategorizePositions(CurrentChar, chromosomeSeq, CurrentRegion_input->PD_Plus, CurrentRegion_output->PD_Plus, i, 1, read.getTOTAL_SNP_ERROR_CHECKED_Minus());
			CategorizePositions(CurrentCharRC, chromosomeSeq, CurrentRegion_input->PD_Minus, CurrentRegion_output->PD_Minus, i, -1, read.getTOTAL_SNP_ERROR_CHECKED_Minus());
		}
		if (CountElements(CurrentRegion_output, read.getTOTAL_SNP_ERROR_CHECKED())) {
			AllEmpty = false;
			WholeGenomeSearchResult_output.push_back(CurrentRegion_output);
		} else {
			delete CurrentRegion_output;
		}
	}
	/*std::cout << "Matching " << CurrentLength << " length and char " << CurrentChar << ", " << CurrentCharRC << std::endl;
	 for (unsigned IndexOfRegion = 0; IndexOfRegion < WholeGenomeSearchResult_input.size(); IndexOfRegion++) {
	 std::cout << "Region index: " << IndexOfRegion << "\n";
	 for (int CheckedIndex = 0; CheckedIndex < read.getTOTAL_SNP_ERROR_CHECKED(); CheckedIndex++) {
	 std::cout << "R["<< CheckedIndex << "]=" << WholeGenomeSearchResult_input[IndexOfRegion].PD_Plus[CheckedIndex].size() << "/" << WholeGenomeSearchResult_input[IndexOfRegion].PD_Minus[CheckedIndex].size() << std::endl;
	 }
	 }
	 std::cout << "UP-size:" << UP.size() << "\n";*/
// this loop looks familiar; candidate for factoring out mini-function?
	if (AllEmpty == false) {
		const short CurrentLengthOutput = CurrentLength + 1;
		CheckBoth(read, readSeq, WholeGenomeSearchResult_output, minimumLengthToReportMatch, BP_End, CurrentLengthOutput, UP);
	} else {
	} // else-if Sum
	for (unsigned int i = 0; i < WholeGenomeSearchResult_output.size(); i++) {
		delete WholeGenomeSearchResult_output[i];
	}
}

unsigned int minimumNumberOfMismatches(const std::vector<FarEndSearchPerRegion*> &WholeGenomeSearchResult_input, const unsigned int maxNumberMismatches) {
	unsigned int Sum = 0;
	unsigned int numberOfMismatches = 0;
	for (; numberOfMismatches <= maxNumberMismatches; numberOfMismatches++) {
		for (unsigned RegionIndex = 0; RegionIndex < WholeGenomeSearchResult_input.size(); RegionIndex++) {
			Sum += WholeGenomeSearchResult_input[RegionIndex]->PD_Plus[numberOfMismatches].size() + WholeGenomeSearchResult_input[RegionIndex]->PD_Minus[numberOfMismatches].size();
		}
		if (Sum != 0) {
			break;
		}
	}
	return numberOfMismatches;
}

void CheckBothPerfect(SPLIT_READ &read, const std::string &readSeq, const std::vector<FarEndSearchPerRegion*> &WholeGenomeSearchResult_input, const short minimumLengthToReportMatch, const short BP_End, const short CurrentLength, SortedUniquePoints &UP) {
//UserDefinedSettings *userSettings = UserDefinedSettings::Instance();
	unsigned NumberOfMatchPositionsWithLessMismatches = 0;
	int Sum = 0;

	if (CurrentLength >= minimumLengthToReportMatch && CurrentLength <= BP_End) {
		if (minimumNumberOfMismatches(WholeGenomeSearchResult_input, userSettings->ADDITIONAL_MISMATCH + 1) > g_maxMismatch[CurrentLength]) {
			return;
		}
		for (short numberOfMismatches = 0; numberOfMismatches < 1; numberOfMismatches++) {
			if (NumberOfMatchPositionsWithLessMismatches) {
				break;
			}

			Sum = 0;
			for (unsigned RegionIndex = 0; RegionIndex < WholeGenomeSearchResult_input.size(); RegionIndex++) {
				Sum += WholeGenomeSearchResult_input[RegionIndex]->PD_Plus[numberOfMismatches].size() + WholeGenomeSearchResult_input[RegionIndex]->PD_Minus[numberOfMismatches].size();
			}
			NumberOfMatchPositionsWithLessMismatches = Sum;
			if (Sum == 1 && CurrentLength >= minimumLengthToReportMatch + numberOfMismatches) {
				Sum = 0;
				if (userSettings->ADDITIONAL_MISMATCH > 0) { // what if this is ADDITIONAL_MISMATCH is 0? Do you save anything then?

					// what if j +ADD_MISMATCH exceeds the max allowed number of mismatches (the PD_element does not exist?)

					// feeling the need to comment here - so factor out?
					// only report reads if there are no reads with fewer mismatches or only one or two more mismatches
					unsigned int regionWithMatch = 0;
					for (short mismatchCount = 0; mismatchCount <= numberOfMismatches + userSettings->ADDITIONAL_MISMATCH; mismatchCount++) {
						for (unsigned RegionIndex = 0; RegionIndex < WholeGenomeSearchResult_input.size(); RegionIndex++) {
							unsigned int hitsInRegion = WholeGenomeSearchResult_input[RegionIndex]->PD_Plus[mismatchCount].size() + WholeGenomeSearchResult_input[RegionIndex]->PD_Minus[mismatchCount].size();
							Sum += hitsInRegion;
							if (hitsInRegion > 0) {
								regionWithMatch = RegionIndex;
							}
						}
					}
					/*if (read.Name=="@read_6990/2" ) {
					 std::cout << "In CFE: CurrentLength = " << CurrentLength << ", mismatch count = " << numberOfMismatches << ", maxMismatch = " << g_maxMismatch[CurrentLength] << std::endl;
					 for (unsigned RegionIndex = 0; RegionIndex < WholeGenomeSearchResult_input.size(); RegionIndex++) {
					 std::cout << "Region " << RegionIndex << " is " <<  WholeGenomeSearchResult_input[ RegionIndex ].CurrentChromosome->getName() << ":" << WholeGenomeSearchResult_input[ RegionIndex ].PD_Plus[0].size() << "-" <<
					 WholeGenomeSearchResult_input[ RegionIndex ].PD_Minus[0].size()<< "\n";
					 }
					 for (short k=0;k<=read.getMAX_SNP_ERROR(); k++) {
					 std::cout << k << "\t" << WholeGenomeSearchResult_input[0].PD_Plus[k].size() + WholeGenomeSearchResult_input[0].PD_Minus[k].size() << "\n";
					 }}*/
					if (Sum == 1 && (unsigned) numberOfMismatches <= g_maxMismatch[CurrentLength]) {
						// why I love constructors
						UniquePoint MatchPosition;

						const FarEndSearchPerRegion *hitRegion = WholeGenomeSearchResult_input[regionWithMatch];
						if (WholeGenomeSearchResult_input[regionWithMatch]->PD_Plus[numberOfMismatches].size() == 1) {
							UniquePoint PlusMatch(hitRegion->CurrentChromosome, CurrentLength, hitRegion->PD_Plus[numberOfMismatches][0], FORWARD, SENSE, numberOfMismatches);
							MatchPosition = PlusMatch;
						} else {
							UniquePoint MinMatch(hitRegion->CurrentChromosome, CurrentLength, hitRegion->PD_Minus[numberOfMismatches][0], BACKWARD, ANTISENSE, numberOfMismatches);
							MatchPosition = MinMatch;
						}

						if (CheckMismatches(WholeGenomeSearchResult_input[regionWithMatch]->CurrentChromosome->getSeq(), readSeq, MatchPosition, read.FarEndMismatch)) {
							UP.push_back(MatchPosition);
							break;
						} // if CheckMismatches
					} // if Sum==1
				} // if AdditionalMismatches
			} // if sumsize ==1
		} // for-loop
	} // if length of match is sufficient to be reportable

	if (CurrentLength < BP_End) {
		ExtendMatchPerfect(read, readSeq, WholeGenomeSearchResult_input, minimumLengthToReportMatch, BP_End, CurrentLength, UP);
	}
}

//CheckBoth(Temp_One_Read, Temp_One_Read.getUnmatchedSeq(), WholeGenomeSearchResult, BP_Start, BP_End, 1, UP);
void CheckBoth(SPLIT_READ &read, const std::string &readSeq, const std::vector<FarEndSearchPerRegion*> &WholeGenomeSearchResult_input, const short minimumLengthToReportMatch, const short BP_End, const short CurrentLength, SortedUniquePoints &UP) {
//UserDefinedSettings *userSettings = UserDefinedSettings::Instance();
	unsigned NumberOfMatchPositionsWithLessMismatches = 0;
	int Sum = 0;

	if (CurrentLength >= minimumLengthToReportMatch && CurrentLength <= BP_End) {
		if (minimumNumberOfMismatches(WholeGenomeSearchResult_input, read.getMAX_SNP_ERROR()) > g_maxMismatch[CurrentLength]) {
			return;
		}
		for (short numberOfMismatches = 0; numberOfMismatches <= read.getMAX_SNP_ERROR(); numberOfMismatches++) {
			if (NumberOfMatchPositionsWithLessMismatches) {
				break;
			}

			Sum = 0;
			for (unsigned RegionIndex = 0; RegionIndex < WholeGenomeSearchResult_input.size(); RegionIndex++) {
				Sum += WholeGenomeSearchResult_input[RegionIndex]->PD_Plus[numberOfMismatches].size() + WholeGenomeSearchResult_input[RegionIndex]->PD_Minus[numberOfMismatches].size();
			}
			NumberOfMatchPositionsWithLessMismatches = Sum;
			if (Sum == 1 && CurrentLength >= minimumLengthToReportMatch + numberOfMismatches) {
				Sum = 0;
				if (userSettings->ADDITIONAL_MISMATCH > 0) { // what if this is ADDITIONAL_MISMATCH is 0? Do you save anything then?

					// what if j +ADD_MISMATCH exceeds the max allowed number of mismatches (the PD_element does not exist?)

					// feeling the need to comment here - so factor out?
					// only report reads if there are no reads with fewer mismatches or only one or two more mismatches
					unsigned int regionWithMatch = 0;
					for (short mismatchCount = 0; mismatchCount <= numberOfMismatches + userSettings->ADDITIONAL_MISMATCH; mismatchCount++) {
						for (unsigned RegionIndex = 0; RegionIndex < WholeGenomeSearchResult_input.size(); RegionIndex++) {
							unsigned int hitsInRegion = WholeGenomeSearchResult_input[RegionIndex]->PD_Plus[mismatchCount].size() + WholeGenomeSearchResult_input[RegionIndex]->PD_Minus[mismatchCount].size();
							Sum += hitsInRegion;
							if (hitsInRegion > 0) {
								regionWithMatch = RegionIndex;
							}
						}
					}
					/*if (read.Name=="@read_6990/2" ) {
					 std::cout << "In CFE: CurrentLength = " << CurrentLength << ", mismatch count = " << numberOfMismatches << ", maxMismatch = " << g_maxMismatch[CurrentLength] << std::endl;
					 for (unsigned RegionIndex = 0; RegionIndex < WholeGenomeSearchResult_input.size(); RegionIndex++) {
					 std::cout << "Region " << RegionIndex << " is " <<  WholeGenomeSearchResult_input[ RegionIndex ].CurrentChromosome->getName() << ":" << WholeGenomeSearchResult_input[ RegionIndex ].PD_Plus[0].size() << "-" <<
					 WholeGenomeSearchResult_input[ RegionIndex ].PD_Minus[0].size()<< "\n";
					 }
					 for (short k=0;k<=read.getMAX_SNP_ERROR(); k++) {
					 std::cout << k << "\t" << WholeGenomeSearchResult_input[0].PD_Plus[k].size() + WholeGenomeSearchResult_input[0].PD_Minus[k].size() << "\n";
					 }}*/
					if (Sum == 1 && (unsigned) numberOfMismatches <= g_maxMismatch[CurrentLength]) {
						// why I love constructors
						UniquePoint MatchPosition;

						const FarEndSearchPerRegion *hitRegion = WholeGenomeSearchResult_input[regionWithMatch];
						if (WholeGenomeSearchResult_input[regionWithMatch]->PD_Plus[numberOfMismatches].size() == 1) {
							UniquePoint PlusMatch(hitRegion->CurrentChromosome, CurrentLength, hitRegion->PD_Plus[numberOfMismatches][0], FORWARD, SENSE, numberOfMismatches);
							MatchPosition = PlusMatch;
						} else {
							UniquePoint MinMatch(hitRegion->CurrentChromosome, CurrentLength, hitRegion->PD_Minus[numberOfMismatches][0], BACKWARD, ANTISENSE, numberOfMismatches);
							MatchPosition = MinMatch;
						}

						if (CheckMismatches(WholeGenomeSearchResult_input[regionWithMatch]->CurrentChromosome->getSeq(), readSeq, MatchPosition, read.FarEndMismatch)) {
							UP.push_back(MatchPosition);
							break;
						} // if CheckMismatches
					} // if Sum==1
				} // if AdditionalMismatches
			} // if sumsize ==1
		} // for-loop
	} // if length of match is sufficient to be reportable

	if (CurrentLength < BP_End) {
		ExtendMatch(read, readSeq, WholeGenomeSearchResult_input, minimumLengthToReportMatch, BP_End, CurrentLength, UP);
	}
}

void CleanUniquePoints(SortedUniquePoints &Input_UP) {
	SortedUniquePoints TempUP; //vector <UniquePoint> UP_Close; UP_Far
	UniquePoint LastUP = Input_UP[Input_UP.size() - 1];
	char LastDirection = LastUP.Direction;
	char LastStrand = LastUP.Strand;
	const Chromosome *LastChr = LastUP.chromosome_p;
//LastChr == .chromosome_p
	unsigned int Terminal;

	if (LastDirection == FORWARD) {
		Terminal = LastUP.AbsLoc - LastUP.LengthStr;
		for (unsigned i = 0; i < Input_UP.size(); i++) {
			if (LastChr != Input_UP[i].chromosome_p) {
				continue;
			}
			if (Input_UP[i].Direction == LastDirection && Input_UP[i].Strand == LastStrand) {
				if (Terminal == Input_UP[i].AbsLoc - Input_UP[i].LengthStr) {
					TempUP.push_back(Input_UP[i]);
				}
			}
		}
	} else if (LastDirection == BACKWARD) {
		Terminal = LastUP.AbsLoc + LastUP.LengthStr;
		for (unsigned i = 0; i < Input_UP.size(); i++) {
			if (LastChr != Input_UP[i].chromosome_p) {
				continue;
			}
			if (Input_UP[i].Direction == LastDirection && Input_UP[i].Strand == LastStrand) {
				if (Terminal == Input_UP[i].AbsLoc + Input_UP[i].LengthStr) {
					TempUP.push_back(Input_UP[i]);
				}
			}
		}
	}
	Input_UP.clear();
	Input_UP = TempUP;
}

