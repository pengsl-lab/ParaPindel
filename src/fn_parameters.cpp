/*** fn_parameters.cpp contains functions related to the parameters, which have been refactored out of the pindel.main */

#include <iostream>
#include <sstream>

#include "logstream.h"
#include "user_defined_settings.h"
#include "fn_parameters.h"
#include "logdef.h"
#include "paraPindel.h" // for logStream

/* 'defineParameters' defines the parameters to be used by Pindel. Takes the variables from the calling function as argument for those variables which
 do not need to be stored in the par structure. */
// _________________________
// ___D_FG__K___O_____UVWX_Z
//g_MinClose
void defineParameters(std::vector<Parameter *>& parameters)
{
   //UserDefinedSettings* userSettings = UserDefinedSettings::Instance();
   parameters.push_back(
      new StringParameter(&userSettings->referenceFilename, "-f", "--fasta",
                          "the reference genome sequences in fasta format", true, ""));
   parameters.push_back(
      new StringParameter(
         &(userSettings->pindelFilename),
         "-p",
         "--pindel-file",
         "the Pindel input file; either this, a pindel configuration file (consisting of multiple pindel filenames) or a bam configuration file is required",
         false, ""));
   parameters.push_back(
      new StringParameter(
         &userSettings->bamConfigFilename,
         "-i",
         "--config-file",
         "the bam config file; either this, a pindel input file, or a pindel config file is required. Per line: path and file name of bam, insert size and sample tag.     For example: /data/tumour.bam  400  tumour",
         false, ""));
   parameters.push_back(
      new StringParameter(
         &userSettings->pindelConfigFilename,
         "-P",
         "--pindel-config-file",
         "the pindel config file, containing the names of all Pindel files that need to be sampled; either this, a bam config file or a pindel input file is required. Per line: path and file name of pindel input. Example: /data/tumour.txt",
         false, ""));
   parameters.push_back(
      new StringParameter(&userSettings->outputFilename, "-o", "--output-prefix", "Output prefix;", true, ""));
   parameters.push_back(
      new StringParameter(
         &userSettings->userDefinedRegion,
         "-c",
         "--chromosome",
         "Which chr/fragment. Pindel will process reads for one chromosome each time. ChrName must be the same as in reference sequence and in read file. '-c ALL' will make Pindel loop over all chromosomes. The search for indels and SVs can also be limited to a specific region; -c 20:10,000,000 will only look for indels and SVs after position 10,000,000 = [10M, end], -c 20:5,000,000-15,000,000 will report indels in the range between and including the bases at position 5,000,000 and 15,000,000 = [5M, 15M]. (default ALL)",
         false, "ALL"));

   parameters.push_back(
      new BoolParameter(&userSettings->showHelp, "-h", "--help",
                        "show the command line options of Pindel", false, false));
   parameters.push_back(
      new BoolParameter(&userSettings->SearchDiscordantReadPair, "-R", "--RP",
                        "search for discordant read-pair to improve sensitivity (default true)", false, true));

   parameters.push_back(
      new IntParameter(
         &userSettings->minClose,
         "-H",
         "--min_distance_to_the_end",
         "the minimum number of bases required to match reference (default 8).",
         false, 8));

   parameters.push_back(
      new IntParameter(
         &userSettings->numThreads,
         "-T",
         "--number_of_threads",
         "the number of threads Pindel will use (default 1).",
         false, 1));

   std::stringstream ss;
   ss << "the maximum size of structural variations to be detected; the higher this number, the greater "
      "the number of SVs reported, but the computational cost and memory requirements increase, as "
      "does the rate of false "
      "positives. 1=128, 2=512, 3=2,048, 4=8,092, 5=32,368, 6=129,472, 7=517,888, 8=2,071,552, 9=8,286,208. "
      "(maximum "<< g_MAX_RANGE_INDEX <<", default 2)";
   parameters.push_back(
      new IntParameter(
         &userSettings->MaxRangeIndex,
         "-x",
         "--max_range_index",
         ss.str(), false, 2));
   parameters.push_back(
      new FloatParameter(
         &userSettings->FLOAT_WINDOW_SIZE,
         "-w",
         "--window_size",
         "for saving RAM, divides the reference in bins of X million bases and only analyzes the reads that belong in the current bin, (default 5 (=5 million))",
		 false, 5.0));

   parameters.push_back(
            new FloatParameter(
               &userSettings->FLOAT_CROSS_WINDOW_SIZE,
               "-W",
               "--cross_window_size",
               "In a bin, some reads may exceed the boundary of the current bin, so the range of the current bin is expanded. The expanded sequence is a sequence of length cross_window_size in the next window in sequence. (default 5.0 (=5k)).",
               false, 5.0));


   parameters.push_back(
      new FloatParameter(&userSettings->Seq_Error_Rate, "-e",
                         "--sequencing_error_rate",
                         "the expected fraction of sequencing errors " // re-explain this or split this
                         "(default 0.01)", false, 0.01));

   parameters.push_back(
      new FloatParameter(&userSettings->sensitivity, "-E",
                         "--sensitivity",
                         "Pindel only reports reads if they can be fit around an event within a certain number of mismatches. If the fraction of sequencing errors is 0.01, (so we'd expect a total "
                         "error rate of 0.011 since on average 1 in 1000 bases is a SNP) and "
                         "pindel calls a deletion, but there are 4 mismatched bases in the new fit of the pindel read (100 bases) to the reference genome, Pindel would calculate that with an error rate of "
                         "0.01 (=0.011 including SNPs) the chance that there are 0, 1 or 2 mismatched bases in the reference genome is 90%. Setting -E to .90 (=90%) will therefore"
                         "throw away all reads with 3 or more mismatches, even though that means that you throw away 1 in 10 valid reads. Increasing this parameter to say 0.99 will increase "
                         "the sensitivity of pindel though you may get more false positives, decreasing the parameter ensures you only get very good matches but pindel may not find as many events. "
                         "(default 0.95)", false, 0.95));

   parameters.push_back(
      new FloatParameter(
         &userSettings->MaximumAllowedMismatchRate,
         "-u",
         "--maximum_allowed_mismatch_rate",
         "Only reads with more than this fraction of mismatches than the reference genome will be considered as harboring potential SVs. "
         "(default 0.02)", false, 0.02));

   parameters.push_back(
      new IntParameter(&userSettings->NM, "-n", "--NM",
                       "the minimum number of edit distance between reads and reference genome (default 2)." " reads at least NM edit distance (>= NM) will be realigned",
                       false, 2));

   parameters.push_back(
      new BoolParameter(&userSettings->Analyze_INV, "-r", "--report_inversions",
                        "report inversions " "(default true)", false, true));
   parameters.push_back(
      new BoolParameter(&userSettings->Analyze_TD, "-t", "--report_duplications",
                        "report tandem duplications " "(default true)", false, true));

   parameters.push_back(
      new BoolParameter(
         &userSettings->Analyze_LI,
         "-l",
         "--report_long_insertions",
         "report insertions of which the full sequence cannot be deduced because of their length "
         "(default false)", false, false));
   parameters.push_back(
      new BoolParameter(&userSettings->Analyze_BP, "-k", "--report_breakpoints",
                        "report breakpoints " "(default false)", false, false));

   parameters.push_back(
      new BoolParameter(
         &userSettings->ReportCloseMappedRead,
         "-s",
         "--report_close_mapped_reads",
         "report reads of which only one end (the one closest to the mapped read of the paired-end read) "
         "could be mapped. " "(default false)", false, false));

   parameters.push_back(
      new BoolParameter(
         &userSettings->reportOnlyCloseMappedReads,
         "-S",
         "--report_only_close_mapped_reads",
         "do not search for SVs, only report reads of which only one end (the one closest to the mapped read of the paired-end read) "
         "could be mapped (the output file can then be used as an input file for another run of pindel, which may save size if you need to transfer files). " "(default false)", false, false));
   parameters.push_back(
      new BoolParameter(
         &userSettings->reportInterchromosomalEvents,
         "-I",
         "--report_interchromosomal_events",
         "search for interchromosomal events. Note: will require the computer to have at least 4 GB of memory (default false)", false, false));

   parameters.push_back(
      new BoolParameter(
         &userSettings->IndelCorrection,
         "-C",
         "--IndelCorrection",
         "search for consensus indels to corret contigs (default false)", false, false));

   parameters.push_back(
      new BoolParameter(
         &userSettings->NormalSamples,
         "-N",
         "--NormalSamples",
         "Turn on germline filtering, less sensistive and you may miss somatic calls (default false)", false, false));

   parameters.push_back(
      new StringParameter(
         &userSettings->breakdancerFilename,
         "-b",
         "--breakdancer",
         "Pindel is able to use calls from other SV methods such as BreakDancer to further increase sensitivity and specificity.                    BreakDancer result or calls from any methods must in the format:   ChrA LocA stringA ChrB LocB stringB other",
         false, ""));

   parameters.push_back(
      new StringParameter(
         &userSettings->inf_InclusiveBedFileName,
         "-j",
         "--include",
         "If you want Pindel to process a set of regions, please provide a bed file here: chr start end",
         false, ""));

   parameters.push_back(
      new StringParameter(
         &userSettings->inf_ExclusiveBedFileName,
         "-J",
         "--exclude",
         "If you want Pindel to skip a set of regions, please provide a bed file here: chr start end",
         false, ""));

   parameters.push_back(
      new IntParameter(
         &userSettings->ADDITIONAL_MISMATCH,
         "-a",
         "--additional_mismatch",
         "Pindel will only map part of a read to the reference genome if there are no other candidate positions with no more than the specified number of mismatches position. The bigger the value, the more accurate but less sensitive. (minimum value 1, default value 1)",
         false, 1));

   parameters.push_back(
      new IntParameter(
         &userSettings->Min_Perfect_Match_Around_BP,
         "-m",
         "--min_perfect_match_around_BP",
         "at the point where the read is split into two, there should at least be "
         "this number of perfectly matching bases between read and reference (default value 3)",
         false, 3));
   /*parameters.push_back(
       new IntParameter(&userSettings->MIN_IndelSize_NT, "-n", "--min_NT_size",
                        "only report inserted (NT) sequences in deletions greater than this size "
                        "(default 50)", false, 50));*/
   // TODO: Make sure MIN_IndelSize_NT is > 0, make into unsigned.
   parameters.push_back(
      new IntParameter(&userSettings->MIN_IndelSize_Inversion, "-v",
                       "--min_inversion_size",
                       "only report inversions greater than this number of bases "
                       "(default 50)", false, 50));
   parameters.push_back(
      new IntParameter(
         &userSettings->Min_Num_Matched_Bases,
         "-d",
         "--min_num_matched_bases",
         "only consider reads as evidence if they map with more than X bases to the reference. "
         "(default 30)", false, 30));
   parameters.push_back(
      new UIntParameter(
         &userSettings->BalanceCutoff,
         "-B",
         "--balance_cutoff",
         "the number of bases of a SV above which a more stringent filter is applied which demands "
         "that both sides of the SV are mapped with sufficiently long strings of bases "
         "(default 100)", false, 100));
   parameters.push_back(
      new UIntParameter(
         &userSettings->minimalAnchorQuality,
         "-A",
         "--anchor_quality",
         "the minimal mapping quality of the reads Pindel uses as anchor "
         "If you only need high confident calls, set to 30 or higher"
         "(default 0)", false, 0));
   parameters.push_back(
      new UIntParameter(
         &userSettings->NumRead2ReportCutOff,
         "-M",
         "--minimum_support_for_event",
         "Pindel only calls events which have this number or more supporting reads "
         "(default 1)", false, 1));
   parameters.push_back(
      new StringParameter(
         &userSettings->inf_AssemblyInputFilename,
         "-z",
         "--input_SV_Calls_for_assembly",
         "A filename of a list of SV calls for assembling breakpoints \n"
         "Types: DEL, INS, DUP, INV, CTX and ITX \n"
         "File format: Type chrA posA Confidence_Range_A chrB posB Confidence_Range_B \n"
         "Example: DEL chr1 10000 50 chr2 20000 100 "
         "", false, ""));
   parameters.push_back(
      new StringParameter(
         &userSettings->inf_GenotypingInputFilename,
         "-g",
         "--genotyping",
         "gentype variants if -i is also used."
         "", false, ""));
   parameters.push_back(
      new StringParameter(
         &userSettings->breakdancerOutputFilename,
         "-Q",
         "--output_of_breakdancer_events",
         "If breakdancer input is used, you can specify a filename here to write the confirmed breakdancer "
         "events with their exact breakpoints to " "The list of BreakDancer calls with Pindel support information. Format: chr   Loc_left   Loc_right   size   type   index. " "            For example, \"1	72766323 	72811840 	45516	D	11970\" means the deletion event chr1:72766323-72811840 of size 45516 is reported as an event with index 11970 in Pindel report of deletion. "
         "", false, ""));
   parameters.push_back(
      new StringParameter(
         &userSettings->logFilename,
         "-L",
         "--name_of_logfile",
         "Specifies a file to write Pindel's log to (default: no logfile, log is written to the screen/stdout)",
         false, ""));
   parameters.push_back(
      new StringParameter(
         &userSettings->PloidyFileName,
         "-Y",
         "--Ploidy",
         "a file with Ploidy information per chr for genotype."
         " per line: ChrName Ploidy."
         " For example, chr1 2",
         false, ""));

   // Dispersed duplication detection
   parameters.push_back(new BoolParameter(&userSettings->Analyze_DD, "-q", "--detect_DD",
                                          "Flag indicating whether to detect dispersed duplications."
                                          " (default: false)",
                                          false, false));

   // Dispersed duplication detection: maximum distance between breakpoints estimated
   // from opposite strands to assume that they refer to the same event.
   parameters.push_back(new IntParameter(&userSettings->MAX_DD_BREAKPOINT_DISTANCE,
                                         "", "--MAX_DD_BREAKPOINT_DISTANCE",
                                         "Maximum distance between dispersed duplication breakpoints to assume"
                                         " they refer to the same event. (default: 350)",
                                         false, 350));

   // Dispersed duplication detection: maximum distance between reads for them to
   // provide evidence for a single breakpoint.
   parameters.push_back(new IntParameter(&userSettings->MAX_DISTANCE_CLUSTER_READS,
                                         "", "--MAX_DISTANCE_CLUSTER_READS",
                                         "Maximum distance between reads for them to provide evidence for a single"
                                         " breakpoint for dispersed duplications. (default: 100)",
                                         false, 100));
   // Dispersed duplication detection: minimum number of reads needed for calling
   // a breakpoint.
   parameters.push_back(new IntParameter(&userSettings->MIN_DD_CLUSTER_SIZE,
                                         "", "--MIN_DD_CLUSTER_SIZE",
                                         "Minimum number of reads needed for calling a breakpoint for dispersed"
                                         " duplications. (default: 3)",
                                         false, 3));
   // Dispersed duplication detection: minimum number of split reads for calling
   // exact breakpoint location.
   parameters.push_back(new IntParameter(&userSettings->MIN_DD_BREAKPOINT_SUPPORT,
                                         "", "--MIN_DD_BREAKPOINT_SUPPORT",
                                         "Minimum number of split reads for calling an exact breakpoint for"
                                         " dispersed duplications. (default: 3)",
                                         false, 3));
   // Dispersed duplication detection: minimum distance dispersed duplication events.
   parameters.push_back(new IntParameter(&userSettings->MIN_DD_MAP_DISTANCE,
                                         "", "--MIN_DD_MAP_DISTANCE",
                                         "Minimum mapping distance of read pairs for them to be considered"
                                         " discordant. (default: 8000)",
                                         false, 8000));
   // Dispersed duplication detection: whether or not to report discordant mates.
   parameters.push_back(new BoolParameter(&userSettings->DD_REPORT_DUPLICATION_READS,
                                          "", "--DD_REPORT_DUPLICATION_READS",
                                          "Report discordant sequences and positions for mates of reads mapping"
                                          " inside dispersed duplications. (default: false)",
                                          false, false));
}

/* 'findParameter' returns the index of the parameter with name 'name'; -1 if not found.*/
unsigned int findParameter(const std::string& name, const std::vector<Parameter *>& parameters)
{
   for (unsigned int parameterCounter = 0; parameterCounter < parameters.size(); parameterCounter++) {
      if (parameters[parameterCounter]->hasName(name)) {
         return parameterCounter;
      }
   }
   LOG_DEBUG(std::cout << "Result of FindParameter is -1\n");
   return -1;
}

/* 'readParameters' reads the parameters as entered in the command line. */
void readParameters(int argc, char *argv[], std::vector<Parameter *>& parameters)
{

   for (int argumentIndex = 1; argumentIndex < argc; argumentIndex++) {
      std::string currentArgument = argv[argumentIndex];

      //find argument in parameterlist
      int parameterIndex = findParameter(currentArgument, parameters);
      if (parameterIndex == -1) {
         LOG_ERROR(*logStream << "unknown argument: " << currentArgument << std::endl);
         return;
      }

      if (parameters[parameterIndex]->isUnary()) {
         parameters[parameterIndex]->setValue(true); // default
         if ((argumentIndex + 1 < argc) && (argv[argumentIndex + 1][0]
                                            != '-')) { // so there are more arguments, and next one isn't regular -x
            if (tolower(argv[argumentIndex + 1][0]) == 'f'
                  || (argv[argumentIndex + 1][0] == '0')) {
               parameters[parameterIndex]->setValue(false);
            }
            argumentIndex++; // in any case increase the argument index
         }
      } else { // argument needs a parameter
         argumentIndex++; // move on to next argument in the list
         if (argumentIndex >= argc) {
            LOG_ERROR(*logStream << "argument of " << currentArgument << " lacking.\n");
            return;
         }
         if (argv[argumentIndex][0] == '-') {
            LOG_ERROR(*logStream << "argument of " << currentArgument
                      << " seems erroneous.\n");
            return;
         }
         // but if everything is allright,
         LOG_DEBUG(*logStream << "Giving " << currentArgument << " the value " << argv[ argumentIndex ] << std::endl);
         parameters[parameterIndex]->setValue(
            std::string(argv[argumentIndex]));
      }
   }
}

/* isReadsFileParam returns whether the parameter points to a read-containing file, and is therefore required,
 even though not both are required. */
bool isReadsFileParam(Parameter * param)
{
   return (param->hasName("-i") || param->hasName("-p"));
}

/* 'printHelp' prints all parameters available. */
void printHelp(const std::vector<Parameter *>& parameters)
{
   *logStream << std::endl
              << "Program:   ParaPindel (multi-processes parallel implementation of Pindel)"
              << std::endl;
   *logStream << "Multi-processes parallel implementation is based on " << Pindel_Version_str << std::endl;
   *logStream << "Contact:   Kai Ye <kaiye@xjtu.edu.cn> or Yaning Yang <yangyn@hnu.edu.cn>" << std::endl << std::endl;

   *logStream << "Usage:     paraPindel -f <reference.fa> -p <pindel_input>"
              << std::endl;
   *logStream << "           [and/or -i bam_configuration_file]" << std::endl;
   *logStream << "           -c <chromosome_name> -o <prefix_for_output_file>"
              << std::endl << std::endl;

   *logStream << "Required parameters:" << std::endl;

   for (unsigned int i = 0; i < parameters.size(); i++) {
      if (parameters[i]->isRequired() || isReadsFileParam(parameters[i])) {
         parameters[i]->describe();
      }
   }
   *logStream << "\nOptional parameters:" << std::endl;

   for (unsigned int parameterIndex = 0; parameterIndex < parameters.size(); parameterIndex++) {
      if (!parameters[parameterIndex]->isRequired() && !isReadsFileParam(
               parameters[parameterIndex])) {
         parameters[parameterIndex]->describe();
      }
   }
}

/* 'checkParameters' checks whether all required parameters have been set. */
bool checkParameters(const std::vector<Parameter *>& parameters)
{
   if (userSettings->showHelp) {
      printHelp( parameters );
      return false;
   }

   for (unsigned int parameterIndex = 0; parameterIndex < parameters.size(); parameterIndex++) {
      if (parameters[parameterIndex]->isRequired()
            && !parameters[parameterIndex]->isSet()) {
         LOG_ERROR(*logStream << "Required parameter "
                   << parameters[parameterIndex]-> getShortName() << "/"
                   << parameters[parameterIndex]-> getLongName() << " "
                   << parameters[parameterIndex]-> getDescription()
                   << " needs to be set." << std::endl);
         return false;
      } //if
   }
   // here handle the tricky fact that at least one of -i or -p needs be set; both are not required.
   bool hasBam = parameters[findParameter("-i",parameters)]->isSet();
   bool hasPin = parameters[findParameter("-p",parameters)]->isSet();
   bool hasPinConfig = parameters[findParameter("-P",parameters)]->isSet();
   if (!hasBam && !hasPin && !hasPinConfig) {
      LOG_ERROR(*logStream
                << "Bam and/or pindel input file required, use -p, -P and/or -i to designate input file(s)."
                << std::endl);
      return false;
   }
   return true;
}

