/*
 *   Options.h
 *
 *   Author: jtr
 */

#ifndef FLEXBAR_OPTIONS_H
#define FLEXBAR_OPTIONS_H

#include <string>
#include <iostream>

#include <tbb/concurrent_vector.h>

#include <seqan/basic.h>
#include <seqan/arg_parse.h>

#include "Enums.h"
#include "SeqRead.h"


typedef std::pair< SeqRead<seqan::CharString, seqan::CharString>*,
                   std::pair< tbb::atomic<unsigned long>, tbb::atomic<unsigned long> > > TAdapter;

struct Options{
	
	std::string readsFile, readsFile2, barReadsFile;
	std::string barcodeFile, adapterFile, barcode2File, adapter2File;
	std::string adapterSeq, targetName, logLevelStr, outCompression;
	
	bool isPaired, useAdapterFile, useNumberTag, useRemovalTag, randTag;
	bool switch2Fasta, writeUnassigned, writeSingleReads, writeLengthDist;
	bool useStdin, useStdout, relaxRegion, revCompAdapter, qtrimPostRm;
	
	int cutLen_begin, cutLen_end, cutLen_read, a_tail_len, b_tail_len;
	int qtrimThresh, qtrimWinSize;
	int maxUncalled, min_readLen, a_min_overlap, b_min_overlap, nThreads;
	int match, mismatch, gapCost, b_match, b_mismatch, b_gapCost;
	
	float a_threshold, b_threshold;
	
	flexbar::TrimEnd         end, b_end;
	flexbar::FileFormat      format;
	flexbar::QualityType     qual;
	flexbar::QualTrimType    qTrim;
	flexbar::LogLevel        logLevel;
	flexbar::CompressionType cmprsType;
	flexbar::RunType         runType;
	flexbar::BarcodeDetect   barDetect;
	flexbar::AdapterRemoval  adapRm;
	
	tbb::concurrent_vector<TAdapter> barcodes, adapters, barcodes2, adapters2;
	
	std::ostream *out;
	std::fstream fstrmOut;
	
	Options(){
		readsFile      = "";
		readsFile2     = "";
		barReadsFile   = "";
		barcodeFile    = "";
		adapterFile    = "";
		barcode2File   = "";
		adapter2File   = "";
		outCompression = "";
		
		isPaired         = false;
		useAdapterFile   = false;
		useNumberTag     = false;
		useRemovalTag    = false;
		writeUnassigned  = false;
		writeSingleReads = false;
		writeLengthDist  = false;
		switch2Fasta     = false;
		randTag          = false;
		useStdin         = false;
		useStdout        = false;
		relaxRegion      = false;
		revCompAdapter   = false;
		qtrimPostRm      = false;
		
		cutLen_begin  = 0;
		cutLen_end    = 0;
		cutLen_read   = 0;
		qtrimThresh   = 0;
		qtrimWinSize  = 0;
		a_tail_len    = 0;
		b_tail_len    = 0;
		b_min_overlap = 0;
		
		format    = flexbar::FASTA;
		qual      = flexbar::SANGER;
		qTrim     = flexbar::QOFF;
		logLevel  = flexbar::NONE;
		cmprsType = flexbar::UNCOMPRESSED;
		barDetect = flexbar::BOFF;
		adapRm    = flexbar::AOFF;
    }
};


const std::string getFlexbarBanner(const seqan::CharString version){
	
	std::string banner = "";
	
	banner += "               ________          __              \n";
	banner += "              / ____/ /__  _  __/ /_  ____ ______\n";
	banner += "             / /_  / / _ \\| |/ / __ \\/ __ `/ ___/\n";
	banner += "            / __/ / /  __/>  </ /_/ / /_/ / /    \n";
	banner += "           /_/   /_/\\___/_/|_/_.___/\\__,_/_/     \n\n";
	
	banner += "Flexbar - flexible barcode and adapter removal, version ";
	
	append(banner, version);
	
	banner += "\n";
	
	return banner;
}


const std::string getFlexbarCitation(){
	return "\nMatthias Dodt, Johannes T. Roehr, Rina Ahmed, Christoph Dieterich:\nFlexbar - flexible barcode and adapter processing for next-generation\nsequencing platforms. Biology 2012, 1(3):895-905.\n";
}


void defineOptionsAndHelp(seqan::ArgumentParser &parser, const std::string version, const std::string date){
	
	using namespace seqan;
	
	typedef ArgParseArgument ARG;
	
	setVersion(parser, version);
	setDate(parser, date);
	
	setShortDescription(parser, "flexible barcode and adapter removal");
	
	addUsageLine(parser, "\\fB-r\\fP reads [\\fB-t\\fP target] [\\fB-b\\fP barcodes] [\\fB-a\\fP adapters] [options]");
	
	// addOption(parser, ArgParseOption("v", "version", "Display program version."));
	addOption(parser, ArgParseOption("H", "advanced", "Print help with advanced options."));
	addOption(parser, ArgParseOption("M", "man", "Print advanced options as man document."));
	addOption(parser, ArgParseOption("c", "cite", "Show program reference for citation."));
	
	// OUTPUTPREFIX
	addSection(parser, "Basic options");
	addOption(parser, ArgParseOption("n", "threads", "Number of threads to employ.", ARG::INTEGER));
	addOption(parser, ArgParseOption("t", "target", "Prefix for output file names or paths.", ARG::STRING));
	addOption(parser, ArgParseOption("r", "reads", "Fasta/q file or stdin (-) with reads that may contain barcodes.", ARG::INPUTFILE));
	addOption(parser, ArgParseOption("p", "reads2", "Second input file of paired reads, gz and bz2 files supported.", ARG::INPUTFILE));
	
	addSection(parser, "Barcode detection");
	addOption(parser, ArgParseOption("b",  "barcodes", "Fasta file with barcodes for demultiplexing, may contain N.", ARG::INPUTFILE));
	addOption(parser, ArgParseOption("b2", "barcodes2", "Additional barcodes file for second read set in paired mode.", ARG::INPUTFILE));
	addOption(parser, ArgParseOption("br", "barcode-reads", "Fasta/q file containing separate barcode reads for detection.", ARG::INPUTFILE));
	addOption(parser, ArgParseOption("be", "barcode-trim-end", "Type of detection, see section trim-end modes.", ARG::STRING));
	addOption(parser, ArgParseOption("bn", "barcode-tail-length", "Region size in tail trim-end modes. Default: barcode length.", ARG::INTEGER));
	addOption(parser, ArgParseOption("bo", "barcode-min-overlap", "Minimum overlap of barcode and read. Default: barcode length.", ARG::INTEGER));
	addOption(parser, ArgParseOption("bt", "barcode-threshold", "Allowed mismatches and gaps per overlap of 10.", ARG::DOUBLE));
	addOption(parser, ArgParseOption("bk", "barcode-keep", "Keep barcodes within reads instead of removal."));
	addOption(parser, ArgParseOption("bu", "barcode-unassigned", "Include unassigned reads in output generation."));
	addOption(parser, ArgParseOption("bm", "barcode-match", "Alignment match score.", ARG::INTEGER));
	addOption(parser, ArgParseOption("bi", "barcode-mismatch", "Alignment mismatch score.", ARG::INTEGER));
	addOption(parser, ArgParseOption("bg", "barcode-gap", "Alignment gap score.", ARG::INTEGER));
	
	addSection(parser, "Adapter removal");
	addOption(parser, ArgParseOption("a",  "adapters", "Fasta file with adapters for removal that may contain N.", ARG::INPUTFILE));
	addOption(parser, ArgParseOption("a2", "adapters2", "File with extra adapters for second read set in paired mode.", ARG::INPUTFILE));
	addOption(parser, ArgParseOption("as", "adapter-seq", "Single adapter sequence as alternative to adapters option.", ARG::STRING));
	addOption(parser, ArgParseOption("ar", "adapter-read-set", "Consider only single read set for adapters.", ARG::STRING));
	addOption(parser, ArgParseOption("ac", "adapter-revcomp", "Consider also reverse complement of each adapter in search."));
	addOption(parser, ArgParseOption("ae", "adapter-trim-end", "Type of removal, see section trim-end modes.", ARG::STRING));
	addOption(parser, ArgParseOption("an", "adapter-tail-length", "Region size for tail trim-end modes. Default: adapter length.", ARG::INTEGER));
	addOption(parser, ArgParseOption("ad", "adapter-relaxed", "Skip restriction to pass read ends in right and left modes."));
	addOption(parser, ArgParseOption("ao", "adapter-min-overlap", "Minimum overlap of adapter and read for removal.", ARG::INTEGER));
	addOption(parser, ArgParseOption("at", "adapter-threshold", "Allowed mismatches and gaps per overlap of 10.", ARG::DOUBLE));
	addOption(parser, ArgParseOption("am", "adapter-match", "Alignment match score.", ARG::INTEGER));
	addOption(parser, ArgParseOption("ai", "adapter-mismatch", "Alignment mismatch score.", ARG::INTEGER));
	addOption(parser, ArgParseOption("ag", "adapter-gap", "Alignment gap score.", ARG::INTEGER));
	
	// addSection(parser, "Joining paired reads");
	// addOption(parser, ArgParseOption("j",  "join", "Fasta file with adapters for removal that may contain N.", ARG::STRING));
	// addOption(parser, ArgParseOption("jo", "join-min-overlap", "Minimum overlap of adapter and read sequence.", ARG::INTEGER));
	// addOption(parser, ArgParseOption("jn", "join-max-overlap", "Region size for tail trim-end modes. Default: adapter length.", ARG::INTEGER));
	// addOption(parser, ArgParseOption("jt", "join-threshold", "Allowed mismatches and gaps per 10 bases overlap.", ARG::DOUBLE));
	// addOption(parser, ArgParseOption("jm", "join-match", "Alignment match score.", ARG::INTEGER));
	// addOption(parser, ArgParseOption("ji", "join-mismatch", "Alignment mismatch score.", ARG::INTEGER));
	// addOption(parser, ArgParseOption("jg", "join-gap", "Alignment gap score.", ARG::INTEGER));
	
	addSection(parser, "Filtering and trimming");
	addOption(parser, ArgParseOption("u", "max-uncalled", "Allowed uncalled bases (N or .) for each read.", ARG::INTEGER));
	addOption(parser, ArgParseOption("x", "pre-trim-left", "Trim given number of bases on 5' read end before detection.", ARG::INTEGER));
	addOption(parser, ArgParseOption("y", "pre-trim-right", "Trim specified number of bases on 3' end prior to detection.", ARG::INTEGER));
	addOption(parser, ArgParseOption("k", "post-trim-length", "Trim to specified read length from 3' end after removal.", ARG::INTEGER));
	addOption(parser, ArgParseOption("m", "min-read-length", "Minimum read length to remain after removal.", ARG::INTEGER));
	
	addSection(parser, "Quality-based trimming");
	addOption(parser, ArgParseOption("q",  "qtrim", "Quality-based trimming mode.", ARG::STRING));
	addOption(parser, ArgParseOption("qf", "qtrim-format", "Quality format: sanger, solexa, i1.3, i1.5, i1.8 (illumina).", ARG::STRING));
	addOption(parser, ArgParseOption("qt", "qtrim-threshold", "Minimum quality as threshold for trimming.", ARG::INTEGER));
	// addOption(parser, ArgParseOption("qm", "qtrim-win-mean", "Different threshold for min mean quality of window.", ARG::INTEGER));
	addOption(parser, ArgParseOption("qw", "qtrim-win-size", "Region size for sliding window approach.", ARG::INTEGER));
	addOption(parser, ArgParseOption("qa", "qtrim-post-removal", "Perform quality-based trimming after removal steps."));
	
	addSection(parser, "Output selection");
	addOption(parser, ArgParseOption("o", "fasta-output", "Prefer non-quality format fasta for output."));
	addOption(parser, ArgParseOption("z", "zip-output", "Direct compression of output files.", ARG::STRING));
	addOption(parser, ArgParseOption("1", "stdout-reads", "Write reads to stdout, tagged and interleaved if needed."));
	addOption(parser, ArgParseOption("j", "length-dist", "Generate length distribution for read output files."));
	addOption(parser, ArgParseOption("s", "single-reads", "Write single paired reads for too short counterparts."));
	
	addSection(parser, "Logging and tagging");
	addOption(parser, ArgParseOption("l", "log-level", "Print chosen read alignments.", ARG::STRING));
	addOption(parser, ArgParseOption("g", "removal-tags", "Tag reads that are subject to adapter or barcode removal."));
	addOption(parser, ArgParseOption("e", "number-tags", "Replace read tags by ascending number to save space."));
	addOption(parser, ArgParseOption("d", "random-tags", "Capture read sequence at barcode or adapter N positions."));
	
	addSection(parser, "Trim-end modes");
	addText(parser._toolDoc, "ANY: longer side of read remains after removal of overlap", false);
	addText(parser._toolDoc, "LEFT: right side remains after removal, align <= read end", false);
	addText(parser._toolDoc, "RIGHT: left part remains after removal, align >= read start", false);
	addText(parser._toolDoc, "LEFT_TAIL: consider first n bases of reads in alignment", false);
	addText(parser._toolDoc, "RIGHT_TAIL: use only last n bases, see tail-length options", false);
	
	
	hideOption(parser, "barcodes2");
	hideOption(parser, "barcode-tail-length");
	hideOption(parser, "barcode-keep");
	hideOption(parser, "barcode-match");
	hideOption(parser, "barcode-mismatch");
	hideOption(parser, "barcode-gap");
	
	hideOption(parser, "adapters2");
	hideOption(parser, "adapter-revcomp");
	hideOption(parser, "adapter-tail-length");
	hideOption(parser, "adapter-relaxed");
	hideOption(parser, "adapter-read-set");
	hideOption(parser, "adapter-match");
	hideOption(parser, "adapter-mismatch");
	hideOption(parser, "adapter-gap");
	
	hideOption(parser, "qtrim-win-size");
	hideOption(parser, "qtrim-post-removal");
	
	hideOption(parser, "man");
	hideOption(parser, "version");
	hideOption(parser, "stdout-reads");
	hideOption(parser, "length-dist");
	hideOption(parser, "number-tags");
	hideOption(parser, "random-tags");
	
	
	setCategory(parser, "Trimming");
	// setRequired(parser, "reads");
	// setMinValue(parser, "threads", "1");
	// setValidValues(parser, "qtrim-format", "sanger solexa i1.3 i1.5 i1.8");
	
	// setValidValues(parser, "target", "fasta fa fastq fq");
	// setValidValues(parser, "reads", "fasta fa fastq fq");
	// setValidValues(parser, "reads2", "fasta fa fastq fq");
	// setValidValues(parser, "barcode-reads", "fasta fa fastq fq");
	
	// setValidValues(parser, "barcodes", "fasta fa");
	// setValidValues(parser, "barcodes2", "fasta fa");
	// setValidValues(parser, "adapters", "fasta fa");
	// setValidValues(parser, "adapters2", "fasta fa");
	
	// setValidValues(parser, "adapter-trim-end", "ANY LEFT RIGHT LEFT_TAIL RIGHT_TAIL");
	// setMinValue(parser, "adapter-tail-length", "1");
	// setMinValue(parser, "adapter-min-overlap", "1");
	// setMinValue(parser, "adapter-threshold",   "0");
	// setMaxValue(parser, "adapter-threshold",   "10");
	// 
	// setValidValues(parser, "barcode-trim-end", "ANY LEFT RIGHT LEFT_TAIL RIGHT_TAIL");
	// setMinValue(parser, "barcode-tail-length", "1");
	// setMinValue(parser, "barcode-min-overlap", "1");
	// setMinValue(parser, "barcode-threshold",   "0");
	// setMaxValue(parser, "barcode-threshold",   "10");
	// 
	// setMinValue(parser, "max-uncalled",     "0");
	// setMinValue(parser, "pre-trim-left",    "1");
	// setMinValue(parser, "pre-trim-right",   "1");
	// setMinValue(parser, "post-trim-length", "1");
	// setMinValue(parser, "min-read-length",  "1");
	// 
	// setMinValue(parser, "qtrim-threshold",  "0");
	
	
	setValidValues(parser, "qtrim", "TAIL WIN BWA");
	setValidValues(parser, "log-level", "ALL MOD TAB");
	setValidValues(parser, "zip-output", "GZ BZ2");
	setValidValues(parser, "adapter-read-set", "1 2");
	
	setDefaultValue(parser, "target",          "flexbar");
	setDefaultValue(parser, "threads",         "1");
	setDefaultValue(parser, "max-uncalled",    "0");
	setDefaultValue(parser, "min-read-length", "18");
	
	setDefaultValue(parser, "barcode-trim-end",  "ANY");
	setDefaultValue(parser, "barcode-threshold", "1.0");
	setDefaultValue(parser, "barcode-match",     "1");
	setDefaultValue(parser, "barcode-mismatch", "-1");
	setDefaultValue(parser, "barcode-gap",      "-9");
	
	setDefaultValue(parser, "adapter-trim-end",    "RIGHT");
	setDefaultValue(parser, "adapter-min-overlap", "3");
	setDefaultValue(parser, "adapter-threshold",   "3.0");
	setDefaultValue(parser, "adapter-match",       "1");
	setDefaultValue(parser, "adapter-mismatch",   "-1");
	setDefaultValue(parser, "adapter-gap",        "-6");
	
	setDefaultValue(parser, "qtrim-threshold",    "20");
	setDefaultValue(parser, "qtrim-win-size",     "5");
	
	addTextSection(parser, "EXAMPLES");
	addText(parser._toolDoc, "\\fBflexbar\\fP \\fB-r\\fP reads.fq \\fB-t\\fP target \\fB-b\\fP brc.fa \\fB-be\\fP LEFT_TAIL \\fB-a\\fP adp.fa", false);
	addText(parser._toolDoc, "\\fBflexbar\\fP \\fB-r\\fP reads.fq.gz \\fB-q\\fP TAIL \\fB-qf\\fP i1.8 \\fB-a\\fP adp.fa \\fB-ao\\fP 5 \\fB-at\\fP 4");
}


void printLocalTime(Options &o){
	time_t t_current;
	time(&t_current);
	*o.out << "Local time:            " << asctime(localtime(&t_current)) << "\n";
}


void parseCommandLine(seqan::ArgumentParser &parser, std::string version, int argc, char const ** argv){
	
	using namespace std;
	
	using seqan::ArgumentParser;
	
	
	bool useStdout = false;
	
	for (int i=0; i<argc; i++){
		if(strncmp(argv[i], "-1", 2) == 0 || strncmp(argv[i], "--stdout-reads", 14) == 0)
			useStdout = true;
	}
	if(! useStdout) cout << endl;
	
	
	ArgumentParser::ParseResult res = parse(parser, argc, argv);
	
	if(res != ArgumentParser::PARSE_OK){
		
		if(isSet(parser, "help")){
			cerr << "\nAdvanced options: flexbar -H\n";
		}
		cerr << endl;
		
		exit(res == ArgumentParser::PARSE_ERROR);
	}
	
	if(isSet(parser, "cite")){
		cout << getFlexbarBanner(version) << getFlexbarCitation() << endl;
		exit(0);
	}
	
	if(isSet(parser, "advanced") || isSet(parser, "man")){
		
		hideOption(parser, "barcodes2",           false);
		hideOption(parser, "barcode-tail-length", false);
		hideOption(parser, "barcode-keep",        false);
		hideOption(parser, "barcode-match",       false);
		hideOption(parser, "barcode-mismatch",    false);
		hideOption(parser, "barcode-gap",         false);
		
		hideOption(parser, "adapter-revcomp",     false);
		hideOption(parser, "adapter-tail-length", false);
		hideOption(parser, "adapter-relaxed",     false);
		hideOption(parser, "adapter-read-set",    false);
		hideOption(parser, "adapter-match",       false);
		hideOption(parser, "adapter-mismatch",    false);
		hideOption(parser, "adapter-gap",         false);
		
		hideOption(parser, "qtrim-win-size",      false);
		hideOption(parser, "qtrim-post-removal",  false);
		
		hideOption(parser, "adapters2",    false);
		hideOption(parser, "version",      false);
		hideOption(parser, "stdout-reads", false);
		hideOption(parser, "length-dist",  false);
		hideOption(parser, "number-tags",  false);
		hideOption(parser, "random-tags",  false);
		
		if(isSet(parser, "man")) printHelp(parser, cerr, "man");
		else{
			printHelp(parser);
			cerr << "\nFurther documentation on github.com/seqan/flexbar\n" << endl;
		}
		exit(0);
	}
	
	if(! isSet(parser, "reads")){
		printShortHelp(parser);
		cerr << "\nPlease set required input file.\n" << endl;
		exit(1);
	}
}


void loadProgramOptions(Options &o, seqan::ArgumentParser &parser){
	
	using namespace std;
	using namespace flexbar;
	
	ostream *out = o.out;
	
	*out << getFlexbarBanner(getVersion(parser)) << "\n" << endl;
	printLocalTime(o);
	
	
	// basic options
	
	getOptionValue(o.targetName, parser, "target");
	*out << "Target name:           " << o.targetName << endl;
	
	*out << "File type:             ";
	     if(o.format == FASTA)   *out << "fasta";
	else if(o.format == FASTQ)   *out << "fastq";
	*out << endl;
	
	
	getOptionValue(o.readsFile, parser, "reads");
	*out << "Reads file:            ";
	
	if(o.readsFile == "-"){
		*out << "stdin" << endl;
		o.useStdin = true;
	}
	else *out << o.readsFile << endl;
	
	o.runType = SINGLE;
	
	if(isSet(parser, "reads2")){
		getOptionValue(o.readsFile2, parser, "reads2");
		*out << "Reads file 2:          " << o.readsFile2 << "   (paired run)" << endl;
		o.runType  = PAIRED;
		o.isPaired = true;
	}
	
	
	// barcode and adapter file options
	
	if(isSet(parser, "barcodes")){
		
		if(isSet(parser, "barcode-reads")){
			getOptionValue(o.barReadsFile, parser, "barcode-reads");
			*out << "Barcode reads file:    " << o.barReadsFile << endl;
			
			o.barDetect = BARCODE_READ;
		}
		else o.barDetect = WITHIN_READ_REMOVAL;
		
		getOptionValue(o.barcodeFile, parser, "barcodes");
		*out << "Barcode file:          " << o.barcodeFile << endl;
		
		     if(o.runType == SINGLE) o.runType = SINGLE_BARCODED;
		else if(o.runType == PAIRED) o.runType = PAIRED_BARCODED;
		
		if(o.barDetect == WITHIN_READ_REMOVAL && isSet(parser, "barcode-keep")){
			o.barDetect = WITHIN_READ;
		}
		
		if(isSet(parser, "barcodes2") && o.barDetect != BARCODE_READ && o.isPaired){
			getOptionValue(o.barcode2File, parser, "barcodes2");
			*out << "Barcode file 2:        " << o.barcode2File << endl;
			
			if(o.barDetect == WITHIN_READ_REMOVAL) o.barDetect = WITHIN_READ_REMOVAL2;
			else if(o.barDetect == WITHIN_READ)    o.barDetect = WITHIN_READ2;
		}
	}
	
	if(isSet(parser, "adapters")){
		getOptionValue(o.adapterFile, parser, "adapters");
		*out << "Adapter file:          " << o.adapterFile << endl;
		o.adapRm = NORMAL;
		o.useAdapterFile = true;
	}
	else if(isSet(parser, "adapter-seq")){
		getOptionValue(o.adapterSeq, parser, "adapter-seq");
		o.adapRm = NORMAL;
	}
	
	if(isSet(parser, "adapters2") && o.adapRm == NORMAL && o.isPaired){
		getOptionValue(o.adapter2File, parser, "adapters2");
		*out << "Adapter file 2:        " << o.adapter2File << endl;
		o.adapRm = NORMAL2;
	}
	*out << endl;
	
	
	// filtering and trimming options
	
	getOptionValue(o.nThreads, parser, "threads");
	*out << "threads:               " << o.nThreads << endl;
	
	getOptionValue(o.maxUncalled, parser, "max-uncalled");
	*out << "max-uncalled:          " << o.maxUncalled << endl;
	
	if(isSet(parser, "pre-trim-left")){
		getOptionValue(o.cutLen_begin, parser, "pre-trim-left");
		*out << "pre-trim-left:         " << o.cutLen_begin << endl;
	}
	
	if(isSet(parser, "pre-trim-right")){
		getOptionValue(o.cutLen_end, parser, "pre-trim-right");
		*out << "pre-trim-right:        " << o.cutLen_end << endl;
	}
	
	if(isSet(parser, "post-trim-length")){
		getOptionValue(o.cutLen_read, parser, "post-trim-length");
		*out << "post-trim-length:      " << o.cutLen_read << endl;
	}
	
	getOptionValue(o.min_readLen, parser, "min-read-length");
	*out << "min-read-length:       " << o.min_readLen << endl;
	
	
	// quality-based trimming
	
	if(isSet(parser, "qtrim") && o.format == FASTQ){
		
		string qt;
		getOptionValue(qt, parser, "qtrim");
		
		     if(qt == "TAIL") o.qTrim = TAIL;
 		else if(qt == "WIN")  o.qTrim = WIN;
 		else if(qt == "BWA")  o.qTrim = BWA;
		else{
			cerr << "\n\n" << "Specified qtrim mode is unknown!\n" << endl;
			exit(1);
		}
		*out << "qtrim:                 " << qt << endl;
		
		if(isSet(parser, "qtrim-format")){
			
			string quality;
			getOptionValue(quality, parser, "qtrim-format");
			
			     if(quality == "sanger") o.qual = SANGER;
			else if(quality == "solexa") o.qual = SOLEXA;
			else if(quality == "i1.3")   o.qual = ILLUMINA13;
			else if(quality == "i1.5")   o.qual = ILLUMINA13;
			else if(quality == "i1.8")   o.qual = SANGER;
			else{
				cerr << "\n\n" << "Specified quality format is unknown!\n" << endl;
				exit(1);
			}
			*out << "qtrim-format:          " << quality << endl;
		}
		else{
			cerr << "\n\n" << "Specify qtrim-format for quality-based trimming.\n" << endl;
			exit(1);
		}
		
		getOptionValue(o.qtrimThresh, parser, "qtrim-threshold");
		
		if(o.qtrimThresh > 0){
			*out << "qtrim-threshold:       " << o.qtrimThresh;
			
			switch(o.qual){
				case SANGER:      o.qtrimThresh += 33;
					break;
				case SOLEXA:      o.qtrimThresh += 59;
					break;
				case ILLUMINA13:  o.qtrimThresh += 64;
			}
			*out << "  (" << o.qtrimThresh << ")" << endl;
		}
		
		if(o.qTrim == WIN || o.qTrim == WINTAIL){
			
			// if(isSet(parser, "qtrim-win-mean")){
			// 	getOptionValue(o.qtrimWinMean, parser, "qtrim-win-mean");
			// 	*out << "qtrim-win-mean:        " << o.qtrimWinMean << endl;
			// }
			
			getOptionValue(o.qtrimWinSize, parser, "qtrim-win-size");
			*out << "qtrim-win-size:        " << o.qtrimWinSize << endl;
		}
		
		if(isSet(parser, "qtrim-post-removal")) o.qtrimPostRm = true;
	}
	
	
	// output, logging and tagging options
	
	if(isSet(parser, "log-level")){
		getOptionValue(o.logLevelStr, parser, "log-level");
		
		     if(o.logLevelStr == "ALL") o.logLevel = ALL;
		else if(o.logLevelStr == "TAB") o.logLevel = TAB;
		else if(o.logLevelStr == "MOD") o.logLevel = MOD;
	}
	
	if(isSet(parser, "zip-output")){
		getOptionValue(o.outCompression, parser, "zip-output");
		
		if(o.outCompression == "GZ"){
			o.cmprsType = GZ;
			o.outCompression = ".gz";
		}
		else if(o.outCompression == "BZ2"){
			o.cmprsType = BZ2;
			o.outCompression = ".bz2";
		}
	}
	
	if(isSet(parser, "fasta-output")){
		if(o.format == FASTQ){
			o.format = FASTA;
			o.switch2Fasta = true;
		}
	}
	
	if(isSet(parser, "single-reads")) o.writeSingleReads = true;
	if(isSet(parser, "length-dist"))  o.writeLengthDist  = true;
	if(isSet(parser, "number-tags"))  o.useNumberTag     = true;
	if(isSet(parser, "removal-tags")) o.useRemovalTag    = true;
	if(isSet(parser, "random-tags"))  o.randTag          = true;
	
	*out << endl;
	
	
	// barcode options
	
	if(o.barDetect != BOFF){
		
		string b_trim_end;
		getOptionValue(b_trim_end, parser, "barcode-trim-end");
		
		     if(b_trim_end == "LEFT")        o.b_end = LEFT;
		else if(b_trim_end == "RIGHT")       o.b_end = RIGHT;
		else if(b_trim_end == "ANY")         o.b_end = ANY;
		else if(b_trim_end == "LEFT_TAIL")   o.b_end = LEFT_TAIL;
		else if(b_trim_end == "RIGHT_TAIL")  o.b_end = RIGHT_TAIL;
		else{
			cerr << "Specified barcode trim-end is unknown!\n" << endl;
			exit(1);
		}
		*out << "barcode-trim-end:      " << b_trim_end << endl;
		
		
		if(isSet(parser, "barcode-tail-length")){
			getOptionValue(o.b_tail_len, parser, "barcode-tail-length");
			*out << "barcode-tail-length:   " << o.b_tail_len << endl;
		}
		
		if(isSet(parser, "barcode-min-overlap")){
			getOptionValue(o.b_min_overlap, parser, "barcode-min-overlap");
			*out << "barcode-min-overlap:   " << o.b_min_overlap << endl;
		}
		
		getOptionValue(o.b_threshold, parser, "barcode-threshold");
		*out << "barcode-threshold:     " << o.b_threshold << endl;
		
		if(isSet(parser, "barcode-unassigned")) o.writeUnassigned = true;
		
		getOptionValue(o.b_match, parser, "barcode-match");
		getOptionValue(o.b_mismatch, parser, "barcode-mismatch");
		getOptionValue(o.b_gapCost, parser, "barcode-gap");
		
		*out << "barcode-match:        ";
		if(o.b_match >= 0) *out << " ";
		*out << o.b_match << endl;
		
		*out << "barcode-mismatch:     ";
		if(o.b_mismatch >= 0) *out << " ";
		*out << o.b_mismatch << endl;
		
		*out << "barcode-gap:          ";
		if(o.b_gapCost >= 0) *out << " ";
		*out << o.b_gapCost << "\n" << endl;
	}
	
	
	// adapter options
	
	if(o.adapRm != AOFF){
		
		string a_trim_end;
		getOptionValue(a_trim_end, parser, "adapter-trim-end");
		
		if     (a_trim_end == "LEFT")        o.end = LEFT;
		else if(a_trim_end == "RIGHT")       o.end = RIGHT;
		else if(a_trim_end == "ANY")         o.end = ANY;
		else if(a_trim_end == "LEFT_TAIL")   o.end = LEFT_TAIL;
		else if(a_trim_end == "RIGHT_TAIL")  o.end = RIGHT_TAIL;
		else {
			cerr << "Specified adapter trim-end is unknown!\n" << endl;
			exit(1);
		}
		*out << "adapter-trim-end:      " << a_trim_end << endl;
		
		
		if(isSet(parser, "adapter-tail-length")){
			getOptionValue(o.a_tail_len, parser, "adapter-tail-length");
			*out << "adapter-tail-length:   " << o.a_tail_len << endl;
		}
		
		if(isSet(parser, "adapter-revcomp")){
			*out << "adapter-revcomp:       yes" << endl;
			o.revCompAdapter = true;
		}
		
		if(isSet(parser, "adapter-relaxed")){
			*out << "adapter-relaxed:       yes" << endl;
			o.relaxRegion = true;
		}
		
		if(isSet(parser, "adapter-read-set") && o.isPaired && o.adapRm != NORMAL2){
			string a_read_set;
			getOptionValue(a_read_set, parser, "adapter-read-set");
			*out << "adapter-read-set:      " << a_read_set << endl;
			
			     if(a_read_set == "1") o.adapRm = AONE;
			else if(a_read_set == "2") o.adapRm = ATWO;
		}
		
		getOptionValue(o.a_min_overlap, parser, "adapter-min-overlap");
		*out << "adapter-min-overlap:   " << o.a_min_overlap << endl;
		
		getOptionValue(o.a_threshold, parser, "adapter-threshold");
		*out << "adapter-threshold:     " << o.a_threshold << endl;
		
		
		getOptionValue(o.match, parser, "adapter-match");
		getOptionValue(o.mismatch, parser, "adapter-mismatch");
		getOptionValue(o.gapCost, parser, "adapter-gap");
		
		*out << "adapter-match:        ";
		if(o.match >= 0) *out << " ";
		*out << o.match << endl;
		
		*out << "adapter-mismatch:     ";
		if(o.mismatch >= 0) *out << " ";
		*out << o.mismatch << endl;
		
		*out << "adapter-gap:          ";
		if(o.gapCost >= 0) *out << " ";
		*out << o.gapCost << "\n" << endl;
	}
	
	
	// option compatibility tests
	
	if(o.cutLen_read != 0 && o.cutLen_read < o.min_readLen){
		o.cutLen_read = 0;
		cerr << "\nOption post-trim-length omitted, as it is shorter than min read length.\n" << endl;
	}
	
}


#endif
