/*
 *   Flexbar.h
 *
 *   Author: jtr
 */

#ifndef FLEXBAR_FLEXBAR_H
#define FLEXBAR_FLEXBAR_H

#include <string>
#include <sstream>
#include <iostream>
#include <vector>

#include <tbb/pipeline.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/concurrent_vector.h>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/align.h>
// #include <seqan/find.h>
// #include <seqan/index.h>
#include <seqan/arg_parse.h>

#include "FlexbarTypes.h"
#include "Options.h"
#include "FlexbarIO.h"
#include "LoadFasta.h"
#include "SeqInput.h"
#include "PairedInput.h"
#include "PairedOutput.h"
#include "PairedAlign.h"


template <typename TSeqStr, typename TString>
void loadBarcodes(Options &o, const bool secondSet){
	
	using namespace std;
	using namespace flexbar;
	
	string barFile = secondSet ? o.barcode2File : o.barcodeFile;
	
	LoadFasta<TSeqStr, TString> lf(o, false);
	
	lf.loadSequences(barFile);
	
	if(secondSet){
		o.barcodes2 = lf.getBars();
		lf.printBars("Barcode2");
		
		if(o.barcodes2.size() == 0){
			cerr << "\nERROR: No barcodes found in file.\n" << endl;
			exit(1);
		}
	}
	else{
		o.barcodes = lf.getBars();
		lf.printBars("Barcode");

		if(o.barcodes.size() == 0){
			cerr << "\nERROR: No barcodes found in file.\n" << endl;
			exit(1);
		}
	}
}


template <typename TSeqStr, typename TString>
void loadAdapters(Options &o, const bool secondSet, const bool useAdapterFile){

	using namespace std;
	using namespace flexbar;
	
	LoadFasta<TSeqStr, TString> lf(o, true);
	
	if(useAdapterFile){
		
		string adapFile = secondSet ? o.adapter2File : o.adapterFile;
		
		lf.loadSequences(adapFile);
		
		if(secondSet){
			o.adapters2 = lf.getBars();
			
			if(o.adapters2.size() == 0){
				cerr << "\nERROR: No adapters found in file.\n" << endl;
				exit(1);
			}
		}
		else{
			o.adapters = lf.getBars();
			
			if(o.adapters.size() == 0){
				cerr << "\nERROR: No adapters found in file.\n" << endl;
				exit(1);
			}
		}
	}
	else{
		TBar bar;
		bar.id  = "cmdline";
		bar.seq = o.adapterSeq;
		o.adapters.push_back(bar);
		
		if(o.revCompAdapter){
			TSeqStr adapterSeqRC = o.adapterSeq;
			seqan::reverseComplement(adapterSeqRC);
			
			TBar barRC;
			barRC.id  = "cmdline revcomp";
			barRC.seq = adapterSeqRC;
			o.adapters.push_back(barRC);
		}
		lf.setBars(o.adapters);
	}
	
	if(secondSet) lf.printBars("Adapter2");
	else          lf.printBars("Adapter");
}


template <typename TSeqStr, typename TString>
void loadBarcodesAndAdapters(Options &o){
	
	using namespace std;
	using namespace flexbar;
	
	if(o.barDetect != BOFF){
		loadBarcodes<TSeqStr, TString>(o, false);
		
		if(o.barDetect == WITHIN_READ2 || o.barDetect == WITHIN_READ_REMOVAL2)
		loadBarcodes<TSeqStr, TString>(o, true);
	}
	
	if(o.adapRm != AOFF){
		loadAdapters<TSeqStr, TString>(o, false, o.useAdapterFile);
		
		if(o.adapRm == NORMAL2)
		loadAdapters<TSeqStr, TString>(o, true, true);
	}
}


void printComputationTime(Options &o, const time_t start){
	
	using namespace std;
	
	time_t end;
	time(&end);
	
	int totalTime = int(difftime(end, start));
	
	int hours   = div(totalTime, 3600).quot;
	int rest    = div(totalTime, 3600).rem;
	int minutes = div(rest, 60).quot;
	int seconds = div(rest, 60).rem;
	
	ostream *out = o.out;
	
	*out << "Computation time:  ";
	if(hours > 0)                               *out << hours     << " h ";
	if(hours > 0 || minutes > 0)                *out << minutes   << " min ";
	if(hours > 0 || minutes > 0 || seconds > 0) *out << seconds   << " sec\n\n\n";
	else                                        *out              << "< 1 sec\n\n\n";
}


std::string alignValue(const int refLength, const unsigned long value){
	using namespace std;
	
	stringstream s; s << value;
	
	int wSpaceLen = refLength - s.str().length();
	if(wSpaceLen < 0) wSpaceLen = 0;
	
	return string(wSpaceLen, ' ') + s.str();
}


void printMessage(Options &o){
	
	using namespace std;
	using namespace flexbar;
	
	string s = "Flexbar completed ";
	
	if(o.barDetect != BOFF)                      s += "barcode";
	if(o.barDetect == WITHIN_READ_REMOVAL)       s += " removal within reads";
	if(o.barDetect == WITHIN_READ)               s += " detection within reads";
	if(o.barDetect == BARCODE_READ)              s += " detection with separate reads";
	if(o.barDetect != BOFF && o.adapRm != AOFF)  s += " and ";
	if(o.barDetect == BOFF && o.adapRm == AOFF)  s += "basic processing";
	if(o.adapRm    != AOFF)                      s += "adapter removal";
	
	*o.out << s << ".\n" << endl;
	
	if(! o.logStdout) closeFile(o.fstrmOut);
}


template <typename TSeqStr, typename TString>
void startProcessing(Options &o){
	
	using namespace std;
	using namespace flexbar;
	
	time_t start;
	time(&start);
	
	ostream *out = o.out;
	
	*out << "\nProcessing reads ..." << flush;
	
	if(o.logAlign != NONE) *out << "\n\nAlignment " << o.logAlignStr << " logging:\n\n" << endl;
	
	PairedInput<TSeqStr, TString>  inputFilter(o);
	PairedAlign<TSeqStr, TString>  alignFilter(o);
	PairedOutput<TSeqStr, TString> outputFilter(o);
	
	tbb::task_scheduler_init init_serial(o.nThreads);
	tbb::pipeline pipe;
	
	pipe.add_filter(inputFilter);
	pipe.add_filter(alignFilter);
	pipe.add_filter(outputFilter);
	pipe.run(o.nThreads);
	
	if(o.logAlign == TAB) *out << "\n";
	*out << "done.\n" << endl;
	
	printComputationTime(o, start);
	
	
	// barcode and adapter removal statistics
	
	if(o.writeLengthDist) outputFilter.writeLengthDist();
	
	if(o.adapRm != AOFF){
		outputFilter.printAdapterRemovalStats();
		alignFilter.printAdapterOverlapStats();
		
		if(o.adapRm == NORMAL2){
			outputFilter.printAdapterRemovalStats2();
			alignFilter.printAdapterOverlapStats2();
		}
	}
	
	outputFilter.printFileSummary();
	
	
	// summary statistics of filtering
	
	const unsigned long nReads   = inputFilter.getNrProcessedReads();
	const unsigned long nChars   = inputFilter.getNrProcessedChars();
	const unsigned long uncalled = inputFilter.getNrUncalledReads();
	const unsigned long uPairs   = inputFilter.getNrUncalledPairedReads();
	
	unsigned long nGoodReads     = outputFilter.getNrGoodReads();
	unsigned long nGoodChars     = outputFilter.getNrGoodChars();
	
	if(o.isPaired && o.writeSingleReadsP){
		nGoodReads -= outputFilter.getNrSingleReads();
		nGoodChars -= outputFilter.getNrSingleReads();
	}
	
	stringstream s; s << nReads;
	int len = s.str().length();
	
	*out << "Filtering statistics\n";
	*out << "====================\n";
	*out << "Processed reads                   " << nReads << endl;
	*out << "  skipped due to uncalled bases   ";
	
	if(o.isPaired){
		*out << alignValue(len, 2 * uPairs);
		
		if(uncalled > 0)
		*out << "   (" << uncalled << " uncalled in " << uPairs << " pairs)";
		*out << endl;
	}
	else *out << alignValue(len, uncalled) << endl;
	
	if(o.qTrim != QOFF && ! o.qtrimPostRm)
	*out << "  trimmed due to low quality      " << alignValue(len, inputFilter.getNrLowPhredReads()) << endl;
	
	if(o.barDetect != BOFF && ! o.writeUnassigned)
	*out << "  skipped unassigned reads        " << alignValue(len, alignFilter.getNrUnassignedReads()) << endl;
	
	if(o.adapRm != AOFF)
	*out << "  short prior to adapter removal  " << alignValue(len, alignFilter.getNrPreShortReads()) << endl;
	
	if(o.qTrim != QOFF && o.qtrimPostRm)
	*out << "  trimmed due to low quality      " << alignValue(len, outputFilter.getNrLowPhredReads()) << endl;
	
	*out << "  finally skipped short reads     " << alignValue(len, outputFilter.getNrShortReads()) << endl;
	
	if(o.isPaired && ! o.writeSingleReads && ! o.writeSingleReadsP)
	*out << "  skipped paired single reads     " << alignValue(len, outputFilter.getNrSingleReads()) << endl;
	
	*out << "Discarded reads overall           " << alignValue(len, nReads - nGoodReads) << endl;
	*out << "Remaining reads                   " << alignValue(len, nGoodReads);
	
	if(nReads > 0)
	*out << "   (" << fixed << setprecision(2) << 100 * nGoodReads / nReads << "% of input)";
	
	stringstream schar; schar << inputFilter.getNrProcessedChars();
	int clen = schar.str().length();
	
	*out << "\n" << endl;
	
	*out << "Processed bases   " << alignValue(clen, nChars) << endl;
	*out << "Remaining bases   " << alignValue(clen, nGoodChars);
	
	if(nChars > 0)
	*out << "   (" << fixed << setprecision(2) << 100 * nGoodChars / nChars << "% of input)";
	
	*out << "\n\n" << endl;
	
	printMessage(o);
}


void performTest(){
	
	using namespace std;
	using namespace seqan;
	
}


void startComputation(Options &o){
	
	// performTest();
	
	using namespace std;
	using namespace flexbar;
	
	loadBarcodesAndAdapters<FSeqStr, FString>(o);
	
	if(o.cmprsType == GZ){
		
		#if SEQAN_HAS_ZLIB
			startProcessing<FSeqStr, FString>(o);
		#else
			o.outCompression = "";
			o.cmprsType = UNCOMPRESSED;
			cerr << "Output file compression inactive.\n"
			     << "This build does not support zlib!\n" << endl;
		#endif
	}
	
	else if(o.cmprsType == BZ2){
		
		#if SEQAN_HAS_BZIP2
			startProcessing<FSeqStr, FString>(o);
		#else
			o.outCompression = "";
			o.cmprsType = UNCOMPRESSED;
			cerr << "Output file compression inactive.\n"
			     << "This build does not support bzip2!\n" << endl;
		#endif
	}
	
	if(o.cmprsType == UNCOMPRESSED){
		startProcessing<FSeqStr, FString>(o);
	}
}


#endif
