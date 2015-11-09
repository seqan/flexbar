/*
 *   Flexbar.h
 *
 *   Author: jtr
 */

#ifndef FLEXBAR_FLEXBAR_H
#define FLEXBAR_FLEXBAR_H

#include <string>
#include <fstream>
#include <iostream>

#include <tbb/pipeline.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/concurrent_vector.h>

#include <seqan/basic.h>
#include <seqan/stream.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>

#include "Enums.h"
#include "Options.h"
#include "FlexbarIO.h"
#include "AdapterLoader.h"
#include "SeqRead.h"
#include "SeqInputFilter.h"
#include "PairedInputFilter.h"
#include "PairedOutputFilter.h"
#include "PairedAlignmentFilter.h"


void loadBarcodes(Options &o, const bool secondSet){
	
	using namespace std;
	using namespace flexbar;
	
	using seqan::CharString;
	
	if(o.barDetect != BOFF){
		tbb::task_scheduler_init init_serial(1);
		tbb::pipeline bpipeline;
		
		string barFile = secondSet ? o.barcode2File : o.barcodeFile;
		
		SeqInputFilter<CharString, CharString, fstream> adapter_filter(o, barFile, true, false, false);
		bpipeline.add_filter(adapter_filter);
		
		AdapterLoader<CharString, CharString> adapterLoader(o, false);
		bpipeline.add_filter(adapterLoader);
		bpipeline.run(1);
		
		if(secondSet){
			o.barcodes2 = adapterLoader.getAdapters();
			adapterLoader.printAdapters("Barcode2");

			if(o.barcodes2.size() == 0){
				cerr << "No barcodes found in file!\n" << endl;
				exit(1);
			}
		}
		else{
			o.barcodes = adapterLoader.getAdapters();
			adapterLoader.printAdapters("Barcode");

			if(o.barcodes.size() == 0){
				cerr << "No barcodes found in file!\n" << endl;
				exit(1);
			}
		}
	}
}

	
void loadAdapters(Options &o, const bool secondSet, const bool useAdapterFile){
	
	using namespace std;
	using namespace flexbar;
	
	using seqan::CharString;
	
	if(o.adapRm != AOFF){
		
		AdapterLoader<CharString, CharString> adapterLoader(o, true);
		
		if(useAdapterFile){
			tbb::task_scheduler_init init_serial(1);
			tbb::pipeline prepipe;
			
			string adapFile = secondSet ? o.adapter2File : o.adapterFile;
			
			SeqInputFilter<CharString, CharString, fstream> adapter_filter(o, adapFile, true, false, false);
			prepipe.add_filter(adapter_filter);
			prepipe.add_filter(adapterLoader);
			prepipe.run(1);
			
			if(secondSet){
				o.adapters2 = adapterLoader.getAdapters();

				if(o.adapters2.size() == 0){
					cerr << "No adapters found in file!\n" << endl;
					exit(1);
				}
			}
			else{
				o.adapters = adapterLoader.getAdapters();

				if(o.adapters.size() == 0){
					cerr << "No adapters found in file!\n" << endl;
					exit(1);
				}
			}
		}
		else{
			CharString adapterSeq = o.adapterSeq;
			
			SeqRead<CharString, CharString> *myRead;
			myRead = new SeqRead<CharString, CharString>(adapterSeq, "cmdline");
			
			TAdapter adap;
			adap.first = myRead;
			o.adapters.push_back(adap);
			
			if(o.revCompAdapter){
				CharString adapterSeqRC = o.adapterSeq;
				seqan::reverseComplement(adapterSeqRC);
				
				SeqRead<CharString, CharString> *myReadRC;
				myReadRC = new SeqRead<CharString, CharString>(adapterSeqRC, "cmdline revcomp");
				
				TAdapter adapRC;
				adapRC.first = myReadRC;
				o.adapters.push_back(adapRC);
			}
			
			adapterLoader.setAdapters(o.adapters);
		}
		
		if(secondSet) adapterLoader.printAdapters("Adapter2");
		else          adapterLoader.printAdapters("Adapter");
	}
}


void loadBarcodesAndAdapters(Options &o){
	
	using namespace std;
	using namespace flexbar;
	
	loadBarcodes(o, false);
	
	if(o.barDetect == WITHIN_READ2 || o.barDetect == WITHIN_READ_REMOVAL2)
	loadBarcodes(o, true);
	
	loadAdapters(o, false, o.useAdapterFile);
	
	if(o.adapRm == NORMAL2)
	loadAdapters(o, true, true);
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


void printCompletedMessage(Options &o){
	
	using namespace std;
	using namespace flexbar;
	
	stringstream s;
	s << "Flexbar completed ";
	
	if(o.barDetect != BOFF)                      s << "barcode";
	if(o.barDetect == WITHIN_READ_REMOVAL)       s << " removal within reads";
	if(o.barDetect == WITHIN_READ)               s << " detection within reads";
	if(o.barDetect == BARCODE_READ)              s << " detection with separate reads";
	if(o.barDetect != BOFF && o.adapRm != AOFF)  s << " and ";
	if(o.barDetect == BOFF && o.adapRm == AOFF)  s << "basic processing";
	if(o.adapRm    != AOFF)                      s << "adapter removal";
	
	*o.out << s.str() << ".\n" << endl;
	
	if(o.useStdout) closeFile(o.fstrmOut);
}


template <typename TStreamR, typename TStreamP, typename TStreamB, typename TStreamOut>
void startProcessing(Options &o){
	
	using namespace std;
	using namespace flexbar;
	
	// Check TSeqStr everywhere, e.g. qual
	// typedef seqan::Dna5String TSeqStr;
	
	typedef seqan::CharString TSeqStr;
	typedef seqan::CharString TString;
	
	time_t start;
	time(&start);
	
	ostream *out = o.out;
	
	*out << "\nProcessing reads ..." << flush;
	
	if(o.logLevel != NONE) *out << "\n\nLog level " << o.logLevelStr << " output generation:\n\n" << endl;
	
	PairedInputFilter<TSeqStr, TString, TStreamR, TStreamP, TStreamB> inputFilter(o);
	PairedAlignmentFilter<TSeqStr, TString> alignFilter(o);
	PairedOutputFilter<TSeqStr, TString, TStreamOut> outputFilter(o);
	
	tbb::task_scheduler_init init_serial(o.nThreads);
	tbb::pipeline pipe;
	
	pipe.add_filter(inputFilter);
	pipe.add_filter(alignFilter);
	pipe.add_filter(outputFilter);
	pipe.run(o.nThreads);
	
	if(o.logLevel == TAB) *out << "\n";
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
	
	
	const unsigned long nReads     = inputFilter.getNrProcessedReads();
	const unsigned long nChars     = inputFilter.getNrProcessedChars();
	const unsigned long uncalled   = inputFilter.getNrUncalledReads();
	const unsigned long uPairs     = inputFilter.getNrUncalledPairedReads();
	const unsigned long nGoodReads = outputFilter.getNrGoodReads();
	const unsigned long nGoodChars = outputFilter.getNrGoodChars();
	
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
	*out << "  short prior adapter removal     " << alignValue(len, alignFilter.getNrPreShortReads()) << endl;
	
	if(o.qTrim != QOFF && o.qtrimPostRm)
	*out << "  trimmed due to low quality      " << alignValue(len, outputFilter.getNrLowPhredReads()) << endl;
	
	*out << "  finally skipped short reads     " << alignValue(len, outputFilter.getNrShortReads()) << endl;
	
	if(o.isPaired && ! o.writeSingleReads)
	*out << "  skipped single paired reads     " << alignValue(len, outputFilter.getNrSingleReads()) << endl;
	
	*out << "Discarded reads overall           " << alignValue(len, nReads - nGoodReads) << endl;
	*out << "Remaining reads                   " << alignValue(len, nGoodReads);
	
	if(nReads > 0)
	*out << "   (" << fixed << setprecision(2) << 100 * nGoodReads / nReads << "% of input)";
	
	stringstream schar; schar << inputFilter.getNrProcessedChars();
	int clen = schar.str().length();
	
	*out << "\n" << endl;
	
	*out << "Processed bases:   " << alignValue(clen, nChars) << endl;
	*out << "Remaining bases:   " << alignValue(clen, nGoodChars);
	
	if(nChars > 0)
	*out << "   (" << fixed << setprecision(2) << 100 * nGoodChars / nChars << "% of input)";
	
	*out << "\n\n" << endl;
}


template <typename TStreamR, typename TStreamP, typename TStreamB>
void startProcessing(Options &o){
	
	using namespace std;
	using namespace flexbar;
	
	if(o.cmprsType == GZ){
		
		#if SEQAN_HAS_ZLIB
			startProcessing<TStreamR, TStreamP, TStreamB, seqan::Stream<seqan::GZFile> >(o);
		#else
			o.outCompression = "";
			o.cmprsType = UNCOMPRESSED;
			cerr << "Output file compression inactive.\n"
			     << "This build does not support zlib!\n" << endl;
		#endif
	}
	
	else if(o.cmprsType == BZ2){
		
		#if SEQAN_HAS_BZIP2
			startProcessing<TStreamR, TStreamP, TStreamB, seqan::Stream<seqan::BZ2File> >(o);
		#else
			o.outCompression = "";
			o.cmprsType = UNCOMPRESSED;
			cerr << "Output file compression inactive.\n"
			     << "This build does not support bzip2!\n" << endl;
		#endif
	}
	
	if(o.cmprsType == UNCOMPRESSED){
		startProcessing<TStreamR, TStreamP, TStreamB, fstream>(o);
	}
}


template <typename TStreamR, typename TStreamP>
void startProcessing(Options &o){
	
	using namespace flexbar;
	
	CompressionType cmprsType;
	checkFileCompression(o.barReadsFile, cmprsType);
	
	#if SEQAN_HAS_ZLIB
		if(cmprsType == GZ){
			startProcessing<TStreamR, TStreamP, seqan::Stream<seqan::GZFile> >(o);
		}
	#endif
	
	#if SEQAN_HAS_BZIP2
		if(cmprsType == BZ2){
			startProcessing<TStreamR, TStreamP, seqan::Stream<seqan::BZ2File> >(o);
		}
	#endif
	
	if(cmprsType == UNCOMPRESSED){
		startProcessing<TStreamR, TStreamP, std::fstream>(o);
	}
}


template <typename TStreamR>
void startProcessing(Options &o){
	
	using namespace flexbar;
	
	CompressionType cmprsType;
	checkFileCompression(o.readsFile2, cmprsType);
	
	#if SEQAN_HAS_ZLIB
		if(cmprsType == GZ){
			startProcessing<TStreamR, seqan::Stream<seqan::GZFile> >(o);
		}
	#endif
	
	#if SEQAN_HAS_BZIP2
		if(cmprsType == BZ2){
			startProcessing<TStreamR, seqan::Stream<seqan::BZ2File> >(o);
		}
	#endif
	
	if(cmprsType == UNCOMPRESSED){
		startProcessing<TStreamR, std::fstream>(o);
	}
}


void startProcessing(Options &o, const bool start){
	
	using namespace flexbar;
	
	CompressionType cmprsType;
	checkFileCompression(o.readsFile, cmprsType);
	
	#if SEQAN_HAS_ZLIB
		if(cmprsType == GZ){
			if(start) startProcessing<seqan::Stream<seqan::GZFile> >(o);
			else      checkInputType<seqan::Stream<seqan::GZFile> >(o.readsFile, o.format);
		}
	#endif
	
	#if SEQAN_HAS_BZIP2
		if(cmprsType == BZ2){
			if(start) startProcessing<seqan::Stream<seqan::BZ2File> >(o);
		  	else      checkInputType<seqan::Stream<seqan::BZ2File> >(o.readsFile, o.format);
		}
	#endif
	
	if(cmprsType == UNCOMPRESSED){
		if(start) startProcessing<std::fstream>(o);
		else      checkInputType<std::fstream>(o.readsFile, o.format);
	}
}


void initOptions(Options &o, seqan::ArgumentParser &parser){
	
	using namespace std;
	
	if(isSet(parser, "stdout-reads")){
		
		string s;
		getOptionValue(s, parser, "target");
		openOutputFile(o.fstrmOut, s + ".out");
		
		o.out = &o.fstrmOut;
		o.useStdout = true;
		*o.out << endl;
	}
	else{
		o.out = &cout;
	}
	
	getOptionValue(o.readsFile, parser, "reads");
	startProcessing(o, false);
}


// #include <seqan/find.h>
// #include <seqan/align.h>

void performTest(){
	
	using namespace std;
	using namespace flexbar;
	
	using seqan::CharString;
	
	// typedef seqan::String<char, seqan::MMap<> > TMMapString;
	// TMMapString mmapStr;
	// 
	// if(! open(mmapStr, "test/test.fasta", seqan::OPEN_RDONLY)){
	// 	cout << "Error opening File." << std::endl;
	// 	exit(1);
	// }
	// seqan::RecordReader<TMMapString, seqan::SinglePass<seqan::StringReader> > mmReader(mmapStr);
	// string text2 = "";
	// readLine(text2, mmReader);
	// cout << text2 << endl;
	
	// CharString haystack = "ATGGATTGCG", needle = "ATGCAT";
	// 
	// seqan::Finder<CharString> finder(haystack);
	// seqan::Pattern<CharString, seqan::DPSearch<seqan::SimpleScore, seqan::FindInfix> > pattern(needle, seqan::SimpleScore(0, -1, -7));
	// 
	// while (find(finder, pattern, -2)){
	//     while (findBegin(finder, pattern, getScore(pattern))){
	//         cout << '[' << beginPosition(finder) << ',' << endPosition(finder) << ")\t" << infix(finder) << endl;
	//        // cout << end(finder) << endl; //',' << position(pattern) << endl;
	// 	} }
	// clear(finder);
	// seqan::Pattern<CharString, seqan::AbndmAlgo > pattern2(needle, -2);
	// 
	// //seqan::Score<int> sc(0,-3,-2);  // = scoringScheme(pattern2);
	// //setScoringScheme(pattern2, sc);
	// 
	// while (find(finder, pattern2)){
	//     while (findBegin(finder, pattern2, getScore(pattern2))){
	//         cout << '[' << beginPosition(finder) << ',' << endPosition(finder) << ")\t" << infix(finder) << endl;
	// 	}
	// }
}


void startComputation(Options &o){
	
	using namespace std;
	
	// performTest();
	
	startProcessing(o, true);
}


#endif
