/*
 *   FlexbarIO.h
 *
 *   Author: jtr
 */

#ifndef FLEXBAR_FLEXBARIO_H
#define FLEXBAR_FLEXBARIO_H

#include <string>
#include <fstream>
#include <iostream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include "Enums.h"

#if SEQAN_HAS_ZLIB
#include <zlib.h>
#endif

#if SEQAN_HAS_BZIP2
#include <bzlib.h>
#endif


void openInputFile(std::fstream &strm, std::string path){
	using namespace std;
	
	strm.open(path.c_str(), ios::in | ios::binary);
	
	if(! strm.good()){
		cerr << "Error opening file: " << path << "\n" << endl;
		exit(1);
	}
}

void openOutputFile(std::fstream &strm, std::string path){
	using namespace std;
	
	strm.open(path.c_str(), ios::out | ios::binary);
	
	if(! strm.good()){
		cerr << "Error opening file: " << path << "\n" << endl;
		exit(1);
	}
}

void closeFile(std::fstream &strm){
	strm.close();
}


void checkFileCompression(const std::string path){
	
	using namespace std;
	using namespace flexbar;
	
	using seqan::CharString;
	using seqan::suffix;
	using seqan::length;
	
	CompressionType cmprsType = UNCOMPRESSED;
	
	if(length(path) > 3){
		CharString ending = suffix(path, length(path) - 3);
		
		if(ending == ".gz"){
			
			#if SEQAN_HAS_ZLIB
				cmprsType = GZ;
			#else
				cerr << "Input file decompression canceled.\n";
				cerr << "This build does not support zlib!\n" << endl;
				exit(1);
			#endif
		}
		else if(length(path) > 4){
			ending = suffix(path, length(path) - 4);
			
			if(ending == ".bz2"){
				
				#if SEQAN_HAS_BZIP2
					cmprsType = BZ2;
				#else
					cerr << "Input file decompression canceled.\n";
		   			cerr << "This build does not support bzip2!\n" << endl;
					exit(1);
				#endif
			}
		}
	}
}


void checkInputType(const std::string path, flexbar::FileFormat &format, const bool isReadsFile){
	
	using namespace std;
	using namespace flexbar;
	
	checkFileCompression(path);
	
	if(path == "-" && isReadsFile){
		
		char c;
		
		if(cin) c = cin.peek();
		else{
			cerr << "Standard input reading error.\n" << endl;
			exit(1);
		}
		
		     if(c == '>') format = FASTA;
		else if(c == '@') format = FASTQ;
		else{
			cerr << "Reads file type not conform.\n";
			cerr << "Uncompressed fasta or fastq for stdin.\n" << endl;
			exit(1);
		}
	}
	else{
		
		seqan::SeqFileIn seqFileIn;
		
		if(!open(seqFileIn, path.c_str())){
			cerr << "ERROR: Could not open file: " << path << "\n" << endl;
			exit(1);
		}
		
		if(! atEnd(seqFileIn)){
			
			try{
				seqan::Dna5String rseq;
				seqan::CharString tag, qual;
				
				readRecord(tag, rseq, qual, seqFileIn);
				
				if(qual == ""){
					format = FASTA;
				}
				else{
					format = FASTQ;
				}
			}
			catch(seqan::Exception const &e){
				cerr << "\n\n" << "ERROR: " << e.what() << "\nProgram execution aborted.\n" << endl;
				close(seqFileIn);
				exit(1);
			}
		}
		else{
			cerr << "Reads file seems to be empty.\n\n" << endl;
			close(seqFileIn);
			exit(1);
		}
		
		close(seqFileIn);
	}
}


std::string toFormatString(const flexbar::FileFormat format){
	
	using namespace flexbar;
	
	switch(format){
		case FASTA:   return ".fasta";
		case FASTQ:   return ".fastq";
	}
	return ".unknown";
}


void runQualityCheck(std::string path){
	
	using namespace std;
	
	if(! system(NULL)) exit(EXIT_FAILURE);
	
	string call = "qcCommand " + path + " &> qc.out";
	
	if(system(call.c_str()) != 0){
		cerr << "Error in quality control.\n" << endl;
	}
}


#endif
