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


void checkFileCompression(std::string path, flexbar::CompressionType &cmprsType){
	
	using namespace std;
	using namespace flexbar;
	
	using seqan::CharString;
	using seqan::suffix;
	using seqan::length;
	
	cmprsType = UNCOMPRESSED;
	
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


void checkInputType(std::string path, flexbar::FileFormat &format){
	
	using namespace std;
	using namespace flexbar;
	
	char c;
	
	if(path == "-"){
		if(cin) c = cin.peek();
		else{
			cerr << "Standard input reading error.\n" << endl;
			exit(1);
		}
	}
	else{
		fstream fstrm;
		
		openInputFile(fstrm, path);
		
		if(fstrm.good()){
			fstrm >> c;
			closeFile(fstrm);
		}
		else{
			cerr << "ERROR: Reads file seems to be empty.\n" << endl;
			closeFile(fstrm);
			exit(1);
		}
	}
	
	     if(c == '>') format = FASTA;
	else if(c == '@') format = FASTQ;
	else{
		cerr << "Reads file type not conform.\n";
		cerr << "Neither fasta nor fastq header.\n" << endl;
		exit(1);
	}
}


std::string toFormatString(flexbar::FileFormat format){
	
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
