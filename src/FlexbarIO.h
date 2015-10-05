/*
 *   FlexbarIO.h
 *
 *   Author: jtr
 */

#ifndef FLEXBAR_FLEXBARIO_H_
#define FLEXBAR_FLEXBARIO_H_

#include <string>
#include <fstream>
#include <iostream>

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/stream.h>
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


// void openInputFile(std::istream &strm, std::string path){}
// void closeFile(std::istream &strm){}


#if SEQAN_HAS_ZLIB

void openInputFile(seqan::Stream<seqan::GZFile> &strm, std::string path){
	using namespace std;
	
	// if(path == "-.gz") path = "-";
	
	if(! open(strm, path.c_str(), "rb")){
		cerr << "Error opening gzip file: " << path << "\n" << endl;
		exit(1);
	}
}

void openOutputFile(seqan::Stream<seqan::GZFile> &strm, std::string path){
	using namespace std;
	
	// bool ok;
	// 
	// if(path == "-") ok = open(strm, path.c_str(), "w");
	// else            ok = open(strm, path.c_str(), "wb");
	
	if(! open(strm, path.c_str(), "wb")){
		cerr << "Error opening gzip file: " << path << "\n" << endl;
		exit(1);
	}
}

void closeFile(seqan::Stream<seqan::GZFile> &strm){}

#endif


#if SEQAN_HAS_BZIP2

void openInputFile(seqan::Stream<seqan::BZ2File> &strm, std::string path){
	using namespace std;
	
	// if(path == "-.bz2") path = "-";
	
	if(! open(strm, path.c_str(), "rb")){
		cerr << "Error opening bz2 file: " << path << "\n" << endl;
		exit(1);
	}
}

void openOutputFile(seqan::Stream<seqan::BZ2File> &strm, std::string path){
	using namespace std;
	
	if(! open(strm, path.c_str(), "wb")){
		cerr << "Error opening bz2 file: " << path << "\n" << endl;
		exit(1);
	}
}

void closeFile(seqan::Stream<seqan::BZ2File> &strm){}

#endif


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


template <typename TStream>
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
		TStream fstrm;
		openInputFile(fstrm, path);
		// streamPeek(c, fstrm);
		
		seqan::RecordReader<TStream, seqan::SinglePass<> > reader(fstrm);
		
		if(! atEnd(reader)){
			c = value(reader);
			
			// seqan::CharString text;
			// if(readDigits(text, reader) != 0){
			// 	cerr << "File reading error occured.\n" << endl;
			// 	exit(1); }; cout << text << endl;
		}
		else{
			cerr << "Reads file seems to be empty.\n" << endl;
			exit(1);
		}
		closeFile(fstrm);
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


#endif /* FLEXBAR_FLEXBARIO_H_ */
