/*
 *   FlexbarIO.h
 *
 *   Author: jtr
 */

#ifndef FLEXBAR_FLEXBARIO_H
#define FLEXBAR_FLEXBARIO_H

#include <fstream>
#include <seqan/seq_io.h>
#include <seqan/stream.h>

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
		cerr << "\nERROR: Could not open file " << path << "\n" << endl;
		exit(1);
	}
}

void openOutputFile(std::fstream &strm, std::string path){
	using namespace std;
	
	strm.open(path.c_str(), ios::out | ios::binary);
	
	if(! strm.good()){
		cerr << "\nERROR: Could not open file " << path << "\n" << endl;
		exit(1);
	}
}

void closeFile(std::fstream &strm){
	strm.close();
}



namespace seqan{
	
	// fasta file with dat ending
	
	// Your custom file format
	struct DatFastaAdaptor_;
	using DatFastaAdaptor = Tag<DatFastaAdaptor_>;
	
	// Specilaize sequence input file with custom tag
	using DatFastaSeqFileIn = FormattedFile<Fastq, Input, DatFastaAdaptor>;
	
	// Your custom format tag
	struct DatFastaSeqFormat_;
	using DatFastaSeqFormat = Tag<DatFastaSeqFormat_>;
	
	// The extended TagList containing our custom format
	using DatFastaSeqInFormats = TagList<DatFastaSeqFormat, SeqInFormats>;
	
	// Overloaded file format metafunction
	template <>
	struct FileFormat<FormattedFile<Fastq, Input, DatFastaAdaptor> >{
	    using Type = TagSelector<DatFastaSeqInFormats>;
	};
	
	// Set magic header
	template <typename T>
	struct MagicHeader<DatFastaSeqFormat, T> : public MagicHeader<Fasta, T>{};
	
	// Specify the valid ending for your fasta adaptor
	template <typename T>
	struct FileExtensions<DatFastaSeqFormat, T>{
	    static char const * VALUE[1];
	};
	
	template <typename T>
	char const * FileExtensions<DatFastaSeqFormat, T>::VALUE[1] = { ".dat" };
	
	// Overload an inner readRecord function to delegate to the actual fasta parser
	template <typename TIdString, typename TSeqString, typename TSpec>
	inline void
	readRecord(TIdString & meta, TSeqString & seq, FormattedFile<Fastq, Input, TSpec> & file, DatFastaSeqFormat){
	    readRecord(meta, seq, file.iter, Fasta());  // Just delegate to Fasta parser
	}
	
	// DatFastaSeqFileIn seqFile(path.c_str());
	
	
	// fastq file with dat ending
	
	// Your custom file format
	struct DatFastqAdaptor_;
	using DatFastqAdaptor = Tag<DatFastqAdaptor_>;
	
	// Specilaize sequence input file with custom tag
	using DatFastqSeqFileIn = FormattedFile<Fastq, Input, DatFastqAdaptor>;
	
	// Your custom format tag
	struct DatFastqSeqFormat_;
	using DatFastqSeqFormat = Tag<DatFastqSeqFormat_>;
	
	// The extended TagList containing our custom format
	using DatFastqSeqInFormats = TagList<DatFastqSeqFormat, SeqInFormats>;
	
	// Overloaded file format metafunction
	template <>
	struct FileFormat<FormattedFile<Fastq, Input, DatFastqAdaptor> >{
	    using Type = TagSelector<DatFastqSeqInFormats>;
	};
	
	// Set magic header
	template <typename T>
	struct MagicHeader<DatFastqSeqFormat, T> : public MagicHeader<Fastq, T>{};
	
	// Specify the valid ending for your fasta adaptor
	template <typename T>
	struct FileExtensions<DatFastqSeqFormat, T>{
	    static char const * VALUE[1];
	};
	
	template <typename T>
	char const * FileExtensions<DatFastqSeqFormat, T>::VALUE[1] = { ".dat" };
	
	// Overload an inner readRecord function to delegate to the actual fasta parser
	template <typename TIdString, typename TSeqString, typename TSpec>
	inline void
	readRecord(TIdString & meta, TSeqString & seq, TIdString & qual, FormattedFile<Fastq, Input, TSpec> & file, DatFastqSeqFormat){
	    readRecord(meta, seq, qual, file.iter, Fastq());  // Just delegate to Fastq parser
	}
	
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
				cerr << "\nInput file decompression canceled.\n";
				cerr << "This build does not support zlib.\n" << endl;
				exit(1);
			#endif
		}
		else if(length(path) > 4){
			ending = suffix(path, length(path) - 4);
			
			if(ending == ".bz2"){
				
				#if SEQAN_HAS_BZIP2
					cmprsType = BZ2;
				#else
					cerr << "\nInput file decompression canceled.\n";
		   			cerr << "This build does not support bzip2.\n" << endl;
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
			cerr << "\nERROR: Could not read from standard input stream.\n" << endl;
			exit(1);
		}
		
		     if(c == '>') format = FASTA;
		else if(c == '@') format = FASTQ;
		else{
			cerr << "\nERROR: Reads file type not conform.\n";
			cerr << "Use uncompressed fasta or fastq for stdin.\n" << endl;
			exit(1);
		}
	}
	else{
		
		seqan::SeqFileIn seqFileIn;
		
		if(!open(seqFileIn, path.c_str())){
			cerr << "\nERROR: Could not open file " << path << "\n" << endl;
			exit(1);
		}
		
		if(! atEnd(seqFileIn)){
			
			try{
				FSeqStr rseq;
				FString tag, qual;
				
				readRecord(tag, rseq, qual, seqFileIn);
				
				if(qual == "") format = FASTA;
				else           format = FASTQ;
			}
			catch(seqan::Exception const &e){
				cerr << "\nERROR: " << e.what() << "\nProgram execution aborted.\n" << endl;
				close(seqFileIn);
				exit(1);
			}
		}
		else{
			cerr << "\nReads file seems to be empty.\n\n" << endl;
			close(seqFileIn);
			exit(1);
		}
		
		close(seqFileIn);
	}
}


std::string getExtension(const flexbar::FileFormat format){
	
	using namespace flexbar;
	
	if(format == FASTA) return ".fasta";
	else                return ".fastq";
}


// void runQualityCheck(std::string path){
//
// 	using namespace std;
//
// 	if(! system(NULL)) exit(EXIT_FAILURE);
//
// 	string call = "qcCommand " + path + " &> qc.out";
//
// 	if(system(call.c_str()) != 0){
// 		cerr << "\nERROR: quality control program execution.\n" << endl;
// 	}
// }


#endif
