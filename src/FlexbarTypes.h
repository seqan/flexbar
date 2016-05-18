/*
 *   FlexbarTypes.h
 *
 */

#ifndef FLEXBAR_FLEXBARTYPES_H
#define FLEXBAR_FLEXBARTYPES_H

#include <vector>
#include <seqan/align.h>


template <typename TSeqStr, typename TString>
class SeqRead {
	
	public:
	TSeqStr seq;
	TString id, qual;
	
	SeqRead(const TSeqStr& sequence, const TString& seqID) :
		seq(sequence),
		id(seqID){
	}
	
	SeqRead(const TSeqStr& sequence, const TString& seqID, const TString& quality) :
		seq(sequence),
		id(seqID),
		qual(quality){
	}
};


template <typename TSeqStr, typename TString>
class PairedRead {
	
	typedef SeqRead<TSeqStr, TString> TSeqRead;
	
	public:
	TSeqRead *r1, *r2, *b;
	unsigned int barID, barID2;
	
	PairedRead(TSeqRead *p_r1, TSeqRead *p_r2, TSeqRead *p_b) :
		r1(p_r1),
		r2(p_r2),
		b(p_b),
		barID(0),
		barID2(0){
	}
	
	virtual ~PairedRead(){
		delete r1;
		delete r2;
		delete b;
	}
};


template <typename TSeqStr>
struct AlignResults{
	
	int score, mismatches, gapsR, gapsA;
	int startPos, startPosA, startPosS;
	int endPos, endPosS, endPosA;
	int overlapLength, queryLength, tailLength;
	
	float allowedErrors;
	TSeqStr randTag;
	seqan::CharString alString;
	
	AlignResults(){
	}
};


namespace flexbar{
	
	const unsigned int MAX_READLENGTH = 2048;
	
	typedef seqan::Dna5String FSeqStr;
	typedef seqan::CharString FString;
	
	typedef seqan::StringSet<FSeqStr> TSeqStrs;
	typedef seqan::StringSet<FString> TStrings;
	
	typedef SeqRead<FSeqStr, FString>    TSeqRead;
	typedef PairedRead<FSeqStr, FString> TPairedRead;
	typedef std::vector<TSeqRead* >      TSeqReads;
	
	typedef seqan::Align<FSeqStr, seqan::ArrayGaps> TAlign;
	typedef seqan::StringSet<TAlign>                TAlignSet;
	typedef seqan::String<int>                      TAlignScores;
	typedef std::pair<TAlignSet, TAlignScores>      TAlignments;
	
	// typedef seqan::StringSet<TAlign, seqan::Dependent<> >             TAlignSet;
	// typedef seqan::StringSet<TAlign, seqan::Dependent<seqan::Tight> > TAlignSet;
	
	typedef std::vector<TAlignments>                    TAlignBundle;
	typedef std::vector<PairedRead<FSeqStr, FString>* > TPairedReadBundle;
	
	
	class TBar {
		
		public:
		FString id;
		FSeqStr seq;
		
		tbb::atomic<unsigned long> rmOverlap, rmFull;
		
		TBar() :
			rmOverlap(0),
			rmFull(0){
	    }
	};
	
	
   	enum ComputeCycle {
   		PRELOAD,
   		COMPUTE,
   		RESULTS
   	};
	
	enum LogAlign {
		NONE,
		ALL,
		TAB,
		MOD
	};
	
	enum CompressionType {
		UNCOMPRESSED,
		GZ,
		BZ2
	};
	
	enum TrimEnd {
		ANY,
		LEFT,
		RIGHT,
		LEFT_TAIL,
		RIGHT_TAIL
	};
	
	enum FileFormat {
		FASTA,
		FASTQ
	};
	
	enum QualityType {
		SANGER,
		SOLEXA,
		ILLUMINA13
	};
	
	enum QualTrimType {
		QOFF,
		TAIL,
		WIN,
		WINTAIL,
		BWA
	};
	
	enum BarcodeDetect {
		BARCODE_READ,
		WITHIN_READ,
		WITHIN_READ_REMOVAL,
		WITHIN_READ2,
		WITHIN_READ_REMOVAL2,
		BOFF
	};
	
	enum AdapterRemoval {
		NORMAL,
		NORMAL2,
		AONE,
		ATWO,
		AOFF
	};
	
	enum RunType {
		SINGLE,
		PAIRED,
		SINGLE_BARCODED,
		PAIRED_BARCODED
	};
}

#endif
