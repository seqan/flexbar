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
	TString tag, qual;
	
	SeqRead(const TSeqStr& sequence, const TString& seqTag) :
		seq(sequence),
		tag(seqTag){
	}
	
	SeqRead(const TSeqStr& sequence, const TString& seqTag, const TString& quality) :
		seq(sequence),
		tag(seqTag),
		qual(quality){
	}
	
	virtual ~SeqRead(){}
};


template <typename TSeqStr, typename TString>
class PairedRead {
	
	public:
	
	typedef SeqRead<TSeqStr, TString> TSeqRead;
	
	TSeqRead *r1, *r2, *b;
	int barID, barID2;
	
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


template <typename TSeqStr, typename TString>
class Bar {
	
	public:
	
	typedef SeqRead<TSeqStr, TString> TSeqRead;
	
	TSeqRead seqRead;
	unsigned long rmOvl, rmFull;
	
	Bar() :
		rmOvl(0),
		rmFull(0){
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
	
	typedef seqan::Align<FSeqStr, seqan::ArrayGaps> TAlign;
	typedef seqan::StringSet<TAlign>                TAlignSet;
	typedef seqan::String<int>                      TAlignScores;
	typedef std::pair<TAlignSet, TAlignScores>      TAlignments;
	
	// typedef seqan::StringSet<TAlign, seqan::Dependent<> >              TAlignSet;
	// typedef seqan::StringSet<TAlign, seqan::Dependent<seqan::Tight> >  TAlignSet;
	
	typedef std::vector<TAlignments>                    TAlignBundle;
	typedef std::vector<PairedRead<FSeqStr, FString>* > TPairedReadBundle;
	
	
	typedef std::pair< SeqRead<FSeqStr, FString>*,
	                   std::pair< tbb::atomic<unsigned long>, tbb::atomic<unsigned long> > > TBar;
	
	
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
