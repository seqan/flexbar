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
	
	SeqRead(const TSeqStr& sequence, const TString& seqTag)
		: seq(sequence),
		  tag(seqTag){
	}
	
	SeqRead(const TSeqStr& sequence, const TString& seqTag, const TString& quality)
		: seq(sequence),
		  tag(seqTag),
		  qual(quality){
	}
	
	virtual ~SeqRead(){}
};


template <typename TSeqStr, typename TString>
class PairedRead {

	public:
	
	typedef SeqRead<TSeqStr, TString> TSeqRead;
	
	TSeqRead *m_r1, *m_r2, *m_b;
	
	TSeqStr m_randTag;
	int m_barcode_id, m_barcode_id2;
	
	PairedRead(TSeqRead *r1, TSeqRead *r2, TSeqRead *b) :
		m_r1(r1),
		m_r2(r2),
		m_b(b),
		m_barcode_id(0),
		m_barcode_id2(0),
		m_randTag(""){
	}
	
	virtual ~PairedRead(){
		delete m_r1;
		delete m_r2;
		delete m_b;
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
	
	typedef std::vector<TAlignments>                    TAlignBundle;
	typedef std::vector<PairedRead<FSeqStr, FString>* > TPairedReadBundle;
	
	
	typedef std::pair< SeqRead<FSeqStr, FString>*,
	                   std::pair< tbb::atomic<unsigned long>, tbb::atomic<unsigned long> > > TAdapter;
	
	
   	enum ComputeCycle {
   		PRECYCLE,
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
