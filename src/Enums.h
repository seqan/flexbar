/*
 *   Enums.h
 *
 *   Authors: mat and jtr
 */

#ifndef FLEXBAR_ENUMS_H_
#define FLEXBAR_ENUMS_H_


// These enums are used by almost every class.

namespace flexbar{
	
	const unsigned int MAX_READLENGTH = 2048;
	
	enum LogLevel {
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
		FASTQ,
		CSFASTA,
		CSFASTQ
	};
	
	enum QualityType {
		SANGER,
		SOLEXA,
		ILLUMINA13
	};
	
	enum BarcodeDetect {
		BARCODE_READ,
		WITHIN_READ,
		WITHIN_READ_REMOVAL,
		BOFF
	};
	
	enum AdapterRemoval {
		NORMAL,
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

#endif /* FLEXBAR_ENUMS_H_ */
