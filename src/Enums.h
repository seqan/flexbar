/*
 *   Enums.h
 *
 *   Authors: mat and jtr
 */

#ifndef FLEXBAR_ENUMS_H
#define FLEXBAR_ENUMS_H


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
		FASTQ
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
