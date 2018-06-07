// Adapters.h

// Oligonucleotide sequences © 2018 Illumina, Inc.  All rights reserved.
// Obtained from https://support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html

#ifndef FLEXBAR_ADAPTERS_H
#define FLEXBAR_ADAPTERS_H

namespace flexbar{
	
	Adapters TrueSeq_ltht = Adapters("TruSeq");
	TrueSeq_ltht.seq1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
	TrueSeq_ltht.seq2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";
	TrueSeq_ltht.info = "TruSeq LT and TruSeq HT-based kits";
	
	Adapters TrueSeq_methyl = Adapters("TrueSeq-Methyl");
	TrueSeq_methyl.seq1 = "AGATCGGAAGAGCACACGTCTGAAC";
	TrueSeq_methyl.seq2 = "AGATCGGAAGAGCGTCGTGTAGGGA";
	TrueSeq_methyl.info = "ScriptSeq and TruSeq DNA Methylation";
	
	Adapters TrueSeq_smallRNA = Adapters("TrueSeq-smallRNA");
	TrueSeq_smallRNA.seq1 = "TGGAATTCTCGGGTGCCAAGG";
	TrueSeq_smallRNA.info = "TruSeq Small RNA";
	
	Adapters TrueSeq_ribo = Adapters("TrueSeq-Ribo");
	TrueSeq_ribo.seq1 = "AGATCGGAAGAGCACACGTCT";
	TrueSeq_ribo.info = "TruSeq Ribo Profile";
	
	Adapters Nextera_TruSight = Adapters("Nextera-TruSight");
	Nextera_TruSight.seq1 = "CTGTCTCTTATACACATCT";
	Nextera_TruSight.seq2 = "CTGTCTCTTATACACATCT";
	Nextera_TruSight.info = "AmpliSeq, Nextera, Nextera DNA Flex, Nextera DNA, Nextera XT, Nextera Enrichment, Nextera Rapid Capture Enrichment, TruSight Enrichment, TruSight Rapid Capture Enrichment, TruSight HLA";
	
	Adapters Nextera_matepair = Adapters("Nextera-Matepair");
	Nextera_matepair.seq1 = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
	Nextera_matepair.seq2 = "GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";
	Nextera_matepair.seqc = "CTGTCTCTTATACACATCT";
	Nextera_matepair.info = "Nextera Mate Pair";
	
}

#endif
