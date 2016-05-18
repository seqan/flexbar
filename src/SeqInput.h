/*
 *   SeqInput.h
 *
 *   Authors: mat and jtr
 */

#ifndef FLEXBAR_SEQINPUT_H
#define FLEXBAR_SEQINPUT_H

#include <seqan/seq_io.h>
#include "QualTrimming.h"


template <typename TSeqStr, typename TString>
class SeqInput {

private:
	
	seqan::SeqFileIn seqFileIn;
	
	const flexbar::QualTrimType m_qtrim;
	flexbar::FileFormat m_format;
	
	// typedef seqan::String<char, seqan::MMap<> > TMMapString;
	
	const bool m_preProcess, m_useStdin, m_qtrimPostRm;
	const int m_maxUncalled, m_preTrimBegin, m_preTrimEnd, m_qtrimThresh, m_qtrimWinSize;
	
	tbb::atomic<unsigned long> m_nrReads, m_nrChars, m_nLowPhred;
	
public:
	
	SeqInput(const Options &o, const std::string filePath, const bool preProcess, const bool useStdin) :
		
		m_preProcess(preProcess),
		m_useStdin(useStdin),
		m_maxUncalled(o.maxUncalled),
		m_preTrimBegin(o.cutLen_begin),
		m_preTrimEnd(o.cutLen_end),
		m_qtrim(o.qTrim),
		m_qtrimThresh(o.qtrimThresh),
		m_qtrimWinSize(o.qtrimWinSize),
		m_qtrimPostRm(o.qtrimPostRm),
		m_format(o.format){
		
		m_nrReads   = 0;
		m_nrChars   = 0;
		m_nLowPhred = 0;
		
		using namespace std;
		
		if(m_useStdin){
			if(! open(seqFileIn, cin)){
				cerr << "ERROR: Could not open input stream.\n" << endl;
				exit(1);
			}
		}
		else{
			if(! open(seqFileIn, filePath.c_str())){
				cerr << "ERROR: Could not open file: " << filePath << "\n" << endl;
				exit(1);
			}
		}
	};
	
	virtual ~SeqInput(){
		close(seqFileIn);
	};
	
	
	// returns number of read SeqReads
	unsigned int getSeqReads(seqan::StringSet<bool> &uncalled, flexbar::TStrings &ids, flexbar::TSeqStrs &seqs, flexbar::TStrings &quals, const unsigned int nReads){
		
		using namespace std;
		using namespace flexbar;
		
		using seqan::prefix;
		using seqan::suffix;
		using seqan::length;
		
		try{
			if(! atEnd(seqFileIn)){
				
				reserve(ids,      nReads);
				reserve(seqs,     nReads);
				reserve(uncalled, nReads);
				
				if(m_format == FASTA){
					readRecords(ids, seqs, seqFileIn, nReads);
				}
				else{
					reserve(quals, nReads);
					readRecords(ids, seqs, quals, seqFileIn, nReads);
				}
				
				for(unsigned int i = 0; i < length(ids); ++i){
					
					TString &id  =  ids[i];
					TSeqStr &seq = seqs[i];
					
					if(length(id) < 1){
						cerr << "\n\n" << "ERROR: read without name in input.\n" << endl;
						close(seqFileIn);
						exit(1);
					}
					if(length(seq) < 1){
						cerr << "\n\n" << "ERROR: read without sequence in input.\n" << endl;
						close(seqFileIn);
						exit(1);
					}
					
					m_nrChars += length(seq);
					
					appendValue(uncalled, isUncalledSequence(seq));
					
					if(m_preProcess){
						
						if(m_preTrimBegin > 0 && length(seq) > 1){
							
							int idx = m_preTrimBegin;
							if(idx >= length(seq)) idx = length(seq) - 1;
							
							erase(seq, 0, idx);
							
							if(m_format == FASTQ)
							erase(quals[i], 0, idx);
						}
						if(m_preTrimEnd > 0 && length(seq) > 1){
							
							int idx = m_preTrimEnd;
							if(idx >= length(seq)) idx = length(seq) - 1;
							
							seq = prefix(seq, length(seq) - idx);
							
							if(m_format == FASTQ)
							quals[i] = prefix(quals[i], length(quals[i]) - idx);
						}
						if(m_format == FASTQ && m_qtrim != QOFF && ! m_qtrimPostRm){
							if(qualTrim(seq, quals[i], m_qtrim, m_qtrimThresh, m_qtrimWinSize)) ++m_nLowPhred;
						}
					}
				}
				m_nrReads += length(ids);
				
				return length(ids);
			}
			
			else return 0;  // end of file
		}
		catch(seqan::Exception const &e){
			cerr << "\n\n" << "ERROR: " << e.what() << "\nProgram execution aborted.\n" << endl;
			close(seqFileIn);
			exit(1);
		}
	}
	
	
	// returns TRUE if read contains too many uncalled bases
	bool isUncalledSequence(TSeqStr &seq){
		
		using namespace seqan;
		
		typename Iterator<TSeqStr>::Type it, itEnd;
		
		it    = begin(seq);
		itEnd = end(seq);
		int n = 0;
		
		while(it != itEnd){
			 if(*it == 'N') n++;
			 ++it;
		}
		return(n > m_maxUncalled);
	}
	
	
	unsigned long getNrLowPhredReads() const {
		return m_nLowPhred;
	}
	
	unsigned long getNrProcessedReads() const {
		return m_nrReads;
	}
	
	unsigned long getNrProcessedChars() const {
		return m_nrChars;
	}
	
};

#endif
