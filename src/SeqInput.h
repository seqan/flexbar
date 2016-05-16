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
	
	const bool m_switch2Fasta, m_preProcess, m_useStdin, m_qtrimPostRm;
	const int m_maxUncalled, m_preTrimBegin, m_preTrimEnd, m_qtrimThresh, m_qtrimWinSize;
	
	tbb::atomic<unsigned long> m_nrReads, m_nrChars, m_nLowPhred;
	
public:
	
	SeqInput(const Options &o, const std::string filePath, const bool fastaFormat, const bool preProcess, const bool useStdin) :
		
		m_preProcess(preProcess),
		m_useStdin(useStdin),
		m_switch2Fasta(o.switch2Fasta),
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
		
		if(fastaFormat){
			m_format = flexbar::FASTA;
		}
		else if(m_switch2Fasta){
			m_format = flexbar::FASTQ;
		}
		
		if(m_useStdin){
			if(!open(seqFileIn, cin)){
				cerr << "ERROR: Could not open input stream.\n" << endl;
				exit(1);
			}
		}
		else{
			if(!open(seqFileIn, filePath.c_str())){
				cerr << "ERROR: Could not open file: " << filePath << "\n" << endl;
				exit(1);
			}
		}
	};
	
	
	virtual ~SeqInput(){
		close(seqFileIn);
	};
	
	
	unsigned long getNrLowPhredReads() const {
		return m_nLowPhred;
	}
	
	
	unsigned long getNrProcessedReads() const {
		return m_nrReads;
	}
	
	
	unsigned long getNrProcessedChars() const {
		return m_nrChars;
	}
	
	
	// returns SeqRead or NULL if end of file
	void* getRead(bool &isUncalled){
		
		using namespace std;
		using namespace flexbar;
		
		using seqan::prefix;
		using seqan::suffix;
		using seqan::length;
		
		SeqRead<TSeqStr, TString> *seqRead = NULL;
		
		if(! atEnd(seqFileIn)){
			
			isUncalled = false;
			
			try{
				if(m_format == FASTA){
					
					TSeqStr rseq;
					TString tag;
					
					readRecord(tag, rseq, seqFileIn);
					
					if(length(tag) < 1){
						cerr << "\n\n" << "ERROR: Read without name in input.\n" << endl;
						close(seqFileIn);
						exit(1);
					}
					if(length(rseq) < 1){
						cerr << "\n\n" << "ERROR: Read without sequence in input.\n" << endl;
						close(seqFileIn);
						exit(1);
					}
					
					m_nrChars += length(rseq);
					
					
					if(m_preProcess){
						isUncalled = isUncalledSequence(rseq);
						
						if(m_preTrimBegin > 0 && length(rseq) > 1){
							
							int idx = m_preTrimBegin;
							if(idx >= length(rseq)) idx = length(rseq) - 1;
							
							erase(rseq, 0, idx);
						}
						
						if(m_preTrimEnd > 0 && length(rseq) > 1){
							
							int idx = m_preTrimEnd;
							if(idx >= length(rseq)) idx = length(rseq) - 1;
							
							rseq = prefix(rseq, length(rseq) - idx);
						}
					}
					
					seqRead = new SeqRead<TSeqStr, TString>(rseq, tag);
					
					++m_nrReads;
				}
				
				else{  // fastq
					
					TSeqStr rseq;
					TString qual, tag;
					
					readRecord(tag, rseq, qual, seqFileIn);
					
					if(length(tag) < 1){
						cerr << "\n\n" << "ERROR: Read without name in input.\n" << endl;
						close(seqFileIn);
						exit(1);
					}
					if(length(rseq) < 1){
						cerr << "\n\n" << "ERROR: Read without sequence in input.\n" << endl;
						close(seqFileIn);
						exit(1);
					}
					
					m_nrChars += length(rseq);
					
					
					if(m_preProcess){
						isUncalled = isUncalledSequence(rseq);
						
						if(m_preTrimBegin > 0 && length(rseq) > 1){
							
							int idx = m_preTrimBegin;
							if(idx >= length(rseq)) idx = length(rseq) - 1;
							
							erase(rseq, 0, idx);
							erase(qual, 0, idx);
						}
						
						if(m_preTrimEnd > 0 && length(rseq) > 1){
							
							int idx = m_preTrimEnd;
							if(idx >= length(rseq)) idx = length(rseq) - 1;
							
							rseq = prefix(rseq, length(rseq) - idx);
							qual = prefix(qual, length(qual) - idx);
						}
						
						if(m_qtrim != QOFF && ! m_qtrimPostRm){
							
							if(qualTrim(rseq, qual, m_qtrim, m_qtrimThresh, m_qtrimWinSize)) ++m_nLowPhred;
						}
					}
					
					if(m_switch2Fasta) seqRead = new SeqRead<TSeqStr, TString>(rseq, tag);
					else               seqRead = new SeqRead<TSeqStr, TString>(rseq, tag, qual);
					
					++m_nrReads;
				}
				
				return seqRead;
			}
			catch(seqan::Exception const &e){
				cerr << "\n\n" << "ERROR: " << e.what() << "\nProgram execution aborted.\n" << endl;
				close(seqFileIn);
				exit(1);
			}
		}
		
		// end of stream
		else return NULL;
	}
	
	
	// returns TRUE if read contains too many uncalled bases
	bool isUncalledSequence(TSeqStr &seq){
		int n = 0;
		
		using namespace seqan;
		
		typename Iterator<TSeqStr >::Type it, itEnd;
		
		it    = begin(seq);
		itEnd = end(seq);
		
		while(it != itEnd){
			 if(*it == 'N') n++;
			 ++it;
		}
		
		return(n > m_maxUncalled);
	}
 	
};

#endif
