// PairedInput.h

#ifndef FLEXBAR_PAIREDINPUT_H
#define FLEXBAR_PAIREDINPUT_H

#include "SeqInput.h"


template <typename TSeqStr, typename TString>
class PairedInput : public tbb::filter {

private:
	
	const flexbar::FileFormat m_format;
	const bool m_isPaired, m_useBarRead, m_useNumberTag;
	const unsigned int m_bundleSize;
	
	tbb::atomic<unsigned long> m_uncalled, m_uncalledPairs, m_tagCounter;
	SeqInput<TSeqStr, TString> *m_f1, *m_f2, *m_b;
	
public:
	
	PairedInput(const Options &o) :
		
		filter(serial_in_order),
		m_format(o.format),
		m_useNumberTag(o.useNumberTag),
		m_isPaired(o.isPaired),
		m_useBarRead(o.barDetect == flexbar::BARCODE_READ),
		m_bundleSize(o.bundleSize),
		m_tagCounter(0),
		m_uncalled(0),
		m_uncalledPairs(0){
		
		m_f1 = new SeqInput<TSeqStr, TString>(o, o.readsFile, true, o.useStdin);
		
		m_f2 = NULL;
		m_b  = NULL;
		
		if(m_isPaired)
		m_f2 = new SeqInput<TSeqStr, TString>(o, o.readsFile2, true, false);
		
		if(m_useBarRead)
		m_b = new SeqInput<TSeqStr, TString>(o, o.barReadsFile, false, false);
	}
	
	virtual ~PairedInput(){
		delete m_f1;
		delete m_f2;
		delete m_b;
	}
	
	
	void* loadPairedReadBundle(){
		
		using namespace std;
		using namespace flexbar;
		
		PairedReadBundle *b = new PairedReadBundle();
		
		unsigned int nReads = m_f1->loadSeqReads(b->srd.uncalled, b->srd.ids, b->srd.seqs, b->srd.quals, m_bundleSize);
		
		if(m_isPaired){
			unsigned int nReads2 = m_f2->loadSeqReads(b->srd2.uncalled, b->srd2.ids, b->srd2.seqs, b->srd2.quals, m_bundleSize);
			
			if(nReads != nReads2){
				cerr << "\nERROR: Single read in paired input mode.\n" << endl;
				exit(1);
			}
		}
		if(m_useBarRead){
			unsigned int nBarReads = m_b->loadSeqReads(b->srdBR.uncalled, b->srdBR.ids, b->srdBR.seqs, b->srdBR.quals, m_bundleSize);
			
			if(nReads < nBarReads){
				cerr << "\nERROR: Barcode read without read in input.\n" << endl;
				exit(1);
			}
			else if(nReads > nBarReads){
				cerr << "\nERROR: Read without barcode read in input.\n" << endl;
				exit(1);
			}
		}
		
		if(nReads == 0){
			delete b;
			return NULL;
		}
		
		unsigned int nEntries = 0;
		
		for(unsigned int i = 0; i < length(b->srd.ids); ++i){
			
			if(b->srd.uncalled[i] || (m_isPaired && b->srd2.uncalled[i])){
				
				if(b->srd.uncalled[i])                ++m_uncalled;
				if(m_isPaired && b->srd2.uncalled[i]) ++m_uncalled;
				if(m_isPaired)                        ++m_uncalledPairs;
			}
			// else if(m_useBarRead && uncalledBR[i]){
			//
			// 	// to be handled
			// }
			else{
				++nEntries;
				
				if(m_useNumberTag){
					stringstream converter;
					converter << ++m_tagCounter;
					TString tagCount = converter.str();
					
					                 b->srd.ids[i]   = tagCount;
					if(m_isPaired)   b->srd2.ids[i]  = tagCount;
					if(m_useBarRead) b->srdBR.ids[i] = tagCount;
				}
				
				TSeqRead *read1 = NULL, *read2 = NULL, *barRead = NULL;
				
				if(m_format == FASTA){
					                 read1   = new TSeqRead(b->srd.seqs[i],   b->srd.ids[i],   b->srd.ids[i]);
					if(m_isPaired)   read2   = new TSeqRead(b->srd2.seqs[i],  b->srd2.ids[i],  b->srd2.ids[i]);
					if(m_useBarRead) barRead = new TSeqRead(b->srdBR.seqs[i], b->srdBR.ids[i], b->srdBR.ids[i]);
				}
				else{
					                 read1   = new TSeqRead(b->srd.seqs[i],   b->srd.ids[i],   b->srd.quals[i]);
					if(m_isPaired)   read2   = new TSeqRead(b->srd2.seqs[i],  b->srd2.ids[i],  b->srd2.quals[i]);
					if(m_useBarRead) barRead = new TSeqRead(b->srdBR.seqs[i], b->srdBR.ids[i], b->srdBR.quals[i]);
				}
				
				b->pReads.push_back(new TPairedRead(read1, read2, barRead));
			}
		}
		
		if(nEntries == 0){
			delete b;
			b = NULL;
			return loadPairedReadBundle();
		}
		
		return b;
	}
	
	
	// tbb filter operator
	void* operator()(void*){
		
		using namespace flexbar;
		
		PairedReadBundle *prBundle = NULL;
		
		prBundle = static_cast< PairedReadBundle* >(loadPairedReadBundle());
		
		return prBundle;
	}
	
	// virtual
	void finalize(void* item){
		
	}
	
	
	unsigned long getNrUncalledReads() const{
		return m_uncalled;
	}
	
	
	unsigned long getNrUncalledPairedReads() const{
		return m_uncalledPairs;
	}
	
	
	unsigned long getNrProcessedReads() const{
		if(m_isPaired) return m_f1->getNrProcessedReads() + m_f2->getNrProcessedReads();
		else           return m_f1->getNrProcessedReads();
	}
	
	
	unsigned long getNrProcessedChars() const{
		if(m_isPaired) return m_f1->getNrProcessedChars() + m_f2->getNrProcessedChars();
		else           return m_f1->getNrProcessedChars();
	}
	
	
	unsigned long getNrLowPhredReads() const {
		if(m_isPaired) return m_f1->getNrLowPhredReads() + m_f2->getNrLowPhredReads();
		else           return m_f1->getNrLowPhredReads();
	}
	
};

#endif
