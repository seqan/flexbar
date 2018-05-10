// PairedInput.h

#ifndef FLEXBAR_PAIREDINPUT_H
#define FLEXBAR_PAIREDINPUT_H

#include "SeqInput.h"


template <typename TSeqStr, typename TString>
class PairedInput : public tbb::filter {

private:
	
	const flexbar::FileFormat m_format;
	const bool m_isPaired, m_useBarRead, m_useNumberTag, m_interleaved;
	const unsigned int m_bundleSize;
	
	tbb::atomic<unsigned long> m_uncalled, m_uncalledPairs, m_tagCounter;
	SeqInput<TSeqStr, TString> *m_f1, *m_f2, *m_b;
	
public:
	
	PairedInput(const Options &o) :
		
		filter(serial_in_order),
		m_format(o.format),
		m_useNumberTag(o.useNumberTag),
		m_interleaved(o.interleavedInput),
		m_isPaired(o.isPaired),
		m_useBarRead(o.barDetect == flexbar::BARCODE_READ),
		m_bundleSize(o.bundleSize),
		m_tagCounter(0),
		m_uncalled(0),
		m_uncalledPairs(0){
		
		m_f1 = new SeqInput<TSeqStr, TString>(o, o.readsFile, true, o.useStdin);
		
		m_f2 = NULL;
		m_b  = NULL;
		
		if(m_isPaired && ! m_interleaved)
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
		
		TSeqStrs seqs,     seqs2,     seqsBR;
		TStrings ids,      ids2,      idsBR;
		TStrings quals,    quals2,    qualsBR;
		TBools   uncalled, uncalled2, uncalledBR;
		
		unsigned int nReads = m_f1->loadSeqReads(uncalled, ids, seqs, quals, m_bundleSize);
		
		if(m_interleaved && nReads % 2 == 1){
			cerr << "\nERROR: Interleaved reads input does not contain even number of reads.\n" << endl;
			exit(1);
		}
		
		if(m_isPaired && ! m_interleaved){
			unsigned int nReads2 = m_f2->loadSeqReads(uncalled2, ids2, seqs2, quals2, m_bundleSize);
			
			if(nReads != nReads2){
				cerr << "\nERROR: Read without counterpart in paired input mode.\n" << endl;
				exit(1);
			}
		}
		
		if(m_useBarRead){
			
			unsigned int bundleSize = m_bundleSize;
			unsigned int multi      = 1;
			
			if(m_interleaved){
				bundleSize = (unsigned int) (m_bundleSize / 2);
				multi      = 2;
			}
			
			unsigned int nBarReads = m_b->loadSeqReads(uncalledBR, idsBR, seqsBR, qualsBR, bundleSize);
			
			if(nReads > nBarReads * multi){
				cerr << "\nERROR: Read without barcode read in input.\n" << endl;
				exit(1);
			}
			else if(nReads < nBarReads * multi){
				cerr << "\nERROR: Barcode read without read in input.\n" << endl;
				exit(1);
			}
		}
		
		if(nReads == 0) return NULL;
		
		
		TPairedReadBundle *prBundle = new TPairedReadBundle();
		
		if(! m_interleaved){
			
			for(unsigned int i = 0; i < length(ids); ++i){
				
				if(uncalled[i] || (m_isPaired && uncalled2[i])){
					
					if(uncalled[i])                ++m_uncalled;
					if(m_isPaired && uncalled2[i]) ++m_uncalled;
					if(m_isPaired)                 ++m_uncalledPairs;
				}
				// else if(m_useBarRead && uncalledBR[i]){
				//
				// 	// to be handled
				// }
				else{
					
					if(m_useNumberTag){
						stringstream converter;
						converter << ++m_tagCounter;
						TString tagCount = converter.str();
						
						                 ids[i]   = tagCount;
						if(m_isPaired)   ids2[i]  = tagCount;
						if(m_useBarRead) idsBR[i] = tagCount;
					}
					
					TSeqRead *read1 = NULL, *read2 = NULL, *barRead = NULL;
					
					if(m_format == FASTA){
						                 read1   = new TSeqRead(seqs[i],   ids[i]);
						if(m_isPaired)   read2   = new TSeqRead(seqs2[i],  ids2[i]);
						if(m_useBarRead) barRead = new TSeqRead(seqsBR[i], idsBR[i]);
					}
					else{
						                 read1   = new TSeqRead(seqs[i],   ids[i],   quals[i]);
						if(m_isPaired)   read2   = new TSeqRead(seqs2[i],  ids2[i],  quals2[i]);
						if(m_useBarRead) barRead = new TSeqRead(seqsBR[i], idsBR[i], qualsBR[i]);
					}
					
					prBundle->push_back(new TPairedRead(read1, read2, barRead));
				}
			}
		}
		else{  // interleaved paired input
			
			unsigned int nEntries = (unsigned int) (length(ids) / 2);
			
			for(unsigned int i = 0; i < nEntries; ++i){
				
				unsigned int r = (i * 2);
				unsigned int p = (i * 2) + 1;
				
				if(uncalled[r] || uncalled[p]){
					
					if(uncalled[r]) ++m_uncalled;
					if(uncalled[p]) ++m_uncalled;
					
					++m_uncalledPairs;
				}
				// else if(m_useBarRead && uncalledBR[i]){
				//
				// 	// to be handled
				// }
				else{
					
					if(m_useNumberTag){
						stringstream converter;
						converter << ++m_tagCounter;
						TString tagCount = converter.str();
						
						                 ids[r]   = tagCount;
						                 ids[p]   = tagCount;
						if(m_useBarRead) idsBR[i] = tagCount;
					}
					
					TSeqRead *read1 = NULL, *read2 = NULL, *barRead = NULL;
					
					if(m_format == FASTA){
						                 read1   = new TSeqRead(seqs[r],   ids[r]);
						                 read2   = new TSeqRead(seqs[p],   ids[p]);
						if(m_useBarRead) barRead = new TSeqRead(seqsBR[i], idsBR[i]);
					}
					else{
						                 read1   = new TSeqRead(seqs[r],   ids[r],   quals[r]);
						                 read2   = new TSeqRead(seqs[p],   ids[p],   quals[p]);
						if(m_useBarRead) barRead = new TSeqRead(seqsBR[i], idsBR[i], qualsBR[i]);
					}
					
					prBundle->push_back(new TPairedRead(read1, read2, barRead));
				}
			}
		}
		
		return prBundle;
	}
	
	
	// tbb filter operator
	void* operator()(void*){
		
		using namespace flexbar;
		
		TPairedReadBundle *prBundle = NULL;
		
		prBundle = static_cast< TPairedReadBundle* >(loadPairedReadBundle());
		
		if(prBundle != NULL){
			
			while(prBundle->size() == 0){
				delete prBundle;
				prBundle = NULL;
				
				prBundle = static_cast< TPairedReadBundle* >(loadPairedReadBundle());
				
				if(prBundle == NULL) return prBundle;
			}
		}
		
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
