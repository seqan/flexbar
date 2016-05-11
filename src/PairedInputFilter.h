/*
 *   PairedInputFilter.h
 *
 *   Authors: mat and jtr
 */

#ifndef FLEXBAR_PAIREDINPUTFILTER_H
#define FLEXBAR_PAIREDINPUTFILTER_H

#include "SeqInput.h"


template <typename TSeqStr, typename TString>
class PairedInputFilter : public tbb::filter {

private:
	
	const bool m_isPaired, m_useBarcodeRead, m_useNumberTag;
	const unsigned int m_bundleSize;
	
	tbb::atomic<unsigned long> m_uncalled, m_uncalledPairs, m_tagCounter;
	
	SeqInput<TSeqStr, TString> *m_f1, *m_f2, *m_b;
	
public:
	
	PairedInputFilter(const Options &o) :
		
		filter(serial_in_order),
		m_useNumberTag(o.useNumberTag),
		m_isPaired(o.isPaired),
		m_useBarcodeRead(o.barDetect == flexbar::BARCODE_READ),
		m_bundleSize(o.bundleSize){
		
		m_tagCounter    = 0;
		m_uncalled      = 0;
		m_uncalledPairs = 0;
		
		m_f1 = new SeqInput<TSeqStr, TString>(o, o.readsFile, false, true, o.useStdin);
		
		m_f2 = NULL;
		m_b  = NULL;
		
		if(m_isPaired){
			m_f2 = new SeqInput<TSeqStr, TString>(o, o.readsFile2, false, true, false);
		}
		
		if(m_useBarcodeRead){
			m_b = new SeqInput<TSeqStr, TString>(o, o.barReadsFile, false, false, false);
		}
	}
	
	
	virtual ~PairedInputFilter(){
		delete m_f1;
		delete m_f2;
		delete m_b;
	}
	
	
	void* getPairedRead(){
		
		using namespace std;
		
		SeqRead<TSeqStr, TString> *read1 = NULL, *read2 = NULL, *barRead = NULL;
		
		bool uncalled = true, uncalled2 = true, uBR = true;
		
		if(! m_isPaired){
			
			while(uncalled){
				read1 = static_cast< SeqRead<TSeqStr, TString>* >(m_f1->getRead(uncalled));
				
				if(m_useBarcodeRead) barRead = static_cast< SeqRead<TSeqStr, TString>* >(m_b->getRead(uBR));
				
				if(read1 == NULL && m_useBarcodeRead && barRead != NULL){
					cerr << "Error: barcode read without read, or file reading error.\n" << endl;
					exit(1);
				}
				else if(read1 == NULL){
					return NULL;
				}
				else if(m_useBarcodeRead && barRead == NULL){
					cerr << "Error: read without barcode read, or file reading error.\n" << endl;
					exit(1);
				}
				
				if(uncalled){
					++m_uncalled;
					delete read1;
					delete barRead;
				}
			}
		}
		
		// paired read input
		else{
			
			while(uncalled || uncalled2){
				
				read1 = static_cast< SeqRead<TSeqStr, TString>* >(m_f1->getRead(uncalled));
				read2 = static_cast< SeqRead<TSeqStr, TString>* >(m_f2->getRead(uncalled2));
				
				if(m_useBarcodeRead) barRead = static_cast< SeqRead<TSeqStr, TString>* >(m_b->getRead(uBR));
				
				
				if((read1 == NULL || read2 == NULL) && m_useBarcodeRead && barRead != NULL){
					cerr << "Error: barcode read without read-pair or file reading error.\n" << endl;
					exit(1);
				}
				else if(read1 == NULL && read2 == NULL){
					return NULL;
				}
				else if(read1 == NULL || read2 == NULL){
					cerr << "Error: single read in paired mode, or file reading error.\n" << endl;
					exit(1);
				}
				else if(m_useBarcodeRead && barRead == NULL){
					cerr << "Error: reads without barcode read or file reading error.\n" << endl;
					exit(1);
				}
				
				if(uncalled || uncalled2){
					++m_uncalledPairs;
					if(uncalled)  ++m_uncalled;
					if(uncalled2) ++m_uncalled;
					
					delete read1;
					delete read2;
					delete barRead;
				}
			}
		}
		
		if(m_useNumberTag){
			stringstream converter;
			converter << ++m_tagCounter;
			TString tagCount = converter.str();
			
			read1->setSequenceTag(tagCount);
			if(m_isPaired) read2->setSequenceTag(tagCount);
			if(m_useBarcodeRead) barRead->setSequenceTag(tagCount);
		}
		
		return new PairedRead<TSeqStr, TString>(read1, read2, barRead);
	}
	
	
	// tbb filter operator
	void* operator()(void*){
		
		using namespace flexbar;
		
		TPairedReadBundle *prBundle = new TPairedReadBundle();
		
		prBundle->reserve(m_bundleSize);
		
		for(unsigned int i = 0; i < m_bundleSize; ++i){
			
			PairedRead<TSeqStr, TString> *pRead = NULL;
			
			pRead = static_cast< PairedRead<TSeqStr, TString>* >(getPairedRead());
			
			     if(pRead == NULL && i == 0) return NULL;
			else if(pRead == NULL)           break;
			else                             prBundle->push_back(pRead);
		}
		return prBundle;
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
