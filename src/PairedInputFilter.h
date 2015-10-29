/*
 *   PairedInputFilter.h
 *
 *   Authors: mat and jtr
 */

#ifndef FLEXBAR_PAIREDINPUTFILTER_H
#define FLEXBAR_PAIREDINPUTFILTER_H

#include <tbb/pipeline.h>

#include "Options.h"
#include "PairedRead.h"
#include "SeqInputFilter.h"


template <typename TSeqStr, typename TString, typename TStreamR, typename TStreamP, typename TStreamB>
class PairedInputFilter : public tbb::filter {

private:
	
	const bool m_isPaired, m_useBarcodeRead, m_useNumberTag;
	tbb::atomic<unsigned long> m_uncalled, m_uncalledPairs, m_tagCounter;
	
	SeqInputFilter<TSeqStr, TSeqStr, TStreamR> *m_f1;
	SeqInputFilter<TSeqStr, TSeqStr, TStreamP> *m_f2;
	SeqInputFilter<TSeqStr, TSeqStr, TStreamB> *m_b;
	
public:
	
	PairedInputFilter(const Options &o) :
		
		filter(serial_in_order),
		m_useNumberTag(o.useNumberTag),
		m_isPaired(o.isPaired),
		m_useBarcodeRead(o.barDetect == flexbar::BARCODE_READ){
		
		m_tagCounter    = 0;
		m_uncalled      = 0;
		m_uncalledPairs = 0;
		
		m_f1 = new SeqInputFilter<TSeqStr, TSeqStr, TStreamR>(o, o.readsFile, false, true, o.useStdin);
		
		m_f2 = NULL;
		m_b  = NULL;
		
		if(m_isPaired){
			m_f2 = new SeqInputFilter<TSeqStr, TSeqStr, TStreamP>(o, o.readsFile2, false, true, false);
		}
		
		if(m_useBarcodeRead){
			m_b = new SeqInputFilter<TSeqStr, TSeqStr, TStreamB>(o, o.barReadsFile, false, false, false);
		}
	}
	
	
	virtual ~PairedInputFilter(){
		delete m_f1;
		delete m_f2;
		delete m_b;
	}
	
	
	void* operator()(void*){
		
		using namespace std;
		
		SeqRead<TSeqStr, TString> *myRead1 = NULL, *myRead2 = NULL, *myBarcode = NULL;
		
		bool uncalled = true, uncalled2 = true, uBR = true;
		
		if(! m_isPaired){
			
			while(uncalled){
				myRead1 = static_cast< SeqRead<TSeqStr, TString>* >(m_f1->getRead(uncalled));
				
				if(m_useBarcodeRead) myBarcode = static_cast< SeqRead<TSeqStr, TString>* >(m_b->getRead(uBR));
				
				if(myRead1 == NULL) return NULL;
				
				else if(m_useBarcodeRead && myBarcode == NULL){
					cerr << "Error: read without barcode read, or file reading error!\n" << endl;
					exit(1);
				}
				
				if(uncalled){
					++m_uncalled;
					delete myRead1;
					delete myBarcode;
				}
			}
		}
		
		// paired read input
		else{
			
			while(uncalled || uncalled2){
				
				myRead1 = static_cast< SeqRead<TSeqStr, TString>* >(m_f1->getRead(uncalled));
				myRead2 = static_cast< SeqRead<TSeqStr, TString>* >(m_f2->getRead(uncalled2));
				
				if(m_useBarcodeRead) myBarcode = static_cast< SeqRead<TSeqStr, TString>* >(m_b->getRead(uBR));
				
				// end of files reached
				if(myRead1 == NULL && myRead2 == NULL) return NULL;
				
				else if(myRead1 == NULL || myRead2 == NULL){
					cerr << "Error: single read in paired mode, or file reading error!\n" << endl;
					exit(1);
				}
				else if(m_useBarcodeRead && myBarcode == NULL){
					cerr << "Error: reads without barcode read or file reading error!\n" << endl;
					exit(1);
				}
				
				if(uncalled || uncalled2){
					++m_uncalledPairs;
					if(uncalled)  ++m_uncalled;
					if(uncalled2) ++m_uncalled;
					
					delete myRead1;
					delete myRead2;
					delete myBarcode;
				}
			}
		}
		
		if(m_useNumberTag){
			stringstream converter;
			converter << ++m_tagCounter;
			TSeqStr tagCount = converter.str();
			
			myRead1->setSequenceTag(tagCount);
			if(m_isPaired) myRead2->setSequenceTag(tagCount);
			if(m_useBarcodeRead) myBarcode->setSequenceTag(tagCount);
		}
		
		return new PairedRead<TSeqStr, TString>(myRead1, myRead2, myBarcode);
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
