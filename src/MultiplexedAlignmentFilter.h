/*
 *   MultiplexedAlignmentFilter.h
 *
 *   Authors: mat and jtr
 */

#ifndef FLEXBAR_MULTIPLEXEDALIGNMENTFILTER_H_
#define FLEXBAR_MULTIPLEXEDALIGNMENTFILTER_H_

#include <string>

#include <tbb/pipeline.h>
#include <tbb/concurrent_vector.h>

#include "Enums.h"
#include "Options.h"
#include "MultiplexedRead.h"
#include "AlignmentFilter.h"
#include "AlignmentAlgorithm.h"
#include "AdapterLoader.h"


// Processes MultiplexedRead and assigns barcode to read or removes adapters.

template <typename TString, typename TIDString>
class MultiplexedAlignmentFilter : public tbb::filter {

private:
	
	const bool m_writeUnassigned;
	
	const flexbar::LogLevel       m_verb;
	const flexbar::RunType        m_runType;
	const flexbar::BarcodeDetect  m_barType;
	const flexbar::AdapterRemoval m_adapRem;
	
	tbb::atomic<unsigned long> m_unassigned;
	
	tbb::concurrent_vector<TAdapter> *m_adapters;
	tbb::concurrent_vector<TAdapter> *m_barcodes;
	
	typedef AlignmentFilter<TString, TIDString, AlignmentAlgorithm<TString> > AliFilter;
	AliFilter *m_afilter, *m_bfilter;
	
	std::ostream *out;
	
public:
	
	MultiplexedAlignmentFilter(Options &o) :
		
		filter(parallel),
		m_verb(o.logLevel),
		m_writeUnassigned(o.writeUnassigned),
		m_runType(o.runType),
		m_barType(o.barDetect),
		m_adapRem(o.adapRm),
		out(o.out){
		
		m_barcodes = &o.barcodes;
		m_adapters = &o.adapters;
		
		m_unassigned = 0;
		
		m_afilter = new AliFilter(m_adapters, o, o.a_min_overlap, o.a_threshold, o.a_tail_len, o.match, o.mismatch, o.gapCost, o.end, false);
		m_bfilter = new AliFilter(m_barcodes, o, o.b_min_overlap, o.b_threshold, o.b_tail_len, o.b_match, o.b_mismatch, o.b_gapCost, o.b_end, true);
		
		if(m_verb == flexbar::TAB)
		*out << "ReadTag\tQueryTag\tQueryStart\tQueryEnd\tOverlapLength\tMismatches\tIndels\tAllowedErrors" << std::endl;
	}
	
	
	virtual ~MultiplexedAlignmentFilter(){
		delete m_afilter;
		delete m_bfilter;
	};
	
	
	void* operator()(void* item){
		
		using namespace flexbar;
		
		if(item != NULL){
			MultiplexedRead<TString, TIDString> *myRead = static_cast< MultiplexedRead<TString, TIDString>* >(item);
			
			bool skipAdapRem = false;
			
			// barcode detection
			if(m_barType != BOFF){
				switch(m_barType){
					case BARCODE_READ:         myRead->m_barcode_id = m_bfilter->align(myRead->m_b,  false); break;
					case WITHIN_READ_REMOVAL:  myRead->m_barcode_id = m_bfilter->align(myRead->m_r1, true);  break;
					case WITHIN_READ:          myRead->m_barcode_id = m_bfilter->align(myRead->m_r1, false); break;
					case BOFF:                                                                               break;
				}
				
				if(myRead->m_barcode_id == 0){
					m_unassigned++;
					
					if(! m_writeUnassigned) skipAdapRem = true;
				}
			}
			
			// adapter removal
			if(m_adapRem != AOFF && ! skipAdapRem){
				if(m_adapRem != ATWO)
				m_afilter->align(myRead->m_r1, true);
				
				if(myRead->m_r2 != NULL && m_adapRem != AONE)
				m_afilter->align(myRead->m_r2, true);
			}
			return myRead;
		}
		else return NULL;
	}
	
	
	unsigned long getNrUnassignedReads() const {
		if(m_runType == flexbar::PAIRED_BARCODED) return m_unassigned * 2;
		else                                      return m_unassigned;
	}
	
	
	unsigned long getNrPreShortReads() const {
		return m_afilter->getNrPreShortReads();
	}
	
	
	void printAdapterOverlapStats(){
		
		if(m_afilter->getNrModifiedReads() > 0){
			*out << m_afilter->getOverlapStatsString() << "\n\n";
		}
		*out << std::endl;
	}
	
};

#endif /* FLEXBAR_MULTIPLEXEDALIGNMENTFILTER_H_ */
