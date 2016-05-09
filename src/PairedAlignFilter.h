/*
 *   PairedAlignFilter.h
 *
 *   Authors: mat and jtr
 */

#ifndef FLEXBAR_PAIREDALIGNFILTER_H
#define FLEXBAR_PAIREDALIGNFILTER_H

#include "SeqAlignFilter.h"
#include "SeqAlignAlgorithm.h"


template <typename TSeqStr, typename TString>
class PairedAlignFilter : public tbb::filter {

private:
	
	const bool m_writeUnassigned, m_twoBarcodes;
	
	const flexbar::LogLevel       m_verb;
	const flexbar::RunType        m_runType;
	const flexbar::BarcodeDetect  m_barType;
	const flexbar::AdapterRemoval m_adapRem;
	
	tbb::atomic<unsigned long> m_unassigned;
	
	tbb::concurrent_vector<flexbar::TAdapter> *m_adapters, *m_adapters2;
	tbb::concurrent_vector<flexbar::TAdapter> *m_barcodes, *m_barcodes2;
	
	typedef SeqAlignFilter<TSeqStr, TString, SeqAlignAlgorithm<TSeqStr> > AlignFilter;
	AlignFilter *m_afilter, *m_bfilter, *m_a2filter, *m_b2filter;
	
	std::ostream *out;
	
public:
	
	PairedAlignFilter(Options &o) :
		
		filter(parallel),
		m_verb(o.logLevel),
		m_runType(o.runType),
		m_barType(o.barDetect),
		m_adapRem(o.adapRm),
		m_writeUnassigned(o.writeUnassigned),
		m_twoBarcodes(o.barDetect == flexbar::WITHIN_READ_REMOVAL2 || o.barDetect == flexbar::WITHIN_READ2),
		out(o.out){
		
		m_unassigned = 0;
		
		m_barcodes  = &o.barcodes;
		m_adapters  = &o.adapters;
		m_barcodes2 = &o.barcodes2;
		m_adapters2 = &o.adapters2;
		
		m_bfilter = new AlignFilter(m_barcodes, o, o.b_min_overlap, o.b_threshold, o.b_tail_len, o.b_match, o.b_mismatch, o.b_gapCost, o.b_end, true);
		m_afilter = new AlignFilter(m_adapters, o, o.a_min_overlap, o.a_threshold, o.a_tail_len, o.match, o.mismatch, o.gapCost, o.end, false);
		
		m_b2filter = new AlignFilter(m_barcodes2, o, o.b_min_overlap, o.b_threshold, o.b_tail_len, o.b_match, o.b_mismatch, o.b_gapCost, o.b_end, true);
		m_a2filter = new AlignFilter(m_adapters2, o, o.a_min_overlap, o.a_threshold, o.a_tail_len, o.match, o.mismatch, o.gapCost, o.end, false);
		
		if(m_verb == flexbar::TAB)
		*out << "ReadTag\tQueryTag\tQueryStart\tQueryEnd\tOverlapLength\tMismatches\tIndels\tAllowedErrors" << std::endl;
	}
	
	
	virtual ~PairedAlignFilter(){
		delete m_bfilter;
		delete m_afilter;
		delete m_b2filter;
		delete m_a2filter;
	};
	
	
	void alignPairedRead(void* item, flexbar::TAlignBundle &alignBundle, const bool preCycle){
		
		using namespace flexbar;
		
		if(item != NULL){
			PairedRead<TSeqStr, TString> *pRead = static_cast< PairedRead<TSeqStr, TString>* >(item);
			
			bool skipAdapRem = false;
			
			// barcode detection
			if(m_barType != BOFF){
				switch(m_barType){
					case BARCODE_READ:         pRead->m_barcode_id  = m_bfilter->align(pRead->m_b,   false, alignBundle.at(0), preCycle); break;
					case WITHIN_READ_REMOVAL2: pRead->m_barcode_id2 = m_b2filter->align(pRead->m_r2, true,  alignBundle.at(2), preCycle);
					case WITHIN_READ_REMOVAL:  pRead->m_barcode_id  = m_bfilter->align(pRead->m_r1,  true,  alignBundle.at(1), preCycle); break;
					case WITHIN_READ2:         pRead->m_barcode_id2 = m_b2filter->align(pRead->m_r2, false, alignBundle.at(2), preCycle);
					case WITHIN_READ:          pRead->m_barcode_id  = m_bfilter->align(pRead->m_r1,  false, alignBundle.at(1), preCycle); break;
					case BOFF:                                                                                                            break;
				}
				
				if(pRead->m_barcode_id == 0 || (m_twoBarcodes && pRead->m_barcode_id2 == 0)){
					
					if(! preCycle) m_unassigned++;
					
					if(! m_writeUnassigned) skipAdapRem = true;
				}
			}
			
			// adapter removal
			if(m_adapRem != AOFF && ! skipAdapRem){
				if(m_adapRem != ATWO)
				m_afilter->align(pRead->m_r1, true, alignBundle.at(3), preCycle);
				
				if(pRead->m_r2 != NULL && m_adapRem != AONE){
					if(m_adapRem != NORMAL2) m_afilter->align(pRead->m_r2,  true, alignBundle.at(4), preCycle);
					else                     m_a2filter->align(pRead->m_r2, true, alignBundle.at(4), preCycle);
				}
			}
		}
	}
	
	
	// tbb filter operator
	void* operator()(void* item){
		
		using namespace flexbar;
		
		if(item != NULL){
			
			TPairedReadBundle *prBundle = static_cast< TPairedReadBundle* >(item);
			
			TAlignBundle alignBundle;
			alignBundle.reserve(5);
			
			TAlignments r1AlignmentsB, r2AlignmentsB, bAlignmentsB;
			TAlignments r1AlignmentsA, r2AlignmentsA;
			
			alignBundle.push_back(bAlignmentsB);
			alignBundle.push_back(r1AlignmentsB);
			alignBundle.push_back(r2AlignmentsB);
			alignBundle.push_back(r1AlignmentsA);
			alignBundle.push_back(r2AlignmentsA);
			
			for(unsigned int i = 0; i < prBundle->size(); ++i){
				
				// alignPairedRead(prBundle->at(i), alignBundle, true);
				// alignPairedRead(prBundle->at(i), alignBundle, false);
				
				alignPairedRead(prBundle->at(i), alignBundle, false);
			}
			return prBundle;
		}
		else return NULL;
	}
	
	
	unsigned long getNrUnassignedReads() const {
		
		using namespace flexbar;
		
		if(m_runType == PAIRED_BARCODED) return m_unassigned * 2;
		else                             return m_unassigned;
	}
	
	
	unsigned long getNrPreShortReads() const {
		
		using namespace flexbar;
		
		if(m_adapRem != NORMAL2) return m_afilter->getNrPreShortReads();
		else return m_afilter->getNrPreShortReads() + m_a2filter->getNrPreShortReads();
	}
	
	
	void printAdapterOverlapStats(){
		
		using namespace flexbar;
		
		if(m_afilter->getNrModifiedReads() > 0){
			*out << m_afilter->getOverlapStatsString() << "\n\n";
		}
		
		if(m_adapRem != NORMAL2) *out << std::endl;
	}
	
	
	void printAdapterOverlapStats2(){
		
		if(m_a2filter->getNrModifiedReads() > 0){
			*out << m_a2filter->getOverlapStatsString() << "\n\n";
		}
		*out << std::endl;
	}
	
};

#endif
