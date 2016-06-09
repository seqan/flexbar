/*
 *   PairedAlign.h
 *
 *   Authors: mat and jtr
 */

#ifndef FLEXBAR_PAIREDALIGN_H
#define FLEXBAR_PAIREDALIGN_H

#include "SeqAlign.h"
#include "SeqAlignAlgo.h"


template <typename TSeqStr, typename TString>
class PairedAlign : public tbb::filter {

private:
	
	const bool m_writeUnassigned, m_twoBarcodes;
	
	const flexbar::LogAlign       m_log;
	const flexbar::RunType        m_runType;
	const flexbar::BarcodeDetect  m_barType;
	const flexbar::AdapterRemoval m_adapRem;
	
	tbb::atomic<unsigned long> m_unassigned;
	
	tbb::concurrent_vector<flexbar::TBar> *m_adapters, *m_adapters2;
	tbb::concurrent_vector<flexbar::TBar> *m_barcodes, *m_barcodes2;
	
	typedef SeqAlign<TSeqStr, TString, SeqAlignAlgo<TSeqStr> > TAlignFilter;
	TAlignFilter *m_afilter, *m_bfilter, *m_a2filter, *m_b2filter;
	
	std::ostream *out;
	
public:
	
	PairedAlign(Options &o) :
		
		filter(parallel),
		m_log(o.logAlign),
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
		
		m_bfilter = new TAlignFilter(m_barcodes, o, o.b_min_overlap, o.b_errorRate, o.b_tail_len, o.b_match, o.b_mismatch, o.b_gapCost, o.b_end, true);
		m_afilter = new TAlignFilter(m_adapters, o, o.a_min_overlap, o.a_errorRate, o.a_tail_len, o.match, o.mismatch, o.gapCost, o.end, false);
		
		m_b2filter = new TAlignFilter(m_barcodes2, o, o.b_min_overlap, o.b_errorRate, o.b_tail_len, o.b_match, o.b_mismatch, o.b_gapCost, o.b_end, true);
		m_a2filter = new TAlignFilter(m_adapters2, o, o.a_min_overlap, o.a_errorRate, o.a_tail_len, o.match, o.mismatch, o.gapCost, o.end, false);
		
		if(m_log == flexbar::TAB)
		*out << "ReadTag\tQueryTag\tQueryStart\tQueryEnd\tOverlapLength\tMismatches\tIndels\tAllowedErrors" << std::endl;
	}
	
	
	virtual ~PairedAlign(){
		delete m_bfilter;
		delete m_afilter;
		delete m_b2filter;
		delete m_a2filter;
	};
	
	
	void alignPairedRead(flexbar::TPairedRead* pRead, flexbar::TAlignBundle &alBundle, flexbar::ComputeCycle cycle, unsigned int &idxAl){
		
		using namespace flexbar;
		
		bool skipAdapRem = false;
		
		// barcode detection
		if(m_barType != BOFF){
			switch(m_barType){
				case BARCODE_READ:         pRead->barID  = m_bfilter->alignSeqRead(pRead->b,   false, alBundle[0], cycle, idxAl); break;
				
				case WITHIN_READ_REMOVAL2: pRead->barID2 = m_b2filter->alignSeqRead(pRead->r2, true,  alBundle[2], cycle, idxAl);
				case WITHIN_READ_REMOVAL:  pRead->barID  = m_bfilter->alignSeqRead(pRead->r1,  true,  alBundle[1], cycle, idxAl); break;
				
				case WITHIN_READ2:         pRead->barID2 = m_b2filter->alignSeqRead(pRead->r2, false, alBundle[2], cycle, idxAl);
				case WITHIN_READ:          pRead->barID  = m_bfilter->alignSeqRead(pRead->r1,  false, alBundle[1], cycle, idxAl); break;
				
				case BOFF: break;
			}
			
			if(pRead->barID == 0 || (m_twoBarcodes && pRead->barID2 == 0)){
				
				if(cycle != PRELOAD)    m_unassigned++;
				if(! m_writeUnassigned) skipAdapRem = true;
			}
		}
		
		// adapter removal
		if(m_adapRem != AOFF && ! skipAdapRem){
			if(m_adapRem != ATWO)
			m_afilter->alignSeqRead(pRead->r1, true, alBundle[3], cycle, idxAl);
			
			if(pRead->r2 != NULL && m_adapRem != AONE){
				if(m_adapRem != NORMAL2) m_afilter->alignSeqRead(pRead->r2,  true, alBundle[4], cycle, idxAl);
				else                     m_a2filter->alignSeqRead(pRead->r2, true, alBundle[4], cycle, idxAl);
			}
		}
	}
	
	
	// tbb filter operator
	void* operator()(void* item){
		
		using namespace flexbar;
		
		if(item != NULL){
			
			PairedReadBundle *prBundle = static_cast<PairedReadBundle* >(item);
			
			TAlignBundle alBundle;
			alBundle.reserve(5);
			
			TAlignments r1AlignmentsB, r2AlignmentsB, bAlignmentsB;
			TAlignments r1AlignmentsA, r2AlignmentsA;
			
			alBundle.push_back(bAlignmentsB);
			alBundle.push_back(r1AlignmentsB);
			alBundle.push_back(r2AlignmentsB);
			alBundle.push_back(r1AlignmentsA);
			alBundle.push_back(r2AlignmentsA);
			
			unsigned int idxAl = 0;
			ComputeCycle cycle = PRELOAD;
			
			// TPairedReadBundleVec prbv;
			//
			// for(unsigned int i = 0; i < length(prBundle->srd->ids); ++i){
			//
			// 	if(! srd->uncalled[i] && ! (m_isPaired && srd2->uncalled[i])){
			//
			//
			// 	}
			// }
			
			for(unsigned int i = 0; i < prBundle->prbv.size(); ++i)
				alignPairedRead(prBundle->prbv.at(i), alBundle, cycle, idxAl);
			
			idxAl = 0;
			cycle = COMPUTE;
			
			for(unsigned int i = 0; i < prBundle->prbv.size(); ++i){
				
				if(i > 0) cycle = RESULTS;
				
				alignPairedRead(prBundle->prbv.at(i), alBundle, cycle, idxAl);
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
