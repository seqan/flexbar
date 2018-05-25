// PairedAlign.h

#ifndef FLEXBAR_PAIREDALIGN_H
#define FLEXBAR_PAIREDALIGN_H

#include "SeqAlign.h"
#include "SeqAlignPair.h"
#include "SeqAlignAlgo.h"


template <typename TSeqStr, typename TString>
class PairedAlign : public tbb::filter {

private:
	
	const bool m_writeUnassigned, m_twoBarcodes, m_umiTags, m_useRcTrimEnd;
	const bool m_htrim, m_htrimAdapterRm, m_htrimMaxFirstOnly;
	
	const std::string m_htrimLeft, m_htrimRight;
	
	const int m_htrimMinLength, m_htrimMaxLength;
	const float m_htrimErrorRate;
	
	const unsigned int m_arTimes;
	
	const flexbar::FileFormat     m_format;
	const flexbar::LogAlign       m_log;
	const flexbar::RunType        m_runType;
	const flexbar::BarcodeDetect  m_barType;
	const flexbar::AdapterRemoval m_adapRem;
	const flexbar::TrimEnd        m_aTrimEnd, m_arcTrimEnd, m_bTrimEnd;
	const flexbar::PairOverlap    m_poMode;
	
	tbb::atomic<unsigned long> m_unassigned;
	tbb::concurrent_vector<flexbar::TBar> *m_adapters, *m_adapters2;
	tbb::concurrent_vector<flexbar::TBar> *m_barcodes, *m_barcodes2;
	
	typedef SeqAlign<TSeqStr, TString, SeqAlignAlgo<TSeqStr> > TSeqAlign;
	TSeqAlign *m_a1, *m_b1, *m_a2, *m_b2;
	
	typedef SeqAlignPair<TSeqStr, TString, SeqAlignAlgo<TSeqStr> > TSeqAlignPair;
	TSeqAlignPair *m_p;
	
	std::ostream *out;
	
public:
	
	PairedAlign(Options &o) :
		
		filter(parallel),
		m_format(o.format),
		m_log(o.logAlign),
		m_runType(o.runType),
		m_barType(o.barDetect),
		m_adapRem(o.adapRm),
		m_poMode(o.poMode),
		m_aTrimEnd(o.a_end),
		m_arcTrimEnd(o.arc_end),
		m_bTrimEnd(o.b_end),
		m_arTimes(o.a_cycles),
		m_umiTags(o.umiTags),
		m_useRcTrimEnd(o.useRcTrimEnd),
		m_writeUnassigned(o.writeUnassigned),
		m_htrimLeft(o.htrimLeft),
		m_htrimRight(o.htrimRight),
		m_htrimMinLength(o.htrimMinLength),
		m_htrimMaxLength(o.htrimMaxLength),
		m_htrimMaxFirstOnly(o.htrimMaxFirstOnly),
		m_htrimErrorRate(o.h_errorRate),
		m_htrimAdapterRm(o.htrimAdapterRm),
		m_htrim(o.htrimLeft != "" || o.htrimRight != ""),
		m_twoBarcodes(o.barDetect == flexbar::WITHIN_READ_REMOVAL2 || o.barDetect == flexbar::WITHIN_READ2),
		out(o.out),
		m_unassigned(0){
		
		m_barcodes  = &o.barcodes;
		m_adapters  = &o.adapters;
		m_barcodes2 = &o.barcodes2;
		m_adapters2 = &o.adapters2;
		
		m_b1 = new TSeqAlign(m_barcodes,  o, o.b_min_overlap, o.b_errorRate, o.b_tail_len, o.b_match, o.b_mismatch, o.b_gapCost, true);
		m_b2 = new TSeqAlign(m_barcodes2, o, o.b_min_overlap, o.b_errorRate, o.b_tail_len, o.b_match, o.b_mismatch, o.b_gapCost, true);
		
		m_a1 = new TSeqAlign(m_adapters,  o, o.a_min_overlap, o.a_errorRate, o.a_tail_len, o.match, o.mismatch, o.gapCost, false);
		m_a2 = new TSeqAlign(m_adapters2, o, o.a_min_overlap, o.a_errorRate, o.a_tail_len, o.match, o.mismatch, o.gapCost, false);
		
		m_p  = new TSeqAlignPair(o, o.p_min_overlap, o.a_errorRate, o.match, o.mismatch, o.gapCost);
		
		if(m_log == flexbar::TAB)
		*out << "ReadTag\tQueryTag\tQueryStart\tQueryEnd\tOverlapLength\tMismatches\tIndels\tAllowedErrors" << std::endl;
	}
	
	
	virtual ~PairedAlign(){
		delete m_b1;
		delete m_b2;
		delete m_a1;
		delete m_a2;
		delete m_p;
	};
	
	
	void alignPairedReadToBarcodes(flexbar::TPairedRead* pRead, flexbar::TAlignBundle &alBundle, std::vector<flexbar::ComputeCycle> &cycle, std::vector<unsigned int> &idxAl, const flexbar::AlignmentMode &alMode){
		
		using namespace flexbar;
		
		switch(m_barType){
			case BARCODE_READ:         pRead->barID  = m_b1->alignSeqRead(pRead->b,  false, alBundle[0], cycle[0], idxAl[0], alMode, m_bTrimEnd); break;
			case WITHIN_READ_REMOVAL2: pRead->barID2 = m_b2->alignSeqRead(pRead->r2, true,  alBundle[2], cycle[2], idxAl[2], alMode, m_bTrimEnd);
			case WITHIN_READ_REMOVAL:  pRead->barID  = m_b1->alignSeqRead(pRead->r1, true,  alBundle[1], cycle[1], idxAl[1], alMode, m_bTrimEnd); break;
			case WITHIN_READ2:         pRead->barID2 = m_b2->alignSeqRead(pRead->r2, false, alBundle[2], cycle[2], idxAl[2], alMode, m_bTrimEnd);
			case WITHIN_READ:          pRead->barID  = m_b1->alignSeqRead(pRead->r1, false, alBundle[1], cycle[1], idxAl[1], alMode, m_bTrimEnd); break;
			case BOFF: break;
		}
		
		if(pRead->barID == 0 || (m_twoBarcodes && pRead->barID2 == 0)){
			
			if(cycle[0] != PRELOAD) m_unassigned++;
		}
	}
	
	
	void alignPairedReadToAdapters(flexbar::TPairedRead* pRead, flexbar::TAlignBundle &alBundle, std::vector<flexbar::ComputeCycle> &cycle, std::vector<unsigned int> &idxAl, const flexbar::AlignmentMode &alMode, const flexbar::TrimEnd trimEnd){
		
		using namespace flexbar;
		
		if(m_adapRem != ATWO)
			m_a1->alignSeqRead(pRead->r1, true, alBundle[0], cycle[0], idxAl[0], alMode, trimEnd);
		
		if(pRead->r2 != NULL && m_adapRem != AONE){
			if(m_adapRem != NORMAL2) m_a1->alignSeqRead(pRead->r2, true, alBundle[1], cycle[1], idxAl[1], alMode, trimEnd);
			else                     m_a2->alignSeqRead(pRead->r2, true, alBundle[1], cycle[1], idxAl[1], alMode, trimEnd);
		}
	}
	
	
	void trimLeftHPS(flexbar::TSeqRead* seqRead){
		
		using namespace std;
		using namespace flexbar;
		
		if(m_htrimAdapterRm && m_useRcTrimEnd){
			if     (seqRead->rmAdapter   && (m_aTrimEnd   == RIGHT || m_aTrimEnd   == RTAIL)) return;
			else if(seqRead->rmAdapterRC && (m_arcTrimEnd == RIGHT || m_arcTrimEnd == RTAIL)) return;
		}
		else if(m_htrimAdapterRm && ! m_useRcTrimEnd){
			if(m_aTrimEnd == RIGHT || m_aTrimEnd == RTAIL) return;
		}
		
		if(! m_htrimAdapterRm || seqRead->rmAdapter || seqRead->rmAdapterRC){
			
			for(unsigned int s = 0; s < m_htrimLeft.length(); ++s){
				
				char nuc = m_htrimLeft[s];
				
				unsigned int cutPos = 0;
				unsigned int notNuc = 0;
				
				for(unsigned int i = 0; i < length(seqRead->seq); ++i){
					
					if(seqRead->seq[i] != nuc){
						notNuc++;
					}
					else if(notNuc <= m_htrimErrorRate * (i+1)){
						
						if(m_htrimMaxLength != 0 && i+1 > m_htrimMaxLength && (! m_htrimMaxFirstOnly || s == 0)) break;
						
						cutPos = i+1;
					}
				}
				
				if(cutPos > 0 && cutPos >= m_htrimMinLength){
					erase(seqRead->seq, 0, cutPos);
					
					if(m_format == FASTQ){
						erase(seqRead->qual, 0, cutPos);
					}
				}
			}
		}
	}
	
	
	void trimRightHPS(flexbar::TSeqRead* seqRead){
		
		using namespace std;
		using namespace flexbar;
		
		if(m_htrimAdapterRm && m_useRcTrimEnd){
			if     (seqRead->rmAdapter   && (m_aTrimEnd   == LEFT || m_aTrimEnd   == LTAIL)) return;
			else if(seqRead->rmAdapterRC && (m_arcTrimEnd == LEFT || m_arcTrimEnd == LTAIL)) return;
		}
		else if(m_htrimAdapterRm && ! m_useRcTrimEnd){
			if(m_aTrimEnd == LEFT || m_aTrimEnd == LTAIL) return;
		}
		
		if(! m_htrimAdapterRm || seqRead->rmAdapter || seqRead->rmAdapterRC){
			
			for(unsigned int s = 0; s < m_htrimRight.length(); ++s){
				
				char nuc = m_htrimRight[s];
				
				unsigned int seqLen = length(seqRead->seq);
				unsigned int cutPos = seqLen;
				unsigned int notNuc = 0;
				
				for(int i = seqLen - 1; i >= 0; --i){
					
					if(seqRead->seq[i] != nuc){
						notNuc++;
					}
					else if(notNuc <= m_htrimErrorRate * (seqLen - i)){
						
						if(m_htrimMaxLength != 0 && i < seqLen - m_htrimMaxLength && (! m_htrimMaxFirstOnly || s == 0)) break;
						
						cutPos = i;
					}
				}
				
				if(cutPos < seqLen && cutPos <= seqLen - m_htrimMinLength){
					erase(seqRead->seq, cutPos, length(seqRead->seq));
					
					if(m_format == FASTQ){
						erase(seqRead->qual, cutPos, length(seqRead->qual));
					}
				}
			}
		}
	}
	
	
	// tbb filter operator
	void* operator()(void* item){
		
		using namespace flexbar;
		
		if(item != NULL){
			
			TPairedReadBundle *prBundle = static_cast<TPairedReadBundle* >(item);
			
			if(m_umiTags){
				for(unsigned int i = 0; i < prBundle->size(); ++i){
					prBundle->at(i)->r1->umi = "";
					
					if(prBundle->at(i)->r2 != NULL)
					prBundle->at(i)->r2->umi = "";
				}
			}
			
			AlignmentMode alMode = ALIGNALL;
			
			// barcode detection
			
			if(m_barType != BOFF){
				
				TAlignBundle alBundle;
				Alignments r1AlignmentsB, r2AlignmentsB, bAlignmentsB;
				
				alBundle.push_back(bAlignmentsB);
				alBundle.push_back(r1AlignmentsB);
				alBundle.push_back(r2AlignmentsB);
				
				std::vector<unsigned int> idxAl;
				std::vector<ComputeCycle> cycle;
				
				for(unsigned int i = 0; i < 3; ++i){
					idxAl.push_back(0);
					cycle.push_back(PRELOAD);
				}
				
				for(unsigned int i = 0; i < prBundle->size(); ++i){
					alignPairedReadToBarcodes(prBundle->at(i), alBundle, cycle, idxAl, alMode);
				}
				
				for(unsigned int i = 0; i < 3; ++i){
					idxAl[i] = 0;
					cycle[i] = COMPUTE;
				}
				
				for(unsigned int i = 0; i < prBundle->size(); ++i){
					alignPairedReadToBarcodes(prBundle->at(i), alBundle, cycle, idxAl, alMode);
				}
			}
			
			// adapter removal
			
			if(m_poMode != POFF){
				
				Alignments alignments;
				
				unsigned int idxAl = 0;
				ComputeCycle cycle = PRELOAD;
				
				for(unsigned int i = 0; i < prBundle->size(); ++i){
					m_p->alignSeqReadPair(prBundle->at(i)->r1, prBundle->at(i)->r2, alignments, cycle, idxAl);
				}
				
				idxAl = 0;
				cycle = COMPUTE;
				
				for(unsigned int i = 0; i < prBundle->size(); ++i){
					m_p->alignSeqReadPair(prBundle->at(i)->r1, prBundle->at(i)->r2, alignments, cycle, idxAl);
				}
			}
			
			if(m_adapRem != AOFF){
				
				for(unsigned int c = 0; c < m_arTimes; ++c){
					
					flexbar::TrimEnd trimEnd = m_aTrimEnd;
					
					unsigned int rc = 1;
					
					if(m_useRcTrimEnd){
						alMode = ALIGNRCOFF;
						rc = 2;
					}
					
					for(unsigned int r = 0; r < rc; ++r){
						
						if(m_useRcTrimEnd && r == 1){
							alMode  = ALIGNRC;
							trimEnd = m_arcTrimEnd;
						}
						
						TAlignBundle alBundle;
						Alignments r1AlignmentsA, r2AlignmentsA;
						
						alBundle.push_back(r1AlignmentsA);
						alBundle.push_back(r2AlignmentsA);
						
						std::vector<unsigned int> idxAl;
						std::vector<ComputeCycle> cycle;
						
						for(unsigned int i = 0; i < 2; ++i){
							idxAl.push_back(0);
							cycle.push_back(PRELOAD);
						}
						
						for(unsigned int i = 0; i < prBundle->size(); ++i){
							alignPairedReadToAdapters(prBundle->at(i), alBundle, cycle, idxAl, alMode, trimEnd);
						}
						
						for(unsigned int i = 0; i < 2; ++i){
							idxAl[i] = 0;
							cycle[i] = COMPUTE;
						}
						
						for(unsigned int i = 0; i < prBundle->size(); ++i){
							alignPairedReadToAdapters(prBundle->at(i), alBundle, cycle, idxAl, alMode, trimEnd);
						}
					}
				}
			}
			
			if(m_umiTags){
				for(unsigned int i = 0; i < prBundle->size(); ++i){
					
					append(prBundle->at(i)->r1->id, prBundle->at(i)->r1->umi);
					
					if(prBundle->at(i)->r2 != NULL){
						append(prBundle->at(i)->r1->id, prBundle->at(i)->r2->umi);
						
						append(prBundle->at(i)->r2->id, prBundle->at(i)->r1->umi);
						append(prBundle->at(i)->r2->id, prBundle->at(i)->r2->umi);
					}
				}
			}
			
			if(m_htrim){
				if(m_htrimLeft != ""){
					
					for(unsigned int i = 0; i < prBundle->size(); ++i){
						trimLeftHPS(prBundle->at(i)->r1);
						
						if(prBundle->at(i)->r2 != NULL)
						trimLeftHPS(prBundle->at(i)->r2);
					}
				}
				if(m_htrimRight != ""){
					
					for(unsigned int i = 0; i < prBundle->size(); ++i){
						trimRightHPS(prBundle->at(i)->r1);
						
						if(prBundle->at(i)->r2 != NULL)
						trimRightHPS(prBundle->at(i)->r2);
					}
				}
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
		
		if     (m_poMode  != POFF)    return m_p->getNrPreShortReads();
		else if(m_adapRem != NORMAL2) return m_a1->getNrPreShortReads();
		else                          return m_a1->getNrPreShortReads() + m_a2->getNrPreShortReads();
	}
	
	
	void printPairOverlapStats(){
		
		using namespace flexbar;
		
		if(m_p->getNrModifiedReads() > 0)
			*out << m_p->getOverlapStatsString() << "\n\n";
		
		if(m_adapRem == AOFF) *out << std::endl;
	}
	
	
	void printAdapterOverlapStats(){
		
		using namespace flexbar;
		
		if(m_a1->getNrModifiedReads() > 0)
			*out << m_a1->getOverlapStatsString() << "\n\n";
		
		if(m_adapRem != NORMAL2) *out << std::endl;
	}
	
	
	void printAdapterOverlapStats2(){
		
		if(m_a2->getNrModifiedReads() > 0)
			*out << m_a2->getOverlapStatsString() << "\n\n";
		
		*out << std::endl;
	}
	
};

#endif
