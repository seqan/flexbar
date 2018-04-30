// PairedOutput.h

#ifndef FLEXBAR_PAIREDOUTPUT_H
#define FLEXBAR_PAIREDOUTPUT_H

#include "SeqOutput.h"
#include "SeqOutputFiles.h"
#include "QualTrimming.h"


template <typename TSeqStr, typename TString>
class PairedOutput : public tbb::filter {

private:
	
	int m_mapsize;
	const int m_minLength, m_qtrimThresh, m_qtrimWinSize;
	const bool m_isPaired, m_writeUnassigned, m_writeSingleReads, m_writeSingleReadsP;
	const bool m_twoBarcodes, m_qtrimPostRm, m_htrim, m_htrimAdapterRm;
	
	tbb::atomic<unsigned long> m_nSingleReads, m_nLowPhred;
	
	const std::string m_target;
	const std::string m_htrimLeft, m_htrimRight;
	
	const int m_htrimMinLength, m_htrimMaxLength;
	const float m_htrimErrorRate;
	
	const flexbar::FileFormat     m_format;
	const flexbar::RunType        m_runType;
	const flexbar::BarcodeDetect  m_barDetect;
	const flexbar::QualTrimType   m_qtrim;
	
	typedef SeqOutput<TSeqStr, TString>      TSeqOutput;
	typedef SeqOutputFiles<TSeqStr, TString> TOutFiles;
	
	TOutFiles *m_outMap;
	std::ostream *out;
	
	tbb::concurrent_vector<flexbar::TBar> *m_adapters,  *m_barcodes;
	tbb::concurrent_vector<flexbar::TBar> *m_adapters2, *m_barcodes2;
	
public:
	
	PairedOutput(Options &o) :
		
		filter(serial_in_order),
		m_target(o.targetName),
		m_format(o.format),
		m_runType(o.runType),
		m_barDetect(o.barDetect),
		m_minLength(o.min_readLen),
		m_htrimLeft(o.htrimLeft),
		m_htrimRight(o.htrimRight),
		m_htrimMinLength(o.htrimMinLength),
		m_htrimMaxLength(o.htrimMaxLength),
		m_htrimErrorRate(o.h_errorRate),
		m_htrimAdapterRm(o.htrimAdapterRm),
		m_htrim(o.htrimLeft != "" || o.htrimRight != ""),
		m_qtrim(o.qTrim),
		m_qtrimThresh(o.qtrimThresh),
		m_qtrimWinSize(o.qtrimWinSize),
		m_qtrimPostRm(o.qtrimPostRm),
		m_isPaired(o.isPaired),
		m_writeUnassigned(o.writeUnassigned),
		m_writeSingleReads(o.writeSingleReads),
		m_writeSingleReadsP(o.writeSingleReadsP),
		m_twoBarcodes(o.barDetect == flexbar::WITHIN_READ_REMOVAL2 || o.barDetect == flexbar::WITHIN_READ2),
		out(o.out){
		
		using namespace std;
		using namespace flexbar;
		
		m_barcodes  = &o.barcodes;
		m_barcodes2 = &o.barcodes2;
		m_adapters  = &o.adapters;
		m_adapters2 = &o.adapters2;
		
		m_mapsize      = 0;
		m_nSingleReads = 0;
		m_nLowPhred    = 0;
		
		switch(m_runType){
			
			case PAIRED_BARCODED:{
				
				int nBarcodes = m_barcodes->size();
				if(m_twoBarcodes) nBarcodes *= m_barcodes2->size();
				
				m_mapsize = nBarcodes + 1;
				m_outMap  = new TOutFiles[m_mapsize];
				
				for(int i = 0; i < nBarcodes; ++i){
					
					int idxB1 = i % m_barcodes->size();
					int idxB2 = div(i, m_barcodes->size()).quot;
					
					TString barcode = m_barcodes->at(idxB1).id;
					
					if(m_twoBarcodes){
						append(barcode, "-");
						append(barcode, m_barcodes2->at(idxB2).id);
					}
					
					TString barcode1 = barcode;
					TString barcode2 = barcode;
					
					append(barcode1, "_1");
					append(barcode2, "_2");
					
					stringstream b1, b2;
					b1 << barcode1;
					b2 << barcode2;
					
					string s = m_target + "_barcode_" + b1.str();
					TSeqOutput *of1 = new TSeqOutput(s, barcode1, false, o);
					
					s = m_target + "_barcode_" + b2.str();
					TSeqOutput *of2 = new TSeqOutput(s, barcode2, false, o);
					
					TOutFiles& f = m_outMap[i + 1];
					f.f1 = of1;
					f.f2 = of2;
					
					if(m_writeSingleReads){
						s = m_target + "_barcode_" + b1.str() + "_single";
						TSeqOutput *osingle1 = new TSeqOutput(s, "", true, o);
						
						s = m_target + "_barcode_" + b2.str() + "_single";
						TSeqOutput *osingle2 = new TSeqOutput(s, "", true, o);
						
						f.single1 = osingle1;
						f.single2 = osingle2;
					}
				}
				
				if(m_writeUnassigned){
					string s = m_target + "_barcode_unassigned_1";
					TSeqOutput *of1 = new TSeqOutput(s, "unassigned_1", false, o);
					
					s = m_target + "_barcode_unassigned_2";
					TSeqOutput *of2 = new TSeqOutput(s, "unassigned_2", false, o);
					
					TOutFiles& f = m_outMap[0];
					f.f1 = of1;
					f.f2 = of2;
					
					if(m_writeSingleReads){
						s = m_target + "_barcode_unassigned_1_single";
						TSeqOutput *osingle1 = new TSeqOutput(s, "", true, o);
						
						s = m_target + "_barcode_unassigned_2_single";
						TSeqOutput *osingle2 = new TSeqOutput(s, "", true, o);
						
						f.single1 = osingle1;
						f.single2 = osingle2;
					}
				}
				break;
			}
			
			case PAIRED:{
				
				m_mapsize = 1;
				m_outMap  = new TOutFiles[m_mapsize];
				
				string s = m_target + "_1";
				TSeqOutput *of1 = new TSeqOutput(s, "1", false, o);
				
				s = m_target + "_2";
				TSeqOutput *of2 = new TSeqOutput(s, "2", false, o);
				
				TOutFiles& f = m_outMap[0];
				f.f1 = of1;
				f.f2 = of2;
				
				if(m_writeSingleReads){
					s = m_target + "_1_single";
					TSeqOutput *osingle1 = new TSeqOutput(s, "", true, o);
					
					s = m_target + "_2_single";
					TSeqOutput *osingle2 = new TSeqOutput(s, "", true, o);
					
					f.single1 = osingle1;
					f.single2 = osingle2;
				}
				break;
			}
			
			case SINGLE:{
				
				m_mapsize = 1;
				m_outMap  = new TOutFiles[m_mapsize];
				
				string s = m_target;
				TSeqOutput *of1 = new TSeqOutput(s, "", false, o);
				
				TOutFiles& f = m_outMap[0];
				f.f1 = of1;
				
				break;
			}
			
			case SINGLE_BARCODED:{
				
				m_mapsize = m_barcodes->size() + 1;
				m_outMap  = new TOutFiles[m_mapsize];
				
				for(int i = 0; i < m_barcodes->size(); ++i){
					
					TString barcode = m_barcodes->at(i).id;
					
					stringstream b;
					b << barcode;
					
					string s = m_target + "_barcode_" + b.str();
					TSeqOutput *of1 = new TSeqOutput(s, barcode, false, o);
					
					TOutFiles& f = m_outMap[i + 1];
					f.f1 = of1;
				}
				
				if(m_writeUnassigned){
					string s = m_target + "_barcode_unassigned";
					TSeqOutput *of1 = new TSeqOutput(s, "unassigned", false, o);
					
					TOutFiles& f = m_outMap[0];
					f.f1 = of1;
				}
			}
		}
	}
	
	
	virtual ~PairedOutput(){
		delete[] m_outMap;
	};
	
	
	void trimLeftHPS(flexbar::TSeqRead* seqRead){
		
		using namespace std;
		using namespace flexbar;
		
		if(! m_htrimAdapterRm || seqRead->adapterRemoved){
			
			for(unsigned int s = 0; s < m_htrimLeft.length(); ++s){
				
				char nuc = m_htrimLeft[s];
				
				unsigned int cutPos = 0;
				unsigned int notNuc = 0;
				
				for(unsigned int i = 0; i < length(seqRead->seq); ++i){
					
					if(seqRead->seq[i] != nuc){
						notNuc++;
					}
					else if(notNuc <= m_htrimErrorRate * (i+1)){
						if(i+1 <= m_htrimMaxLength) cutPos = i+1;
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
		
		if(! m_htrimAdapterRm || seqRead->adapterRemoved){
			
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
						if(i >= seqLen - m_htrimMaxLength) cutPos = i;
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
	
	
	void writePairedRead(flexbar::TPairedRead* pRead){
		
		using namespace flexbar;
		
		bool l1ok = false, l2ok = false;
		
		switch(m_runType){
			
			case SINGLE:
			case SINGLE_BARCODED:{
				
				if(pRead->r1 != NULL){
					if(m_runType == SINGLE || m_writeUnassigned || pRead->barID > 0){
						
						if(m_htrim){
							if(m_htrimLeft  != "") trimLeftHPS(pRead->r1);
							if(m_htrimRight != "") trimRightHPS(pRead->r1);
						}
						
						if(m_qtrim != QOFF && m_qtrimPostRm){
							if(qualTrim(pRead->r1, m_qtrim, m_qtrimThresh, m_qtrimWinSize)) ++m_nLowPhred;
						}
						
						if(length(pRead->r1->seq) >= m_minLength){
							
							m_outMap[pRead->barID].f1->writeRead(pRead->r1);
						}
						else m_outMap[pRead->barID].m_nShort_1++;
					}
				}
				break;
			}
			
			case PAIRED:
			case PAIRED_BARCODED:{
				
				if(pRead->r1 != NULL && pRead->r2 != NULL){
					
					int outIdx = pRead->barID;
					
					if(m_twoBarcodes){
						if(outIdx == 0 || pRead->barID2 == 0){
							outIdx = 0;
						}
						else outIdx += (pRead->barID2 - 1) * m_barcodes->size();
					}
					
					if(m_runType == PAIRED || m_writeUnassigned || outIdx > 0){
						
						if(m_htrim){
							if(m_htrimLeft != ""){
								trimLeftHPS(pRead->r1);
								trimLeftHPS(pRead->r2);
							}
							if(m_htrimRight != ""){
								trimRightHPS(pRead->r1);
								trimRightHPS(pRead->r2);
							}
						}
						
						if(m_qtrim != QOFF && m_qtrimPostRm){
							if(qualTrim(pRead->r1, m_qtrim, m_qtrimThresh, m_qtrimWinSize)) ++m_nLowPhred;
							if(qualTrim(pRead->r2, m_qtrim, m_qtrimThresh, m_qtrimWinSize)) ++m_nLowPhred;
						}
						
						if(length(pRead->r1->seq) >= m_minLength) l1ok = true;
						if(length(pRead->r2->seq) >= m_minLength) l2ok = true;
						
						if(l1ok && l2ok){
							m_outMap[outIdx].f1->writeRead(pRead->r1);
							m_outMap[outIdx].f2->writeRead(pRead->r2);
						}
						else if(l1ok && ! l2ok){
							m_nSingleReads++;
							
							if(m_writeSingleReads){
								m_outMap[outIdx].single1->writeRead(pRead->r1);
							}
							else if(m_writeSingleReadsP){
								
								pRead->r2->seq = "N";
								
								if(m_format == FASTQ)
								pRead->r2->qual = prefix(pRead->r1->qual, 1);
								
								m_outMap[outIdx].f1->writeRead(pRead->r1);
								m_outMap[outIdx].f2->writeRead(pRead->r2);
							}
						}
						else if(! l1ok && l2ok){
							m_nSingleReads++;
							
							if(m_writeSingleReads){
								m_outMap[outIdx].single2->writeRead(pRead->r2);
							}
							else if(m_writeSingleReadsP){
								
								pRead->r1->seq = "N";
								
								if(m_format == FASTQ)
								pRead->r1->qual = prefix(pRead->r2->qual, 1);
								
								m_outMap[outIdx].f1->writeRead(pRead->r1);
								m_outMap[outIdx].f2->writeRead(pRead->r2);
							}
						}
						
						if(! l1ok) m_outMap[outIdx].m_nShort_1++;
						if(! l2ok) m_outMap[outIdx].m_nShort_2++;
					}
				}
			}
		}
	}
	
	
	// tbb filter operator
	void* operator()(void* item){
		
		using namespace flexbar;
		
		if(item != NULL){
			
			TPairedReadBundle *prBundle = static_cast< TPairedReadBundle* >(item);
			
			for(unsigned int i = 0; i < prBundle->size(); ++i){
				
				writePairedRead(prBundle->at(i));
				delete prBundle->at(i);
			}
			delete prBundle;
		}
		
		return NULL;
	}
	
	
	void writeLengthDist(){
		
		for(unsigned int i = 0; i < m_mapsize; i++){
			
			m_outMap[i].f1->writeLengthDist();
			
			if(m_outMap[i].f2 != NULL)
			m_outMap[i].f2->writeLengthDist();
		}
	}
	
	
	unsigned long getNrSingleReads() const {
		return m_nSingleReads;
	}
	
	
	unsigned long getNrLowPhredReads() const {
		return m_nLowPhred;
	}
	
	
	unsigned long getNrGoodReads(){
		using namespace flexbar;
		
		unsigned long nGood = 0;
		
		for(unsigned int i = 0; i < m_mapsize; i++){
			if(m_barDetect == BOFF || m_writeUnassigned || i > 0){
				
				nGood += m_outMap[i].f1->getNrGoodReads();
				
				if(m_outMap[i].f2 != NULL){
					nGood += m_outMap[i].f2->getNrGoodReads();
					
					if(m_writeSingleReads){
						nGood += m_outMap[i].single1->getNrGoodReads();
						nGood += m_outMap[i].single2->getNrGoodReads();
					}
				}
			}
		}
		return nGood;
	}
	
	
	unsigned long getNrGoodChars(){
		using namespace flexbar;
		
		unsigned long nGood = 0;
		
		for(unsigned int i = 0; i < m_mapsize; i++){
			if(m_barDetect == BOFF || m_writeUnassigned || i > 0){
				
				nGood += m_outMap[i].f1->getNrGoodChars();
				
				if(m_outMap[i].f2 != NULL){
					nGood += m_outMap[i].f2->getNrGoodChars();
					
					if(m_writeSingleReads){
						nGood += m_outMap[i].single1->getNrGoodChars();
						nGood += m_outMap[i].single2->getNrGoodChars();
					}
				}
			}
		}
		return nGood;
	}
	
	
	unsigned long getNrShortReads(){
		using namespace flexbar;
		
		unsigned long nShort = 0;
		
		for(unsigned int i = 0; i < m_mapsize; i++){
			if(m_barDetect == BOFF || m_writeUnassigned || i > 0){
				
				nShort += m_outMap[i].m_nShort_1;
				if(m_isPaired)
				nShort += m_outMap[i].m_nShort_2;
			}
		}
		return nShort;
	}
	
	
	void printAdapterRemovalStats(const bool secondSet){
		
		using namespace std;
		
		tbb::concurrent_vector<flexbar::TBar> *adapters;
		const unsigned int maxSpaceLen = 20;
		
		int startLen = 8;
		
		if(secondSet){
			adapters = m_adapters2;
			*out << "Adapter2";
			startLen++;
		}
		else{
			adapters = m_adapters;
			*out << "Adapter removal statistics\n";
			*out << "==========================\n";
			*out << "Adapter";
		}
		
		*out << ":" << string(maxSpaceLen -  startLen, ' ') << "Overlap removal:"
		            << string(maxSpaceLen -        16, ' ') << "Full length:\n";
		
		for(unsigned int i = 0; i < adapters->size(); i++){
			
			TString seqTag = adapters->at(i).id;
			
			int wsLen = maxSpaceLen - length(seqTag);
			if(wsLen < 2) wsLen = 2;
			string whiteSpace = string(wsLen, ' ');
			
			unsigned long nAdapOvl  = adapters->at(i).rmOverlap;
			unsigned long nAdapFull = adapters->at(i).rmFull;
			
			stringstream s; s << nAdapOvl;
			
			int wsLen2 = maxSpaceLen - s.str().length();
			if(wsLen2 < 2) wsLen2 = 2;
			string whiteSpace2 = string(wsLen2, ' ');
			
			*out << seqTag << whiteSpace << nAdapOvl << whiteSpace2 << nAdapFull << "\n";
		}
		*out << endl;
	}
	
	
	void printAdapterRemovalStats(){
		printAdapterRemovalStats(false);
	}
	
	
	void printAdapterRemovalStats2(){
		printAdapterRemovalStats(true);
	}
	
	
	void printFileSummary(){
		
		using namespace std;
		using namespace flexbar;
		
		*out << "Output file statistics\n";
		*out << "======================\n";
		
		for(unsigned int i = 0; i < m_mapsize; i++){
			
			if(m_barDetect == BOFF || m_writeUnassigned || i > 0){
				*out << "Read file:               " << m_outMap[i].f1->getFileName()    << "\n";
				*out << "  written reads          " << m_outMap[i].f1->getNrGoodReads() << "\n";
				*out << "  short reads            " << m_outMap[i].m_nShort_1           << "\n";
				
				if(m_isPaired){
					*out << "Read file 2:             " << m_outMap[i].f2->getFileName()    << "\n";
					*out << "  written reads          " << m_outMap[i].f2->getNrGoodReads() << "\n";
					*out << "  short reads            " << m_outMap[i].m_nShort_2           << "\n";
					
					if(m_writeSingleReads){
						*out << "Single read file:        " << m_outMap[i].single1->getFileName()    << "\n";
						*out << "  written reads          " << m_outMap[i].single1->getNrGoodReads() << "\n";
						*out << "Single read file 2:      " << m_outMap[i].single2->getFileName()    << "\n";
						*out << "  written reads          " << m_outMap[i].single2->getNrGoodReads() << "\n";
					}
				}
				*out << endl;
			}
		}
		*out << endl;
	}
	
};

#endif
