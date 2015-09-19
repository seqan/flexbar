/*
 *   MultiplexedOutputFilter.h
 *
 *   Authors: mat and jtr
 */

#ifndef FLEXBAR_MULTIPLEXEDOUTPUTFILTER_H_
#define FLEXBAR_MULTIPLEXEDOUTPUTFILTER_H_

#include <tbb/pipeline.h>
#include <tbb/concurrent_vector.h>

#include <seqan/basic.h>

#include "Enums.h"
#include "Options.h"
#include "FlexbarIO.h"
#include "MultiplexedRead.h"
#include "SequenceOutputFilter.h"
#include "OutputFileStruct.h"
#include "AdapterLoader.h"


/* This class will process a MultiplexedRead and write it to a file
   depending on the runtype: single-end, paired-end and/or barcoded. */

template <typename TString, typename TIDString, typename TStream>
class MultiplexedOutputFilter : public tbb::filter {

private:
	
	int m_mapsize;
	const int m_minLength, m_cutLen_read;
	const bool m_isPaired, m_writeUnassigned, m_writeSingleReads;
	
	tbb::atomic<unsigned long> m_nSingleReads;
	
	const std::string m_target;
	
	const flexbar::FileFormat     m_format;
	const flexbar::RunType        m_runType;
	const flexbar::BarcodeDetect  m_barDetect;
	
	typedef SequenceOutputFilter<TString, TIDString, TStream> TOutputFilter;
	typedef OutputFileStruct<TString, TIDString, TStream> filters;
	
	filters *m_outputMap;
	std::ostream *out;
	
	tbb::concurrent_vector<TAdapter> *m_adapters, *m_barcodes;
	
public:
	
	MultiplexedOutputFilter(Options &o) :
		
		filter(serial_in_order),
		m_target(o.targetName),
		m_format(o.format),
		m_runType(o.runType),
		m_barDetect(o.barDetect),
		m_minLength(o.min_readLen),
		m_cutLen_read(o.cutLen_read),
		m_isPaired(o.isPaired),
		m_writeUnassigned(o.writeUnassigned),
		m_writeSingleReads(o.writeSingleReads),
		out(o.out){
		
		using namespace std;
		using namespace flexbar;
		
		m_adapters = &o.adapters;
		m_barcodes = &o.barcodes;
		
		m_nSingleReads = 0;
		
		m_mapsize = 0;
		
		switch(m_runType){
			
			case PAIRED_BARCODED:{
				
				m_mapsize = m_barcodes->size() + 1;
				m_outputMap = new filters[m_mapsize];
				
				stringstream ss;
				
				for(unsigned int i = 0; i < m_barcodes->size(); ++i){
					
					ss << m_target << "_barcode_" << m_barcodes->at(i).first->getSequenceTag() << "_1" << toFormatString(m_format);
					TOutputFilter *of1 = new TOutputFilter(ss.str(), "", o);
					ss.str("");
					ss.clear();
					
					ss << m_target << "_barcode_" << m_barcodes->at(i).first->getSequenceTag() << "_2"<< toFormatString(m_format);
					TOutputFilter *of2 = new TOutputFilter(ss.str(), "", o);
					ss.str("");
					ss.clear();
					
					filters& f = m_outputMap[i + 1];
					f.f1       = of1;
					f.f2       = of2;
					
					if(m_writeSingleReads){
						ss << m_target << "_barcode_" << m_barcodes->at(i).first->getSequenceTag() << "_1_single" << toFormatString(m_format);
						TOutputFilter *osingle1 = new TOutputFilter(ss.str(), "", o);
						ss.str("");
						ss.clear();
						
						ss << m_target << "_barcode_" << m_barcodes->at(i).first->getSequenceTag() << "_2_single"<< toFormatString(m_format);
						TOutputFilter *osingle2 = new TOutputFilter(ss.str(), "", o);
						ss.str("");
						ss.clear();
						
						f.single1 = osingle1;
						f.single2 = osingle2;
					}
				}
				
				if(m_writeUnassigned){
					ss << m_target << "_barcode_unassigned_1"<< toFormatString(m_format);
					TOutputFilter *of1 = new TOutputFilter(ss.str(), "", o);
					ss.str("");
					ss.clear();
					
					ss << m_target << "_barcode_unassigned_2"<< toFormatString(m_format);
					TOutputFilter *of2 = new TOutputFilter(ss.str(), "", o);
					ss.str("");
					ss.clear();
					
					filters& f = m_outputMap[0];
					f.f1       = of1;
					f.f2       = of2;
					
					if(m_writeSingleReads){
						ss << m_target << "_barcode_unassigned_1_single"<< toFormatString(m_format);
						TOutputFilter *osingle1 = new TOutputFilter(ss.str(), "", o);
						ss.str("");
						ss.clear();
						
						ss << m_target << "_barcode_unassigned_2_single"<< toFormatString(m_format);
						TOutputFilter *osingle2 = new TOutputFilter(ss.str(), "", o);
						
						f.single1 = osingle1;
						f.single2 = osingle2;
					}
				}
				break;
			}
			
			case PAIRED:{
				
				m_mapsize = 1;
				m_outputMap = new filters[m_mapsize];
				stringstream ss;
				
				ss << m_target << "_1"<< toFormatString(m_format);
				TOutputFilter *of1 = new TOutputFilter(ss.str(), "", o);
				ss.str("");
				ss.clear();
				
				ss << m_target << "_2"<< toFormatString(m_format);
				TOutputFilter *of2 = new TOutputFilter(ss.str(), "", o);
				ss.str("");
				ss.clear();
				
				filters& f = m_outputMap[0];
				f.f1       = of1;
				f.f2       = of2;
				
				if(m_writeSingleReads){
					ss << m_target << "_1_single" << toFormatString(m_format);
					TOutputFilter *osingle1 = new TOutputFilter(ss.str(), "", o);
					ss.str("");
					ss.clear();
					
					ss << m_target << "_2_single"<< toFormatString(m_format);
					TOutputFilter *osingle2 = new TOutputFilter(ss.str(), "", o);
					
					f.single1 = osingle1;
					f.single2 = osingle2;
				}
				break;
			}
			
			case SINGLE:{
				
				m_mapsize = 1;
				m_outputMap = new filters[m_mapsize];
				
				stringstream ss;
				ss << m_target << toFormatString(m_format);
				TOutputFilter *of1 = new TOutputFilter(ss.str(), "", o);
				
				filters& f = m_outputMap[0];
				f.f1 = of1;
				
				break;
			}
			
			case SINGLE_BARCODED:{
				
				m_mapsize = m_barcodes->size() + 1;
				m_outputMap = new filters[m_mapsize];
				
				for(unsigned int i=0; i < m_barcodes->size(); ++i){
					
					TIDString barcode = m_barcodes->at(i).first->getSequenceTag();
					
					stringstream ss;
					ss << m_target << "_barcode_" << barcode << toFormatString(m_format);
					TOutputFilter *of1 = new TOutputFilter(ss.str(), barcode, o);
					
					filters& f = m_outputMap[i + 1];
					f.f1 = of1;
				}
				
				if(m_writeUnassigned){
					stringstream ss;
					ss << m_target << "_barcode_unassigned" << toFormatString(m_format);
					TOutputFilter *of1 = new TOutputFilter(ss.str(), "unassigned", o);
					
					filters& f = m_outputMap[0];
					f.f1 = of1;
				}
			}
		}
	}
	
	
	virtual ~MultiplexedOutputFilter(){
		delete[] m_outputMap;
	};
	
	
	void* operator()(void* item) {
		
		using namespace flexbar;
		
		MultiplexedRead<TString, TIDString> *myRead = static_cast< MultiplexedRead<TString, TIDString>* >(item);
		
		bool l1ok = false, l2ok = false;
		
		switch(m_runType){
			
			case SINGLE:
			case SINGLE_BARCODED:{
				
				if(myRead->m_r1 != NULL){
					if(m_runType == SINGLE || m_writeUnassigned || myRead->m_barcode_id > 0){
						
						if(length(myRead->m_r1->getSequence()) >= m_minLength){
							
							m_outputMap[myRead->m_barcode_id].f1->writeRead(myRead->m_r1);
						}
						else m_outputMap[myRead->m_barcode_id].m_nShort_1++;
					}
				}
				break;
			}
			
			case PAIRED:
			case PAIRED_BARCODED:{
				
				if(myRead->m_r1 != NULL && myRead->m_r2 != NULL){
					if(m_runType == PAIRED || m_writeUnassigned || myRead->m_barcode_id > 0){
						
						// now check if both reads have min length
						if(length(myRead->m_r1->getSequence()) >= m_minLength) l1ok = true;
						if(length(myRead->m_r2->getSequence()) >= m_minLength) l2ok = true;
						
						if(l1ok && l2ok){
							m_outputMap[myRead->m_barcode_id].f1->writeRead(myRead->m_r1);
							m_outputMap[myRead->m_barcode_id].f2->writeRead(myRead->m_r2);
						}
						else if(l1ok && ! l2ok){
							m_nSingleReads++;
							
							if(m_writeSingleReads){
								m_outputMap[myRead->m_barcode_id].single1->writeRead(myRead->m_r1);
							}
						}
						else if(! l1ok && l2ok){
							m_nSingleReads++;
							
							if(m_writeSingleReads){
								m_outputMap[myRead->m_barcode_id].single2->writeRead(myRead->m_r2);
							}
						}
						
						if(! l1ok) m_outputMap[myRead->m_barcode_id].m_nShort_1++;
						if(! l2ok) m_outputMap[myRead->m_barcode_id].m_nShort_2++;
					}
				}
			}
		}
		
		delete myRead;
		
		return NULL;
	}
	
	
	void writeLengthDist(){
		
		for(unsigned int i = 0; i < m_mapsize; i++){
			m_outputMap[i].f1->writeLengthDist();
			if(m_outputMap[i].f2 != NULL)
			m_outputMap[i].f2->writeLengthDist();
		}
	}
	
	
	unsigned long getNrSingleReads() const {
		return m_nSingleReads;
	}
	
	
	unsigned long getNrGoodReads(){
		using namespace flexbar;
		
		unsigned long nGood = 0;
		
		for(unsigned int i = 0; i < m_mapsize; i++){
			if(m_barDetect == BOFF || m_writeUnassigned || i > 0){
				
				nGood += m_outputMap[i].f1->getNrGoodReads();
				
				if(m_outputMap[i].f2 != NULL){
					nGood += m_outputMap[i].f2->getNrGoodReads();
					
					if(m_writeSingleReads){
						nGood += m_outputMap[i].single1->getNrGoodReads();
						nGood += m_outputMap[i].single2->getNrGoodReads();
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
				
				nShort += m_outputMap[i].m_nShort_1;
				if(m_isPaired)
				nShort += m_outputMap[i].m_nShort_2;
			}
		}
		return nShort;
	}
	
	
	void printAdapterRemovalStats(){
		using namespace std;
		
		*out << "Adapter removal statistics\n";
		*out << "==========================\n";
		
		const unsigned int maxSpaceLen = 20;
		
		*out << "Adapter:" << string(maxSpaceLen -  8, ' ') << "Overlap removal:"
		                   << string(maxSpaceLen - 16, ' ') << "Full length:" << "\n";
		
		for(unsigned int i = 0; i < m_adapters->size(); i++){
			
			seqan::CharString seqTag = m_adapters->at(i).first->getSequenceTag();
			
			int wsLen = maxSpaceLen - length(seqTag);
			if(wsLen < 2) wsLen = 2;
			string whiteSpace = string(wsLen, ' ');
			
			unsigned long nAdapOvl  = m_adapters->at(i).second.first;
			unsigned long nAdapFull = m_adapters->at(i).second.second;
			
			stringstream ss;  ss << nAdapOvl;
			
			int wsLen2 = maxSpaceLen - ss.str().length();
			if(wsLen2 < 2) wsLen2 = 2;
			string whiteSpace2 = string(wsLen2, ' ');
			
			*out << seqTag << whiteSpace << nAdapOvl << whiteSpace2 << nAdapFull << "\n";
		}
		*out << endl;
	}
	
	
	void printFileSummary(){
		
		using namespace std;
		using namespace flexbar;
		
		*out << "Output file statistics\n";
		*out << "======================\n";
		
		for(unsigned int i = 0; i < m_mapsize; i++){
			
			if(m_barDetect == BOFF || m_writeUnassigned || i > 0){
				*out << "Read file:               " << m_outputMap[i].f1->getFileName()    << "\n";
				*out << "  written reads          " << m_outputMap[i].f1->getNrGoodReads() << "\n";
				*out << "  skipped short reads    " << m_outputMap[i].m_nShort_1           << "\n";
				
				if(m_isPaired){
					*out << "Read file 2:             " << m_outputMap[i].f2->getFileName()    << "\n";
					*out << "  written reads          " << m_outputMap[i].f2->getNrGoodReads() << "\n";
					*out << "  too short reads        " << m_outputMap[i].m_nShort_2           << "\n";
					
					if(m_writeSingleReads){
						*out << "Single read file:        " << m_outputMap[i].single1->getFileName()    << "\n";
						*out << "  written reads          " << m_outputMap[i].single1->getNrGoodReads() << "\n";
						*out << "Single read file 2:      " << m_outputMap[i].single2->getFileName()    << "\n";
						*out << "  written reads          " << m_outputMap[i].single2->getNrGoodReads() << "\n";
					}
				}
				*out << endl;
			}
		}
		*out << endl;
	}
	
};

#endif /* FLEXBAR_MULTIPLEXEDOUTPUTFILTER_H_ */
