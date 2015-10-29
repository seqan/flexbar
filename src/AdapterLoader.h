/*
 *   AdapterLoader.h
 *
 *   Authors: mat and jtr
 */

#ifndef FLEXBAR_ADAPTERLOADER_H
#define FLEXBAR_ADAPTERLOADER_H

#include <sstream>
#include <string>

#include <tbb/pipeline.h>
#include <tbb/concurrent_vector.h>

#include <seqan/basic.h>

#include "Enums.h"
#include "Options.h"
#include "SeqRead.h"


template <typename TSeqStr, typename TString>
class AdapterLoader : public tbb::filter{

private:
	
	std::ostream *out;
	flexbar::FileFormat m_format;
	tbb::concurrent_vector<TAdapter> adapters;
	
	bool m_revComp, m_isAdapter;
	
public:
	
	AdapterLoader(const Options &o, const bool isAdapter) :
		
		filter(serial),
		out(o.out),
		m_format(o.format),
		m_isAdapter(isAdapter){
			
			m_revComp = o.revCompAdapter && isAdapter;
	};
	
	
	virtual ~AdapterLoader(){};
	
	
	void* operator()( void* item ){
		
		using namespace std;
		using namespace flexbar;
		
		SeqRead<TSeqStr, TString> *myRead = static_cast< SeqRead<TSeqStr, TString>* >(item);
		SeqRead<TSeqStr, TString> *myReadRC;
		
		TString tag = myRead->getSequenceTag();
		
		if(adapters.size() < 1000){
			for(int i = 0; i < adapters.size(); ++i){
				
				if(tag == adapters.at(i).first->getSequenceTag()){
					cerr << "Two ";
					
					if(m_isAdapter) cerr << "adapters";
					else            cerr << "barcodes";
					
					cerr << " have the same name.\n";
					cerr << "Please use unique names and restart.\n" << endl;
					
					exit(1);
				}
			}
		}
		
		if(m_revComp){
			TSeqStr seq = myRead->getSequence();
			seqan::reverseComplement(seq);
			
			append(tag, " revcomp");
			
			myReadRC = new SeqRead<TSeqStr, TString>(seq, tag);
		}
		
		TAdapter adap;
		adap.first = myRead;
		adapters.push_back(adap);
		
		if(m_revComp){
			TAdapter adapRC;
			adapRC.first = myReadRC;
			adapters.push_back(adapRC);
		}
		
		return NULL;
	};
	
	
	tbb::concurrent_vector<TAdapter> getAdapters(){
		return adapters;
	}
	
	
	void setAdapters(tbb::concurrent_vector<TAdapter> &adapterVec){
		adapters = adapterVec;
	}
	
	
	void printAdapters(std::string adapterName) const {
		using namespace std;
		
		const unsigned int maxSpaceLen = 23;
		
		stringstream s; s << adapterName;
		int len = s.str().length() + 1;
		
		if(len + 2 > maxSpaceLen) len = maxSpaceLen - 2;
		
		*out << adapterName << ":" << string(maxSpaceLen - len, ' ') << "Sequence:" << "\n";
		
		for(unsigned int i=0; i < adapters.size(); ++i){
			TSeqStr seqTag = adapters.at(i).first->getSequenceTag();
			
			int whiteSpaceLen = maxSpaceLen - length(seqTag);
			if(whiteSpaceLen < 2) whiteSpaceLen = 2;
			
			string whiteSpace = string(whiteSpaceLen, ' ');
			
			*out << seqTag << whiteSpace << adapters.at(i).first->getSequence() << "\n";
		}
		*out << endl;
	}
	
};

#endif
