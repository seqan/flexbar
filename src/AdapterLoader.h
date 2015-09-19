/*
 *   AdapterLoader.h
 *
 *   Authors: mat and jtr
 */

#ifndef FLEXBAR_ADAPTERLOADER_H_
#define FLEXBAR_ADAPTERLOADER_H_

#include <sstream>
#include <string>

#include <tbb/pipeline.h>
#include <tbb/concurrent_vector.h>

#include <seqan/basic.h>

#include "Enums.h"
#include "Options.h"
#include "SequencingRead.h"
#include "SequenceConverter.h"


// This class will store each processed read plus it's ID in a vector.

template <typename TString, typename TIDString>
class AdapterLoader : public tbb::filter{

private:
	
	std::ostream *out;
	flexbar::FileFormat m_format;
	tbb::concurrent_vector<TAdapter> adapters;
	
public:
	
	AdapterLoader(const Options &o) :
		
		filter(serial),
		m_format(o.format),
		out(o.out){
	};
	
	
	virtual ~AdapterLoader(){};
	
	
	void* operator()( void* item ){
		
		SequencingRead<TString, TIDString> *myRead = static_cast< SequencingRead<TString, TIDString>* >(item);
		
		if(m_format == flexbar::CSFASTA || m_format == flexbar::CSFASTQ){
			TString csRead = SequenceConverter<TString>::getInstance()->bpToColorSpace(myRead->getSequence());
			myRead->setSequence(csRead);
		}
		
		TAdapter adap;
		adap.first = myRead;
		adapters.push_back(adap);
		
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
		
		*out << adapterName << ":" << string(maxSpaceLen - 8, ' ') << "Sequence:" << "\n";
		
		for(unsigned int i=0; i < adapters.size(); ++i){
			TString seqTag = adapters.at(i).first->getSequenceTag();
			
			int whiteSpaceLen = maxSpaceLen - length(seqTag);
			if(whiteSpaceLen < 2) whiteSpaceLen = 2;
			
			string whiteSpace = string(whiteSpaceLen, ' ');
			
			*out << seqTag << whiteSpace << adapters.at(i).first->getSequence() << "\n";
		}
		*out << endl;
	}
	
};

#endif /* FLEXBAR_ADAPTERLOADER_H_ */
