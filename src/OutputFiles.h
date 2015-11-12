/*
 *   OutputFileStruct.h
 *
 *   Author: mat and jtr
 */

#ifndef FLEXBAR_OUTPUTFILESTRUCT_H
#define FLEXBAR_OUTPUTFILESTRUCT_H

#include "SeqOutputFilter.h"


template <typename TSeqStr, typename TString>
class OutputFileStruct {
	
public:
	
	typedef SeqOutputFilter<TSeqStr, TString> TOutputFilter;
	
	TOutputFilter *f1, *f2, *single1, *single2;
	
	tbb::atomic<unsigned long> m_nShort_1, m_nShort_2;
	
	OutputFileStruct() :
		f1(0),
		f2(0),
		single1(0),
		single2(0){
		
		m_nShort_1 = 0;
		m_nShort_2 = 0;
	};
	
	
	virtual ~OutputFileStruct(){
    	delete f1;
    	delete f2;
    	delete single1;
    	delete single2;
	};
	

private:
	
	// forbid copying this object to call destructor only once (pointing to unique objects)
	OutputFileStruct(OutputFileStruct&);
	
	OutputFileStruct& operator =(const OutputFileStruct& rhs);
	
};

#endif
