/*
 *   SeqOutputFiles.h
 *
 *   Author: mat and jtr
 */

#ifndef FLEXBAR_SEQOUTPUTFILES_H
#define FLEXBAR_SEQOUTPUTFILES_H

#include "SeqOutput.h"


template <typename TSeqStr, typename TString>
class SeqOutputFiles {
	
public:
	
	typedef SeqOutput<TSeqStr, TString> TSeqOutput;
	
	TSeqOutput *f1, *f2, *single1, *single2;
	tbb::atomic<unsigned long> m_nShort_1, m_nShort_2;
	
	SeqOutputFiles() :
		f1(0),
		f2(0),
		single1(0),
		single2(0){
		
		m_nShort_1 = 0;
		m_nShort_2 = 0;
	};
	
	
	virtual ~SeqOutputFiles(){
    	delete f1;
    	delete f2;
    	delete single1;
    	delete single2;
	};
	

private:
	
	// forbid copying this object to call destructor only once
	// (pointing to unique objects)
	SeqOutputFiles(SeqOutputFiles&);
	
	SeqOutputFiles& operator =(const SeqOutputFiles& rhs);
	
};

#endif
