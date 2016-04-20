/*
 *   PairedReadBundle.h
 *
 *   Author: jtr
 */

#ifndef FLEXBAR_PAIREDREADBUNDLE_H
#define FLEXBAR_PAIREDREADBUNDLE_H


template <typename TSeqStr, typename TString>
class PairedReadBundle {

public:
	
	typedef std::vector<PairedRead<TSeqStr, TString>* > TPairedReadBundle;
	
	TPairedReadBundle *m_bundle;
	
	PairedReadBundle(TPairedReadBundle *bundle) :
		m_bundle(bundle){
	};
	
	virtual ~PairedReadBundle(){
		
		for(unsigned int i = 0; i < m_bundle->size(); ++i){
			delete m_bundle->at(i);
		}
		delete m_bundle;
	};
	
};

#endif
