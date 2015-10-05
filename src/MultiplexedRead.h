/*
 *   MultiplexedRead.h
 *
 *   Author: mat and jtr
 */

#ifndef FLEXBAR_MULTIPLEXEDREAD_H_
#define FLEXBAR_MULTIPLEXEDREAD_H_


template <typename TString, typename TIDString>
class MultiplexedRead {

public:
	
	typedef SequencingRead<TString, TIDString> TSequencingRead;
	
	TSequencingRead *m_r1;
	TSequencingRead *m_r2;
	TSequencingRead *m_b;
	
	TString m_randTag;
	int m_barcode_id, m_barcode_id2;
	
	MultiplexedRead(TSequencingRead *r1, TSequencingRead *r2, TSequencingRead *b) :
		m_r1(r1),
		m_r2(r2),
		m_b(b),
		m_barcode_id(0),
		m_barcode_id2(0),
		m_randTag(""){
	};
	
	virtual ~MultiplexedRead(){
		delete m_r1;
		delete m_r2;
		delete m_b;
	};
	
};

#endif /* FLEXBAR_MULTIPLEXEDREAD_H_ */
