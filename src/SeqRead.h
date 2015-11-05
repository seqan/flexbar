/*
 *   SeqRead.h
 *
 *   Author: mat and jtr
 */

#ifndef FLEXBAR_SEQREAD_H
#define FLEXBAR_SEQREAD_H


template <typename TSeqStr, typename TString>
class SeqRead {

private:
	
	TSeqStr m_seq;
	TString m_tag, m_qual;
	
public:
	
	SeqRead()
    	: m_tag(),
      	  m_seq(){
	}
	
	SeqRead(const TSeqStr& source, const TString& sequence_tag)
		: m_tag(sequence_tag),
		  m_seq(source){
	}
	
	SeqRead(const TSeqStr& source, const TString& sequence_tag, const TString& qual)
		: m_tag(sequence_tag),
		  m_seq(source),
		  m_qual(qual){
	}
	
	
	void setSequenceTag(const TString& tag){
		m_tag = tag;
	}
	
	void setSequence(const TSeqStr& seq){
		m_seq = seq;
	}
	
	void setQuality(const TString& qual){
		m_qual = qual;
	}
	
	
	const TString& getSequenceTag() const {
		return m_tag;
	}
	
	const TSeqStr& getSequence() const {
		return m_seq;
	}
	
	const TString& getQuality() const{
		return m_qual;
	}
	
	
	virtual ~SeqRead(){};
	
};

#endif
