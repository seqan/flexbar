/*
 *   SeqRead.h
 *
 *   Author: mat and jtr
 */

#ifndef FLEXBAR_SEQREAD_H
#define FLEXBAR_SEQREAD_H


template <typename TString, typename TIDString>
class SeqRead {

private:
	
	TString m_seq;
	TIDString m_tag, m_qual;
	
public:
	
	SeqRead()
    	: m_tag(),
      	  m_seq(){
	}
	
	SeqRead(const TString& source, const TIDString& sequence_tag)
		: m_tag(sequence_tag),
		  m_seq(source){
	}
	
	SeqRead(const TString& source, const TIDString& sequence_tag, const TIDString& qual)
		: m_tag(sequence_tag),
		  m_seq(source),
		  m_qual(qual){
	}
	
	
	void setSequenceTag(const TString& tag){
		m_tag = tag;
	}
	
	void setSequence(const TString& seq){
		m_seq = seq;
	}
	
	void setQuality(const TString& qual){
		m_qual = qual;
	}
	
	
	const TIDString& getSequenceTag() const {
		return m_tag;
	}
	
	const TString& getSequence() const {
		return m_seq;
	}
	
	const TIDString& getQuality() const{
		return m_qual;
	}
	
	
	virtual ~SeqRead(){};
	
};

#endif
