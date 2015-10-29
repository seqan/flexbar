/*
 *   SeqRead.h
 *
 *   Author: mat and jtr
 */

#ifndef FLEXBAR_SEQREAD_H
#define FLEXBAR_SEQREAD_H


template <typename TSeqStr, typename TIDString>
class SeqRead {

private:
	
	TSeqStr m_seq;
	TIDString m_tag, m_qual;
	
public:
	
	SeqRead()
    	: m_tag(),
      	  m_seq(){
	}
	
	SeqRead(const TSeqStr& source, const TIDString& sequence_tag)
		: m_tag(sequence_tag),
		  m_seq(source){
	}
	
	SeqRead(const TSeqStr& source, const TIDString& sequence_tag, const TIDString& qual)
		: m_tag(sequence_tag),
		  m_seq(source),
		  m_qual(qual){
	}
	
	
	void setSequenceTag(const TSeqStr& tag){
		m_tag = tag;
	}
	
	void setSequence(const TSeqStr& seq){
		m_seq = seq;
	}
	
	void setQuality(const TSeqStr& qual){
		m_qual = qual;
	}
	
	
	const TIDString& getSequenceTag() const {
		return m_tag;
	}
	
	const TSeqStr& getSequence() const {
		return m_seq;
	}
	
	const TIDString& getQuality() const{
		return m_qual;
	}
	
	
	virtual ~SeqRead(){};
	
};

#endif
