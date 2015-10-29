/*
 *   SequenceOutputFilter.h
 *
 *   Authors: mat and jtr
 */

#ifndef FLEXBAR_SEQUENCEOUTPUTFILTER_H
#define FLEXBAR_SEQUENCEOUTPUTFILTER_H

#include <fstream>

#include <tbb/concurrent_vector.h>

#include "Enums.h"
#include "FlexbarIO.h"
#include "SequencingRead.h"


template <typename TString, typename TIDString, typename TStream>
class SequenceOutputFilter {

private:
	
	TStream m_targetStream;
	
	const bool m_writeLenDist, m_useStdout;
	const unsigned int m_minLength, m_cutLen_read;
	
	const std::string m_filePath;
	const TIDString m_tagStr;
	
	const flexbar::FileFormat m_format;
	const flexbar::CompressionType m_cmprsType;
	
	tbb::atomic<unsigned long> m_countGood, m_countGoodChars;
	
	tbb::concurrent_vector<unsigned long> *m_lengthDist;
	
public:
	
	SequenceOutputFilter(const std::string &filePath, const TIDString tagStr, const bool alwaysFile, const Options &o) :
		m_format(o.format),
		m_tagStr(tagStr),
		m_minLength(o.min_readLen),
		m_cutLen_read(o.cutLen_read),
		m_writeLenDist(o.writeLengthDist),
		m_useStdout(o.useStdout && ! alwaysFile),
		m_cmprsType(o.cmprsType),
		m_filePath(filePath + o.outCompression){
		
		using namespace flexbar;
		
		m_countGood      = 0;
		m_countGoodChars = 0;
		
		m_lengthDist = new tbb::concurrent_vector<unsigned long>(MAX_READLENGTH + 1, 0);
		
		if(! m_useStdout) openOutputFile(m_targetStream, m_filePath);
		
		// if(m_useStdout && m_cmprsType != UNCOMPRESSED) openOutputFile(m_targetStream, "-");
		// else if(! m_useStdout)                         openOutputFile(m_targetStream, m_filePath);
	};
	
	
	virtual ~SequenceOutputFilter(){
		if(! m_useStdout) closeFile(m_targetStream);
		delete m_lengthDist;
	};
	
	
	const std::string getFileName() const {
		if(! m_useStdout) return m_filePath;
		else              return "stdout";
	}
	
	
	void writeLengthDist() const {
		using namespace std;
		
		string fname = m_filePath + ".lengthdist";
		fstream lstream;
		
		lstream.open(fname.c_str(), ios::out | ios::binary);
		
		if(! lstream.is_open()){
			cerr << "Error opening File: " << fname << "\n";
		}
		else{
			lstream << "Readlength\tCount" << "\n";
			
			for (int i = 0; i <= flexbar::MAX_READLENGTH; ++i){
				if(m_lengthDist->at(i) > 0)
					lstream << i << "\t" << m_lengthDist->at(i) << "\n";
			}
			lstream.close();
		}
	}
	
	
	void writeFastString(const SequencingRead<TString, TIDString>& myRead){
		
		using namespace std;
		using namespace flexbar;
		
		seqan::CharString s = "";
		
		switch(m_format){
			case FASTQ:
				append(s, "@");
				append(s, myRead.getSequenceTag());
				
				if(m_useStdout && m_tagStr != ""){
					append(s, "_");
					append(s, m_tagStr);
				}
				append(s, "\n");
				
				append(s, myRead.getSequence());
				append(s, "\n+\n");
				append(s, myRead.getQuality());
				append(s, "\n");
			break;
			
			case FASTA:
				append(s, ">");
				append(s, myRead.getSequenceTag());
				
				if(m_useStdout && m_tagStr != ""){
					append(s, "_");
					append(s, m_tagStr);
				}
				append(s, "\n");
				
				append(s, myRead.getSequence());
				append(s, "\n");
		}
		
		// if(m_useStdout && m_cmprsType == UNCOMPRESSED) cout << s;
		if(m_useStdout) cout << s;
		else{
			if(streamPut(m_targetStream, s) != 0){
				cerr << "File writing error occured!\n" << endl;
				exit(1);
			}
		}
	}
	
	
	unsigned long getNrGoodReads() const {
		return m_countGood;
	}
	
	
	unsigned long getNrGoodChars() const {
		return m_countGoodChars;
	}
	
	
	void *writeRead(void *item){
		
		using namespace std;
		using namespace flexbar;
		
		if(item){
			SequencingRead<TString, TIDString> *myRead = static_cast< SequencingRead<TString, TIDString>* >(item);
			
			unsigned int readLength = length(myRead->getSequence());
			
			if(m_cutLen_read > 1 && m_cutLen_read >= m_minLength && m_cutLen_read < readLength){
				
				myRead->setSequence(prefix(myRead->getSequence(), m_cutLen_read));
				
				if(m_format == FASTQ){
					myRead->setQuality(prefix(myRead->getQuality(), m_cutLen_read));
				}
				
				readLength = m_cutLen_read;
			}
			
			m_countGoodChars += readLength;
			
			++m_countGood;
			
			// store read length distribution
			if(m_writeLenDist && readLength <= MAX_READLENGTH)
				m_lengthDist->at(readLength)++;
			else if(m_writeLenDist)
				cerr << "\nCompile Flexbar with larger max read length to get correct length dist.\n" << endl;
			
			writeFastString(*myRead);
		}
		
		return NULL;
	}
	
};

#endif
