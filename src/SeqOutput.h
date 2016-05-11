/*
 *   SeqOutput.h
 *
 *   Authors: mat and jtr
 */

#ifndef FLEXBAR_SEQOUTPUT_H
#define FLEXBAR_SEQOUTPUT_H

#include <fstream>


template <typename TSeqStr, typename TString>
class SeqOutput {

private:
	
	seqan::SeqFileOut seqFileOut;
	
	const bool m_writeLenDist, m_useStdout;
	const unsigned int m_minLength, m_cutLen_read;
	
	const std::string m_filePath;
	const TString m_tagStr;
	
	const flexbar::FileFormat m_format;
	const flexbar::CompressionType m_cmprsType;
	
	tbb::atomic<unsigned long> m_countGood, m_countGoodChars;
	
	tbb::concurrent_vector<unsigned long> *m_lengthDist;
	
public:
	
	SeqOutput(const std::string &filePath, const TString tagStr, const bool alwaysFile, const Options &o) :
		m_format(o.format),
		m_tagStr(tagStr),
		m_minLength(o.min_readLen),
		m_cutLen_read(o.cutLen_read),
		m_writeLenDist(o.writeLengthDist),
		m_useStdout(o.useStdout && ! alwaysFile),
		m_cmprsType(o.cmprsType),
		m_filePath(filePath + o.outCompression){
		
		using namespace std;
		using namespace flexbar;
		
		m_countGood      = 0;
		m_countGoodChars = 0;
		
		m_lengthDist = new tbb::concurrent_vector<unsigned long>(MAX_READLENGTH + 1, 0);
		
		if(m_useStdout){
			
			if(m_format == FASTA){
				setFormat(seqFileOut, seqan::Fasta());
			}
			else if(m_format == FASTQ){
				setFormat(seqFileOut, seqan::Fastq());
			}
			
			if(! open(seqFileOut, cout)){
				cerr << "ERROR: Could not open output stream." << "\n" << endl;
				exit(1);
			}
		}
		else{
			if(! open(seqFileOut, m_filePath.c_str())){
				cerr << "ERROR: Could not open file: " << m_filePath << "\n" << endl;
				exit(1);
			}
		}
	};
	
	
	virtual ~SeqOutput(){
		delete m_lengthDist;
		if(! m_useStdout) close(seqFileOut);
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
			cerr << "ERROR: Could not open file: " << fname << "\n";
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
	
	
	void writeSeqRead(const SeqRead<TSeqStr, TString>& myRead){
		
		using namespace std;
		using namespace flexbar;
		
		TString tag = myRead.getSequenceTag();
		
		if(m_useStdout && m_tagStr != ""){
			append(tag, "_");
			append(tag, m_tagStr);
		}
		
		try{
			if(m_format == FASTA){
				writeRecord(seqFileOut, tag, myRead.getSequence());
			}
			else if(m_format == FASTQ){
				writeRecord(seqFileOut, tag, myRead.getSequence(), myRead.getQuality());
			}
		}
		catch(seqan::Exception const &e){
			cerr << "\n\n" << "ERROR: " << e.what() << "\nProgram execution aborted.\n" << endl;
			
			close(seqFileOut);
			delete m_lengthDist;
			exit(1);
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
			SeqRead<TSeqStr, TString> *myRead = static_cast< SeqRead<TSeqStr, TString>* >(item);
			
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
			
			writeSeqRead(*myRead);
		}
		
		return NULL;
	}
	
};

#endif
