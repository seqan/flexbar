// SeqOutput.h

#ifndef FLEXBAR_SEQOUTPUT_H
#define FLEXBAR_SEQOUTPUT_H


template <typename TSeqStr, typename TString>
class SeqOutput {

private:
	
	seqan::FlexbarReadsSeqFileOut seqFileOut;
	std::string m_filePath;
	
	const TString m_tagStr;
	const flexbar::FileFormat m_format;
	const flexbar::CompressionType m_cmprsType;
	const bool m_switch2Fasta, m_writeLenDist, m_useStdout;
	const unsigned int m_minLength, m_cutLen_read;
	
	tbb::atomic<unsigned long> m_countGood, m_countGoodChars;
	tbb::concurrent_vector<unsigned long> m_lengthDist;
	
public:
	
	SeqOutput(const std::string &filePath, const TString tagStr, const bool alwaysFile, const Options &o) :
		m_format(o.format),
		m_switch2Fasta(o.switch2Fasta),
		m_tagStr(tagStr),
		m_minLength(o.min_readLen),
		m_cutLen_read(o.cutLen_read),
		m_writeLenDist(o.writeLengthDist),
		m_useStdout(o.useStdout && ! alwaysFile),
		m_cmprsType(o.cmprsType),
		m_countGood(0),
		m_countGoodChars(0){
		
		using namespace std;
		using namespace flexbar;
		
		m_filePath = filePath;
		
		if(filePath != o.outReadsFile && filePath != o.outReadsFile2){
			
			if(m_format == FASTA || m_switch2Fasta)
			     m_filePath += getExtension(FASTA);
			else m_filePath += getExtension(FASTQ);
		}
		m_filePath += o.outCompression;
		
		m_lengthDist = tbb::concurrent_vector<unsigned long>(MAX_READLENGTH + 1, 0);
		
		if(m_useStdout){
			
			if(m_format == FASTA || m_switch2Fasta)
			     setFormat(seqFileOut, seqan::Fasta());
			else setFormat(seqFileOut, seqan::Fastq());
			
			if(! open(seqFileOut, cout)){
				cerr << "\nERROR: Could not open output stream." << "\n" << endl;
				exit(1);
			}
		}
		else{
			if(! open(seqFileOut, m_filePath.c_str())){
				cerr << "\nERROR: Could not open file " << m_filePath << "\n" << endl;
				exit(1);
			}
		}
	};
	
	
	virtual ~SeqOutput(){
		if(! m_useStdout) close(seqFileOut);
	};
	
	
	const std::string getFileName(){
		if(! m_useStdout) return m_filePath;
		else              return "stdout";
	}
	
	
	void writeLengthDist(){
		
		using namespace std;
		
		string fname = m_filePath + ".lengthdist";
		fstream lstream;
		
		lstream.open(fname.c_str(), ios::out | ios::binary);
		
		if(! lstream.is_open()){
			cerr << "\nERROR: Could not open file " << fname << "\n";
		}
		else{
			lstream << "Readlength\tCount" << "\n";
			
			for (int i = 0; i <= flexbar::MAX_READLENGTH; ++i){
				if(m_lengthDist.at(i) > 0)
					lstream << i << "\t" << m_lengthDist.at(i) << "\n";
			}
			lstream.close();
		}
	}
	
	
	void writeSeqRead(flexbar::TSeqRead &seqRead){
		
		using namespace std;
		using namespace flexbar;
		
		if(m_useStdout && m_tagStr != ""){
			append(seqRead.id, "_");
			append(seqRead.id, m_tagStr);
		}
		
		try{
			if(m_format == FASTA || m_switch2Fasta){
				writeRecord(seqFileOut, seqRead.id, seqRead.seq);
			}
			else{
				writeRecord(seqFileOut, seqRead.id, seqRead.seq, seqRead.qual);
			}
		}
		catch(seqan::Exception const &e){
			cerr << "\nERROR: " << e.what() << "\nProgram execution aborted.\n" << endl;
			close(seqFileOut);
			exit(1);
		}
	}
	
	
	unsigned long getNrGoodReads() const {
		return m_countGood;
	}
	
	
	unsigned long getNrGoodChars() const {
		return m_countGoodChars;
	}
	
	
	void* writeRead(void* item){
		
		using namespace std;
		using namespace flexbar;
		
		if(item){
			SeqRead<TSeqStr, TString> *seqRead = static_cast< SeqRead<TSeqStr, TString>* >(item);
			
			unsigned int readLength = length(seqRead->seq);
			
			if(m_cutLen_read > 1 && m_cutLen_read >= m_minLength && m_cutLen_read < readLength){
				
				seqRead->seq = prefix(seqRead->seq, m_cutLen_read);
				
				if(m_format == FASTQ)
				seqRead->qual = prefix(seqRead->qual, m_cutLen_read);
				
				readLength = m_cutLen_read;
			}
			
			m_countGoodChars += readLength;
			
			++m_countGood;
			
			// store read length distribution
			
			if(m_writeLenDist && readLength <= MAX_READLENGTH)
				m_lengthDist.at(readLength)++;
			else if(m_writeLenDist)
				cerr << "\nCompile Flexbar with larger max read length to get correct length dist.\n" << endl;
			
			writeSeqRead(*seqRead);
		}
		
		return NULL;
	}
	
};

#endif
