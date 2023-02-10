// LoadAdapters.h

#ifndef FLEXBAR_LOADADAPTERS_H
#define FLEXBAR_LOADADAPTERS_H


template <typename TSeqStr, typename TString>
class LoadAdapters {

private:
	
	std::ostream *out;
	oneapi::tbb::concurrent_vector<flexbar::TBar> adapters;
	
	flexbar::Adapters a;
	
	const flexbar::AdapterPreset m_aPreset;
	const flexbar::RevCompMode   m_rcMode;
	
public:
	
	LoadAdapters(const Options &o) :
		
		out(o.out),
		m_aPreset(o.aPreset),
		m_rcMode(o.rcMode){
		
		using namespace flexbar;
		
		// Illumina sequencing adapters
		// Oligonucleotide sequences © 2018 Illumina, Inc.  All rights reserved.
		// Obtained from https://support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html
		
		if(m_aPreset == TRUSEQ){
			a.id   = "TruSeq";
			a.seq1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
			a.seq2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";
			a.info = "TruSeq LT and TruSeq HT-based kits";
		}
		else if(m_aPreset == METHYL){
			a.id   = "TrueSeq-Methyl";
			a.seq1 = "AGATCGGAAGAGCACACGTCTGAAC";
			a.seq2 = "AGATCGGAAGAGCGTCGTGTAGGGA";
			a.info = "ScriptSeq and TruSeq DNA Methylation";
		}
		else if(m_aPreset == SMALLRNA){
			a.id   = "TrueSeq-smallRNA";
			a.seq1 = "TGGAATTCTCGGGTGCCAAGG";
			a.info = "TruSeq Small RNA";
		}
		else if(m_aPreset == RIBO){
			a.id   = "TrueSeq-Ribo";
			a.seq1 = "AGATCGGAAGAGCACACGTCT";
			a.info = "TruSeq Ribo Profile";
		}
		else if(m_aPreset == NEXTERA){
			a.id   = "Nextera-TruSight";
			a.seq1 = "CTGTCTCTTATACACATCT";
			a.info = "AmpliSeq, Nextera, Nextera DNA Flex, Nextera DNA, Nextera XT, Nextera Enrichment, Nextera Rapid Capture Enrichment, TruSight Enrichment, TruSight Rapid Capture Enrichment, TruSight HLA";
		}
		else if(m_aPreset == NEXTERAMP){
			a.id   = "Nextera-Matepair";
			a.seq1 = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
			a.seq2 = "GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";
			a.seqc = "CTGTCTCTTATACACATCT";
			a.info = "Nextera Mate Pair";
		}
		
		// IonTorrent sequencing adapters
		
		Adapters IonTorrent;
		IonTorrent.id   = "IonTorrent";
		IonTorrent.seq1 = "ATCACCGACTGCCCATAGAGAGGCTGAGAC";
		IonTorrent.seq1 = "CCATCTCATCCCTGCGTGTCTCCGACTCAG";
		IonTorrent.seq2 = "CCTCTCTATGGGCAGTCGGTGAT";
		IonTorrent.info = "IonTorrent";
	};
	
	
	virtual ~LoadAdapters(){};
	
	
	void loadSequences(const bool secondSet){
		
		using namespace std;
		using namespace flexbar;
		
		TString id = a.id;
		TSeqStr seq;
		
		if(! secondSet) seq = a.seq1;
		else            seq = a.seq2;
		
		if(m_rcMode == RCOFF || m_rcMode == RCON){
			TBar adapter;
			adapter.id  = id;
			adapter.seq = seq;
			adapters.push_back(adapter);
		}
		
		if(m_rcMode == RCON || m_rcMode == RCONLY){
			TString  idRC = id;
			TSeqStr seqRC = seq;
			
			append(idRC, "_rc");
			seqan::reverseComplement(seqRC);
			
			TBar adapterRC;
			adapterRC.id        = idRC;
			adapterRC.seq       = seqRC;
			adapterRC.rcAdapter = true;
			adapters.push_back(adapterRC);
		}
		
		if(m_aPreset == NEXTERAMP){
			TString  idc = id;
			TSeqStr seqc = a.seqc;
			
			append(idc, "_circ");
			
			TBar adapter;
			adapter.id  = idc;
			adapter.seq = seqc;
			adapters.push_back(adapter);
			
			append(idc, "_rc");
			seqan::reverseComplement(seqc);
			
			TBar adapterRC;
			adapterRC.id  = idc;
			adapterRC.seq = seqc;
			adapters.push_back(adapterRC);
		}
	};
	
	
	oneapi::tbb::concurrent_vector<flexbar::TBar> getAdapters(){
		return adapters;
	}
	
	
	void printAdapters(std::string adapterName) const {
		
		using namespace std;
		
		const unsigned int maxSpaceLen = 23;
		
		stringstream s; s << adapterName;
		int len = s.str().length() + 1;
		
		if(len + 2 > maxSpaceLen) len = maxSpaceLen - 2;
		
		*out << adapterName << ":" << string(maxSpaceLen - len, ' ') << "Sequence:" << "\n";
		
		for(unsigned int i=0; i < adapters.size(); ++i){
			TString seqTag = adapters.at(i).id;
			
			int whiteSpaceLen = maxSpaceLen - length(seqTag);
			if(whiteSpaceLen < 2) whiteSpaceLen = 2;
			
			string whiteSpace = string(whiteSpaceLen, ' ');
			
			*out << seqTag << whiteSpace << adapters.at(i).seq << "\n";
		}
		*out << endl;
	}
	
};

#endif
