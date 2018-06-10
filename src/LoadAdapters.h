// LoadAdapters.h

#ifndef FLEXBAR_LOADADAPTERS_H
#define FLEXBAR_LOADADAPTERS_H


template <typename TSeqStr, typename TString>
class LoadAdapters {

private:
	
	std::ostream *out;
	tbb::concurrent_vector<flexbar::TBar> adapters;
	
	flexbar::Adapters a;
	
	const flexbar::AdapterPreset m_aPreset;
	const flexbar::RevCompMode   m_rcMode;
	
public:
	
	LoadAdapters(const Options &o) :
		
		out(o.out),
		m_aPreset(o.aPreset),
		m_rcMode(o.rcMode){
		
		using namespace flexbar;
		
		// Oligonucleotide sequences © 2018 Illumina, Inc.  All rights reserved.
		// Obtained from https://support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html
		
		Adapters TrueSeq;
		TrueSeq.id   = "TruSeq";
		TrueSeq.seq1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
		TrueSeq.seq2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";
		TrueSeq.info = "TruSeq LT and TruSeq HT-based kits";
		
		Adapters TrueSeq_methyl;
		TrueSeq_methyl.id   = "TrueSeq-Methyl";
		TrueSeq_methyl.seq1 = "AGATCGGAAGAGCACACGTCTGAAC";
		TrueSeq_methyl.seq2 = "AGATCGGAAGAGCGTCGTGTAGGGA";
		TrueSeq_methyl.info = "ScriptSeq and TruSeq DNA Methylation";
		
		Adapters TrueSeq_smallRNA;
		TrueSeq_smallRNA.id   = "TrueSeq-smallRNA";
		TrueSeq_smallRNA.seq1 = "TGGAATTCTCGGGTGCCAAGG";
		TrueSeq_smallRNA.info = "TruSeq Small RNA";
		
		Adapters TrueSeq_ribo;
		TrueSeq_ribo.id   = "TrueSeq-Ribo";
		TrueSeq_ribo.seq1 = "AGATCGGAAGAGCACACGTCT";
		TrueSeq_ribo.info = "TruSeq Ribo Profile";
		
		Adapters Nextera_TruSight;
		Nextera_TruSight.id   = "Nextera-TruSight";
		Nextera_TruSight.seq1 = "CTGTCTCTTATACACATCT";
		Nextera_TruSight.info = "AmpliSeq, Nextera, Nextera DNA Flex, Nextera DNA, Nextera XT, Nextera Enrichment, Nextera Rapid Capture Enrichment, TruSight Enrichment, TruSight Rapid Capture Enrichment, TruSight HLA";
		
		Adapters Nextera_matepair;
		Nextera_matepair.id   = "Nextera-Matepair";
		Nextera_matepair.seq1 = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
		Nextera_matepair.seq2 = "GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";
		Nextera_matepair.seqc = "CTGTCTCTTATACACATCT";
		Nextera_matepair.info = "Nextera Mate Pair";
		
		Adapters IonTorrent;
		IonTorrent.id   = "IonTorrent";
		IonTorrent.seq1 = "ATCACCGACTGCCCATAGAGAGGCTGAGAC";
		IonTorrent.seq1 = "CCATCTCATCCCTGCGTGTCTCCGACTCAG";
		IonTorrent.seq2 = "CCTCTCTATGGGCAGTCGGTGAT";
		IonTorrent.info = "IonTorrent";
		
		     if(m_aPreset == TRUSEQ)    a = TrueSeq;
		else if(m_aPreset == SMALLRNA)  a = TrueSeq_smallRNA;
		else if(m_aPreset == METHYL)    a = TrueSeq_methyl;
		else if(m_aPreset == RIBO)      a = TrueSeq_ribo;
		else if(m_aPreset == NEXTERA)   a = Nextera_TruSight;
		else if(m_aPreset == NEXTERAMP) a = Nextera_matepair;
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
	
	
	tbb::concurrent_vector<flexbar::TBar> getAdapters(){
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
