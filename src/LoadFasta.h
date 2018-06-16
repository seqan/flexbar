// LoadFasta.h

#ifndef FLEXBAR_LOADFASTA_H
#define FLEXBAR_LOADFASTA_H


template <typename TSeqStr, typename TString>
class LoadFasta {

private:
	
	std::ostream *out;
	tbb::concurrent_vector<flexbar::TBar> bars;
	
	const bool m_isAdapter;
	const flexbar::RevCompMode m_rcMode;
	const flexbar::InputFileEnd m_inFileEnd;
	
public:
	
	LoadFasta(const Options &o, const bool isAdapter) :
		
		out(o.out),
		m_rcMode(o.rcMode),
		m_inFileEnd(o.inFileEnd),
		m_isAdapter(isAdapter){
	};
	
	
	virtual ~LoadFasta(){};
	
	
	void loadSequences(const std::string filePath){
		
		using namespace std;
		using namespace flexbar;
		
		seqan::SeqFileIn         seqFileIn;
		seqan::DatFastaSeqFileIn datSeqFileIn(filePath.c_str());
		
		setFormat(seqFileIn,    seqan::Fasta());
		setFormat(datSeqFileIn, seqan::Fasta());
		
		if(m_inFileEnd == IOFF || m_inFileEnd == TXT){
			if(! open(seqFileIn, filePath.c_str())){
				cerr << "\nERROR: Could not open file " << filePath << "\n" << endl;
				exit(1);
			}
		}
		// else if(m_inFileEnd == DAT){
		// 	if(! open(datSeqFileIn, filePath.c_str())){
		// 		cerr << "\nERROR: Could not open file " << filePath << "\n" << endl;
		// 		exit(1);
		// 	}
		// }
		
		TSeqStrs seqs;
		TStrings ids;
		
		try{
			if(m_inFileEnd == IOFF || m_inFileEnd == TXT){
				readRecords(ids, seqs, seqFileIn);
			}
			else if(m_inFileEnd == DAT){
				readRecords(ids, seqs, datSeqFileIn);
			}
			
			map<TString, short> idMap;
			
			for(unsigned int i = 0; i < length(ids); ++i){
				
				if(idMap.count(ids[i]) == 1){
					cerr << "Two ";
					
					if(m_isAdapter) cerr << "adapters";
					else            cerr << "barcodes";
					
					cerr << " have the same name.\n";
					cerr << "Please use unique names and restart.\n" << endl;
					exit(1);
				}
				else idMap[ids[i]] = 1;
				
				if(! m_isAdapter || m_rcMode == RCOFF || m_rcMode == RCON){
					TBar bar;
					bar.id  =  ids[i];
					bar.seq = seqs[i];
					bars.push_back(bar);
				}
				
				if(m_isAdapter && (m_rcMode == RCON || m_rcMode == RCONLY)){
					TString  id =  ids[i];
					TSeqStr seq = seqs[i];
					
					append(id, "_rc");
					seqan::reverseComplement(seq);
					
					TBar barRC;
					barRC.id        = id;
					barRC.seq       = seq;
					barRC.rcAdapter = true;
					bars.push_back(barRC);
				}
			}
		}
		catch(seqan::Exception const &e){
			cerr << "\nERROR: " << e.what() << "\nProgram execution aborted.\n" << endl;
			
			if(m_inFileEnd == IOFF || m_inFileEnd == TXT){
				close(seqFileIn);
			}
			else if(m_inFileEnd == DAT){
				close(datSeqFileIn);
			}
			exit(1);
		}
		
		if(m_inFileEnd == IOFF || m_inFileEnd == TXT){
			close(seqFileIn);
		}
		else if(m_inFileEnd == DAT){
			close(datSeqFileIn);
		}
	};
	
	
	tbb::concurrent_vector<flexbar::TBar> getBars(){
		return bars;
	}
	
	
	void setBars(tbb::concurrent_vector<flexbar::TBar> &newBars){
		bars = newBars;
	}
	
	
	void printBars(std::string adapterName) const {
		
		using namespace std;
		
		const unsigned int maxSpaceLen = 23;
		
		stringstream s; s << adapterName;
		int len = s.str().length() + 1;
		
		if(len + 2 > maxSpaceLen) len = maxSpaceLen - 2;
		
		*out << adapterName << ":" << string(maxSpaceLen - len, ' ') << "Sequence:" << "\n";
		
		for(unsigned int i=0; i < bars.size(); ++i){
			TString seqTag = bars.at(i).id;
			
			int whiteSpaceLen = maxSpaceLen - length(seqTag);
			if(whiteSpaceLen < 2) whiteSpaceLen = 2;
			
			string whiteSpace = string(whiteSpaceLen, ' ');
			
			*out << seqTag << whiteSpace << bars.at(i).seq << "\n";
		}
		*out << endl;
	}
	
};

#endif
