/*
 *   LoadFasta.h
 *
 */

#ifndef FLEXBAR_LOADFASTA_H
#define FLEXBAR_LOADFASTA_H


template <typename TSeqStr, typename TString>
class LoadFasta : public tbb::filter{

private:
	
	std::ostream *out;
	tbb::concurrent_vector<flexbar::TAdapter> adapters;
	
	bool m_revComp, m_isAdapter;
	
public:
	
	LoadFasta(const Options &o, const bool isAdapter) :
		
		filter(serial),
		out(o.out),
		m_isAdapter(isAdapter){
			
			m_revComp = o.revCompAdapter && isAdapter;
	};
	
	
	virtual ~LoadFasta(){};
	
	
	void* operator()(void* item){
		
		using namespace std;
		using namespace flexbar;
		
		SeqRead<TSeqStr, TString> *seqRead = static_cast< SeqRead<TSeqStr, TString>* >(item);
		SeqRead<TSeqStr, TString> *seqReadRC;
		
		TString tag = seqRead->tag;
		
		if(adapters.size() < 1000){
			for(int i = 0; i < adapters.size(); ++i){
				
				if(tag == adapters.at(i).first->tag){
					cerr << "Two ";
					
					if(m_isAdapter) cerr << "adapters";
					else            cerr << "barcodes";
					
					cerr << " have the same name.\n";
					cerr << "Please use unique names and restart.\n" << endl;
					
					exit(1);
				}
			}
		}
		
		if(m_revComp){
			TSeqStr seq = seqRead->seq;
			seqan::reverseComplement(seq);
			
			append(tag, " revcomp");
			
			seqReadRC = new SeqRead<TSeqStr, TString>(seq, tag);
		}
		
		TAdapter adap;
		adap.first = seqRead;
		adapters.push_back(adap);
		
		if(m_revComp){
			TAdapter adapRC;
			adapRC.first = seqReadRC;
			adapters.push_back(adapRC);
		}
		
		return NULL;
	};
	
	
	tbb::concurrent_vector<flexbar::TAdapter> getAdapters(){
		return adapters;
	}
	
	
	void setAdapters(tbb::concurrent_vector<flexbar::TAdapter> &adapterVec){
		adapters = adapterVec;
	}
	
	
	void printAdapters(std::string adapterName) const {
		using namespace std;
		
		const unsigned int maxSpaceLen = 23;
		
		stringstream s; s << adapterName;
		int len = s.str().length() + 1;
		
		if(len + 2 > maxSpaceLen) len = maxSpaceLen - 2;
		
		*out << adapterName << ":" << string(maxSpaceLen - len, ' ') << "Sequence:" << "\n";
		
		for(unsigned int i=0; i < adapters.size(); ++i){
			TString seqTag = adapters.at(i).first->tag;
			
			int whiteSpaceLen = maxSpaceLen - length(seqTag);
			if(whiteSpaceLen < 2) whiteSpaceLen = 2;
			
			string whiteSpace = string(whiteSpaceLen, ' ');
			
			*out << seqTag << whiteSpace << adapters.at(i).first->seq << "\n";
		}
		*out << endl;
	}
	
};

#endif
