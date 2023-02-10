// SeqAlignPair.h

#ifndef FLEXBAR_SEQALIGNPAIR_H
#define FLEXBAR_SEQALIGNPAIR_H


template <typename TSeqStr, typename TString, class TAlgorithm>
class SeqAlignPair {

private:
	
	typedef AlignResults<TSeqStr> TAlignResults;
	
	const flexbar::LogAlign    m_log;
	const flexbar::FileFormat  m_format;
	const flexbar::PairOverlap m_poMode;
	
	const bool m_writeTag;
	const int m_minLength, m_minOverlap, m_aMinOverlap;
	const float m_errorRate;
	const unsigned int m_bundleSize;
	
	std::atomic<unsigned long> m_nPreShortReads, m_overlaps, m_modified;
	oneapi::tbb::concurrent_vector<unsigned long> m_overlapLengths;
	
	std::ostream *m_out;
	TAlgorithm m_algo;
	
public:
	
	SeqAlignPair(const Options &o, const int minOverlap, const float errorRate, const int match, const int mismatch, const int gapCost):
			
			m_minOverlap(minOverlap),
			m_aMinOverlap(o.a_min_overlap),
			m_errorRate(errorRate),
			m_minLength(o.min_readLen),
			m_poMode(o.poMode),
			m_log(o.logAlign),
			m_format(o.format),
			m_writeTag(o.useRemovalTag),
			m_bundleSize(o.bundleSize),
			m_out(o.out),
			m_nPreShortReads(0),
			m_overlaps(0),
			m_modified(0),
			m_algo(TAlgorithm(o, match, mismatch, gapCost, true)){
		
		m_overlapLengths = oneapi::tbb::concurrent_vector<unsigned long>(flexbar::MAX_READLENGTH + 1, 0);
	};
	
	
	void alignSeqReadPair(flexbar::TSeqRead* sr, flexbar::TSeqRead* sr2, flexbar::Alignments &alignments, flexbar::ComputeCycle &cycle, unsigned int &idxAl){
		
		using namespace std;
		using namespace flexbar;
		
		TSeqRead &seqRead  = *sr;
		TSeqRead &seqRead2 = *sr2;
		
		int readLength  = length(seqRead.seq);
		int readLength2 = length(seqRead2.seq);
		
		if(cycle != PRELOAD){
			if(readLength  < m_minLength) ++m_nPreShortReads;
			if(readLength2 < m_minLength) ++m_nPreShortReads;
		}
		
		if(readLength < 1 || readLength2 < 1) return;
		
		
		if(cycle == PRELOAD){
			
			if(idxAl == 0) reserve(alignments.aset, m_bundleSize);
			
			TSeqStr rcSeq2 = seqRead2.seq;
			seqan::reverseComplement(rcSeq2);
			
			TAlign align;
			appendValue(alignments.aset, align);
			resize(rows(alignments.aset[idxAl]), 2);
			
			assignSource(row(alignments.aset[idxAl], 0), seqRead.seq);
			assignSource(row(alignments.aset[idxAl], 1), rcSeq2);
			
			++idxAl;
			return;
		}
		
		TAlignResults a;
		
		m_algo.alignGlobal(a, alignments, cycle, idxAl++, ANY);
		
		a.overlapLength = a.endPos - a.startPos;
		a.allowedErrors = m_errorRate * a.overlapLength;
		
		float madeErrors = static_cast<float>(a.mismatches + a.gapsR + a.gapsA);
		
		stringstream s;
		
		// check if alignment is valid, number of errors and overlap length
		if((a.startPosA < a.startPosS || a.endPosA < a.endPosS) && madeErrors <= a.allowedErrors && a.overlapLength >= m_minOverlap){
			
			if(a.startPosA < a.startPosS){
				
				seqRead2.pairOverlap = true;
				
				if(m_poMode == PONLY || (m_poMode == PSHORT && a.startPosS < m_aMinOverlap)){
					
					unsigned int rCutPos = readLength2 - a.startPosS;
					erase(seqRead2.seq, rCutPos, readLength2);
					
					if(m_format == FASTQ)
					erase(seqRead2.qual, rCutPos, readLength2);
					
					++m_modified;
					
					seqRead2.poRemoval = true;
					
					if(m_writeTag) append(seqRead2.id, "_Flexbar_removal_PO");
				}
			}
			
			if(a.endPosA < a.endPosS){
				
				seqRead.pairOverlap = true;
				
				if(m_poMode == PONLY || (m_poMode == PSHORT && (a.endPosS - a.endPosA) < m_aMinOverlap)){
					
					unsigned int rCutPos = readLength - (a.endPosS - a.endPosA);
					erase(seqRead.seq, rCutPos, readLength);
					
					if(m_format == FASTQ)
					erase(seqRead.qual, rCutPos, readLength);
					
					++m_modified;
					
					seqRead.poRemoval = true;
					
					if(m_writeTag) append(seqRead.id, "_Flexbar_removal_PO");
				}
			}
			
			++m_overlaps;
			
			// store overlap occurrences
			if(a.overlapLength <= MAX_READLENGTH) m_overlapLengths.at(a.overlapLength)++;
			else cerr << "\nCompile Flexbar with larger max read length for correct overlap stats.\n" << endl;
			
			// alignment stats
			if(m_log == ALL || m_log == MOD){
				
				s << "Sequence removal:\n";
				
				s << "  read id          " << seqRead.id                          << "\n"
				  << "  read pos         " << a.startPosS << "-" << a.endPosS     << "\n"
				  << "  read2 id         " << seqRead2.id                         << "\n"
  				  << "  read2 pos        " << a.startPosA << "-" << a.endPosA     << "\n"
				  << "  score            " << a.score                             << "\n"
				  << "  overlap          " << a.overlapLength                     << "\n"
				  << "  errors           " << a.gapsR + a.gapsA + a.mismatches    << "\n"
				  << "  error threshold  " << a.allowedErrors                     << "\n"
				  << "  remaining read   " << seqRead.seq                         << "\n";
				
				if(m_format == FASTQ)
				s << "  remaining qual   " << seqRead.qual  << "\n";
				
				s << "  remaining read2  " << seqRead2.seq  << "\n";
				
				if(m_format == FASTQ)
				s << "  remaining qual2  " << seqRead2.qual << "\n";
				
				s << "\n  Alignment:\n" << endl << a.alString;
			}
			else if(m_log == TAB){
				s << seqRead.id   << "\t" << seqRead2.id       << "\t"
				  << a.startPosA  << "\t" << a.endPosA         << "\t" << a.overlapLength << "\t"
				  << a.mismatches << "\t" << a.gapsR + a.gapsA << "\t" << a.allowedErrors << endl;
			}
		}
		else if(m_log == ALL){
			s << "Unvalid alignment:"        << "\n"
			  << "read id   " << seqRead.id  << "\n"
			  << "read2 id  " << seqRead2.id << "\n\n" << endl;
		}
		
		*m_out << s.str();
		
		return;
	}
	
	
	std::string getOverlapStatsString(){
		
		using namespace std;
		using namespace flexbar;
		
		unsigned long nValues = 0, halfValues = 0, cumValues = 0, lenSum = 0;
		unsigned int max = 0, median = 0, mean = 0;
		
		unsigned int min = numeric_limits<unsigned int>::max();
		
		for(unsigned int i = 0; i <= MAX_READLENGTH; ++i){
			unsigned long lenCount = m_overlapLengths.at(i);
			
			if(lenCount > 0 && i < min) min = i;
			if(lenCount > 0 && i > max) max = i;
			
			nValues += lenCount;
			lenSum  += lenCount * i;
		}
		
		halfValues = nValues / 2;
		
		for(unsigned int i = 0; i <= MAX_READLENGTH; ++i){
			cumValues += m_overlapLengths.at(i);
			
			if(cumValues >= halfValues){
				median = i;
				break;
			}
		}
		
		if(m_overlaps > 0) mean = lenSum / m_overlaps;
		
		stringstream s;
		
		if(m_modified > 0){
			s << "Number of trimmed reads based on pair overlap:     ";
			s << m_modified << "\n";
		}
		
		s << "Min, max, mean and median overlap of paired reads: ";
		s << min << " / " << max << " / " << mean << " / " << median;
		
		return s.str();
	}
	
	
	unsigned long getNrPreShortReads() const {
		return m_nPreShortReads;
	}
	
	
	unsigned long getNrOverlappingReads() const {
		return m_overlaps;
	}
	
};

#endif
