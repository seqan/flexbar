/*
 *   SeqAlign.h
 *
 *   Authors: mat and jtr
 */

#ifndef FLEXBAR_SEQALIGN_H
#define FLEXBAR_SEQALIGN_H


template <typename TSeqStr, typename TString, class TAlgorithm>
class SeqAlign {

private:
	
	typedef AlignResults<TSeqStr> TAlignResults;
	
	const flexbar::TrimEnd m_trimEnd;
	const flexbar::LogAlign m_log;
	const flexbar::FileFormat m_format;
	
	const bool m_isBarcoding, m_writeTag, m_randTag, m_strictRegion;
	const int m_minLength, m_minOverlap, m_tailLength;
	const float m_threshold;
	
	const unsigned int m_bundleSize;
	
	tbb::atomic<unsigned long> m_nPreShortReads, m_modified;
	
	tbb::concurrent_vector<flexbar::TBar> *m_queries;
	tbb::concurrent_vector<unsigned long> *m_rmOverlaps;
	
	std::ostream *m_out;
	TAlgorithm *algo;
	
public:
	
	SeqAlign(tbb::concurrent_vector<flexbar::TBar> *queries, const Options &o, int minOverlap, float threshold, const int tailLength, const int match, const int mismatch, const int gapCost, const flexbar::TrimEnd end, const bool isBarcoding):
			
			m_minOverlap(minOverlap),
			m_threshold(threshold),
			m_tailLength(tailLength),
			m_trimEnd(end),
			m_isBarcoding(isBarcoding),
			m_randTag(o.randTag),
			m_minLength(o.min_readLen),
			m_log(o.logAlign),
			m_format(o.format),
			m_writeTag(o.useRemovalTag),
			m_strictRegion(! o.relaxRegion),
			m_bundleSize(o.bundleSize),
			m_out(o.out){
		
		m_queries = queries;
		
		m_nPreShortReads = 0;
		m_modified       = 0;
		
		algo = new TAlgorithm(o, match, mismatch, gapCost, m_trimEnd);
		
		m_rmOverlaps = new tbb::concurrent_vector<unsigned long>(flexbar::MAX_READLENGTH + 1, 0);
	};
	
	
	virtual ~SeqAlign(){
		delete algo;
		delete m_rmOverlaps;
	};
	
	
	int alignSeqRead(flexbar::TSeqRead* sr, const bool performRemoval, flexbar::TAlignments &alignments, flexbar::ComputeCycle &cycle, unsigned int &idxAl){
		
		using namespace std;
		using namespace flexbar;
		
		using seqan::prefix;
		using seqan::suffix;
		using seqan::infix;
		
		TSeqRead &seqRead = *sr;
		
		int readLength = length(seqRead.seq);
		
		if(! m_isBarcoding && readLength < m_minLength){
			if(cycle != PRELOAD) ++m_nPreShortReads;
			return 0;
		}
		
		
		if(cycle == PRELOAD){
			
			if(idxAl == 0) reserve(alignments.first, m_bundleSize * m_queries->size());
			
			for(unsigned int i = 0; i < m_queries->size(); ++i){
				
				TSeqStr &qseq = m_queries->at(i).seq;
				TSeqStr *rseq;
				
				rseq = &seqRead.seq;
				
				if(m_trimEnd == LEFT_TAIL || m_trimEnd == RIGHT_TAIL){
					int tailLength  = (m_tailLength > 0) ? m_tailLength : length(qseq);
					
					if(tailLength < readLength){
						
						TSeqStr seq;
						if(m_trimEnd == LEFT_TAIL) seq = prefix(seqRead.seq, tailLength);
						else                       seq = suffix(seqRead.seq, readLength - tailLength);
						
						rseq = &seq;
					}
				}
				
				TAlign align;
				
				appendValue(alignments.first, align);
				resize(rows(alignments.first[idxAl]), 2);
				
				assignSource(row(alignments.first[idxAl], 0), *rseq);
				assignSource(row(alignments.first[idxAl], 1),  qseq);
				
				++idxAl;
			}
			return 0;
		}
		
		
		TAlignResults *am;
		
		int qIndex  = -1;
		int amScore = numeric_limits<int>::min();
		
		// align each query sequence and store best one
		for(unsigned int i = 0; i < m_queries->size(); ++i){
			
			if(i > 0) cycle = RESULTS;
			
			TAlignResults a;
			
			// global sequence alignment
			algo->alignGlobal(a, alignments, cycle, idxAl++);
			
			a.queryLength = length(m_queries->at(i).seq);
			a.tailLength  = (m_tailLength > 0) ? m_tailLength : a.queryLength;
			
			a.overlapLength = a.endPos - a.startPos;
			a.allowedErrors = m_threshold * a.overlapLength / 10.0f;
			
			float madeErrors = static_cast<float>(a.mismatches + a.gapsR + a.gapsA);
			int minOverlap   = (m_isBarcoding && m_minOverlap == 0) ? a.queryLength : m_minOverlap;
			
			bool validAl = true;
			
			if(((m_trimEnd == RIGHT_TAIL || m_trimEnd == RIGHT) && a.startPosA < a.startPosS && m_strictRegion) ||
			   ((m_trimEnd == LEFT_TAIL  || m_trimEnd == LEFT)  && a.endPosA   > a.endPosS   && m_strictRegion) ||
			     a.overlapLength < 1){
				
				validAl = false;
			}
			
			// check if alignment is valid, score is max, number of errors and overlap length
			if(validAl && a.score > amScore && madeErrors <= a.allowedErrors && a.overlapLength >= minOverlap){
				
				am      = &a;
				amScore = a.score;
				qIndex  = i;
			}
		}
		
		
		// valid alignment
		if(qIndex >= 0){
			
			TrimEnd trimEnd = m_trimEnd;
			
			// trim read based on alignment
			if(performRemoval){
				
				if(trimEnd == ANY){
					
					if(am->startPosA <= am->startPosS && am->endPosS <= am->endPosA){
						seqRead.seq = "";
						if(m_format == FASTQ) seqRead.qual = "";
					}
					else if(am->startPosA - am->startPosS >= am->endPosS - am->endPosA){
						trimEnd = RIGHT;
					}
					else trimEnd = LEFT;
				}
				
				switch(trimEnd){
					
					int rCutPos;
					
					case LEFT_TAIL:
					case LEFT:
						rCutPos = am->endPos;
						
						// translate alignment end pos to read idx
						if(am->startPosS > 0) rCutPos -= am->startPosS;
						
						// adjust to inner read gaps
						rCutPos -= am->gapsR;
						
						if(rCutPos > readLength) rCutPos = readLength;
						
						erase(seqRead.seq, 0, rCutPos);
						
						if(m_format == FASTQ)
						erase(seqRead.qual, 0, rCutPos);
						
						break;
					
					case RIGHT_TAIL:
						// adjust cut pos to original read length
						am->startPos += readLength - am->tailLength;
					
					case RIGHT:
						rCutPos = am->startPos;
						
						// skipped restriction
						if(rCutPos < 0) rCutPos = 0;
						
						erase(seqRead.seq, rCutPos, readLength);
						
						if(m_format == FASTQ)
						erase(seqRead.qual, rCutPos, readLength);
						
						break;
						
	                case ANY:;
				}
				
				++m_modified;
				
				// count number of removals for each query
				m_queries->at(qIndex).rmOverlap++;
				
				if(am->overlapLength == am->queryLength)
				m_queries->at(qIndex).rmFull++;
				
				if(m_writeTag){
					append(seqRead.id, "_Flexbar_removal");
					
					if(! m_isBarcoding){
						append(seqRead.id, "_");
						append(seqRead.id, m_queries->at(qIndex).id);
					}
				}
				
				// store overlap occurrences
				if(am->overlapLength <= MAX_READLENGTH) m_rmOverlaps->at(am->overlapLength)++;
				else cerr << "\nCompile Flexbar with larger max read length for correct overlap stats.\n" << endl;
			}
			
			// valid alignment, not neccesarily removal
			
			if(m_randTag && am->randTag != ""){
				append(seqRead.id, "_");
				append(seqRead.id, am->randTag);
			}
			
			// alignment stats
			if(m_log == ALL || (m_log == MOD && performRemoval)){
				
				stringstream s;
				
				if(performRemoval){
					s << "Sequence removal:";
					
					     if(trimEnd == LEFT  || trimEnd == LEFT_TAIL)  s << " left side\n";
					else if(trimEnd == RIGHT || trimEnd == RIGHT_TAIL) s << " right side\n";
					else                                               s << " any side\n";
				}
				else s << "Sequence detection, no removal:\n";
				
				s << "  query tag        " << m_queries->at(qIndex).id      << "\n"
				  << "  read tag         " << seqRead.id                            << "\n"
				  << "  read seq         " << seqRead.seq                            << "\n"
				  << "  read pos         " << am->startPosS << "-" << am->endPosS    << "\n"
				  << "  query pos        " << am->startPosA << "-" << am->endPosA    << "\n"
				  << "  score            " << am->score                              << "\n"
				  << "  overlap          " << am->overlapLength                      << "\n"
				  << "  errors           " << am->gapsR + am->gapsA + am->mismatches << "\n"
				  << "  allowed errors   " << am->allowedErrors                      << "\n";
				
				if(performRemoval){
					s << "  remaining read   "  << seqRead.seq << "\n";
					
					if(m_format == FASTQ)
					s << "  remaining qual   " << seqRead.qual << "\n";
				}
				s << "\n  Alignment:\n" << endl << am->alString;
				*m_out << s.str();
			}
			else if(m_log == TAB){
				
				stringstream s;
				s << seqRead.id    << "\t" << m_queries->at(qIndex).id                  << "\t"
				  << am->startPosA  << "\t" << am->endPosA           << "\t" << am->overlapLength << "\t"
				  << am->mismatches << "\t" << am->gapsR + am->gapsA << "\t" << am->allowedErrors << endl;
				*m_out << s.str();
			}
		}
		else if(m_log == ALL){
			
			stringstream s;
			s << "No valid alignment:"       << "\n"
			  << "read tag  " << seqRead.id << "\n"
			  << "read seq  " << seqRead.seq << "\n\n" << endl;
			*m_out << s.str();
		}
		
		return ++qIndex;
	}
	
	
	std::string getOverlapStatsString(){
		
		using namespace std;
		using namespace flexbar;
		
		unsigned long nValues = 0, halfValues = 0, cumValues = 0, lenSum = 0;
		unsigned int max = 0, median = 0, mean = 0;
		
		unsigned int min = numeric_limits<unsigned int>::max();
		
		for(unsigned int i = 0; i <= MAX_READLENGTH; ++i){
			unsigned long lenCount = m_rmOverlaps->at(i);
			
			if(lenCount > 0 && i < min) min = i;
			if(lenCount > 0 && i > max) max = i;
			
			nValues += lenCount;
			lenSum  += lenCount * i;
		}
		
		halfValues = nValues / 2;
		
		for(unsigned int i = 0; i <= MAX_READLENGTH; ++i){
			cumValues += m_rmOverlaps->at(i);
			
			if(cumValues >= halfValues){
				median = i;
				break;
			}
		}
		
		if(m_modified > 0) mean = lenSum / m_modified;
		
		stringstream s;
		
		s << "Min, max, mean, and median overlap: ";
		s << min << " / " << max << " / " << mean << " / " << median;
		
		return s.str();
	}
	
	
	unsigned long getNrPreShortReads() const {
		return m_nPreShortReads;
	}
	
	
	unsigned long getNrModifiedReads() const {
		return m_modified;
	}
	
};

#endif
