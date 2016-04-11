/*
 *   AlignmentFilter.h
 *
 *   Authors: mat and jtr
 */

#ifndef FLEXBAR_ALIGNMENTFILTER_H
#define FLEXBAR_ALIGNMENTFILTER_H

#include <string>
#include <sstream>

#include <tbb/pipeline.h>
#include <tbb/concurrent_vector.h>

#include "Enums.h"


template <typename TSeqStr, typename TString, class TAlgorithm>
class AlignmentFilter {

private:
	
	const flexbar::TrimEnd m_trimEnd;
	const flexbar::LogLevel m_verb;
	const flexbar::FileFormat m_format;
	
	const bool m_isBarcoding, m_writeTag, m_randTag, m_strictRegion;
	const int m_minLength, m_minOverlap, m_tailLength;
	const float m_threshold;
	
	tbb::atomic<unsigned long> m_nPreShortReads, m_modified;
	
	tbb::concurrent_vector<flexbar::TAdapter> *m_queries;
	tbb::concurrent_vector<unsigned long> *m_rmOverlaps;
	
	std::ostream *m_out;
	TAlgorithm *algo;
	
public:
	
	AlignmentFilter(tbb::concurrent_vector<flexbar::TAdapter> *queries, const Options &o, int minOverlap, float threshold, const int tailLength, const int match, const int mismatch, const int gapCost, const flexbar::TrimEnd end, const bool isBarcoding):
			
			m_minOverlap(minOverlap),
			m_threshold(threshold),
			m_tailLength(tailLength),
			m_trimEnd(end),
			m_isBarcoding(isBarcoding),
			m_randTag(o.randTag),
			m_minLength(o.min_readLen),
			m_verb(o.logLevel),
			m_format(o.format),
			m_writeTag(o.useRemovalTag),
			m_strictRegion(! o.relaxRegion),
			m_out(o.out){
		
		m_queries = queries;
		
		m_nPreShortReads = 0;
		m_modified       = 0;
		
		algo = new TAlgorithm(o, match, mismatch, gapCost, m_trimEnd);
		
		m_rmOverlaps = new tbb::concurrent_vector<unsigned long>(flexbar::MAX_READLENGTH + 1, 0);
	};
	
	
	virtual ~AlignmentFilter(){
		delete algo;
		delete m_rmOverlaps;
	};
	
	
	int align(void* item, const bool performRemoval){
		
		using namespace std;
		using namespace flexbar;
		
		using seqan::prefix;
		using seqan::suffix;
		using seqan::infix;
		
		
		SeqRead<TSeqStr, TString> &myRead = *static_cast< SeqRead<TSeqStr, TString>* >(item);
		
		int fmismatches, fgapsR, fgapsA, foverlapLength, fqueryLength, ftailLength;
		int fstartPos, fstartPosA, fstartPosS, fendPos, fendPosS, fendPosA;
		
		int qIndex   = -1;
		int scoreMax = -1000000;
		
		float fallowedErrors;
		
		stringstream ss;
		TSeqStr seqread, finalRandTag;
		
		TString quality, finalAliStr;
		
		TString readTag = myRead.getSequenceTag();
		
		seqread = myRead.getSequence();
		quality = "";
		
		if(m_format == FASTQ) quality = myRead.getQuality();
		
		TSeqStr sequence = seqread;
		int readLength   = length(seqread);
		
		if(! m_isBarcoding && readLength < m_minLength){
			++m_nPreShortReads;
			return ++qIndex;
		}
		
		
		// align each query sequence and keep track of best one
		for(unsigned int i = 0; i < m_queries->size(); ++i){
			
			TSeqStr query = m_queries->at(i).first->getSequence();
			
			int queryLength = length(query);
			int tailLength  = (m_tailLength > 0) ? m_tailLength : queryLength;
			
			
			if(m_trimEnd == LEFT_TAIL || m_trimEnd == RIGHT_TAIL){
				
				if(tailLength < readLength){
					
					if(m_trimEnd == LEFT_TAIL){
						sequence = prefix<TSeqStr>(seqread, tailLength);
					}else{
						sequence = suffix<TSeqStr>(seqread, readLength - tailLength);
					}
					if(m_verb == ALL || m_verb == MOD)
					ss << "Read tail length:  " << tailLength << "\n\n";
				}
			}
			
			
			int startPos = 0, endPos = 0, startPosA = 0, endPosA = 0, startPosS = 0, endPosS = 0;
			int aliScore = 0, mismatches = 0, gapsR = 0, gapsA = 0;
			
			TSeqStr randTag = "";
			stringstream aliString;
			
			// align query with specified algorithm
			algo->align(query, sequence, gapsR, gapsA, mismatches, startPos, endPos,
			            startPosA, endPosA, startPosS, endPosS, aliScore, aliString, randTag);
			
			int overlapLength = endPos - startPos;
			
			float allowedErrors = m_threshold * overlapLength / 10.0f;
			float madeErrors    = static_cast<float>(mismatches + gapsR + gapsA);
			
			int minOverlapValue = (m_isBarcoding && m_minOverlap == 0) ? queryLength : m_minOverlap;
			
			
			bool validAli = true;
			
			if(((m_trimEnd == RIGHT_TAIL || m_trimEnd == RIGHT) && startPosA < startPosS && m_strictRegion) ||
			   ((m_trimEnd == LEFT_TAIL  || m_trimEnd == LEFT)  && endPosA   > endPosS   && m_strictRegion) ||
			     overlapLength < 1){
				
				validAli = false;
			}
			
			
			// check if alignment is valid and score is max as well as if number of errors and overlap length are allowed
			if(validAli && aliScore > scoreMax && madeErrors <= allowedErrors && overlapLength >= minOverlapValue){
				
				qIndex      = i;
				scoreMax    = aliScore;
				fstartPos   = startPos;
				fstartPosA  = startPosA;
				fstartPosS  = startPosS;
				fendPos     = endPos;
				fendPosA    = endPosA;
				fendPosS    = endPosS;
				fgapsR      = gapsR;
				fgapsA      = gapsA;
				
				finalRandTag   = randTag;
				ftailLength    = tailLength;
				foverlapLength = overlapLength;
				fqueryLength   = queryLength;
				
				if(m_verb != NONE){
					fmismatches    = mismatches;
					finalAliStr    = aliString.str();
					fallowedErrors = allowedErrors;
				}
			}
		}
		
		
		// valid alignment
		if(qIndex >= 0){
			
			TrimEnd trimEnd = m_trimEnd;
			
			// cut read according to best alignment
			if(performRemoval){
				
				if(trimEnd == ANY){
					
					if(fstartPosA <= fstartPosS && fendPosS <= fendPosA){
						myRead.setSequence("");
						if(m_format == FASTQ) myRead.setQuality("");
					}
					else if(fstartPosA - fstartPosS >= fendPosS - fendPosA){
						trimEnd = RIGHT;
					}
					else{
						trimEnd = LEFT;
					}
				}
				
				switch(trimEnd){
					
					int rCutPos;
					
					case LEFT_TAIL:
						sequence = seqread;
					
					case LEFT:
						rCutPos = fendPos;
						
						// translate alignment end pos to read idx
						if(fstartPosS > 0) rCutPos -= fstartPosS;
						
						// adjust to inner read gaps
						rCutPos -= fgapsR;
						
						if(rCutPos > readLength) rCutPos = readLength;
						
						erase(sequence, 0, rCutPos);
						myRead.setSequence(sequence);
						
						if(m_format == FASTQ){
							erase(quality, 0, rCutPos);
							myRead.setQuality(quality);
						}
						break;
					
					case RIGHT_TAIL:
						sequence  = seqread;
						// adjust cut pos to original read length
						fstartPos += readLength - ftailLength;
					
					case RIGHT:
						rCutPos = fstartPos;
						
						// skipped restriction
						if(rCutPos < 0) rCutPos = 0;
						
						erase(sequence, rCutPos, readLength);
						myRead.setSequence(sequence);
						
						if(m_format == FASTQ){
							erase(quality, rCutPos, readLength);
							myRead.setQuality(quality);
						}
						break;
						
	                case ANY:;
				}
				
				++m_modified;
				
				
				// count for each query number of removals
				m_queries->at(qIndex).second.first++;
				
				if(foverlapLength == fqueryLength){
					m_queries->at(qIndex).second.second++;
				}
				
				if(m_writeTag){
					TString newTag = myRead.getSequenceTag();
					append(newTag, "_Flexbar_removal");
					
					if(! m_isBarcoding){
						append(newTag, "_");
						append(newTag, m_queries->at(qIndex).first->getSequenceTag());
					}
					
					myRead.setSequenceTag(newTag);
				}
				
				// store overlap occurrences for min, max, mean and median
				if(foverlapLength <= MAX_READLENGTH) m_rmOverlaps->at(foverlapLength)++;
				else cerr << "\nCompile Flexbar with larger max read length to get correct overlap stats.\n" << endl;
			}
			
			
			// valid alignment, not neccesarily removal
			
			if(m_randTag && finalRandTag != ""){
				TString newTag = myRead.getSequenceTag();
				append(newTag, "_");
				append(newTag, finalRandTag);
				myRead.setSequenceTag(newTag);
			}
			
			
			// alignment stats
			TString queryTag = m_queries->at(qIndex).first->getSequenceTag();
			
			if(m_verb == ALL || (m_verb == MOD && performRemoval)){
				
				if(performRemoval){
					ss << "Sequence removal:";
					
					     if(trimEnd == LEFT  || trimEnd == LEFT_TAIL)  ss << "  left side\n";
					else if(trimEnd == RIGHT || trimEnd == RIGHT_TAIL) ss << "  right side\n";
					else                                               ss << "  any side\n";
				}
				else{ ss << "Sequence detection, no removal:\n"; }
				
				ss << "  query tag        " << queryTag                      << "\n"
				   << "  read tag         " << readTag                       << "\n"
				   << "  read             " << seqread                       << "\n"
				   << "  read pos         " << fstartPosS << "-" << fendPosS << "\n"
				   << "  query pos        " << fstartPosA << "-" << fendPosA << "\n"
				   << "  score            " << scoreMax                      << "\n"
				   << "  overlap          " << foverlapLength                << "\n"
				   << "  errors           " << fgapsR + fgapsA + fmismatches << "\n"
				   << "  allowed errors   " << fallowedErrors                << "\n";
				
				if(performRemoval){
					ss << "  remaining read   "  << myRead.getSequence() << "\n";
					
					if(m_format == FASTQ)
					ss << "  remaining qual   " << myRead.getQuality() << "\n";
				}
				
				ss << "\n  Alignment:\n" << endl << finalAliStr;
			}
			else if(m_verb == TAB){
				ss << readTag     << "\t" << queryTag        << "\t"
				   << fstartPosA  << "\t" << fendPosA        << "\t" << foverlapLength << "\t"
				   << fmismatches << "\t" << fgapsR + fgapsA << "\t" << fallowedErrors << endl;
			}
		}
		else if(m_verb == ALL){
			ss << "No valid alignment:"   << "\n"
			   << "read tag  " << readTag << "\n"
			   << "read      " << seqread << "\n\n" << endl;
		}
		
		// bundeled output for multi-threading
		if(m_verb != NONE) *m_out << ss.str();
		
		return ++qIndex;
	}
	
	
	std::string getOverlapStatsString(){
		
		using namespace flexbar;
		
		unsigned long nValues = 0, halfValues = 0, cumValues = 0, lenSum = 0;
		int min = 1000000, max = 0, median = 0, mean = 0;
		
		for (int i = 0; i <= MAX_READLENGTH; ++i){
			unsigned long lenCount = m_rmOverlaps->at(i);
			
			if(lenCount > 0 && i < min) min = i;
			if(lenCount > 0 && i > max) max = i;
			
			nValues += lenCount;
			lenSum  += lenCount * i;
		}
		
		halfValues = nValues / 2;
		
		for (int i = 0; i <= MAX_READLENGTH; ++i){
			cumValues += m_rmOverlaps->at(i);
			
			if(cumValues >= halfValues){
				median = i;
				break;
			}
		}
		
		if(m_modified > 0) mean = lenSum / m_modified;
		
		std::stringstream ss;
		
		ss << "Min, max, mean and median adapter overlap: ";
		ss << min << " / " << max << " / " << mean << " / " << median;
		
		return ss.str();
	}
	
	
	unsigned long getNrPreShortReads() const {
		return m_nPreShortReads;
	}
	
	
	unsigned long getNrModifiedReads() const {
		return m_modified;
	}
	
};

#endif
