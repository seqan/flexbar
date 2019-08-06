// SeqAlign.h

#ifndef FLEXBAR_SEQALIGN_H
#define FLEXBAR_SEQALIGN_H

tbb::mutex ouputMutex;

template <typename TSeqStr, typename TString, class TAlgorithm>
class SeqAlign {

private:

	typedef AlignResults<TSeqStr> TAlignResults;

	const flexbar::LogAlign    m_log;
	const flexbar::FileFormat  m_format;
	const flexbar::PairOverlap m_poMode;

	const bool m_isBarcoding, m_writeTag, m_umiTags, m_strictRegion, m_addBarcodeAdapter, m_logEverything;
	const int m_minLength, m_minOverlap, m_tailLength;
	const float m_errorRate;
	const unsigned int m_bundleSize;

	tbb::atomic<unsigned long> m_nPreShortReads, m_modified;
	tbb::concurrent_vector<flexbar::TBar> *m_queries;
	tbb::concurrent_vector<unsigned long> m_rmOverlaps;

	std::ostream *m_out;
	TAlgorithm m_algo;

public:

	SeqAlign(tbb::concurrent_vector<flexbar::TBar> *queries, const Options &o, int minOverlap, float errorRate, const int tailLength, const int match, const int mismatch, const int gapCost, const bool isBarcoding):

			m_minOverlap(minOverlap),
			m_errorRate(errorRate),
			m_tailLength(tailLength),
			m_isBarcoding(isBarcoding),
			m_umiTags(o.umiTags),
			m_minLength(o.min_readLen),
			m_poMode(o.poMode),
			m_log(o.logAlign),
			m_logEverything(o.logEverything),
			m_format(o.format),
			m_writeTag(o.useRemovalTag),
			m_addBarcodeAdapter(o.addBarcodeAdapter),
			m_strictRegion(! o.relaxRegion),
			m_bundleSize(o.bundleSize),
			m_out(o.out),
			m_nPreShortReads(0),
			m_modified(0),
			m_algo(TAlgorithm(o, match, mismatch, gapCost, ! isBarcoding)){

		m_queries    = queries;
		m_rmOverlaps = tbb::concurrent_vector<unsigned long>(flexbar::MAX_READLENGTH + 1, 0);
	};

	int alignSeqRead(flexbar::TSeqRead* sr, const bool performRemoval, flexbar::Alignments &alignments, flexbar::ComputeCycle &cycle, unsigned int &idxAl, const flexbar::AlignmentMode &alMode, const flexbar::TrimEnd trimEnd, const TSeqStr &addBarcode){

		using namespace std;
		using namespace flexbar;

		using seqan::prefix;
		using seqan::suffix;

		TSeqRead &seqRead = *sr;
		int readLength    = length(seqRead.seq);

		if(! m_isBarcoding && readLength < m_minLength){
			if(cycle != PRELOAD) ++m_nPreShortReads;
			// return 0;
		}

		if(readLength < 1) return 0;

		if(cycle == PRELOAD){

			if(idxAl == 0) reserve(alignments.aset, m_bundleSize * m_queries->size());

			for(unsigned int i = 0; i < m_queries->size(); ++i){

				if     (alMode == ALIGNRCOFF &&   m_queries->at(i).rcAdapter) continue;
				else if(alMode == ALIGNRC    && ! m_queries->at(i).rcAdapter) continue;

				TSeqStr *qseq = &m_queries->at(i).seq;
				TSeqStr *rseq = &seqRead.seq;
				TSeqStr tmp, tmpq;

				if(! m_isBarcoding && m_addBarcodeAdapter && addBarcode != ""){
					tmpq = addBarcode;
					append(tmpq, m_queries->at(i).seq);
					qseq = &tmpq;
				}

				if(trimEnd == LTAIL || trimEnd == RTAIL){
					int tailLength  = (m_tailLength > 0) ? m_tailLength : length(*qseq);

					if(tailLength < readLength){
						if(trimEnd == LTAIL) tmp = prefix(seqRead.seq, tailLength);
						else                 tmp = suffix(seqRead.seq, readLength - tailLength);
						rseq = &tmp;
					}
				}

				TAlign align;
				appendValue(alignments.aset, align);
				resize(rows(alignments.aset[idxAl]), 2);

				assignSource(row(alignments.aset[idxAl], 0), *rseq);
				assignSource(row(alignments.aset[idxAl], 1), *qseq);

				++idxAl;
			}
			return 0;
		}

		TAlignResults am;

		int qIndex  = -1;
		int pos_bestScore = -1;
		int amScore = numeric_limits<int>::min();
		int amScore_bestScore = amScore;

		std::vector<TAlignResults> am_v;
		std::vector<int> qIndex_v;
		std::vector<int> scores;

		// align each query sequence and store best one
		for(unsigned int i = 0; i < m_queries->size(); ++i){

			if     (alMode == ALIGNRCOFF &&   m_queries->at(i).rcAdapter) continue;
			else if(alMode == ALIGNRC    && ! m_queries->at(i).rcAdapter) continue;

			TAlignResults a;

			// global sequence alignment
			m_algo.alignGlobal(a, alignments, cycle, idxAl++, trimEnd);

			a.queryLength = length(m_queries->at(i).seq);

			if(! m_isBarcoding && m_addBarcodeAdapter && addBarcode != ""){
				a.queryLength += length(addBarcode);
			}

			a.tailLength  = (m_tailLength > 0) ? m_tailLength : a.queryLength;

			a.overlapLength = a.endPos - a.startPos;
			a.allowedErrors = m_errorRate * a.overlapLength;

			float madeErrors = static_cast<float>(a.mismatches + a.gapsR + a.gapsA);
			int minOverlap   = (m_isBarcoding && m_minOverlap == 0) ? a.queryLength : m_minOverlap;

			if(! m_isBarcoding && m_poMode == PON && seqRead.pairOverlap &&
				(trimEnd == RIGHT || trimEnd == RTAIL)) minOverlap = 1;

			bool validAl = true;

			if(((trimEnd == RTAIL  || trimEnd == RIGHT) && a.startPosA < a.startPosS && m_strictRegion) ||
			   ((trimEnd == LTAIL  || trimEnd == LEFT)  && a.endPosA   > a.endPosS   && m_strictRegion) ||
			     a.overlapLength < 1){

				validAl = false;
			}

			// check if alignment is valid, score max, number of errors and overlap length
			if(validAl && a.score > amScore && madeErrors <= a.allowedErrors && a.overlapLength >= minOverlap){
				amScore = a.score;
				if(m_logEverything){
					scores.push_back(amScore);
					am_v.push_back(a);
					qIndex_v.push_back(i);
				}

				if(amScore_bestScore < amScore){
					if(!m_logEverything){
						am      = a;
						qIndex  = i;
					}
					amScore_bestScore = amScore;
					pos_bestScore = qIndex_v.size() - 1;
				}
			}
		}

		// If we are only interested in the best alignment put them in the vector after checking all queries (if valid alignment was found)
		if(!m_logEverything && qIndex != -1){
			am_v.push_back(am);
			qIndex_v.push_back(qIndex);
		}
		else if (qIndex_v.size() > 1 && am_v.size() > 1)
		{
//			qIndex_v.push_back(qIndex_v[pos_bestScore]);
//			qIndex_v.erase(qIndex_v.begin() + pos_bestScore);
//			am_v.push_back(am_v[pos_bestScore]);
//			am_v.erase(am_v.begin() + pos_bestScore);
			std::iter_swap(qIndex_v.begin() + pos_bestScore, qIndex_v.begin() + (qIndex_v.size() - 1)); // + (qIndex_v.size() - 1)
			std::iter_swap(am_v.begin() + pos_bestScore, am_v.begin() + (am_v.size() - 1)); 
		}

		int smallest_diff = -1;
		if(m_logEverything){
			sort(scores.rbegin(), scores.rend());
		if(scores.size() > 1)
			smallest_diff = amScore_bestScore - scores[1];
		}

		string smallest_diff_to_best_score;
		if(smallest_diff == -1){
			smallest_diff_to_best_score = "";
		}else{
			smallest_diff_to_best_score = "/" + std::to_string(smallest_diff);
		}

		stringstream s;

//		TAlignResults *am_p;  am_p = am_v.front(); am_p*

		// valid alignment
		if(qIndex_v.size() > 0){
			for(int i = 0; i < qIndex_v.size(); ++i){
				TSeqRead seqReadTmp = seqRead;
				if(!m_logEverything)
				{
					//only do one iteration of the loop
	//				i = qIndex_v.size();
					// use best alignment
					qIndex = qIndex_v.front();
					am = am_v.front();
				}
				else
				{
					if(i < qIndex_v.size() - 1){
						if(m_log == ALL)
							s << "Alternative alignment:\n";
						qIndex = qIndex_v[i];
						am = am_v[i];
					}
					else
					{
						if(m_log == ALL)
							s << "Best alignment (" << qIndex_v.size() << smallest_diff_to_best_score << "):" << "\n";
						qIndex = qIndex_v[i];
						am = am_v[i];
					}
				}
				TrimEnd trEnd = trimEnd;

				// trim read based on alignment
				if(performRemoval){

					if(trEnd == ANY){

						if(am.startPosA <= am.startPosS && am.endPosS <= am.endPosA){
							seqReadTmp.seq = "";
							if(m_format == FASTQ) seqReadTmp.qual = "";
						}
						else if(am.startPosA - am.startPosS >= am.endPosS - am.endPosA){
							trEnd = RIGHT;
						}
						else trEnd = LEFT;
					}

					switch(trEnd){

						int rCutPos;

						case LTAIL:
						case LEFT:
							rCutPos = am.endPos;

							// translate alignment end pos to read idx
							if(am.startPosS > 0) rCutPos -= am.startPosS;

							// adjust to inner read gaps
							rCutPos -= am.gapsR;

							if(rCutPos > readLength) rCutPos = readLength;

							erase(seqReadTmp.seq, 0, rCutPos);

							if(m_format == FASTQ)
								erase(seqReadTmp.qual, 0, rCutPos);

							break;

						case RTAIL:
							// adjust cut pos to original read length
							am.startPos += readLength - am.tailLength;

						case RIGHT:
							rCutPos = am.startPos;

							// skipped restriction
							if(rCutPos < 0) rCutPos = 0;

							erase(seqReadTmp.seq, rCutPos, readLength);

							if(m_format == FASTQ)
								erase(seqReadTmp.qual, rCutPos, readLength);

							break;

			     			case ANY:;
					}

					++m_modified;
					if(! m_isBarcoding){
						if(! m_queries->at(qIndex).rcAdapter){
							seqReadTmp.rmAdapter   = true;
						}else{
		                                  seqReadTmp.rmAdapterRC = true;
						}
					}
					// count number of removals for each query
					m_queries->at(qIndex).rmOverlap++;
					if(am.overlapLength == am.queryLength)
					m_queries->at(qIndex).rmFull++;
					if(m_writeTag){
						append(seqReadTmp.id, "_Flexbar_removal");

						if(! m_isBarcoding){
							append(seqReadTmp.id, "_");
							append(seqReadTmp.id, m_queries->at(qIndex).id);
						}
					}
					// store overlap occurrences
					if(am.overlapLength <= MAX_READLENGTH) m_rmOverlaps.at(am.overlapLength)++;
					else cerr << "\nCompile Flexbar with larger max read length for correct overlap stats.\n" << endl;
				}

				// valid alignment, not neccesarily removal
				if(m_umiTags && am.umiTag != ""){
					append(seqReadTmp.umi, "_");
					append(seqReadTmp.umi, am.umiTag);
				}

				// alignment stats
				if(m_log == ALL || (m_log == MOD && performRemoval)){
					if(performRemoval){
						s << "Sequence removal:";

						     if(trEnd == LEFT  || trEnd == LTAIL) s << " left side\n";
						else if(trEnd == RIGHT || trEnd == RTAIL) s << " right side\n";
						else                                      s << " any side\n";
					}
					else s << "Sequence detection, no removal:\n";

					s << "  query id         " << m_queries->at(qIndex).id            << "\n"
	  				  << "  query pos        " << am.startPosA << "-" << am.endPosA   << "\n"
					  << "  read id          " << seqReadTmp.id                          << "\n"
					  << "  read pos         " << am.startPosS << "-" << am.endPosS   << "\n"
					  << "  score            " << am.score                            << "\n"
					  << "  overlap          " << am.overlapLength                    << "\n"
					  << "  errors           " << am.gapsR + am.gapsA + am.mismatches << "\n"
					  << "  error threshold  " << am.allowedErrors                    << "\n";

				if(performRemoval){
						s << "  remaining read   " << seqReadTmp.seq << "\n";

					if(m_format == FASTQ)
						s << "  remaining qual   " << seqReadTmp.qual << "\n";
				}
				s << "\n  Alignment:\n" << endl << am.alString;
			}
			else if(m_log == TAB){
				s << seqReadTmp.id    << "\t" << m_queries->at(qIndex).id << "\t"
				  << am.startPosA  << "\t" << am.endPosA               << "\t" << am.overlapLength << "\t"
				  << am.mismatches << "\t" << am.gapsR + am.gapsA      << "\t" << am.allowedErrors << am.score;

				if(m_logEverything){
					if(i < qIndex_v.size() - 1)
						s << "\t" << "a" << endl;
					else
						s << "\t" << "b:" << qIndex_v.size() << smallest_diff_to_best_score << endl;
					}
					else
					{
						s << endl;
					}
				}

				if(i == qIndex_v.size() - 1 || !m_logEverything){
					seqRead = seqReadTmp;
				}

			}
		}
		else if(m_log == ALL){
			s << "Unvalid alignment:"        << "\n"
			  << "read id   " << seqRead.id  << "\n"
			  << "read seq  " << seqRead.seq << "\n\n" << endl;
		}
		
		ouputMutex.lock();
		*m_out << s.str();
		ouputMutex.unlock();
		
		return ++qIndex;
	}

	std::string getOverlapStatsString(){

		using namespace std;
		using namespace flexbar;

		unsigned long nValues = 0, halfValues = 0, cumValues = 0, lenSum = 0;
		unsigned int max = 0, median = 0, mean = 0;

		unsigned int min = numeric_limits<unsigned int>::max();

		for(unsigned int i = 0; i <= MAX_READLENGTH; ++i){
			unsigned long lenCount = m_rmOverlaps.at(i);

			if(lenCount > 0 && i < min) min = i;
			if(lenCount > 0 && i > max) max = i;

			nValues += lenCount;
			lenSum  += lenCount * i;
		}

		halfValues = nValues / 2;

		for(unsigned int i = 0; i <= MAX_READLENGTH; ++i){
			cumValues += m_rmOverlaps.at(i);

			if(cumValues >= halfValues){
				median = i;
				break;
			}
		}

		if(m_modified > 0) mean = lenSum / m_modified;

		stringstream s;

		s << "Min, max, mean and median overlap: ";
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
