// SeqAlignAlgo.h

#ifndef FLEXBAR_SEQALIGNALGO_H
#define FLEXBAR_SEQALIGNALGO_H


template <typename TSeqStr>
class SeqAlignAlgo {

private:
	
	typedef typename seqan::Value<TSeqStr>::Type        TChar;
	typedef typename seqan::Row<flexbar::TAlign>::Type  TRow;
	typedef typename seqan::Iterator<TRow>::Type        TRowIterator;
	
	typedef AlignResults<TSeqStr> TAlignResults;
	
	typedef seqan::Score<int, seqan::Simple>                              TScoreSimple;
	typedef seqan::Score<int, seqan::ScoreMatrix<TChar, seqan::Default> > TScoreMatrix;
	
	// TScoreSimple m_score;
	TScoreMatrix m_scoreMatrix;
	
	const bool m_randTag;
	const flexbar::LogAlign m_log;
	const flexbar::TrimEnd m_trimEnd;
	
public:
	
	SeqAlignAlgo(const Options &o, const int match, const int mismatch, const int gapCost, const flexbar::TrimEnd trimEnd):
			m_randTag(o.randTag),
			m_log(o.logAlign),
			m_trimEnd(trimEnd){
		
		using namespace seqan;
		
		// m_score       = Score<int, Simple>(match, mismatch, gapCost);
		m_scoreMatrix = TScoreMatrix(gapCost);
		
		for(unsigned i = 0; i < ValueSize<TChar>::VALUE; ++i){
			for(unsigned j = 0; j < ValueSize<TChar>::VALUE; ++j){
				
				if(i == j || TChar(j) == 'N')
					 setScore(m_scoreMatrix, TChar(i), TChar(j), match);
				else setScore(m_scoreMatrix, TChar(i), TChar(j), mismatch);
			}
		}
		
		// printScoreMatrix(m_scoreMatrix);
	};
	
	
	void alignGlobal(TAlignResults &a, flexbar::Alignments &alignments, flexbar::ComputeCycle &cycle, const unsigned int idxAl){
		
		using namespace std;
		using namespace seqan;
		using namespace flexbar;
		
		// int band1 = overhang;
		// int band2 = readLen - minOvl;
		
		// appendValue(alignments.ascores, 0);
		// AlignConfig<true, true, true, true> ac;
		
		// alignments.ascores[idxAl] = globalAlignment(alignments.aset[idxAl], m_scoreMatrix, ac, band1, band2);
		
		
		if(cycle == COMPUTE){
			
			cycle = RESULTS;
			
			if(m_trimEnd == RIGHT || m_trimEnd == RTAIL){
				
				AlignConfig<true, false, true, true> ac;
				alignments.ascores = globalAlignment(alignments.aset, m_scoreMatrix, ac);
			}
			else if(m_trimEnd == LEFT || m_trimEnd == LTAIL){
				
				AlignConfig<true, true, false, true> ac;
				alignments.ascores = globalAlignment(alignments.aset, m_scoreMatrix, ac);
			}
			else{
				AlignConfig<true, true, true, true> ac;
				alignments.ascores = globalAlignment(alignments.aset, m_scoreMatrix, ac);
			}
		}
		
		TAlign &align = alignments.aset[idxAl];
		a.score       = alignments.ascores[idxAl];
		
		// cout << "Score: " << a.score << endl;
		// cout << "Align: " << align << endl;
		
		
		TRow &row1 = row(align, 0);
		TRow &row2 = row(align, 1);
		
		a.startPosS = toViewPosition(row1, 0);
		a.startPosA = toViewPosition(row2, 0);
		a.endPosS   = toViewPosition(row1, length(source(row1)));
		a.endPosA   = toViewPosition(row2, length(source(row2)));
		
		a.startPos = (a.startPosA > a.startPosS) ? a.startPosA : a.startPosS;
		a.endPos   = (a.endPosA   > a.endPosS)   ? a.endPosS   : a.endPosA;
		
		// cout << startPosS << endl << startPosA << endl;
		// cout << endPosS   << endl << endPosA   << endl;
		
		
		if(m_log != NONE){
			stringstream s;
			s << align;
			a.alString = s.str();
		}
		
		if(m_randTag) a.randTag = "";
		
		
		TRowIterator it1 = begin(row1);
		TRowIterator it2 = begin(row2);
		
		int alPos    = 0;
		a.gapsR      = 0;
		a.gapsA      = 0;
		a.mismatches = 0;
		
		for(; it1 != end(row1); ++it1){
			
			if(a.startPos <= alPos && alPos < a.endPos){
				     if(isGap(it1))                   ++a.gapsR;
				else if(isGap(it2))                   ++a.gapsA;
				else if(*it1 != *it2 && *it2 != 'N')  ++a.mismatches;
				else if(m_randTag    && *it2 == 'N')  append(a.randTag, (TChar) *it1);
			}
			++alPos;
			++it2;
		}
		
		// cout << gapsR << endl << gapsA << endl << mismatches << endl;
	}
	
	
	void printScoreMatrix(TScoreMatrix &scoreMatrix){
		
		using namespace std;
		using namespace seqan;
		
		cout << endl;
		for(unsigned i = 0; i < ValueSize<TChar>::VALUE; ++i)
			cout << "\t" << TChar(i);
		cout << endl;
		
		for(unsigned i = 0; i < ValueSize<TChar>::VALUE; ++i){
			cout << TChar(i);
			for(unsigned j = 0; j < ValueSize<TChar>::VALUE; ++j)
				cout << "\t" << score(scoreMatrix, TChar(i), TChar(j));
			cout << endl;
		}
	}
};


#endif
