/*
 *   SeqAlignAlgorithm.h
 *
 */

#ifndef FLEXBAR_SEQALIGNALGORITHM_H
#define FLEXBAR_SEQALIGNALGORITHM_H


template <typename TSeqStr>
class SeqAlignAlgorithm {

private:
	
	typedef typename seqan::Dna5 TChar;
	typedef typename seqan::Value<TSeqStr>::Type TSeqStrChar;
	
	typedef typename seqan::Row<flexbar::TAlign>::Type TRow;
	typedef typename seqan::Iterator<TRow>::Type TRowIterator;
	
	typedef seqan::Score<int, seqan::ScoreMatrix<TChar, seqan::Default> > TScoreDna5;
	
	TScoreDna5 m_scoreDna5;
	
	const bool m_randTag;
	const flexbar::LogLevel m_verb;
	const flexbar::TrimEnd m_trimEnd;
	
public:
	
	SeqAlignAlgorithm(const Options &o, const int match, const int mismatch, const int gapCost, const flexbar::TrimEnd trimEnd):
			m_randTag(o.randTag),
			m_verb(o.logLevel),
			m_trimEnd(trimEnd){
		
		using namespace std;
		using namespace seqan;
		
		m_scoreDna5 = TScoreDna5(gapCost);
		
		for (unsigned i = 0; i < ValueSize<TChar>::VALUE; ++i){
			for (unsigned j = 0; j < ValueSize<TChar>::VALUE; ++j){
				
				if(i == j || TChar(j) == 'N'){
					setScore(m_scoreDna5, TChar(i), TChar(j), match);
				}
				else{
					setScore(m_scoreDna5, TChar(i), TChar(j), mismatch);
				}
				
				// cout << i << "\t" << TChar(i) << endl;
				// cout << j << "\t" << TChar(j) << endl;
				// cout << ValueSize<TChar>::VALUE << endl << endl;
			}
		}
		
		// cout << endl;
		// for (unsigned i = 0; i < ValueSize<TChar>::VALUE; ++i)
		// 	cout << "\t" << TChar(i);
		// cout << endl;
		// 
		// for (unsigned i = 0; i < ValueSize<TChar>::VALUE; ++i) {
		// 	cout << TChar(i);
		// 	for (unsigned j = 0; j < ValueSize<TChar>::VALUE; ++j){
		// 		cout << "\t" << score(m_scoreDna5, TChar(i), TChar(j));
		// 	}
		// 	cout << endl;
		// }
	};
	
	
	virtual ~SeqAlignAlgorithm(){
	};
	
	
	void alignGlobal(const TSeqStr &querySeq, const TSeqStr &readSeq, int &gapsR, int &gapsA, int &mismatches, int &startPos, int &endPos, int &startPosA, int &endPosA, int &startPosS, int &endPosS, int &aliScore, std::stringstream &aliString, TSeqStr &tagSeq, flexbar::TAlignments &alignments, const flexbar::ComputeCycle cycle){
		
		using namespace std;
		using namespace seqan;
		using namespace flexbar;
		
		if(cycle != PRECYCLE){
			
		}
		
		TAlign align;
		resize(rows(align), 2);
		assignSource(row(align, 0), readSeq);
		assignSource(row(align, 1), querySeq);
		
		
		StringSet<TAlign> aligns;
		
		appendValue(aligns, align);
		appendValue(aligns, align);
		// value(aligns, 0)
		
		AlignConfig<true, false, true, true> ac;
		
		String<int> results = globalAlignment(aligns, m_scoreDna5, ac);
		
		cout << results[0] << endl << results[1] << endl << endl;
		cout << value(aligns, 0) << endl << endl;
		cout << value(aligns, 1) << endl << endl;
		
		
		if(m_trimEnd == RIGHT || m_trimEnd == RIGHT_TAIL){
			
			AlignConfig<true, false, true, true> ac;
			aliScore = globalAlignment(align, m_scoreDna5, ac);
		}
		else if(m_trimEnd == LEFT || m_trimEnd == LEFT_TAIL){
			
			AlignConfig<true, true, false, true> ac;
			aliScore = globalAlignment(align, m_scoreDna5, ac);
		}
		else{
			AlignConfig<true, true, true, true> ac;
			aliScore = globalAlignment(align, m_scoreDna5, ac);
		}
		
		
		TRow &row1 = row(align, 0);
		TRow &row2 = row(align, 1);
		
		startPosS = toViewPosition(row1, 0);
		startPosA = toViewPosition(row2, 0);
		endPosS   = toViewPosition(row1, length(source(row1)));
		endPosA   = toViewPosition(row2, length(source(row2)));
		
		if(startPosA > startPosS) startPos = startPosA;
		else                      startPos = startPosS;
		
		if(endPosA > endPosS) endPos = endPosS;
		else                  endPos = endPosA;
		
		
		// cout << "\n\n" << startPosS << endl << startPosA << endl << endPosS << endl << endPosA;
		// cout << align << endl << aliScore << endl;
		
		if(m_verb != flexbar::NONE) aliString << align;
		
		
		TRowIterator it1 = begin(row1);
		TRowIterator it2 = begin(row2);
		
		int aliPos = 0;
		gapsR      = 0;
		gapsA      = 0;
		mismatches = 0;
		
		for(; it1 != end(row1); ++it1){
			
			if(startPos <= aliPos && aliPos < endPos){
				     if(isGap(it1))                   ++gapsR;
				else if(isGap(it2))                   ++gapsA;
				else if(*it1 != *it2 && *it2 != 'N')  ++mismatches;
				else if(m_randTag    && *it2 == 'N')  append(tagSeq, (TSeqStrChar) *it1);
			}
			++aliPos;
			++it2;
		}
		
		// cout << "\n\n" << gapsR << endl << gapsA << endl << mismatches << endl << align;
	}
};


#endif
