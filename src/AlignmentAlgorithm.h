/*
 *   AlignmentAlgorithm.h
 *
 *   Authors: mat and jtr
 */

#ifndef FLEXBAR_ALIGNMENTALGORITHM_H_
#define FLEXBAR_ALIGNMENTALGORITHM_H_

#include <iostream>
#include <string>

#include <seqan/basic.h>
#include <seqan/align.h>


template <typename TString>
class AlignmentAlgorithm {

private:
	
	typedef typename seqan::Dna5 TChar;
	
	typedef typename seqan::Value<TString>::Type TStringChar;
	// typedef seqan::SimpleType<unsigned char, seqan::Finite<5> > TChar;
	
	typedef seqan::Align<TString, seqan::ArrayGaps> TAlign;
	typedef typename seqan::Row<TAlign>::Type TRow;
	typedef typename seqan::Iterator<TRow>::Type TRowIterator;
	
	typedef seqan::Score<int, seqan::ScoreMatrix<TChar, seqan::Default> > TScoreDna5;
	
	TScoreDna5 m_scoreDna5;
	seqan::Score<int> m_score;
	
	const bool m_isColorSpace, m_randTag;
	const flexbar::LogLevel m_verb;
	const flexbar::TrimEnd m_trimEnd;
	
public:
	
	AlignmentAlgorithm(const Options &o, const int match, const int mismatch, const int gapCost, const flexbar::TrimEnd trimEnd):
			m_randTag(o.randTag),
			m_isColorSpace(o.isColorSpace),
			m_verb(o.logLevel),
			m_trimEnd(trimEnd){
		
		using namespace std;
		using namespace seqan;
		
		m_score = Score<int>(match, mismatch, gapCost);
		
		m_scoreDna5 = TScoreDna5(gapCost);
		
		for (unsigned i = 0; i < ValueSize<TChar>::VALUE; ++i){
			for (unsigned j = 0; j < ValueSize<TChar>::VALUE; ++j){
				
				if(i == j || TChar(i) == 'N' || TChar(j) == 'N'){
					
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
	
	
	virtual ~AlignmentAlgorithm(){
	};
	
	
	void align(const TString &querySeq, const TString &read, int &gapsR, int &gapsA, int &mismatches, int &startPos, int &endPos, int &startPosA, int &endPosA, int &startPosS, int &endPosS, int &aliScore, std::stringstream &aliString, TString &tagSeq){
		
		using namespace std;
		using namespace seqan;
		using namespace flexbar;
		
		TAlign align;
		resize(rows(align), 2);
		assignSource(row(align, 0), read);
		assignSource(row(align, 1), querySeq);
		
		
		if(m_trimEnd == RIGHT || m_trimEnd == RIGHT_TAIL){
			AlignConfig<true, false, true, true> ac;
			if(m_isColorSpace) aliScore = globalAlignment(align, m_score,     ac);
			else               aliScore = globalAlignment(align, m_scoreDna5, ac);
		}
		else if(m_trimEnd == LEFT  || m_trimEnd == LEFT_TAIL){
			AlignConfig<true, true, false, true> ac;
			if(m_isColorSpace) aliScore = globalAlignment(align, m_score,     ac);
			else               aliScore = globalAlignment(align, m_scoreDna5, ac);
		}
		else{
			AlignConfig<true, true, true, true> ac;
			if(m_isColorSpace) aliScore = globalAlignment(align, m_score,     ac);
			else               aliScore = globalAlignment(align, m_scoreDna5, ac);
		}
		
		
		TRow &row1 = row(align, 0);
		TRow &row2 = row(align, 1);
		
		startPosS = toViewPosition(row1, 0);
		startPosA = toViewPosition(row2, 0);
		endPosS   = toViewPosition(row1, length(source(row1)));
		endPosA   = toViewPosition(row2, length(source(row2)));
		
		// calculate overlap start and end
		if(startPosA > startPosS) startPos = startPosA;
		else                      startPos = startPosS;
		
		if(endPosA > endPosS) endPos = endPosS;
		else                  endPos = endPosA;
		
		
		// cout << endl << endl << startPosS << endl << startPosA << endl << endPosS << endl << endPosA;
		
		// int fstartPosS = toViewPosition(row1, 0);
		// int fstartPosA = toViewPosition(row2, 0);
		// int fendPosS   = toViewPosition(row1, length(source(row1)));
		// int fendPosA   = toViewPosition(row2, length(source(row2)));
		// cout << endl << endl << fstartPosS << endl << fstartPosA << endl << fendPosS << endl << fendPosA;
		
		// cout << align << endl << aliScore << endl;
		
		if(m_verb != flexbar::NONE) aliString << align;
		
		
		// compute number of mismatches and gaps
		TRowIterator it1 = begin(row1);
		TRowIterator it2 = begin(row2);
		
		int aliPos = 0;
		gapsR      = 0;
		gapsA      = 0;
		mismatches = 0;
		
		for(; it1 != end(row1); ++it1){
			
			if(startPos <= aliPos && aliPos < endPos){
				     if(isGap(it1))                                  ++gapsR;
				else if(isGap(it2))                                  ++gapsA;
				else if(*it1 != *it2 && *it1 != 'N' && *it2 != 'N')  ++mismatches;
				else if(m_randTag    && *it2 == 'N')                 append(tagSeq, (TStringChar) *it1);
			}
			++aliPos;
			++it2;
		}
		
		// cout << endl << endl << gapsR << endl << gapsA << endl << mismatches << endl << align;
	}
};


#endif /* FLEXBAR_ALIGNMENTALGORITHM_H_ */
