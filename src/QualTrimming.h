// QualTrimming.h

// Authors:  Sebastian Roskosch
//           Benjamin Menkuec
//           Johannes Roehr

#ifndef FLEXBAR_QUALTRIMMING_H
#define FLEXBAR_QUALTRIMMING_H


struct Tail {};

struct BWA {};

struct Window {
	unsigned size;
	Window(unsigned s) : size(s) {}
};


template <typename TString>
inline unsigned getQuality(const TString& qual, unsigned i){
	
	return static_cast<int>(qual[i]);
}


// Tail trimming method
template <typename TString>
unsigned qualTrimming(const TString& qual, unsigned const cutoff, Tail const &){
	
	for (int i = length(qual) - 1; i >= 0; --i){
		
		if(getQuality(qual, i) >= cutoff) return i + 1;
    }
	return 0;
}


// Trim by shifting a window over the seq and cut where avg qual in window turns bad
template <typename TString>
unsigned qualTrimming(const TString& qual, unsigned const _cutoff, Window const & spec){
	
	unsigned window = spec.size;
	unsigned avg = 0, i = 0;
	
	// Absolute cutoff in window to avoid divisions
	unsigned cutoff = _cutoff * window;
	
	// Calculate average quality of initial window
	for (i = 0; i < window; ++i){
		avg += getQuality(qual, i);
    }
	
	// Shift window over read and keep mean quality, update in constant time
	for (i = 0; i < length(qual) && avg >= cutoff; ++i){
		
		// Take care only not to go over the end of the sequence. Shorten window near the end
		avg -= getQuality(qual, i);
		
		if(i + window < length(qual)){
			avg += getQuality(qual, i + window);
		}
		else{
			cutoff = _cutoff * ((length(qual) - 1) - i);
		}
	}
	return i;   // holds start of first window that turned bad
}


// Trimming mechanism using BWA. Trim to argmax_x sum_{i=x+1}^l {cutoff - q_i}
template <typename TString>
unsigned qualTrimming(const TString& qual, unsigned const cutoff, BWA const &){
	
	int max_arg = length(qual) - 1, sum = 0, max = 0;
	
	for (int i = length(qual) - 1; i >= 0; --i){
		
		sum += cutoff - getQuality(qual, i);
		
		if(sum < 0){
			break;
        }
		if(sum > max){
			max = sum;
			max_arg = i;
		}
	}
	return max_arg + 1;
}


template <typename TSeqStr, typename TString>
bool qualTrim(TSeqStr &seq, TString &qual, const flexbar::QualTrimType qtrim, const int cutoff, const int wSize){
	
	using namespace seqan;
	
	unsigned cutPos;
	
	if(qtrim == flexbar::TAIL){
		cutPos = qualTrimming(qual, cutoff, Tail());
	}
	else if(qtrim == flexbar::WIN){
		cutPos = qualTrimming(qual, cutoff, Window(wSize));
	}
	else if(qtrim == flexbar::BWA){
		cutPos = qualTrimming(qual, cutoff, BWA());
	}
	
	if(cutPos < length(qual)){
		
		seq  = prefix(seq,  cutPos);
		qual = prefix(qual, cutPos);
		
		return true;
	}
	else return false;
}


template <typename TSeqStr, typename TString>
bool qualTrim(SeqRead<TSeqStr, TString> *seqRead, const flexbar::QualTrimType qtrim, const int cutoff, const int wSize){
	
	TSeqStr seq  = seqRead->seq;
	TString qual = seqRead->qual;
	
	bool trimmed = qualTrim(seq, qual, qtrim, cutoff, wSize);
	
	if(trimmed){
		seqRead->seq  = seq;
		seqRead->qual = qual;
	}
	
	return trimmed;
}


#endif
