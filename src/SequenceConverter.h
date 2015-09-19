/*
 *   SequenceConverter.h
 *
 *   Authors: mat and jtr
 */

#ifndef FLEXBAR_SEQUENCECONVERTER_H_
#define FLEXBAR_SEQUENCECONVERTER_H_


// This class converts sequences from basepair space to colorspace.

template <typename TString>
class SequenceConverter {

private:
	
	static SequenceConverter<TString>* instance;
	
	SequenceConverter(){};
	
public:
	
	static SequenceConverter<TString>* getInstance(){
		if(instance == NULL) instance = new SequenceConverter();
		return instance;
	}
	
	
	TString bpToColorSpace(TString bpSequence){
		
		TString result = "";
		TString substr = "XX";
		
		for(size_t i = 1; i < length(bpSequence); ++i){
			
			substr[0] = bpSequence[i - 1];
			substr[1] = bpSequence[i];
			
			if(substr=="TT") append(result, "0");
			if(substr=="TG") append(result, "1");
			if(substr=="TC") append(result, "2");
			if(substr=="TA") append(result, "3");
			if(substr=="CC") append(result, "0");
			if(substr=="CA") append(result, "1");
			if(substr=="CT") append(result, "2");
			if(substr=="CG") append(result, "3");
			if(substr=="GG") append(result, "0");
			if(substr=="GT") append(result, "1");
			if(substr=="GA") append(result, "2");
			if(substr=="GC") append(result, "3");
			if(substr=="AA") append(result, "0");
			if(substr=="AC") append(result, "1");
			if(substr=="AG") append(result, "2");
			if(substr=="AT") append(result, "3");
		}
		return result;
	}
	
	
	virtual ~SequenceConverter(){};
	
};

template <typename TString> SequenceConverter<TString>* SequenceConverter<TString>::instance = 0;


#endif /* FLEXBAR_SEQUENCECONVERTER_H_ */
