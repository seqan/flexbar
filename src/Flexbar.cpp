/*==================================================
   
   Flexbar - flexible barcode and adapter removal
   
   Version 3.0.0
   
   uses SeqAn library release 2.1.1
   and TBB library 4.0 or later
   
   
             Developer:  Johannes Roehr
   
   Former contributors:  Matthias Dodt
                         Benjamin Menkuec
                         Sebastian Roskosch
   
   
   https://github.com/seqan/flexbar
   
===================================================*/

#include "Flexbar.h"


int main(int argc, const char* argv[]){
	
	using namespace std;
	using namespace seqan;
	
	const string version = "3.0 dev";
	const string date    = "June 2016";
	
	ArgumentParser parser("flexbar");
	
	defineOptions(parser, version, date);
	parseCmdLine(parser, version, argc, argv);
	
	Options o;
	
	initOptions(o, parser);
	loadOptions(o, parser);
	
	startComputation(o);
	
	return 0;
}
