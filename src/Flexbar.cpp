/*==================================================
   
   Flexbar - flexible barcode and adapter removal
   
   Version 2.7
   
   uses SeqAn library release 2.1.1
   and TBB library 4.0 or later
   
   
           Developer:  Johannes Roehr
   
   Former developers:  Matthias Dodt
                       Benjamin Menkuec
                       Sebastian Roskosch
   
   
   https://github.com/seqan/flexbar
   
===================================================*/

#include "Flexbar.h"


int main(int argc, const char* argv[]){
	
	using namespace std;
	using namespace flexbar;
	
	using seqan::ArgumentParser;
	
	const string version = "2.7";
	const string date    = "April 2016";
	
	ArgumentParser parser("flexbar");
	
	defineOptionsAndHelp(parser, version, date);
	parseCommandLine(parser, version, argc, argv);
	
	Options o;
	initOptions(o, parser);
	
	loadProgramOptions(o, parser);
	loadBarcodesAndAdapters(o);
	
	startComputation(o);
	printCompletedMessage(o);
	
	return 0;
}
