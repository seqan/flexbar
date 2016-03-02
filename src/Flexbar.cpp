/*==================================================
   
   Flexbar - flexible barcode and adapter removal
   
   Version 2.8 alpha
   
   Uses SeqAn library release 2.0.2
    and TBB library 4.0 or later
   
   Authors: Matthias Dodt and Johannes Roehr
   
   https://github.com/seqan/flexbar
   
===================================================*/

#include "Flexbar.h"
#include "Options.h"


int main(int argc, const char* argv[]){
	
	using namespace std;
	using namespace flexbar;
	
	using seqan::ArgumentParser;
	
	const string version = "2.8 alpha";
	const string date    = "February 22, 2016";
	
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
