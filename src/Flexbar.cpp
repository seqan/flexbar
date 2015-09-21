/*==================================================
   
   Flexbar - flexible barcode and adapter removal
   
   Version 2.5  (GPLv3)
   
   Uses seqan library release 1.4 or later
    and tbb library 4.0 or later
   
   Authors: Matthias Dodt and Johannes Roehr
   
===================================================*/


#include "Flexbar.h"
#include "Options.h"
#include "Enums.h"


int main(int argc, const char* argv[]){
	
	using namespace std;
	using namespace flexbar;
	
	using seqan::ArgumentParser;
	
	const string version = "2.5";
	const string date    = "June 30, 2014";
	
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

