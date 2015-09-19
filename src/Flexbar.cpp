/*=======================================================================
   Name:         Flexbar.cpp
   Authors:      Matthias Dodt and Johannes Roehr

   Description:  Flexbar - flexible barcode and adapter removal
   Version:      2.4
   Copyright:    GPL version 3

   SeqAn lib:    release 1.4.1, revision 14262 on July 11, 2013
   TBB   lib:    version 4.0 update 5, stable release June 13, 2012
========================================================================*/

#include "Flexbar.h"
#include "Options.h"
#include "Enums.h"


int main(int argc, const char* argv[]){
	
	using namespace std;
	using namespace flexbar;
	
	using seqan::ArgumentParser;
	
	const string version = "2.4";
	const string date    = "July 29, 2013";
	
	ArgumentParser parser("flexbar");
	
	defineOptionsAndHelp(parser, version, date);
	parseCommandLine(parser, version, argc, argv);
	
	Options o;
	initOptions(o, parser);
	
	loadProgramOptions(o, parser, version);
	loadBarcodesAndAdapters(o);
	
	startComputation(o);
	printCompletedMessage(o);
	
	return 0;
}

