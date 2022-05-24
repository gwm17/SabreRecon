#include "Histogrammer.h"
#include <iostream>
#include <string>


int main(int argc, const char** argv)
{
	if(argc != 2)
	{
		std::cerr<<"SabreRecon requires an input config file! Unable to run."<<std::endl;
		return 1;
	}

	SabreRecon::Histogrammer grammer(argv[1]);

	if(grammer.IsValid())
		grammer.Run();
	else
	{
		std::cerr<<"Config file error! Shutting down."<<std::endl;
		return 1;
	}

	return 0;

}