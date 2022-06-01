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

	std::cout<<"---------- SABRE-SPS Invariant Mass Reconstruction Analysis ----------"<<std::endl;
	std::cout<<"Processing configuration file "<<argv[1]<<std::endl;
	SabreRecon::Histogrammer grammer(argv[1]);

	if(grammer.IsValid())
	{
		std::cout<<"Running analysis..."<<std::endl;
		grammer.Run();
	}
	else
	{
		std::cerr<<"ERR -- Config file error! Shutting down."<<std::endl;
		return 1;
	}
	std::cout<<"----------------------------------------------------------------------"<<std::endl;

	return 0;

}