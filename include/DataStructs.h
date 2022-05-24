#ifndef DATA_STRUCTS_H
#define DATA_STRUCTS_H

#include <vector>

struct SabrePair
{
	int ringch, wedgech; //global
	int detID;
	int local_ring, local_wedge; //local (obv.)
	double ringE, wedgeE;
	double ringT, wedgeT;
};

struct CalEvent
{
	double xavg = -1e6, x1 = -1e6, x2 = -1e6;
	double theta = -1, cathodeE = -1, scintE = -1;
	double anodeFrontE = -1, anodeBackE = -1;
	double scintT = -1;
	std::vector<SabrePair> sabre;
};

bool EnforceDictionaryLinked();

#endif