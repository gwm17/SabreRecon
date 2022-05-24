#include "RandomGenerator.h"

namespace SabreRecon {

	RandomGenerator::RandomGenerator() 
	{
		std::random_device rd;
		rng.seed(rd());
	}

	RandomGenerator::~RandomGenerator() {}
}