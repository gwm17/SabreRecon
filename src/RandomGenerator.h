#ifndef RANDOMGENERATOR_H
#define RANDOMGENERATOR_H

#include <random>
#include <thread>
#include <mutex>

namespace SabreRecon {

	class RandomGenerator {
	public:
		RandomGenerator();
		~RandomGenerator();
		
		inline std::mt19937& GetGenerator() { return rng; }
		inline static RandomGenerator& GetInstance()
		{
			static RandomGenerator s_generator;
			return s_generator;
		}

	private:
		std::mt19937 rng;
	};

}

#endif