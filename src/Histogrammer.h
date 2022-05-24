#ifndef HISTOGRAMMER_H
#define HISTOGRAMMER_H

#include <string>
#include <memory>
#include <unordered_map>
#include <TROOT.h>
#include "CutHandler.h"
#include "Reconstructor.h"

namespace SabreRecon {

	struct Histogram1DParams
	{
		std::string name;
		std::string title;
		int bins;
		double min;
		double max;
	};

	struct Histogram2DParams
	{
		std::string name;
		std::string title;
		int binsX;
		double minX;
		double maxX;
		int binsY;
		double minY;
		double maxY;
	};

	class Histogrammer
	{
	public:
		Histogrammer(const std::string& input);
		~Histogrammer();

		inline const bool IsValid() const { return m_isValid; }
		void Run();

	private:
		void ParseConfig(const std::string& name);
		void FillHistogram1D(const Histogram1DParams& params, double value);
		void FillHistogram2D(const Histogram2DParams& params, double valueX, double valueY);

		std::string m_inputData;
		std::string m_outputData;

		CalEvent* m_eventPtr;

		Reconstructor m_recon;
		CutHandler m_cuts;

		bool m_isValid;

		std::unordered_map<std::string, std::shared_ptr<TObject>> m_histoMap;
	};
}

#endif