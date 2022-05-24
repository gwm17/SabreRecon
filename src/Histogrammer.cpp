#include "Histogrammer.h"
#include "CalDict/DataStructs.h"
#include <fstream>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>


namespace SabreRecon {

	Histogrammer::Histogrammer(const std::string& input) :
		m_inputData(""), m_outputData(""), m_eventPtr(new CalEvent), m_isValid(false)
	{
		TH1::AddDirectory(kFALSE);
		ParseConfig(input);
	}

	Histogrammer::~Histogrammer() {}

	void Histogrammer::ParseConfig(const std::string& name)
	{
		std::ifstream input(name);
		if(!input.is_open())
		{
			m_isValid = false;
			return;
		}

		std::string junk;

		double B, theta;
		std::vector<double> fpCal;

		double thickness;
		std::vector<int> targ_z;
		std::vector<int> targ_s;

		std::vector<ReconCut> cuts;
		ReconCut this_cut;

		input>>junk;
		if(junk == "begin_data")
		{
			input>>junk>>m_inputData;
			input>>junk>>m_outputData;
			input>>junk;
		}
		else
		{
			m_isValid = false;
			return;
		}

		input>>junk;
		if(junk == "begin_reconstructor")
		{
			while(input>>junk)
			{
				if(junk == "begin_focalplane")
				{
					input>>junk>>B;
					input>>junk>>theta;
					input>>junk;
					if(junk == "begin_fpcal")
					{
						while(input>>junk)
						{
							if(junk == "end_fpcal")
								break;
							else
								fpCal.push_back(std::stod(junk));
						}
					}
				}
				else if(junk == "begin_target")
				{
					input>>junk>>thickness;
					input>>junk;
					if(junk == "begin_elements")
					{
						while(input>>junk)
						{
							if(junk == "end_elements")
								break;
							else
							{
								targ_z.push_back(std::stoi(junk));
								input>>junk;
								targ_s.push_back(std::stoi(junk));
							}
						}
					}
				}
				else if(junk == "end_fpcal")
					continue;
				else if(junk == "end_target")
					continue;
				else
					break;
			}
		}

		input>>junk;
		if(junk == "begin_cuts")
		{
			while(input>>junk)
			{
				if(junk == "end_cuts")
					break;
				else
				{
					this_cut.cutname = junk;
					input>>this_cut.filename;
					input>>this_cut.xparam;
					input>>this_cut.yparam;
					cuts.push_back(this_cut);
				}
			}
		}

		//init resources
		Target target(targ_z, targ_s, thickness);
		m_recon.Init(target, theta, B, fpCal);
		m_cuts.InitCuts(cuts);
		m_cuts.InitEvent(m_eventPtr);

	}

	void Histogrammer::FillHistogram1D(const Histogram1DParams& params, double value)
	{
		auto h = std::static_pointer_cast<TH1>(m_histoMap[params.name]);
		if(h)
			h->Fill(value);
		else
		{
			h = std::make_shared<TH1F>(params.name.c_str(), params.title.c_str(), params.bins, params.min, params.max);
			h->Fill(value);
			m_histoMap[params.name] = h;
		}
	}

	void Histogrammer::FillHistogram2D(const Histogram2DParams& params, double valueX, double valueY)
	{
		auto h = std::static_pointer_cast<TH2>(m_histoMap[params.name]);
		if(h)
			h->Fill(valueX, valueY);
		else
		{
			h = std::make_shared<TH2F>(params.name.c_str(), params.title.c_str(), params.binsX, params.minX, params.maxX, params.binsY, params.minY, params.maxY);
			h->Fill(valueX, valueY);
			m_histoMap[params.name] = h;
		}
	}


	void Histogrammer::Run()
	{

	}
}