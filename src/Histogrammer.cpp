#include "Histogrammer.h"
#include "CalDict/DataStructs.h"
#include <iostream>
#include <fstream>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TTree.h>


namespace SabreRecon {

	double Phi360(double phi)
	{
		return phi < 0 ? (2.0*M_PI + phi) : phi;
	}

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
			input>>junk>>m_beamKE;
			input>>junk;
			std::cout<<"Input datafile: "<<m_inputData<<std::endl;
			std::cout<<"Output datafile: "<<m_outputData<<std::endl;
			std::cout<<"Beam Kinetic Energy (MeV): "<<m_beamKE<<std::endl;
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
					std::cout<<"Found Focal Plane Detector with B(kG): "<<B<<" angle(deg): "<<theta<<std::endl;
					if(junk == "begin_fpcal")
					{
						std::cout<<"FP calibration parameters given: ";
						while(input>>junk)
						{
							if(junk == "end_fpcal")
								break;
							else
							{
								fpCal.push_back(std::stod(junk));
								std::cout<<"a"<<fpCal.size()-1<<": "<<fpCal[fpCal.size()-1]<<" ";
							}
						}
						std::cout<<std::endl;
					}
				}
				else if(junk == "begin_target")
				{
					input>>junk>>thickness;
					input>>junk;
					std::cout<<"Found a target with thickness: "<<thickness<<" ug/cm^2"<<std::endl;
					if(junk == "begin_elements")
					{
						std::cout<<"Target elements given: ";
						while(input>>junk)
						{
							if(junk == "end_elements")
								break;
							else
							{
								targ_z.push_back(std::stoi(junk));
								input>>junk;
								targ_s.push_back(std::stoi(junk));
								std::cout<<"e"<<targ_s.size()-1<<": ("<<targ_z[targ_z.size()-1]<<","<<targ_s[targ_z.size()-1]<<") ";
							}
						}
						std::cout<<std::endl;
					}
				}
				else if(junk == "end_focalplane")
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
					std::cout<<"Added cut "<<this_cut.cutname<<" with file "<<this_cut.filename
							 <<" with xpar "<<this_cut.xparam<<" ypar "<<this_cut.yparam<<std::endl;
				}
			}
		}

		//init resources
		std::cout<<"Initializing resources..."<<std::endl;
		Target target(targ_z, targ_s, thickness);
		m_recon.Init(target, theta, B, fpCal);
		m_cuts.InitCuts(cuts);
		m_cuts.InitEvent(m_eventPtr);

		if(m_cuts.IsValid())
			m_isValid = true;
	}

	void Histogrammer::FillHistogram1D(const Histogram1DParams& params, double value)
	{
		std::shared_ptr<TH1> h = std::static_pointer_cast<TH1>(m_histoMap[params.name]);
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
		std::shared_ptr<TH1> h = std::static_pointer_cast<TH2>(m_histoMap[params.name]);
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
		if(!m_isValid)
		{
			std::cerr<<"ERR -- Resources not initialized properly at Histogrammer::Run()."<<std::endl;
			return;
		}

		TFile* input = TFile::Open(m_inputData.c_str(), "READ");
		if(!input->IsOpen())
		{
			std::cerr<<"ERR -- Unable to open input data file "<<m_inputData<<" at Histogrammer::Run()"<<std::endl;
			return;
		}

		TTree* tree = (TTree*) input->Get("CalTree");
		if(tree == nullptr)
		{
			std::cerr<<"ERR -- No tree named CalTree found in input data file "<<m_inputData<<" at Histogrammer::Run()"<<std::endl;
			return;
		}
		tree->SetBranchAddress("event", &m_eventPtr);

		TFile* output = TFile::Open(m_outputData.c_str(), "RECREATE");
		if(!output->IsOpen())
		{
			std::cerr<<"ERR -- Unable to open output data file "<<m_outputData<<" at Histogrammer::Run()"<<std::endl;
			input->Close();
			return;
		}

		uint64_t nevents = tree->GetEntries();
		float flush_frac = 0.01f;
		uint64_t count = 0, flush_count = 0, flush_val = nevents*flush_frac;

		//Analysis results data
		ReconResult recon5Li, recon7Be, recon8Be, recon14N;
		TVector3 sabreCoords;

		//Temp
		TFile* punchCutFile = TFile::Open("/Volumes/Wyndle/10B3He_May2022/cuts/protonPunchGate_strict.root");
		TCutG* protonGate = (TCutG*) punchCutFile->Get("CUTG");
		protonGate->SetName("protonPunchGate");

		for(uint64_t i=0; i<nevents; i++)
		{
			tree->GetEntry(i);
			count++;
			if(count == flush_val)
			{
				count=0;
				flush_count++;
				std::cout<<"\rPercent of data processed: "<<flush_count*flush_frac*100<<"%"<<std::flush;
			}

			//Only analyze data that passes cuts, has sabre, and passes a weak threshold requirement
			if(m_cuts.IsInside())
			{
				FillHistogram1D({"xavg_gated","xavg_gated;xavg;counts",600,-300.0,300.0}, m_eventPtr->xavg);
				if(!m_eventPtr->sabre.empty() && m_eventPtr->sabre[0].ringE > s_weakSabreThreshold)
				{
					auto& biggestSabre = m_eventPtr->sabre[0];
					recon5Li = m_recon.RunSabreExcitation(m_eventPtr->xavg, m_beamKE, biggestSabre, {{5,10},{2,3},{2,4},{2,4}});
					recon8Be = m_recon.RunSabreExcitation(m_eventPtr->xavg, m_beamKE, biggestSabre, {{5,10},{2,3},{2,4},{1,1}});
					recon7Be = m_recon.RunSabreExcitation(m_eventPtr->xavg, m_beamKE, biggestSabre, {{5,10},{2,3},{2,4},{1,2}});
					recon14N = m_recon.RunSabreExcitation(m_eventPtr->xavg, m_beamKE, biggestSabre, {{8,16},{2,3},{2,4},{1,1}});
					sabreCoords = m_recon.GetSabreCoordinates(biggestSabre);

					FillHistogram1D({"xavg_gated_sabre","xavg_gated_sabre;xavg;counts",600,-300.0,300.0}, m_eventPtr->xavg);
					FillHistogram2D({"scintE_cathodeE","scintE_cathodeE;scintE;cathodeE",512,0,4096,512,0,4096},
									m_eventPtr->scintE, m_eventPtr->cathodeE);
					FillHistogram2D({"xavg_theta","xavg_theta;xavg;theta",600,-300.0,300.0,500,0.0,1.5}, m_eventPtr->xavg, m_eventPtr->theta);

					FillHistogram1D({"ex_5Li", "ex_5Li;E_x(MeV);counts",3000,-5.0,25.0}, recon5Li.excitation);
					FillHistogram1D({"ex_7Be", "ex_7Be;E_x(MeV);counts",3000,-20.0,10.0}, recon7Be.excitation);
					FillHistogram1D({"ex_8Be", "ex_8Be;E_x(MeV);counts",3000,-5.0,25.0}, recon8Be.excitation);
					FillHistogram1D({"ex_14N", "ex_14N;E_x(MeV);counts",3000,-20.0,10.0}, recon14N.excitation);
					FillHistogram2D({"ex_14N_7Be","ex_14N_7Be;E_x 14N;E_x 7Be",500,-10.0,10.0,500,-10.,10.0}, recon14N.excitation,
									recon7Be.excitation);
					FillHistogram2D({"sabreTheta_sabreE","sabreTheta_sabreE;#theta (deg); E(MeV)",180,0,180,400,0,20.0},
									sabreCoords.Theta()*s_rad2deg, biggestSabre.ringE);
					FillHistogram2D({"xavg_sabreE","xavg_sabreE;xavg; E(MeV)",600,-300.0,300.0,400,0,20.0},
									m_eventPtr->xavg, biggestSabre.ringE);

					if(m_eventPtr->xavg > -186.0 && m_eventPtr->xavg < -178.0 && (biggestSabre.detID == 2 || biggestSabre.detID == 3))
					{
						FillHistogram2D({"sabreE_sabreTheta_nub","sabreE_sabreTheta_nub;#theta (deg);E(MeV)",
										180,0.0,180.0,400,0.0,20.0},sabreCoords.Theta()*s_rad2deg,biggestSabre.ringE);
						FillHistogram2D({"sabreE_sabrePhi_nub","sabreE_sabreTheta_nub;#phi (deg);E(MeV)",
										360,0.0,360.0,400,0.0,20.0},Phi360(sabreCoords.Phi())*s_rad2deg,biggestSabre.ringE);
					}

					//Gate on reconstr. excitation structures; overlaping cases are possible!
					if(recon5Li.excitation > -1.5 && recon5Li.excitation < 1.5)
					{
						FillHistogram1D({"xavg_gated5Ligs", "xavg_gated5Ligs;xavg;counts",600,-300.0,300.0}, m_eventPtr->xavg);
					}
					if(recon8Be.excitation > -0.1 && recon8Be.excitation < 0.1)
					{
						FillHistogram1D({"xavg_gated8Begs", "xavg_gated8Begs;xavg;counts",600,-300.0,300.0}, m_eventPtr->xavg);
					}
					if(recon7Be.excitation > -0.1 && recon7Be.excitation < 0.15)
					{
						FillHistogram1D({"xavg_gated7Begs", "xavg_gated7Begs;xavg;counts",600,-300.0,300.0}, m_eventPtr->xavg);
						FillHistogram2D({"xavg_sabreE_7Begs","xavg_sabreE_7Begs;xavg;E(MeV)",600,-300.0,300.0,400,0.0,20.0}, m_eventPtr->xavg,
										biggestSabre.ringE);
						if(m_eventPtr->xavg > -186.0 && m_eventPtr->xavg < -178.0)
						{
							FillHistogram2D({"sabreE_sabreTheta_7begs_nub","sabreE_sabreTheta_7begs_nub;#theta (deg);E(MeV)",
											180,0.0,180.0,400,0.0,20.0},sabreCoords.Theta()*s_rad2deg,biggestSabre.ringE);
							FillHistogram2D({"sabreE_sabrePhi_7begs_nub","sabreE_sabreTheta_7begs_nub;#phi (deg);E(MeV)",
											360,0.0,360.0,400,0.0,20.0},Phi360(sabreCoords.Phi())*s_rad2deg,biggestSabre.ringE);
						}
					}
					if(recon14N.excitation > -0.1 && recon14N.excitation < 0.1)
					{
						FillHistogram1D({"xavg_gated14Ngs", "xavg_gated14Ngs;xavg;counts",600,-300.0,300.0}, m_eventPtr->xavg);
					}
					else
					{
						FillHistogram1D({"xavg_notGated14Ngs", "xavg_notGated14Ngs;xavg;counts",600,-300.0,300.0}, m_eventPtr->xavg);
					}

					//Degrader analysis... SABRE detectors 0, 1, 4 are covered with tantalum
					if(biggestSabre.detID == 0 || biggestSabre.detID == 1 || biggestSabre.detID == 4)
					{
						FillHistogram2D({"sabreE_sabreTheta_degraded", "sabreE_sabreTheta_degraded;#theta (deg);E(MeV)",
										180,0.0,180.0,400,0.0,20.0}, sabreCoords.Theta()*s_rad2deg, biggestSabre.ringE);
						FillHistogram2D({"xavg_sabreE_degraded", "xavg_sabreE_degraded;xavg;E(MeV)",
										600,0.-300.0,300.0,400,0.0,20.0}, m_eventPtr->xavg, biggestSabre.ringE);
						if(biggestSabre.local_wedge != 0 && biggestSabre.local_wedge != 7 && biggestSabre.local_ring != 15 && biggestSabre.local_ring != 0) //Edges might not be degraded right
						{
							FillHistogram2D({"sabreE_sabreTheta_degraded_rejectEdge", "sabreE_sabreTheta_degraded;#theta (deg);E(MeV)",
										180,0.0,180.0,400,0.0,20.0}, sabreCoords.Theta()*s_rad2deg, biggestSabre.ringE);
							FillHistogram2D({"xavg_sabreE_degraded_rejectEdge", "xavg_sabreE_degraded;xavg;E(MeV)",
										600,0.-300.0,300.0,400,0.0,20.0}, m_eventPtr->xavg, biggestSabre.ringE);
							FillHistogram1D({"xavg_degraded_rejectEdge","xavg_degraded_rejectEdge;xavg",600,-300.0,300.0}, m_eventPtr->xavg);
							if(protonGate->IsInside(sabreCoords.Theta()*s_rad2deg, biggestSabre.ringE))
							{
								FillHistogram2D({"xavg_sabreE_degraded_rejectEdge_pGate", "xavg_sabreE_degraded;xavg;E(MeV)",
												600,0.-300.0,300.0,400,0.0,20.0}, m_eventPtr->xavg, biggestSabre.ringE);
								FillHistogram1D({"xavg_degraded_rejectEdge_pGate","xavg_degraded_rejectEdge;xavg",600,-300.0,300.0}, m_eventPtr->xavg);
							}
						}
					}
				}
			}	
		}
		std::cout<<std::endl;
		input->Close();
		output->cd();
		for(auto& gram : m_histoMap)
			gram.second->Write(gram.second->GetName(), TObject::kOverwrite);
		protonGate->Write();
		output->Close();
		punchCutFile->Close();
	}
}