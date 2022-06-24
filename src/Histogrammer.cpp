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
		//As of Ubuntu 22.04: linker aggresively strips (particularly for release) to the point where non-header symbols MUST be called AND used
		//within the executable. Use EnforceDictionaryLinked as shown.
		if(EnforceDictionaryLinked())
		{
			std::cout<<"Enforcing dictionary linking"<<std::endl;
		}
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
		std::vector<int> targ_a;
		std::vector<int> targ_s;

		std::vector<ReconCut> cuts;
		ReconCut this_cut;

		std::vector<std::string> ptables;
		std::vector<std::string> etables;

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
								targ_a.push_back(std::stoi(junk));
								input>>junk;
								targ_s.push_back(std::stoi(junk));
								std::cout<<"e"<<targ_s.size()-1<<": ("<<targ_z[targ_z.size()-1]<<","<<targ_a[targ_z.size()-1]<<","<<targ_s[targ_z.size()-1]<<") ";
							}
						}
						std::cout<<std::endl;
					}
				}
				else if(junk == "begin_punchtables")
				{
					std::cout<<"Looking for PunchTables..."<<std::endl;
					while(input>>junk)
					{
						if(junk == "end_punchtables")
							break;
						ptables.push_back(junk);
						std::cout<<"Adding PunchTable: "<<junk<<std::endl;
					}
				}
				else if(junk == "begin_elosstables")
				{
					std::cout<<"Looking for ElossTables..."<<std::endl;
					while(input>>junk)
					{
						if(junk == "end_elosstables")
							break;
						etables.push_back(junk);
						std::cout<<"Adding ElossTable: "<<junk<<std::endl;
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
		Target target(targ_a, targ_z, targ_s, thickness);
		m_recon.Init(target, theta, B, fpCal);
		for(auto& table : ptables)
			m_recon.AddPunchThruTable(table);
		for(auto& table : etables)
			m_recon.AddEnergyLossTable(table);
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
					if(m_eventPtr->sabre[0].detID == 0 || m_eventPtr->sabre[0].detID == 1 || m_eventPtr->sabre[0].detID == 4)
					{
						RunDegradedSabre();
					}
					else
					{
						RunSabre();
					}
				}
			}	
		}
		std::cout<<std::endl;
		input->Close();
		output->cd();
		for(auto& gram : m_histoMap)
			gram.second->Write(gram.second->GetName(), TObject::kOverwrite);
		output->Close();
	}

	void Histogrammer::RunSabre()
	{
		static ReconResult recon5Li, recon7Be, recon8Be, recon14N, recon9B;
		static TVector3 sabreCoords, b9Coords;
		static double relAngle;

		auto& biggestSabre = m_eventPtr->sabre[0];
		recon9B = m_recon.RunFPResidExcitation(m_eventPtr->xavg, m_beamKE, {{5,10},{2,3},{3,4}});
		recon5Li = m_recon.RunSabreExcitation(m_eventPtr->xavg, m_beamKE, biggestSabre, {{5,10},{2,3},{2,4},{2,4}});
		recon8Be = m_recon.RunSabreExcitation(m_eventPtr->xavg, m_beamKE, biggestSabre, {{5,10},{2,3},{2,4},{1,1}});
		recon7Be = m_recon.RunSabreExcitation(m_eventPtr->xavg, m_beamKE, biggestSabre, {{5,10},{2,3},{2,4},{1,2}});
		recon14N = m_recon.RunSabreExcitation(m_eventPtr->xavg, m_beamKE, biggestSabre, {{8,16},{2,3},{2,4},{1,1}});
		sabreCoords = m_recon.GetSabreCoordinates(biggestSabre);
		b9Coords.SetMagThetaPhi(1.0, recon9B.residThetaLab, recon9B.residPhiLab);
		relAngle = std::acos(b9Coords.Dot(sabreCoords)/(sabreCoords.Mag()*b9Coords.Mag()));

		FillHistogram1D({"xavg_gated_sabre","xavg_gated_sabre;xavg;counts",600,-300.0,300.0}, m_eventPtr->xavg);
		FillHistogram2D({"scintE_cathodeE","scintE_cathodeE;scintE;cathodeE",512,0,4096,512,0,4096}, m_eventPtr->scintE, m_eventPtr->cathodeE);
		FillHistogram2D({"xavg_theta","xavg_theta;xavg;theta",600,-300.0,300.0,500,0.0,1.5}, m_eventPtr->xavg, m_eventPtr->theta);

		FillHistogram1D({"ex_5Li", "ex_5Li;E_x(MeV);counts",3000,-5.0,25.0}, recon5Li.excitation);
		FillHistogram1D({"ex_7Be", "ex_7Be;E_x(MeV);counts",3000,-20.0,10.0}, recon7Be.excitation);
		FillHistogram1D({"ex_8Be", "ex_8Be;E_x(MeV);counts",3000,-5.0,25.0}, recon8Be.excitation);
		FillHistogram1D({"ex_14N", "ex_14N;E_x(MeV);counts",3000,-20.0,10.0}, recon14N.excitation);
		FillHistogram2D({"ex_14N_7Be","ex_14N_7Be;E_x 14N;E_x 7Be",500,-10.0,10.0,500,-10.,10.0}, recon14N.excitation, recon7Be.excitation);
		FillHistogram2D({"sabreTheta_sabreE","sabreTheta_sabreE;#theta (deg); E(MeV)",180,0,180,400,0,20.0},sabreCoords.Theta()*s_rad2deg, biggestSabre.ringE);
		FillHistogram2D({"xavg_sabreE","xavg_sabreE;xavg; E(MeV)",600,-300.0,300.0,400,0,20.0},m_eventPtr->xavg, biggestSabre.ringE);
		FillHistogram2D({"9Btheta_sabreTheta","9Btheta_sabreTheta;#theta_{9B};#theta_{SABRE}",180,0.0,180.0,180,0.0,180.0}, recon9B.residThetaLab*s_rad2deg, sabreCoords.Theta()*s_rad2deg);
		FillHistogram2D({"sabreE_relAngle","sabreE_relAngle;#theta_{rel};E(MeV)",180,0.0,180.0,400,0.0,20.0},relAngle*s_rad2deg,biggestSabre.ringE);
		FillHistogram2D({"sabreTheta_5Liex","sabreTheta_5Liex;#theta (deg);E_x (MeV)",180,0.0,180.0,1000,-5.0,25.0},sabreCoords.Theta()*s_rad2deg,recon5Li.excitation);
		FillHistogram2D({"sabreTheta_7Beex","sabreTheta_7Beex;#theta (deg);E_x (MeV)",180,0.0,180.0,1000,-20.0,10.0},sabreCoords.Theta()*s_rad2deg,recon7Be.excitation);
		FillHistogram2D({"sabreTheta_14Nex","sabreTheta_14Nex;#theta (deg);E_x (MeV)",180,0.0,180.0,1000,-20.0,10.0},sabreCoords.Theta()*s_rad2deg,recon14N.excitation);
		FillHistogram2D({"sabrePhi_5Liex","sabrePhi_5Liex;#phi (deg);E_x (MeV)",360,0.0,360.0,1000,-5.0,25.0},Phi360(sabreCoords.Phi())*s_rad2deg,recon5Li.excitation);
		FillHistogram2D({"sabrePhi_7Beex","sabrePhi_7Beex;#phi (deg);E_x (MeV)",360,0.0,360.0,1000,-20.0,10.0},Phi360(sabreCoords.Phi())*s_rad2deg,recon7Be.excitation);
		FillHistogram2D({"sabrePhi_14Nex","sabrePhi_14Nex;#phi (deg);E_x (MeV)",360,0.0,360.0,1000,-20.0,10.0},Phi360(sabreCoords.Phi())*s_rad2deg,recon14N.excitation);

		if(m_eventPtr->xavg > -186.0 && m_eventPtr->xavg < -178.0) //nub
		{
			FillHistogram2D({"sabreE_sabreTheta_nub","sabreE_sabreTheta_nub;#theta (deg);E(MeV)",180,0.0,180.0,400,0.0,20.0},sabreCoords.Theta()*s_rad2deg,biggestSabre.ringE);
			FillHistogram2D({"sabreE_sabrePhi_nub","sabreE_sabreTheta_nub;#phi (deg);E(MeV)",360,0.0,360.0,400,0.0,20.0},Phi360(sabreCoords.Phi())*s_rad2deg,biggestSabre.ringE);
		}
		else if(m_eventPtr->xavg > -195.0 && m_eventPtr->xavg < -185.0) //Nabin peak
		{
			FillHistogram2D({"sabreTheta_5Liex_nabinPeak","sabreTheta_5Liex_nabinPeak;#theta (deg);E_x (MeV)",180,0.0,180.0,1000,-5.0,25.0},sabreCoords.Theta()*s_rad2deg,recon5Li.excitation);
			FillHistogram2D({"sabreTheta_7Beex_nabinPeak","sabreTheta_7Beex_nabinPeak;#theta (deg);E_x (MeV)",180,0.0,180.0,1000,-20.0,10.0},sabreCoords.Theta()*s_rad2deg,recon7Be.excitation);
			FillHistogram2D({"sabrePhi_5Liex_nabinPeak","sabrePhi_5Liex_nabinPeak;#phi (deg);E_x (MeV)",360,0.0,360.0,1000,-5.0,25.0},Phi360(sabreCoords.Phi())*s_rad2deg,recon5Li.excitation);
			FillHistogram2D({"sabrePhi_7Beex_nabinPeak","sabrePhi_7Beex_nabinPeak;#phi (deg);E_x (MeV)",360,0.0,360.0,1000,-20.0,10.0},Phi360(sabreCoords.Phi())*s_rad2deg,recon7Be.excitation);
		}

		//Gate on reconstr. excitation structures; overlaping cases are possible!
		if(recon5Li.excitation > -2.0 && recon5Li.excitation < 2.0)
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
			FillHistogram2D({"xavg_sabreE_7Begs","xavg_sabreE_7Begs;xavg;E(MeV)",600,-300.0,300.0,400,0.0,20.0}, m_eventPtr->xavg, biggestSabre.ringE);
			if(!(recon14N.excitation > -0.1 && recon14N.excitation < 2.0))
				FillHistogram1D({"xavg_gated7Begs_reject14Ngs", "xavg_gated7Begs_reject14Ngs;xavg;counts",600, -300.0, 300.0}, m_eventPtr->xavg);
			if(m_eventPtr->xavg > -186.0 && m_eventPtr->xavg < -178.0)
			{
				FillHistogram2D({"sabreE_sabreTheta_7begs_nub","sabreE_sabreTheta_7begs_nub;#theta (deg);E(MeV)",180,0.0,180.0,400,0.0,20.0},sabreCoords.Theta()*s_rad2deg,biggestSabre.ringE);
				FillHistogram2D({"sabreE_sabrePhi_7begs_nub","sabreE_sabreTheta_7begs_nub;#phi (deg);E(MeV)",360,0.0,360.0,400,0.0,20.0},Phi360(sabreCoords.Phi())*s_rad2deg,biggestSabre.ringE);
			}
		}
		if(recon14N.excitation > -0.1 && recon14N.excitation < 0.2)
		{
			FillHistogram1D({"xavg_gated14Ngs", "xavg_gated14Ngs;xavg;counts",600,-300.0,300.0}, m_eventPtr->xavg);
		}
		if(!(recon14N.excitation > -0.1 && recon14N.excitation < 0.2) && !(recon7Be.excitation > -0.1 && recon7Be.excitation < 0.15)
			&& !(recon8Be.excitation > -0.1 && recon8Be.excitation < 0.1) && !(recon5Li.excitation > -2.0 && recon5Li.excitation < 2.0))
		{
			FillHistogram1D({"xavg_notGatedAllChannels", "xavg_notGatedAllChannels;xavg;counts",600,-300.0,300.0}, m_eventPtr->xavg);
		}
	}

	void Histogrammer::RunDegradedSabre()
	{
		static ReconResult recon8Be, recon8BePunch, recon9Be, recon9BePunch, recon9B;
		static TVector3 sabreCoords, b9Coords;
		static double relAngle;

		auto& biggestSabre = m_eventPtr->sabre[0];

		recon8Be = m_recon.RunSabreExcitation(m_eventPtr->xavg, m_beamKE, biggestSabre, {{5,10},{2,3},{2,4},{1,1}});
		recon8BePunch = m_recon.RunSabreExcitationPunchDegraded(m_eventPtr->xavg, m_beamKE, biggestSabre, {{5,10},{2,3},{2,4},{1,1}});
		recon9Be = m_recon.RunSabreExcitation(m_eventPtr->xavg, m_beamKE, biggestSabre, {{5,11},{2,3},{2,4},{1,1}});
		recon9BePunch = m_recon.RunSabreExcitationPunchDegraded(m_eventPtr->xavg, m_beamKE, biggestSabre, {{5,11},{2,3},{2,4},{1,1}});
		recon9B = m_recon.RunFPResidExcitation(m_eventPtr->xavg, m_beamKE, {{5,10},{2,3},{3,4}});
		sabreCoords = m_recon.GetSabreCoordinates(biggestSabre);
		b9Coords.SetMagThetaPhi(1.0, recon9B.residThetaLab, recon9B.residPhiLab);
		relAngle = std::acos(b9Coords.Dot(sabreCoords)/(sabreCoords.Mag()*b9Coords.Mag()));

		FillHistogram1D({"ex_8be_degraded","ex_8be_degraded; E_x(MeV); counts",3000,-5.0,25.0}, recon8Be.excitation);
		FillHistogram2D({"xavg_ex8be_degraded","xavg_ex8be_degraded;xavg;E_x(MeV)",600,-300.0,300.0,300,-5.0,25.0}, m_eventPtr->xavg, recon8Be.excitation);
		FillHistogram1D({"ex_8be_degradedPunched","ex_8be_degradedPunched; E_x(MeV); counts",3000,-15.0,15.0}, recon8BePunch.excitation);
		FillHistogram2D({"xavg_ex8be_degradedPunched","xavg_ex8be_degradedPunched;xavg;E_x(MeV)",600,-300.0,300.0,300,-15.0,15.0}, m_eventPtr->xavg, recon8BePunch.excitation);
		FillHistogram1D({"ex_9be_degraded","ex_9be_degraded; E_x(MeV); counts",3000,-5.0,25.0}, recon9Be.excitation);
		FillHistogram2D({"xavg_ex9be_degraded","xavg_ex9be_degraded;xavg;E_x(MeV)",600,-300.0,300.0,300,-5.0,25.0}, m_eventPtr->xavg, recon9Be.excitation);
		FillHistogram1D({"ex_9be_degradedPunched","ex_9be_degradedPunched; E_x(MeV); counts",3000,-15.0,15.0}, recon9BePunch.excitation);
		FillHistogram2D({"xavg_ex9be_degradedPunched","xavg_ex9be_degradedPunched;xavg;E_x(MeV)",600,-300.0,300.0,300,-15.0,15.0}, m_eventPtr->xavg, recon9BePunch.excitation);
		FillHistogram2D({"relAngle_recovSabreKE_8bePunchRecon","relAngle_recovSabreKe;#theta_{rel}(deg);Recovered KE (MeV)",180,0.0,180.0,400,0.0,20.0},relAngle*s_rad2deg,recon8BePunch.sabreRxnKE);
		FillHistogram2D({"relAngle_recovSabreKE_9bePunchRecon","relAngle_recovSabreKe;#theta_{rel}(deg);Recovered KE (MeV)",180,0.0,180.0,400,0.0,20.0},relAngle*s_rad2deg,recon9BePunch.sabreRxnKE);
		if(recon8Be.excitation > 2.2 && recon8Be.excitation < 3.8)
		{
			FillHistogram1D({"xavg_gated_8be1ex_degraded", "xavg_gated_8be1ex_degraded;xavg;counts",600,-300.0,300.0}, m_eventPtr->xavg);
		}

		FillHistogram2D({"sabreE_sabreTheta_degraded", "sabreE_sabreTheta_degraded;#theta (deg);E(MeV)",180,0.0,180.0,400,0.0,20.0}, sabreCoords.Theta()*s_rad2deg, biggestSabre.ringE);
		FillHistogram2D({"xavg_sabreE_degraded", "xavg_sabreE_degraded;xavg;E(MeV)", 600,0.-300.0,300.0,400,0.0,20.0}, m_eventPtr->xavg, biggestSabre.ringE);
		FillHistogram2D({"9Btheta_sabreTheta_degraded","9Btheta_sabreTheta_degraded;#theta_{9B};#theta_{SABRE}",180,0.0,180.0,180,0.0,180.0}, recon9B.residThetaLab*s_rad2deg, sabreCoords.Theta()*s_rad2deg);
		FillHistogram2D({"sabreE_relAngle_degraded","sabreE_relAngle_degraded;#theta_{rel};E(MeV)",180,0.0,180.0,400,0.0,20.0},relAngle*s_rad2deg, biggestSabre.ringE);
		if(biggestSabre.local_wedge != 0 && biggestSabre.local_wedge != 7 && biggestSabre.local_ring != 15 && biggestSabre.local_ring != 0) //Edges might not be degraded right
		{
			FillHistogram2D({"sabreE_sabreTheta_degraded_rejectEdge", "sabreE_sabreTheta_degraded;#theta (deg);E(MeV)",180,0.0,180.0,400,0.0,20.0}, sabreCoords.Theta()*s_rad2deg, biggestSabre.ringE);
			FillHistogram2D({"xavg_sabreE_degraded_rejectEdge", "xavg_sabreE_degraded;xavg;E(MeV)",600,0.-300.0,300.0,400,0.0,20.0}, m_eventPtr->xavg, biggestSabre.ringE);
			FillHistogram1D({"xavg_degraded_rejectEdge","xavg_degraded_rejectEdge;xavg",600,-300.0,300.0}, m_eventPtr->xavg);
			FillHistogram2D({"9Btheta_sabreTheta_degraded_rejectEdge","9Btheta_sabreTheta_degraded_rejectEdge;#theta_{9B};#theta_{SABRE}",180,0.0,180.0,180,0.0,180.0}, recon9B.residThetaLab*s_rad2deg,
							  sabreCoords.Theta()*s_rad2deg);
			FillHistogram2D({"sabreE_relAngle_degraded_rejectEdge","sabreE_relAngle_degraded_rejectEdge;#theta_{rel};E(MeV)",180,0.0,180.0,400,0.0,20.0},relAngle*s_rad2deg, biggestSabre.ringE);
		}
	}
}