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
					FillHistogram1D({"sabre_counts_gated","sabre_counts_gated;number per event;counts",10,-1.0, 9.0}, m_eventPtr->sabre.size());
					if(m_eventPtr->sabre[0].detID == 0 || m_eventPtr->sabre[0].detID == 1 || m_eventPtr->sabre[0].detID == 4)
					{
						RunDegradedSabre(m_eventPtr->sabre[0]);
					}
					else
					{
						RunSabre(m_eventPtr->sabre[0]);
					}
					/*
					if(m_eventPtr->sabre.size() > 1 && m_eventPtr->sabre[1].ringE > s_weakSabreThreshold &&
					   (m_eventPtr->sabre[1].detID == 0 || m_eventPtr->sabre[1].detID == 1 || m_eventPtr->sabre[1].detID == 4)
					  )
					{
						//Attemp to find some maybe double stuff?
						RunDegradedSabre(m_eventPtr->sabre[1]);
					}
					*/
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

	void Histogrammer::RunSabre(const SabrePair& pair)
	{
		static ReconResult recon5Li, recon7Be, recon8Be, recon14N, recon9B;
		static TVector3 sabreCoords, b9Coords;
		static double relAngle;

		recon9B = m_recon.RunFPResidExcitation(m_eventPtr->xavg, m_beamKE, {{5,10},{2,3},{3,4}});
		recon5Li = m_recon.RunSabreExcitation(m_eventPtr->xavg, m_beamKE, pair, {{5,10},{2,3},{2,4},{2,4}});
		recon8Be = m_recon.RunSabreExcitation(m_eventPtr->xavg, m_beamKE, pair, {{5,10},{2,3},{2,4},{1,1}});
		recon7Be = m_recon.RunSabreExcitation(m_eventPtr->xavg, m_beamKE, pair, {{5,10},{2,3},{2,4},{1,2}});
		recon14N = m_recon.RunSabreExcitation(m_eventPtr->xavg, m_beamKE, pair, {{8,16},{2,3},{2,4},{1,1}});
		sabreCoords = m_recon.GetSabreCoordinates(pair);
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
		FillHistogram2D({"sabreTheta_sabreE","sabreTheta_sabreE;#theta (deg); E(MeV)",180,0,180,400,0,20.0},sabreCoords.Theta()*s_rad2deg, pair.ringE);
		FillHistogram2D({"xavg_sabreE","xavg_sabreE;xavg; E(MeV)",600,-300.0,300.0,400,0,20.0},m_eventPtr->xavg, pair.ringE);
		FillHistogram2D({"9Btheta_sabreTheta","9Btheta_sabreTheta;#theta_{9B};#theta_{SABRE}",180,0.0,180.0,180,0.0,180.0}, recon9B.residThetaLab*s_rad2deg, sabreCoords.Theta()*s_rad2deg);
		FillHistogram2D({"sabreE_relAngle","sabreE_relAngle;#theta_{rel};E(MeV)",180,0.0,180.0,400,0.0,20.0},relAngle*s_rad2deg,pair.ringE);
		FillHistogram2D({"sabreTheta_5Liex","sabreTheta_5Liex;#theta (deg);E_x (MeV)",180,0.0,180.0,1000,-5.0,25.0},sabreCoords.Theta()*s_rad2deg,recon5Li.excitation);
		FillHistogram2D({"sabreTheta_7Beex","sabreTheta_7Beex;#theta (deg);E_x (MeV)",180,0.0,180.0,1000,-20.0,10.0},sabreCoords.Theta()*s_rad2deg,recon7Be.excitation);
		FillHistogram2D({"sabreTheta_14Nex","sabreTheta_14Nex;#theta (deg);E_x (MeV)",180,0.0,180.0,1000,-20.0,10.0},sabreCoords.Theta()*s_rad2deg,recon14N.excitation);
		FillHistogram2D({"sabrePhi_5Liex","sabrePhi_5Liex;#phi (deg);E_x (MeV)",360,0.0,360.0,1000,-5.0,25.0},Phi360(sabreCoords.Phi())*s_rad2deg,recon5Li.excitation);
		FillHistogram2D({"sabrePhi_7Beex","sabrePhi_7Beex;#phi (deg);E_x (MeV)",360,0.0,360.0,1000,-20.0,10.0},Phi360(sabreCoords.Phi())*s_rad2deg,recon7Be.excitation);
		FillHistogram2D({"sabrePhi_14Nex","sabrePhi_14Nex;#phi (deg);E_x (MeV)",360,0.0,360.0,1000,-20.0,10.0},Phi360(sabreCoords.Phi())*s_rad2deg,recon14N.excitation);

		if(m_eventPtr->xavg > -186.0 && m_eventPtr->xavg < -178.0) //nub
		{
			FillHistogram2D({"sabreE_sabreTheta_nub","sabreE_sabreTheta_nub;#theta (deg);E(MeV)",180,0.0,180.0,400,0.0,20.0},sabreCoords.Theta()*s_rad2deg,pair.ringE);
			FillHistogram2D({"sabreE_sabrePhi_nub","sabreE_sabreTheta_nub;#phi (deg);E(MeV)",360,0.0,360.0,400,0.0,20.0},Phi360(sabreCoords.Phi())*s_rad2deg,pair.ringE);
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
			FillHistogram2D({"xavg_sabreE_7Begs","xavg_sabreE_7Begs;xavg;E(MeV)",600,-300.0,300.0,400,0.0,20.0}, m_eventPtr->xavg, pair.ringE);
			if(!(recon14N.excitation > -0.1 && recon14N.excitation < 2.0))
				FillHistogram1D({"xavg_gated7Begs_reject14Ngs", "xavg_gated7Begs_reject14Ngs;xavg;counts",600, -300.0, 300.0}, m_eventPtr->xavg);
			if(m_eventPtr->xavg > -186.0 && m_eventPtr->xavg < -178.0)
			{
				FillHistogram2D({"sabreE_sabreTheta_7begs_nub","sabreE_sabreTheta_7begs_nub;#theta (deg);E(MeV)",180,0.0,180.0,400,0.0,20.0},sabreCoords.Theta()*s_rad2deg,pair.ringE);
				FillHistogram2D({"sabreE_sabrePhi_7begs_nub","sabreE_sabreTheta_7begs_nub;#phi (deg);E(MeV)",360,0.0,360.0,400,0.0,20.0},Phi360(sabreCoords.Phi())*s_rad2deg,pair.ringE);
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

	void Histogrammer::RunDegradedSabre(const SabrePair& pair)
	{
		static ReconResult recon8Be, recon8BeDegrade, recon8BePunch, recon9B, recon5Li, recon7Be;
		static TVector3 sabreCoords, b9Coords, sabreNorm;
		static double relAngle, incidentAngle;

		FillHistogram1D({"sabre_counts_gated_degraderDets","sabre_counts_gated;number per event;counts",10,-1.0, 9.0}, m_eventPtr->sabre.size());

		recon5Li = m_recon.RunSabreExcitation(m_eventPtr->xavg, m_beamKE, pair, {{5,10},{2,3},{2,4},{2,4}});
		recon7Be = m_recon.RunSabreExcitation(m_eventPtr->xavg, m_beamKE, pair, {{5,10},{2,3},{2,4},{1,2}});
		recon8Be = m_recon.RunSabreExcitation(m_eventPtr->xavg, m_beamKE, pair, {{5,10},{2,3},{2,4},{1,1}});
		recon8BeDegrade = m_recon.RunSabreExcitationDegraded(m_eventPtr->xavg, m_beamKE, pair, {{5,10},{2,3},{2,4},{1,1}});
		recon8BePunch = m_recon.RunSabreExcitationPunchDegraded(m_eventPtr->xavg, m_beamKE, pair, {{5,10},{2,3},{2,4},{1,1}});
		recon9B = m_recon.RunFPResidExcitation(m_eventPtr->xavg, m_beamKE, {{5,10},{2,3},{3,4}});
		sabreCoords = m_recon.GetSabreCoordinates(pair);
		b9Coords.SetMagThetaPhi(1.0, recon9B.residThetaLab, recon9B.residPhiLab);
		sabreNorm = m_recon.GetSabreNorm(pair.detID);
		relAngle = std::acos(b9Coords.Dot(sabreCoords)/(sabreCoords.Mag()*b9Coords.Mag()));
		incidentAngle = std::acos(sabreNorm.Dot(sabreCoords)/(sabreCoords.Mag()*sabreNorm.Mag()));
		if(incidentAngle > M_PI/2.0)
			incidentAngle = M_PI - incidentAngle;

		FillHistogram1D({"incidentAngle","incidentAngle;#theta_inc;counts",180,0.0,180.0}, incidentAngle*s_rad2deg);
		FillHistogram1D({"ex_5Li_degDets", "ex_5Li;E_x(MeV);counts",3000,-5.0,25.0}, recon5Li.excitation);
		FillHistogram1D({"ex_7Be_degDets", "ex_5Li;E_x(MeV);counts",3000,-5.0,25.0}, recon7Be.excitation);
		FillHistogram1D({"ex_8be_degdDets","ex_8be_degDets; E_x(MeV); counts",3000,-5.0,25.0}, recon8Be.excitation);
		FillHistogram2D({"xavg_ex8be_degDets","xavg_ex8be_degDets;xavg;E_x(MeV)",600,-300.0,300.0,300,-5.0,25.0}, m_eventPtr->xavg, recon8Be.excitation);
		FillHistogram1D({"ex_8be_degradedPunched","ex_8be_degradedPunched; E_x(MeV); counts",300,-10.0,20.0}, recon8BePunch.excitation);
		FillHistogram1D({"ex_8be_degradedPunched"+std::to_string(pair.detID),"ex_8be_degradedPunched; E_x(MeV); counts",300,-10.0,20.0}, recon8BePunch.excitation);
		FillHistogram1D({"ex_8be_degraded","ex_8be_degraded; E_x(MeV); counts",300,-10.0,20.0}, recon8BeDegrade.excitation);
		FillHistogram1D({"ex_8be_degraded"+std::to_string(pair.detID),"ex_8be_degraded; E_x(MeV); counts",300,-10.0,20.0}, recon8BeDegrade.excitation);
		if(pair.detID == 0 || pair.detID == 4)
		{
			FillHistogram1D({"ex_8be_degradedPunched04","ex_8be_degradedPunched04; E_x(MeV); counts",300,-10.0,20.0}, recon8BePunch.excitation);
			FillHistogram1D({"ex_8be_degraded04","ex_8be_degraded04; E_x(MeV); counts",300,-10.0,20.0}, recon8BeDegrade.excitation);
		}
		FillHistogram2D({"xavg_ex8be_degradedPunched","xavg_ex8be_degradedPunched;xavg;E_x(MeV)",600,-300.0,300.0,300,-10.0,20.0}, m_eventPtr->xavg, recon8BePunch.excitation);
		FillHistogram2D({"xavg_ex8be_degraded","xavg_ex8be_degraded;xavg;E_x(MeV)",600,-300.0,300.0,300,-10.0,20.0}, m_eventPtr->xavg, recon8BeDegrade.excitation);
		FillHistogram2D({"sabrePhi_5Liex_degDets","sabrePhi_5Liex;#phi (deg);E_x (MeV)",360,0.0,360.0,1000,-5.0,25.0},Phi360(sabreCoords.Phi())*s_rad2deg,recon5Li.excitation);
		FillHistogram2D({"sabreTheta_5Liex_degDets","sabreTheta_5Liex;#theta (deg);E_x (MeV)",360,0.0,360.0,1000,-5.0,25.0},sabreCoords.Theta()*s_rad2deg,recon5Li.excitation);
		FillHistogram2D({"sabrePhi_8Beex_degDets","sabrePhi_8Beex;#phi (deg);E_x (MeV)",360,0.0,360.0,1000,-20.0,10.0},Phi360(sabreCoords.Phi())*s_rad2deg,recon8Be.excitation);
		FillHistogram2D({"sabreTheta_8Beex_degDets","sabreTheta_8Beex;#theta (deg);E_x (MeV)",360,0.0,360.0,1000,-20.0,10.0},sabreCoords.Theta()*s_rad2deg,recon8Be.excitation);
		FillHistogram2D({"sabrePhi_8Beex_degradedPunched","sabrePhi_8Beex;#phi (deg);E_x (MeV)",360,0.0,360.0,1000,-20.0,10.0},Phi360(sabreCoords.Phi())*s_rad2deg,recon8BePunch.excitation);
		FillHistogram2D({"sabreTheta_8Beex_degradedPunched","sabreTheta_8Beex;#theta (deg);E_x (MeV)",360,0.0,360.0,1000,-20.0,10.0},sabreCoords.Theta()*s_rad2deg,recon8BePunch.excitation);
		FillHistogram2D({"sabrePhi_8Beex_degraded","sabrePhi_8Beex;#phi (deg);E_x (MeV)",360,0.0,360.0,1000,-20.0,10.0},Phi360(sabreCoords.Phi())*s_rad2deg,recon8BeDegrade.excitation);
		FillHistogram2D({"sabreTheta_8Beex_degraded","sabreTheta_8Beex;#theta (deg);E_x (MeV)",360,0.0,360.0,1000,-20.0,10.0},sabreCoords.Theta()*s_rad2deg,recon8BeDegrade.excitation);
		FillHistogram2D({"relAngle_recovSabreKE_8bePunchRecon","relAngle_recovSabreKe;#theta_{rel}(deg);Recovered KE (MeV)",180,0.0,180.0,400,0.0,20.0},relAngle*s_rad2deg,recon8BePunch.sabreRxnKE);
		
		//Some KE vs. rel angle plots.
		FillHistogram2D({"relAngle_sabreKE_degDets","relAngle_sabreKE_degDets;#theta_{rel};SABRE E(Mev)",180,0.0,180.0,400,0.0,20.0},relAngle*s_rad2deg,pair.ringE);
		FillHistogram2D({"relAngle_sabreKE_degraded","relAngle_sabreKE_degraded;#theta_{rel};SABRE E(Mev)",180,0.0,180.0,400,0.0,20.0},relAngle*s_rad2deg,recon8BeDegrade.sabreRxnKE);
		FillHistogram2D({"relAngle_sabreKE_degradedPunched","relAngle_sabreKE_degradedPunched;#theta_{rel};SABRE E(Mev)",180,0.0,180.0,400,0.0,20.0},relAngle*s_rad2deg,recon8BePunch.sabreRxnKE);

		if(recon8Be.excitation > 2.2 && recon8Be.excitation < 3.8)
		{
			FillHistogram1D({"xavg_gated_8be1ex_degDets", "xavg_gated_8be1ex_degDets;xavg;counts",600,-300.0,300.0}, m_eventPtr->xavg);
		}
		if(recon7Be.excitation > -0.1 && recon7Be.excitation < 0.15)
		{
			FillHistogram1D({"xavg_gated7Begs_degDets", "xavg_gated7Begs_degDets;xavg;counts",600,-300.0,300.0}, m_eventPtr->xavg);
		}

		//Need to switch between cases, reject looking at data that has already been reconstructed correctly
		if(recon8BeDegrade.excitation > -0.5 && recon8BeDegrade.excitation < 0.5)
		{
			FillHistogram1D({"xavg_gated_8begs_degraded","xavg_gated_8begs_degraded;xavg;counts",600,-300,300}, m_eventPtr->xavg);
			FillHistogram1D({"xavg_gated_8begs_recoveredSum","xavg_gated_8begs_recoveredSum;xavg;counts",600,-300,300}, m_eventPtr->xavg);
		}
		else if(recon8BeDegrade.excitation > 2.0 && recon8BeDegrade.excitation < 4.0)
		{
			FillHistogram1D({"xavg_gated_8be1ex_degraded","xavg_gated_8be1ex_degraded;xavg;counts",600,-300,300}, m_eventPtr->xavg);
			FillHistogram1D({"xavg_gated_8be1ex_recoveredSum","xavg_gated_8be1ex_recoveredSum;xavg;counts",600,-300,300}, m_eventPtr->xavg);
		}
		else if(recon8BePunch.excitation > -1.0 && recon8BePunch.excitation < 1.0)
		{
			FillHistogram1D({"xavg_gated_8begs_degradedPunched","xavg_gated_8begs_degradedPunched;xavg;counts",600,-300,300}, m_eventPtr->xavg);
			FillHistogram1D({"xavg_gated_8begs_recoveredSum","xavg_gated_8begs_recoveredSum;xavg;counts",600,-300,300}, m_eventPtr->xavg);
		}
		else if(recon8BePunch.excitation > 1.0 && recon8BePunch.excitation < 5.0)
		{
			FillHistogram1D({"xavg_gated_8be1ex_degradedPunched","xavg_gated_8be1ex_degradedPunched;xavg;counts",600,-300,300}, m_eventPtr->xavg);
			FillHistogram1D({"xavg_gated_8be1ex_recoveredSum","xavg_gated_8be1ex_recoveredSum;xavg;counts",600,-300,300}, m_eventPtr->xavg);
			if(pair.detID == 0 || pair.detID == 4)
			{
				FillHistogram1D({"xavg_gated_8be1ex_degradedPunched_04","xavg_gated_8be1ex_degradedPunched_04;xavg;counts",600,-300,300}, m_eventPtr->xavg);
			}
			if(pair.local_wedge != 0 && pair.local_wedge != 7 && pair.local_ring != 15 && pair.local_ring != 0) //Edges might not be degraded right
			{
				FillHistogram1D({"xavg_gated_8be1ex_degradedPunched_rejectEdge","xavg_gated_8be1ex_degradedPunched_rejectEdge;xavg;counts",600,-300,300}, m_eventPtr->xavg);
				FillHistogram1D({"xavg_gated_8be1ex_recoveredSum_rejectEdge","xavg_gated_8be1ex_recoveredSum_rejectEdge;xavg;counts",600,-300,300}, m_eventPtr->xavg);
				if(pair.detID == 0 || pair.detID == 4)
					FillHistogram1D({"xavg_gated_8be1ex_degradedPunched_rejectEdge_04","xavg_gated_8be1ex_degradedPunched_rejectEdge_04;xavg;counts",600,-300,300}, m_eventPtr->xavg);
			}
		}
		if(!(recon8BeDegrade.excitation > -1.0 && recon8BeDegrade.excitation < 4.0))
		{
			FillHistogram1D({"ex_8be_degradedPunched_rejectPrev","ex_8be_degradedPunched_rejectPrev;E_x(MeV);counts",300,-10.0,20.0}, recon8BePunch.excitation);
			FillHistogram2D({"sabrePhi_8Beex_degradedPunched_rejectPrev","sabrePhi_8Beex;#phi (deg);E_x (MeV)",360,0.0,360.0,1000,-20.0,10.0},Phi360(sabreCoords.Phi())*s_rad2deg,recon8BePunch.excitation);
			FillHistogram2D({"xavg_ex8be_degradedPunched_rejectPrev","xavg_ex8be_degradedPunched;xavg;E_x(MeV)",600,-300.0,300.0,300,-10.0,20.0}, m_eventPtr->xavg, recon8BePunch.excitation);
			if(recon8BePunch.excitation > -1.0 && recon8BePunch.excitation < 5.0)
				FillHistogram1D({"xavg_gated_8beex_low_degradedPunched","xavg_gated_8beex_low_degradedPunched;xavg;counts",600,-300.0,300.0},m_eventPtr->xavg);
			if(pair.detID == 0 || pair.detID == 4)
			{
				FillHistogram1D({"ex_8be_degradedPunched04_rejectPrev","ex_8be_degradedPunched04_rejectPrev; E_x(MeV); counts",300,-10.0,20.0}, recon8BePunch.excitation);
				FillHistogram1D({"ex_8be_degraded04_rejectPrev","ex_8be_degraded04_rejectPrev; E_x(MeV); counts",300,-10.0,20.0}, recon8BeDegrade.excitation);
			}
		}
		else
		{
			FillHistogram1D({"xavg_gated_8beex_low_degraded","xavg_gated_8beex_low_degraded;xavg;counts",600,-300.0,300.0},m_eventPtr->xavg);
		}

		FillHistogram2D({"sabreE_sabreTheta_degDets", "sabreE_sabreTheta_degDets;#theta (deg);E(MeV)",180,0.0,180.0,400,0.0,20.0}, sabreCoords.Theta()*s_rad2deg, pair.ringE);
		FillHistogram2D({"xavg_sabreE_degDets", "xavg_sabreE_degDets;xavg;E(MeV)", 600,0.-300.0,300.0,400,0.0,20.0}, m_eventPtr->xavg, pair.ringE);
		FillHistogram2D({"9Btheta_sabreTheta_degDets","9Btheta_sabreTheta_degDets;#theta_{9B};#theta_{SABRE}",180,0.0,180.0,180,0.0,180.0}, recon9B.residThetaLab*s_rad2deg, sabreCoords.Theta()*s_rad2deg);
		FillHistogram2D({"sabreE_relAngle_degDets","sabreE_relAngle_degDets;#theta_{rel};E(MeV)",180,0.0,180.0,400,0.0,20.0},relAngle*s_rad2deg, pair.ringE);
		if(pair.local_wedge != 0 && pair.local_wedge != 7 && pair.local_ring != 15 && pair.local_ring != 0) //Edges might not be degraded right
		{
			FillHistogram2D({"sabreE_sabreTheta_degDets_rejectEdge", "sabreE_sabreTheta_degDets;#theta (deg);E(MeV)",180,0.0,180.0,400,0.0,20.0}, sabreCoords.Theta()*s_rad2deg, pair.ringE);
			FillHistogram2D({"xavg_sabreE_degDets_rejectEdge", "xavg_sabreE_degDets;xavg;E(MeV)",600,0.-300.0,300.0,400,0.0,20.0}, m_eventPtr->xavg, pair.ringE);
			FillHistogram1D({"xavg_degDets_rejectEdge","xavg_degDets_rejectEdge;xavg",600,-300.0,300.0}, m_eventPtr->xavg);
			FillHistogram2D({"9Btheta_sabreTheta_degDets_rejectEdge","9Btheta_sabreTheta_degDets_rejectEdge;#theta_{9B};#theta_{SABRE}",180,0.0,180.0,180,0.0,180.0}, recon9B.residThetaLab*s_rad2deg,
							  sabreCoords.Theta()*s_rad2deg);
			FillHistogram2D({"sabreE_relAngle_degDets_rejectEdge","sabreE_relAngle_degDets_rejectEdge;#theta_{rel};E(MeV)",180,0.0,180.0,400,0.0,20.0},relAngle*s_rad2deg, pair.ringE);
			if(!(recon8BeDegrade.excitation > -1.0 && recon8BeDegrade.excitation < 5.0))
			{
				FillHistogram1D({"ex_8be_degradedPunched_rejectPrev_rejectEdge","ex_8be_degradedPunched_rejectPrev;E_x(MeV);counts",300,-10.0,20.0}, recon8BePunch.excitation);
				FillHistogram2D({"sabrePhi_8Beex_degradedPunched_rejectPrev_rejectEdge","sabrePhi_8Beex;#phi (deg);E_x (MeV)",360,0.0,360.0,1000,-20.0,10.0},Phi360(sabreCoords.Phi())*s_rad2deg,recon8BePunch.excitation);
				FillHistogram2D({"xavg_ex8be_degradedPunched_rejectPrev_rejectEdge","xavg_ex8be_degradedPunched;xavg;E_x(MeV)",600,-300.0,300.0,300,-10.0,20.0}, m_eventPtr->xavg, recon8BePunch.excitation);
				if(recon8BePunch.excitation > -1.0 && recon8BePunch.excitation < 5.0)
					FillHistogram1D({"xavg_gated_8beex_low_degradedPunched_rejectEdge","xavg_gated_8beex_low_degradedPunched;xavg;counts",600,-300.0,300.0},m_eventPtr->xavg);
			}
			else
			{
				FillHistogram1D({"xavg_gated_8beex_low_degraded_rejectEdge","xavg_gated_8beex_low_degraded;xavg;counts",600,-300.0,300.0},m_eventPtr->xavg);
			}
		}
	}
}
