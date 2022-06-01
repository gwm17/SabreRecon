#include "Reconstructor.h"
#include "MassLookup.h"
#include <iostream>

namespace SabreRecon {

	constexpr double Reconstructor::s_phiDet[5]; //C++11 weirdness with static constexpr

	Reconstructor::Reconstructor()
	{
	}

	Reconstructor::Reconstructor(const Target& target, double spsTheta, double spsB, const std::vector<double>& spsCal)
	{
		Init(target, spsTheta, spsB, spsCal);
	}

	Reconstructor::~Reconstructor() {}

	void Reconstructor::Init(const Target& target, double spsTheta, double spsB, const std::vector<double>& spsCal)
	{
		m_target = target;
		m_focalPlane.Init({spsB, spsTheta, spsCal});
		for(int i=0; i<5; i++)
			m_sabreArray.emplace_back(SabreDetector::Parameters(s_phiDet[i], s_tiltAngle, s_zOffset, false, i));
	}

	TLorentzVector Reconstructor::GetSabre4Vector(const SabrePair& pair, double mass)
	{
		TVector3 coords;
		TLorentzVector result;
		double p, E;
		if(pair.detID == 4)
			coords = m_sabreArray[4].GetHitCoordinates(15-pair.local_ring, pair.local_wedge);
		else
			coords = m_sabreArray[pair.detID].GetHitCoordinates(pair.local_ring, pair.local_wedge);
		p = std::sqrt(pair.ringE*(pair.ringE + 2.0*mass));
		E = pair.ringE + mass;
		result.SetPxPyPzE(p*std::sin(coords.Theta())*std::cos(coords.Phi()),
						  p*std::sin(coords.Theta())*std::sin(coords.Phi()),
						  p*std::cos(coords.Theta()),
						  E);
		return result;
	}

	TLorentzVector Reconstructor::GetSabre4VectorEloss(const SabrePair& pair, double mass, const NucID& id)
	{
		TVector3 coords;
		TLorentzVector result;
		double p, E, rxnKE;
		if(pair.detID == 4)
			coords = m_sabreArray[4].GetHitCoordinates(15-pair.local_ring, pair.local_wedge);
		else
			coords = m_sabreArray[pair.detID].GetHitCoordinates(pair.local_ring, pair.local_wedge);
		rxnKE = pair.ringE + m_target.GetReverseEnergyLossFractionalDepth(id.Z, id.A, pair.ringE, coords.Theta(), 0.5);
		p = std::sqrt(pair.ringE*(pair.ringE + 2.0*mass));
		E = pair.ringE + mass;
		result.SetPxPyPzE(p*std::sin(coords.Theta())*std::cos(coords.Phi()),
						  p*std::sin(coords.Theta())*std::sin(coords.Phi()),
						  p*std::cos(coords.Theta()),
						  E);
		return result;
	}

	TLorentzVector Reconstructor::GetFP4VectorEloss(double xavg, double mass, const NucID& id)
	{
		TLorentzVector result;
		double p = m_focalPlane.GetP(xavg, id.Z);
		double theta = m_focalPlane.GetFPTheta();
		double KE = std::sqrt(p*p + mass*mass) - mass;
		double rxnKE = KE + m_target.GetReverseEnergyLossFractionalDepth(id.Z, id.A, KE, theta, 0.5);
		double rxnP = sqrt(rxnKE*(rxnKE + 2.0*mass));
		double rxnE = rxnKE + mass;
		result.SetPxPyPzE(rxnP*std::sin(theta), 0.0, rxnP*std::cos(theta), rxnE);
		return result;
	}

	TLorentzVector Reconstructor::GetProj4VectorEloss(double beamKE, double mass, const NucID& id)
	{
		TLorentzVector result;
		double rxnKE = beamKE + m_target.GetReverseEnergyLossFractionalDepth(id.Z, id.A, beamKE, 0.0, 0.5);
		result.SetPxPyPzE(0.0,0.0,std::sqrt(rxnKE*(rxnKE+2.0*mass)),rxnKE+mass);
		return result;
	}

	ReconResult Reconstructor::RunThreeParticleExcitation(const SabrePair& p1, const SabrePair& p2, const SabrePair& p3, const std::vector<NucID>& nuclei)
	{
		ReconResult result;

		NucID parent;
		parent.Z = nuclei[0].Z + nuclei[1].Z + nuclei[2].Z;
		parent.A = nuclei[0].A + nuclei[1].A + nuclei[2].A;

		if(parent.Z > parent.A || parent.A <= 0 || parent.Z < 0)
		{
			std::cerr<<"Invalid parent nucleus at Reconstructor::RunThreeParticleExcitation with Z: "<<parent.Z<<" A: "<<parent.A<<std::endl;
			return result;
		}

		MassLookup& masses = MassLookup::GetInstance();
		double massParent = masses.FindMass(parent.Z, parent.A);
		double massP1 = masses.FindMass(nuclei[0].Z, nuclei[0].A);
		double massP2 = masses.FindMass(nuclei[1].Z, nuclei[1].A);
		double massP3 = masses.FindMass(nuclei[2].Z, nuclei[2].A);

		if(massParent == 0.0 || massP1 == 0.0 || massP2 == 0.0 || massP3 == 0.0)
		{
			std::cerr<<"Invalid nuclei at Reconstructor::RunThreeParticleExcitation by mass!"<<std::endl;
			return result;
		}

		auto p1_vec = GetSabre4VectorEloss(p1, massP1, nuclei[0]);
		auto p2_vec = GetSabre4VectorEloss(p2, massP2, nuclei[1]);
		auto p3_vec = GetSabre4VectorEloss(p3, massP3, nuclei[2]);

		auto parent_vec = p1_vec + p2_vec + p3_vec;

		result.excitation = parent_vec.M() - massParent;
		return result;
	}

	ReconResult Reconstructor::RunTwoParticleExcitation(const SabrePair& p1, const SabrePair& p2, const std::vector<NucID>& nuclei)
	{
		ReconResult result;

		NucID parent;
		parent.Z = nuclei[0].Z + nuclei[1].Z;
		parent.A = nuclei[0].A + nuclei[1].A;

		if(parent.Z > parent.A || parent.A <= 0 || parent.Z < 0)
		{
			std::cerr<<"Invalid parent nucleus at Reconstructor::RunTwoParticleExcitation with Z: "<<parent.Z<<" A: "<<parent.A<<std::endl;
			return result;
		}

		MassLookup& masses = MassLookup::GetInstance();
		double massParent = masses.FindMass(parent.Z, parent.A);
		double massP1 = masses.FindMass(nuclei[0].Z, nuclei[0].A);
		double massP2 = masses.FindMass(nuclei[1].Z, nuclei[1].A);

		if(massParent == 0.0 || massP1 == 0.0 || massP2 == 0.0)
		{
			std::cerr<<"Invalid nuclei at Reconstructor::RunTwoParticleExcitation by mass!"<<std::endl;
			return result;
		}

		auto p1_vec = GetSabre4VectorEloss(p1, massP1, nuclei[0]);
		auto p2_vec = GetSabre4VectorEloss(p2, massP2, nuclei[1]);

		auto parent_vec = p1_vec + p2_vec;

		result.excitation = parent_vec.M() - massParent;
		return result;
	}

	ReconResult Reconstructor::RunFPResidExcitation(double xavg, double beamKE, const std::vector<NucID>& nuclei)
	{
		ReconResult result;

		NucID resid;
		resid.Z = nuclei[0].Z + nuclei[1].Z - nuclei[2].Z;
		resid.A = nuclei[0].A + nuclei[1].A - nuclei[2].A;

		if(resid.Z > resid.A || resid.A <= 0 || resid.Z < 0)
		{
			std::cerr<<"Invalid reisdual nucleus at Reconstructor::RunFPResidExcitation with Z: "<<resid.Z<<" A: "<<resid.A<<std::endl;
			return result;
		}

		MassLookup& masses = MassLookup::GetInstance();
		double massTarg = masses.FindMass(nuclei[0].Z, nuclei[0].A);
		double massProj = masses.FindMass(nuclei[1].Z, nuclei[1].A);
		double massEject = masses.FindMass(nuclei[2].Z, nuclei[2].A);
		double massResid = masses.FindMass(resid.Z, resid.A);

		if(massTarg == 0.0 || massProj == 0.0 || massEject == 0.0 || massResid == 0.0)
		{
			std::cerr<<"Invalid nuclei at Reconstructor::RunFPResidExcitation by mass!"<<std::endl;
			return result;
		}

		TLorentzVector targ_vec;
		targ_vec.SetPxPyPzE(0.0, 0.0, 0.0, massTarg);
		auto proj_vec = GetProj4VectorEloss(beamKE, massProj, nuclei[1]);
		auto eject_vec = GetFP4VectorEloss(xavg, massEject, nuclei[2]);

		auto resid_vec = targ_vec + proj_vec - eject_vec;

		result.excitation = resid_vec.M() - massResid;

		auto parent_vec = targ_vec + proj_vec;
		auto boost = parent_vec.BoostVector();
		eject_vec.Boost(-1.0*boost);
		result.theta_cm = eject_vec.Theta();
		result.phi_cm = eject_vec.Phi();

		return result;
	}

	ReconResult Reconstructor::RunSabreResidExcitationDetEject(double beamKE, const SabrePair& pair, const std::vector<NucID>& nuclei)
	{
		ReconResult result;

		NucID resid;
		resid.Z = nuclei[0].Z + nuclei[1].Z - nuclei[2].Z;
		resid.A = nuclei[0].A + nuclei[1].A - nuclei[2].A;

		if(resid.Z > resid.A || resid.A <= 0 || resid.Z < 0)
		{
			std::cerr<<"Invalid reisdual nucleus at Reconstructor::RunFPResidExcitation with Z: "<<resid.Z<<" A: "<<resid.A<<std::endl;
			return result;
		}

		MassLookup& masses = MassLookup::GetInstance();
		double massTarg = masses.FindMass(nuclei[0].Z, nuclei[0].A);
		double massProj = masses.FindMass(nuclei[1].Z, nuclei[1].A);
		double massEject = masses.FindMass(nuclei[2].Z, nuclei[2].A);
		double massResid = masses.FindMass(resid.Z, resid.A);

		if(massTarg == 0.0 || massProj == 0.0 || massEject == 0.0 || massResid == 0.0)
		{
			std::cerr<<"Invalid nuclei at Reconstructor::RunFPResidExcitation by mass!"<<std::endl;
			return result;
		}

		TLorentzVector targ_vec;
		targ_vec.SetPxPyPzE(0.0, 0.0, 0.0, massTarg);
		auto proj_vec = GetProj4VectorEloss(beamKE, massProj, nuclei[1]);
		auto eject_vec = GetSabre4VectorEloss(pair, massEject, nuclei[2]);

		auto resid_vec = targ_vec + proj_vec - eject_vec;

		result.excitation = resid_vec.M() - massResid;

		auto parent_vec = targ_vec + proj_vec;
		auto boost = parent_vec.BoostVector();
		eject_vec.Boost(-1.0*boost);
		result.theta_cm = eject_vec.Theta();
		result.phi_cm = eject_vec.Phi();

		return result;
	}

	ReconResult Reconstructor::RunSabreExcitation(double xavg, double beamKE, const SabrePair& sabre, const std::vector<NucID>& nuclei)
	{
		ReconResult result;

		NucID decayFrag;
		decayFrag.Z = nuclei[0].Z + nuclei[1].Z - nuclei[2].Z - nuclei[3].Z;
		decayFrag.A = nuclei[0].A + nuclei[1].A - nuclei[2].A - nuclei[3].A;

		if(decayFrag.Z > decayFrag.A || decayFrag.A <= 0 || decayFrag.Z < 0)
		{
			std::cerr<<"Invalid reisdual nucleus at Reconstructor::RunSabreExcitation with Z: "<<decayFrag.Z<<" A: "<<decayFrag.A<<std::endl;
			return result;
		}

		MassLookup& masses = MassLookup::GetInstance();
		double massTarg = masses.FindMass(nuclei[0].Z, nuclei[0].A);
		double massProj = masses.FindMass(nuclei[1].Z, nuclei[1].A);
		double massEject = masses.FindMass(nuclei[2].Z, nuclei[2].A);
		double massDecayBreak = masses.FindMass(nuclei[3].Z, nuclei[3].A);
		double massDecayFrag = masses.FindMass(decayFrag.Z, decayFrag.A);

		if(massTarg == 0.0 || massProj == 0.0 || massEject == 0.0 || massDecayBreak == 0.0 || massDecayFrag == 0.0)
		{
			std::cerr<<"Invalid nuclei at Reconstructor::RunSabreExcitation by mass!"<<std::endl;
			return result;
		}

		TLorentzVector targ_vec;
		targ_vec.SetPxPyPzE(0.0, 0.0, 0.0, massTarg);
		auto proj_vec = GetProj4VectorEloss(beamKE, massProj, nuclei[1]);
		auto eject_vec = GetFP4VectorEloss(xavg, massEject, nuclei[2]);
		auto decayBreak_vec = GetSabre4VectorEloss(sabre, massDecayBreak, nuclei[3]);
		TLorentzVector resid_vec = targ_vec + proj_vec - eject_vec;
		TLorentzVector decayFrag_vec = resid_vec - decayBreak_vec;

		result.excitation = decayFrag_vec.M() - massDecayFrag;
		auto boost = resid_vec.BoostVector();
		decayBreak_vec.Boost(-1.0*boost);
		result.theta_cm = decayBreak_vec.Theta();
		result.phi_cm = decayBreak_vec.Phi();

		return result;
	}

	ReconResult Reconstructor::RunSabreExcitationDetEject(double xavg, double beamKE, const SabrePair& sabre, const std::vector<NucID>& nuclei)
	{
		ReconResult result;

		NucID decayFrag;
		decayFrag.Z = nuclei[0].Z + nuclei[1].Z - nuclei[2].Z - nuclei[3].Z;
		decayFrag.A = nuclei[0].A + nuclei[1].A - nuclei[2].A - nuclei[3].A;

		if(decayFrag.Z > decayFrag.A || decayFrag.A <= 0 || decayFrag.Z < 0)
		{
			std::cerr<<"Invalid reisdual nucleus at Reconstructor::RunSabreExcitation with Z: "<<decayFrag.Z<<" A: "<<decayFrag.A<<std::endl;
			return result;
		}

		MassLookup& masses = MassLookup::GetInstance();
		double massTarg = masses.FindMass(nuclei[0].Z, nuclei[0].A);
		double massProj = masses.FindMass(nuclei[1].Z, nuclei[1].A);
		double massEject = masses.FindMass(nuclei[2].Z, nuclei[2].A);
		double massDecayBreak = masses.FindMass(nuclei[3].Z, nuclei[3].A);
		double massDecayFrag = masses.FindMass(decayFrag.Z, decayFrag.A);

		if(massTarg == 0.0 || massProj == 0.0 || massEject == 0.0 || massDecayBreak == 0.0 || massDecayFrag == 0.0)
		{
			std::cerr<<"Invalid nuclei at Reconstructor::RunSabreExcitation by mass!"<<std::endl;
			return result;
		}

		TLorentzVector targ_vec;
		targ_vec.SetPxPyPzE(0.0, 0.0, 0.0, massTarg);
		auto proj_vec = GetProj4VectorEloss(beamKE, massProj, nuclei[1]);
		auto decayBreak_vec = GetFP4VectorEloss(xavg, massEject, nuclei[2]);
		auto eject_vec = GetSabre4VectorEloss(sabre, massDecayBreak, nuclei[3]);
		TLorentzVector resid_vec = targ_vec + proj_vec - eject_vec;
		TLorentzVector decayFrag_vec = resid_vec - decayBreak_vec;

		result.excitation = decayFrag_vec.M() - massDecayFrag;
		auto boost = resid_vec.BoostVector();
		decayBreak_vec.Boost(-1.0*boost);
		result.theta_cm = decayBreak_vec.Theta();
		result.phi_cm = decayBreak_vec.Phi();

		return result;
	}

	TVector3 Reconstructor::GetSabreCoordinates(const SabrePair& pair)
	{
		if(pair.detID == 4)
			return m_sabreArray[4].GetHitCoordinates(15-pair.local_ring, pair.local_wedge);
		else
			return m_sabreArray[pair.detID].GetHitCoordinates(pair.local_ring, pair.local_wedge);
	}
}