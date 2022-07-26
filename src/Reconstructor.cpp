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

		//Setup intermediate energy loss layers
		m_sabreDeadLayer.SetParameters({28}, {14}, {1}, s_sabreDeadlayerThickness);
	}

	void Reconstructor::AddEnergyLossTable(const std::string& filename)
	{
		m_elossTables.emplace_back(filename);
	}

	void Reconstructor::AddPunchThruTable(const std::string& filename)
	{
		m_punchTables.emplace_back(filename);
	}

	PunchTable::ElossTable* Reconstructor::GetElossTable(const NucID& projectile, const NucID& material)
	{
		MassLookup& masses = MassLookup::GetInstance();
		std::string projString = masses.FindSymbol(projectile.Z, projectile.A);
		std::string matString = masses.FindSymbol(material.Z, material.A) + "1"; //temp

		for(auto& table : m_elossTables)
		{
			if(table.GetProjectile() == projString && table.GetMaterial() == matString)
				return &table;
		}

		return nullptr;
	}

	PunchTable::PunchTable* Reconstructor::GetPunchThruTable(const NucID& projectile, const NucID& material)
	{
		MassLookup& masses = MassLookup::GetInstance();
		std::string projString = masses.FindSymbol(projectile.Z, projectile.A);
		std::string matString = masses.FindSymbol(material.Z, material.A) + "1"; //temp

		for(auto& table : m_punchTables)
		{
			if(table.GetProjectile() == projString && table.GetMaterial() == matString)
				return &table;
		}

		return nullptr;
	}

	TLorentzVector Reconstructor::GetSabre4Vector(const SabrePair& pair, double mass)
	{
		TVector3 coords;
		TLorentzVector result;
		double p, E, theta, phi;
		if(pair.detID == 4)
			coords = m_sabreArray[4].GetHitCoordinates(15-pair.local_ring, pair.local_wedge);
		else
			coords = m_sabreArray[pair.detID].GetHitCoordinates(pair.local_ring, pair.local_wedge);
		p = std::sqrt(pair.ringE*(pair.ringE + 2.0*mass));
		E = pair.ringE + mass;
		theta = coords.Theta();
		phi = coords.Phi();
		result.SetPxPyPzE(p*std::sin(theta)*std::cos(phi),
						  p*std::sin(theta)*std::sin(phi),
						  p*std::cos(theta),
						  E);
		return result;
	}

	TLorentzVector Reconstructor::GetSabre4VectorEloss(const SabrePair& pair, double mass, const NucID& id)
	{
		TVector3 coords, sabreNorm;
		TLorentzVector result;
		double incidentAngle, p, E, rxnKE, theta, phi;

		if(pair.detID == 4)
			coords = m_sabreArray[4].GetHitCoordinates(15-pair.local_ring, pair.local_wedge);
		else
			coords = m_sabreArray[pair.detID].GetHitCoordinates(pair.local_ring, pair.local_wedge);
		sabreNorm = m_sabreArray[pair.detID].GetNormTilted();
		incidentAngle = std::acos(sabreNorm.Dot(coords)/(sabreNorm.Mag()*coords.Mag()));

		rxnKE = pair.ringE + m_sabreDeadLayer.GetReverseEnergyLossTotal(id.Z, id.A, pair.ringE, incidentAngle);
		rxnKE += m_target.GetReverseEnergyLossFractionalDepth(id.Z, id.A, rxnKE, coords.Theta(), 0.5);
		p = std::sqrt(rxnKE*(rxnKE + 2.0*mass));
		E = rxnKE + mass;
		theta = coords.Theta();
		phi = coords.Phi();
		result.SetPxPyPzE(p*std::sin(theta)*std::cos(phi),
						  p*std::sin(theta)*std::sin(phi),
						  p*std::cos(theta),
						  E);
		return result;
	}

	TLorentzVector Reconstructor::GetSabre4VectorElossPunchThru(const SabrePair& pair, double mass, const NucID& id)
	{
		TVector3 coords, sabreNorm;
		TLorentzVector result;
		double incidentAngle, p, E, rxnKE, theta, phi;

		PunchTable::PunchTable* table = GetPunchThruTable(id, {14, 28});
		if(table == nullptr)
			return result;

		if(pair.detID == 4)
			coords = m_sabreArray[4].GetHitCoordinates(15-pair.local_ring, pair.local_wedge);
		else
			coords = m_sabreArray[pair.detID].GetHitCoordinates(pair.local_ring, pair.local_wedge);
		sabreNorm = m_sabreArray[pair.detID].GetNormTilted();
		incidentAngle = std::acos(sabreNorm.Dot(coords)/(sabreNorm.Mag()*coords.Mag()));
		if(incidentAngle > M_PI/2.0)
			incidentAngle = M_PI - incidentAngle;
		rxnKE = table->GetInitialKineticEnergy(incidentAngle, pair.ringE);
		if(rxnKE == 0.0)
			return result;
		rxnKE += m_target.GetReverseEnergyLossFractionalDepth(id.Z, id.A, rxnKE, coords.Theta(), 0.5);
		p = std::sqrt(rxnKE*(rxnKE + 2.0*mass));
		E = rxnKE + mass;
		theta = coords.Theta();
		phi = coords.Phi();
		result.SetPxPyPzE(p*std::sin(theta)*std::cos(phi), p*std::sin(theta)*std::sin(phi), p*std::cos(theta), E);
		return result;
	}

	TLorentzVector Reconstructor::GetSabre4VectorElossPunchThruDegraded(const SabrePair& pair, double mass, const NucID& id)
	{
		TVector3 coords, sabreNorm;
		TLorentzVector result;
		double incidentAngle, p, E, rxnKE, theta, phi;

		PunchTable::PunchTable* ptable = GetPunchThruTable(id, {14, 28});
		PunchTable::ElossTable* etable = GetElossTable(id, {73, 181});
		if(ptable == nullptr || etable == nullptr)
			return result;

		if(pair.detID == 4)
			coords = m_sabreArray[4].GetHitCoordinates(15-pair.local_ring, pair.local_wedge);
		else
			coords = m_sabreArray[pair.detID].GetHitCoordinates(pair.local_ring, pair.local_wedge);
		sabreNorm = m_sabreArray[pair.detID].GetNormTilted();
		incidentAngle = std::acos(sabreNorm.Dot(coords)/(sabreNorm.Mag()*coords.Mag()));
		if(incidentAngle > M_PI/2.0)
			incidentAngle = M_PI - incidentAngle;

		rxnKE = ptable->GetInitialKineticEnergy(incidentAngle, pair.ringE);
		if(rxnKE == pair.ringE)
			return result;
		rxnKE += etable->GetEnergyLoss(incidentAngle, rxnKE);
		if(rxnKE == 0.0)
			return result;
		rxnKE += m_target.GetReverseEnergyLossFractionalDepth(id.Z, id.A, rxnKE, coords.Theta(), 0.5);
		p = std::sqrt(rxnKE*(rxnKE + 2.0*mass));
		E = rxnKE + mass;
		theta = coords.Theta();
		phi = coords.Phi();
		result.SetPxPyPzE(p*std::sin(theta)*std::cos(phi), p*std::sin(theta)*std::sin(phi), p*std::cos(theta), E);
		return result;
	}

	TLorentzVector Reconstructor::GetSabre4VectorElossDegraded(const SabrePair& pair, double mass, const NucID& id)
	{
		TVector3 coords, sabreNorm;
		TLorentzVector result;
		double incidentAngle, p, E, rxnKE, theta, phi;

		PunchTable::ElossTable* etable = GetElossTable(id, {73, 181});
		if(etable == nullptr)
			return result;

		if(pair.detID == 4)
			coords = m_sabreArray[4].GetHitCoordinates(15-pair.local_ring, pair.local_wedge);
		else
			coords = m_sabreArray[pair.detID].GetHitCoordinates(pair.local_ring, pair.local_wedge);
		sabreNorm = m_sabreArray[pair.detID].GetNormTilted();
		incidentAngle = std::acos(sabreNorm.Dot(coords)/(sabreNorm.Mag()*coords.Mag()));
		if(incidentAngle > M_PI/2.0)
			incidentAngle = M_PI - incidentAngle;

		rxnKE = pair.ringE + etable->GetEnergyLoss(incidentAngle, pair.ringE);
		if(rxnKE == 0.0)
			return result;
		rxnKE += m_target.GetReverseEnergyLossFractionalDepth(id.Z, id.A, rxnKE, coords.Theta(), 0.5);
		p = std::sqrt(rxnKE*(rxnKE + 2.0*mass));
		E = rxnKE + mass;
		theta = coords.Theta();
		phi = coords.Phi();
		result.SetPxPyPzE(p*std::sin(theta)*std::cos(phi), p*std::sin(theta)*std::sin(phi), p*std::cos(theta), E);
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
		result.residThetaLab = resid_vec.Theta();
		result.residPhiLab = resid_vec.Phi();

		auto parent_vec = targ_vec + proj_vec;
		auto boost = parent_vec.BoostVector();
		eject_vec.Boost(-1.0*boost);

		result.ejectThetaCM = eject_vec.Theta();
		result.ejectPhiCM = eject_vec.Phi();

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
		result.sabreRxnKE = eject_vec.E() - massEject;
		auto boost = parent_vec.BoostVector();
		eject_vec.Boost(-1.0*boost);
		result.ejectThetaCM = eject_vec.Theta();
		result.ejectPhiCM = eject_vec.Phi();

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
		result.sabreRxnKE = decayBreak_vec.E() - massDecayBreak;
		auto boost = resid_vec.BoostVector();
		decayBreak_vec.Boost(-1.0*boost);
		result.ejectThetaCM = decayBreak_vec.Theta();
		result.ejectPhiCM = decayBreak_vec.Phi();

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
		result.sabreRxnKE = eject_vec.E() - massEject;
		auto boost = resid_vec.BoostVector();
		decayBreak_vec.Boost(-1.0*boost);
		result.ejectThetaCM = decayBreak_vec.Theta();
		result.ejectPhiCM = decayBreak_vec.Phi();

		return result;
	}

	ReconResult Reconstructor::RunSabreExcitationPunch(double xavg, double beamKE, const SabrePair& sabre, const std::vector<NucID>& nuclei)
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
		auto decayBreak_vec = GetSabre4VectorElossPunchThru(sabre, massDecayBreak, nuclei[3]);
		if(decayBreak_vec.E() == 0.0)
		{
			return result;
		}
		TLorentzVector resid_vec = targ_vec + proj_vec - eject_vec;
		TLorentzVector decayFrag_vec = resid_vec - decayBreak_vec;

		result.excitation = decayFrag_vec.M() - massDecayFrag;
		result.sabreRxnKE = decayBreak_vec.E() - massDecayBreak;
		auto boost = resid_vec.BoostVector();
		decayBreak_vec.Boost(-1.0*boost);
		result.ejectThetaCM = decayBreak_vec.Theta();
		result.ejectPhiCM = decayBreak_vec.Phi();

		return result;
	}

	ReconResult Reconstructor::RunSabreExcitationPunchDegraded(double xavg, double beamKE, const SabrePair& sabre, const std::vector<NucID>& nuclei)
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
		auto decayBreak_vec = GetSabre4VectorElossPunchThruDegraded(sabre, massDecayBreak, nuclei[3]);
		if(decayBreak_vec.E() == 0.0)
		{
			return result;
		}
		TLorentzVector resid_vec = targ_vec + proj_vec - eject_vec;
		TLorentzVector decayFrag_vec = resid_vec - decayBreak_vec;

		result.excitation = decayFrag_vec.M() - massDecayFrag;
		result.sabreRxnKE = decayBreak_vec.E() - massDecayBreak;
		auto boost = resid_vec.BoostVector();
		decayBreak_vec.Boost(-1.0*boost);
		result.ejectThetaCM = decayBreak_vec.Theta();
		result.ejectPhiCM = decayBreak_vec.Phi();

		return result;
	}

	ReconResult Reconstructor::RunSabreExcitationDegraded(double xavg, double beamKE, const SabrePair& sabre, const std::vector<NucID>& nuclei)
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
		auto decayBreak_vec = GetSabre4VectorElossDegraded(sabre, massDecayBreak, nuclei[3]);
		if(decayBreak_vec.E() == 0.0)
		{
			return result;
		}
		TLorentzVector resid_vec = targ_vec + proj_vec - eject_vec;
		TLorentzVector decayFrag_vec = resid_vec - decayBreak_vec;

		result.excitation = decayFrag_vec.M() - massDecayFrag;
		result.sabreRxnKE = decayBreak_vec.E() - massDecayBreak;
		auto boost = resid_vec.BoostVector();
		decayBreak_vec.Boost(-1.0*boost);
		result.ejectThetaCM = decayBreak_vec.Theta();
		result.ejectPhiCM = decayBreak_vec.Phi();

		return result;
	}

	TVector3 Reconstructor::GetSabreCoordinates(const SabrePair& pair)
	{
		if(pair.detID == 4)
			return m_sabreArray[4].GetHitCoordinates(15-pair.local_ring, pair.local_wedge);
		else
			return m_sabreArray[pair.detID].GetHitCoordinates(pair.local_ring, pair.local_wedge);
	}

	TVector3 Reconstructor::GetSabreNorm(int detID)
	{
		if(detID >= m_sabreArray.size() || detID < 0)
			return TVector3();
		return m_sabreArray[detID].GetNormTilted();
	}
}