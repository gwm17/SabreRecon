#ifndef RECONSTRUCTOR_H
#define RECONSTRUCTOR_H

#include <string>
#include <vector>
#include "EnergyLoss/Target.h"
#include "EnergyLoss/ElossTable.h"
#include "EnergyLoss/PunchTable.h"
#include "CalDict/DataStructs.h"
#include "TLorentzVector.h"
#include "Detectors/SabreDetector.h"
#include "Detectors/FocalPlaneDetector.h"

namespace SabreRecon {

	struct ReconResult
	{
		double excitation = -100.0;
		double sabreRxnKE = -100.0;
		double ejectThetaCM = -100.0;
		double ejectPhiCM = -100.0;
		double residThetaLab = -100.0;
		double residPhiLab = -100.0;
		double residThetaCM = -100.0;
		double residPhiCM = -100.0;
	};

	struct NucID
	{
		int Z, A;

		NucID() :
		Z(0), A(0)
		{
		}

		NucID(int z, int a) :
		Z(z), A(a)
		{
		}
	};

	class Reconstructor
	{
	public:
		Reconstructor();
		Reconstructor(const Target& target, double spsTheta, double spsB, const std::vector<double>& spsCal);
		~Reconstructor();

		void Init(const Target& target, double spsTheta, double spsB, const std::vector<double>& spsCal);

		void AddEnergyLossTable(const std::string& filename);
		void AddPunchThruTable(const std::string& filename);

		ReconResult RunThreeParticleExcitation(const SabrePair& p1, const SabrePair& p2, const SabrePair& p3, const std::vector<NucID>& nuclei);
		ReconResult RunTwoParticleExcitation(const SabrePair& p1, const SabrePair& p2, const std::vector<NucID>& nuclei);
		//nuclei: target, projectile, ejectile
		ReconResult RunFPResidExcitation(double xavg, double beamKE, const std::vector<NucID>& nuclei);
		//nuclei: target, projectile, ejectile
    	ReconResult RunSabreResidExcitationDetEject(double beamKE, const SabrePair& sabre, const std::vector<NucID>& nuclei);
    	//nuclei: target, projectile, ejectile, decaySabre
    	ReconResult RunSabreExcitation(double xavg, double beamKE, const SabrePair& sabre,  const std::vector<NucID>& nuclei);
    	//nuclei: target, projectile, ejectile, decayFP
    	ReconResult RunSabreExcitationDetEject(double xavg, double beamKE, const SabrePair& sabre, const std::vector<NucID>& nuclei);

		ReconResult RunSabreExcitationPunch(double xavg, double beamKE, const SabrePair& sabre, const std::vector<NucID>& nuclei);
		ReconResult RunSabreExcitationPunchDegraded(double xavg, double beamKE, const SabrePair& sabre, const std::vector<NucID>& nuclei);

		TVector3 GetSabreCoordinates(const SabrePair& pair);
    	

	private:
		TLorentzVector GetSabre4Vector(const SabrePair& pair, double mass);
    	TLorentzVector GetSabre4VectorEloss(const SabrePair& pair, double mass, const NucID& id);
    	TLorentzVector GetSabre4VectorElossPunchThru(const SabrePair& pair, double mass, const NucID& id);
    	TLorentzVector GetSabre4VectorElossPunchThruDegraded(const SabrePair& pair, double mass, const NucID& id);
    	TLorentzVector GetFP4VectorEloss(double xavg, double mass, const NucID& id);
    	TLorentzVector GetProj4VectorEloss(double beamKE, double mass, const NucID& id);

		PunchTable::PunchTable* GetPunchThruTable(const NucID& projectile, const NucID& material);
		PunchTable::ElossTable* GetElossTable(const NucID& projectile, const NucID& material);

    	std::vector<SabreDetector> m_sabreArray;
    	FocalPlaneDetector m_focalPlane;
    	Target m_target;
		Target m_sabreDeadLayer;

		std::vector<PunchTable::PunchTable> m_punchTables;
		std::vector<PunchTable::ElossTable> m_elossTables;

    	//SABRE constants
    	static constexpr double s_phiDet[5] = { 306.0, 18.0, 234.0, 162.0, 90.0 };
    	static constexpr double s_tiltAngle = 40.0;
    	static constexpr double s_zOffset = -0.1245;
		static constexpr double s_sabreDeadlayerThickness = 50.0 * 1.0e-7 * 2.3296 * 1.0e6; // 50 nm deadlayer -> ug/cm^2

    	//Kinematics constants
    	static constexpr double s_deg2rad = M_PI/180.0; //rad/deg
	};
}

#endif