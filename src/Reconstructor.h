#ifndef RECONSTRUCTOR_H
#define RECONSTRUCTOR_H

#include <string>
#include <vector>
#include "EnergyLoss/Target.h"
#include "CalDict/DataStructs.h"
#include "TLorentzVector.h"
#include "Detectors/SabreDetector.h"
#include "Detectors/FocalPlaneDetector.h"

namespace SabreRecon {

	struct ReconResult
	{
		double excitation;
		double ejectThetaCM;
		double ejectPhiCM;
		double residThetaLab;
		double residPhiLab;
		double residThetaCM;
		double residPhiCM;
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

		TVector3 GetSabreCoordinates(const SabrePair& pair);
    	

	private:
		TLorentzVector GetSabre4Vector(const SabrePair& pair, double mass);
    	TLorentzVector GetSabre4VectorEloss(const SabrePair& pair, double mass, const NucID& id);
    	TLorentzVector GetFP4VectorEloss(double xavg, double mass, const NucID& id);
    	TLorentzVector GetProj4VectorEloss(double beamKE, double mass, const NucID& id);

    	std::vector<SabreDetector> m_sabreArray;
    	FocalPlaneDetector m_focalPlane;
    	Target m_target;

    	//SABRE constants
    	static constexpr double s_phiDet[5] = { 306.0, 18.0, 234.0, 162.0, 90.0 };
    	static constexpr double s_tiltAngle = 40.0;
    	static constexpr double s_zOffset = -0.1245;

    	//Kinematics constants
    	static constexpr double s_deg2rad = M_PI/180.0; //rad/deg
	};
}

#endif