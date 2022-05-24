#include "FocalPlaneDetector.h"
#include <cstdlib>

namespace SabreRecon {

	FocalPlaneDetector::FocalPlaneDetector() {}

	FocalPlaneDetector::FocalPlaneDetector(const Parameters& params) :
		m_params(params)
	{
	}

	FocalPlaneDetector::~FocalPlaneDetector() {}

	double FocalPlaneDetector::GetRho(double xavg)
	{
		double rho = 0.0;
		for(size_t i=0; i< m_params.calParams.size(); i++)
		{
			rho += m_params.calParams[i]*std::pow(xavg, i);
		}
		return rho;
	}

	double FocalPlaneDetector::GetP(double xavg, int Z)
	{
		double rho = GetRho(xavg);
		return Z*rho*m_params.B*s_qbrho2p;
	}
}