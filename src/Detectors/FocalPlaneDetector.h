#ifndef FOCAL_PLANE_DETECTOR_H
#define FOCAL_PLANE_DETECTOR_H

#include <vector>
#include <cmath>

namespace SabreRecon {

	class FocalPlaneDetector
	{
	public:
		struct Parameters
		{
			double B;
			double angle;
			std::vector<double> calParams;
		};

		FocalPlaneDetector();
		FocalPlaneDetector(const Parameters& params);
		~FocalPlaneDetector();

		void Init(const Parameters& params) { m_params = params; }
		double GetRho(double xavg);
		double GetP(double xavg, int Z);
		inline double GetFPTheta() const { return m_params.angle*s_deg2rad; }

	private:
		Parameters m_params;

		//Constants
		static constexpr double s_lightspeed = 299792458.0; //1/s
    	static constexpr double s_qbrho2p = 1.0e-9 * s_lightspeed; //MeV/(cm*kG)
    	static constexpr double s_deg2rad = M_PI/180.0; //rad/deg
	};
}

#endif