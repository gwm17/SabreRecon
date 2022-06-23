#ifndef PUNCHTABLE_H
#define PUNCHTABLE_H

#include "CubicSpline.h"
#include <vector>
#include <string>
#include <cmath>

namespace PunchTable {

	class PunchTable {
	public:
		PunchTable();
		PunchTable(const std::string& filename);
		~PunchTable();

		void ReadFile(const std::string& filename);
		std::string GetProjectile() { return m_projectileString; }
        std::string GetMaterial() { return m_materialString; }

		double GetInitialKineticEnergy(double theta_incident, double e_deposited); //radians, MeV
		inline bool IsValid() const { return m_validFlag; }

	private:
		std::string m_projectileString;
        std::string m_materialString;
		std::vector<CubicSpline> m_splines;
		double m_thetaStep, m_thetaMin, m_thetaMax;

		bool m_validFlag;

		static constexpr double s_deg2rad = M_PI/180.0;
	};

}
#endif	