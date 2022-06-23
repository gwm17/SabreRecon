#ifndef ELOSS_TABLE_H
#define ELOSS_TABLE_H

#include "CubicSpline.h"
#include <vector>
#include <string>
#include <cmath>

namespace PunchTable {

    /*
        Class for interpolating a reverse energy loss file. Note that
        this is for REVERSE energy loss, which is typically what is interesting
        for experimental data.
    */
    class ElossTable
    {
    public:
        ElossTable();
        ElossTable(const std::string& filename);
        ~ElossTable();

        void ReadFile(const std::string& filename);
        std::string GetProjectile() { return m_projectileString; }
        std::string GetMaterial() { return m_materialString; }

        double GetEnergyLoss(double thetaIncident, double finalEnergy);
        
        inline const bool IsValid() const { return m_isValid; }

    private:
        std::string m_projectileString;
        std::string m_materialString;
        std::vector<CubicSpline> m_splines;
        double m_thetaStep, m_thetaMin, m_thetaMax;

        bool m_isValid;

        static constexpr double s_deg2rad = M_PI/180.0;
    };
}

#endif