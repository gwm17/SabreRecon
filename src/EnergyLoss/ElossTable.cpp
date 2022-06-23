#include "ElossTable.h"
#include <iostream>
#include <iomanip>
#include <fstream>

namespace PunchTable {

    ElossTable::ElossTable() :
        m_isValid(false)
    {
    }

    ElossTable::ElossTable(const std::string& filename) :
        m_isValid(false)
    {
        ReadFile(filename);
    }

    ElossTable::~ElossTable() {}

    void ElossTable::ReadFile(const std::string& filename)
    {
        std::ifstream input(filename);
		if(!input.is_open())
		{
			std::cerr<<"Unable to open table file named "<<filename<<"! Exiting."<<std::endl;
			m_isValid = false;
			return;
		}

		std::string junk;
		double thickness;
		double theta;
		double value;
		std::vector<double> energyLoss, energyFinal;

		input>>junk>>junk>>m_projectileString;
		input>>junk>>junk;

		while(input>>junk)
		{
			if(junk == "---------------------------------")
				break;
			input>>junk;
			m_materialString += junk;
			input>>junk;
		}
		input>>junk>>m_thetaMin>>junk>>m_thetaMax>>junk>>m_thetaStep;
		std::getline(input, junk);
		std::getline(input, junk);
		std::getline(input, junk);
		std::getline(input, junk);
		while(input>>junk)
		{
			if(junk == "begin_theta")
			{
				energyLoss.clear();
				energyFinal.clear();
				input>>theta;
				while(input>>junk)
				{
					if(junk == "end_theta")
						break;
					energyFinal.push_back(std::stod(junk));
					input>>value;
					energyLoss.push_back(value);
				}

				if(!energyFinal.empty())
					m_splines.emplace_back(energyFinal, energyLoss);
				else
					m_splines.emplace_back();
			}
			else
			{
				std::cerr<<"Unexpected expression found when reading punch table: "<<junk<<std::endl;
				m_isValid = false;
				return;
			}
		}

		m_isValid = true;
    }

    double ElossTable::GetEnergyLoss(double thetaIncident, double finalEnergy)
    {
        thetaIncident /= s_deg2rad;
		if(!m_isValid)
		{
			std::cerr<<"ElossTable not initialized at GetEnergyLoss()"<<std::endl;
			return 0.0;
		}
		else if(thetaIncident < m_thetaMin || thetaIncident > m_thetaMax)
		{
			//std::cerr<<"Theta incident "<<thetaIncident<<" outside of range of calculated values for ElossTable::GetEnergyLoss()"<<std::endl;
			return 0.0;
		}

		int theta_bin = (thetaIncident - m_thetaMin)/m_thetaStep;

		//std::cout<<"theta bin: "<<theta_bin<<" theta_inc: "<<thetaIncident<<std::endl;

		if(m_splines[theta_bin].IsValid())
		{
			double eloss =  m_splines[theta_bin].Evaluate(finalEnergy);
			if(eloss == 0.0) //Not in the spline
			{
				//std::cout<<"here"<<std::endl;
				return 0.0;
			}
			else
				return eloss;
		}
		else
        {
            std::cerr<<"Spline error at ElossTable::GetEnergyLoss()!"<<std::endl;
			return 0.0;
        }
    }
}