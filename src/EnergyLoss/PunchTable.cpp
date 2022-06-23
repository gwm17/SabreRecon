#include "PunchTable.h"
#include <fstream>
#include <iostream>

namespace PunchTable {

	PunchTable::PunchTable() :
		m_validFlag(false)
	{
	}

	PunchTable::PunchTable(const std::string& filename) :
		m_validFlag(false)
	{
		ReadFile(filename);
	}

	PunchTable::~PunchTable() {}

	void PunchTable::ReadFile(const std::string& filename)
	{
		std::ifstream input(filename);
		if(!input.is_open())
		{
			std::cerr<<"Unable to open table file named "<<filename<<"! Exiting."<<std::endl;
			m_validFlag = false;
			return;
		}

		std::string junk;
		double thickness;
		double theta;
		double value;
		std::vector<double> energyIn, energyDep;

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
				energyIn.clear();
				energyDep.clear();
				input>>theta;
				while(input>>junk)
				{
					if(junk == "end_theta")
						break;
					energyDep.push_back(std::stod(junk));
					input>>value;
					energyIn.push_back(value);
				}

				if(!energyDep.empty())
					m_splines.emplace_back(energyDep, energyIn);
				else
					m_splines.emplace_back();
			}
			else
			{
				std::cerr<<"Unexpected expression found when reading punch table: "<<junk<<std::endl;
				m_validFlag = false;
				return;
			}
		}

		m_validFlag = true;

	}

	double PunchTable::GetInitialKineticEnergy(double theta_incident, double e_deposited)
	{
		theta_incident /= s_deg2rad;
		if(!m_validFlag)
		{
			std::cerr<<"PunchTable not initialized at GetInitialKineticEnergy()"<<std::endl;
			return 0.0;
		}
		else if(theta_incident < m_thetaMin || theta_incident > m_thetaMax)
		{
			//std::cerr<<"Theta incident outside of range of calculated values for PunchTable::GetInitialKineticEnergy"<<std::endl;
			return 0.0;
		}

		int theta_bin = (theta_incident - m_thetaMin)/m_thetaStep;

		//std::cout<<"theta bin: "<<theta_bin<<" theta_inc: "<<theta_incident<<std::endl;

		if(m_splines[theta_bin].IsValid())
		{
			double initialE =  m_splines[theta_bin].Evaluate(e_deposited);
			if(initialE == 0.0) //Not in the spline, stopped completely
			{
				//std::cout<<"here"<<std::endl;
				return e_deposited;
			}
			else
				return initialE;
		}
		else
			return e_deposited;
	}

}