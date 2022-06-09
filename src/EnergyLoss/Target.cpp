/*

Target.cpp
A basic target unit for use in the SPANCRedux environment. A target
is defined as a single compound with elements Z,A of a given stoichiometry 
Holds an energy loss class

Based on code by D.W. Visser written at Yale for the original SPANC

Written by G.W. McCann Aug. 2020

*/
#include "Target.h"
#include "EnergyLossConstants.h"
#include "MassLookup.h"
#include <cmath>

namespace SabreRecon {

	Target::Target() :
		m_totalThickness(0.0), m_isValid(false)
	{
	}

	/*Targets must be of known thickness*/
	Target::Target(const std::vector<int>& a, const std::vector<int>& z, const std::vector<int>& stoich, double thick) :
		m_isValid(false)
	{
		SetParameters(a, z, stoich, thick);
	}
	
	Target::~Target()
	{
	}
	
	/*Set target elements of given Z, A, S*/
	void Target::SetParameters(const std::vector<int>& a, const std::vector<int>& z, const std::vector<int>& stoich, double thick)
	{
		m_params.ZT = z;
		double denom = 0;
		for(auto& s : stoich)
			denom += s;
		for(auto& s : stoich)
			m_params.composition.push_back(s/denom);
		m_totalThickness = thick;
		m_totalThickness_gcm2 = m_totalThickness*1.0e-6;

		auto& masses = MassLookup::GetInstance();
		for(size_t i=0; i<z.size(); i++)
		{
			m_material.add_element(masses.FindMassU(z[i], a[i]), z[i], stoich[i]);
		}
		m_isValid = true;
	}
	
	/*Calculates energy loss for travelling all the way through the target*/
	double Target::GetEnergyLossTotal(int zp, int ap, double startEnergy, double theta)
	{
		m_params.ZP = zp;
		m_params.massP = MassLookup::GetInstance().FindMass(zp, ap)*EnergyLoss::mev2u;

		if(theta == M_PI/2.) 
			return startEnergy;
		else if (theta > M_PI/2.) 
			theta = M_PI - theta;

		m_params.energy = startEnergy;
		m_params.thickness = m_totalThickness/(std::fabs(std::cos(theta)));

		m_projectile.A = MassLookup::GetInstance().FindMassU(zp, ap);
		m_projectile.Z = zp;
		m_projectile.Q = zp;
		m_projectile.T = startEnergy/m_projectile.A;
		m_material.thickness(m_totalThickness_gcm2/(std::fabs(std::cos(theta))));

		return catima::integrate_energyloss(m_projectile, m_material);
		//return EnergyLoss::GetEnergyLoss(m_params);
	}

	/*Calculates the energy loss for traveling some fraction through the target*/
	double Target::GetEnergyLossFractionalDepth(int zp, int ap, double startEnergy, double theta, double percent_depth)
	{
		m_params.ZP = zp;
		m_params.massP = MassLookup::GetInstance().FindMass(zp, ap)*EnergyLoss::mev2u;

		if(theta == M_PI/2.)
			return startEnergy;
		else if (theta > M_PI/2.)
			theta = M_PI-theta;

		m_params.energy = startEnergy;
		m_params.thickness = m_totalThickness*percent_depth/(std::fabs(std::cos(theta)));

		m_projectile.A = MassLookup::GetInstance().FindMassU(zp, ap);
		m_projectile.Z = zp;
		m_projectile.Q = zp;
		m_projectile.T = startEnergy/m_projectile.A;
		m_material.thickness(m_totalThickness_gcm2*percent_depth/(std::fabs(std::cos(theta))));

		return catima::integrate_energyloss(m_projectile, m_material);

		//return EnergyLoss::GetEnergyLoss(m_params);
	}
	
	/*Calculates reverse energy loss for travelling all the way through the target*/
	double Target::GetReverseEnergyLossTotal(int zp, int ap, double finalEnergy, double theta)
	{
		m_params.ZP = zp;
		m_params.massP = MassLookup::GetInstance().FindMass(zp, ap)*EnergyLoss::mev2u;

		if(theta == M_PI/2.) 
			return finalEnergy;
		else if (theta > M_PI/2.) 
			theta = M_PI - theta;

		m_params.energy = finalEnergy;
		m_params.thickness = m_totalThickness/(std::fabs(std::cos(theta)));

		m_projectile.A = MassLookup::GetInstance().FindMassU(zp, ap);
		m_projectile.Z = zp;
		m_projectile.Q = zp;
		m_projectile.T = finalEnergy/m_projectile.A;
		m_material.thickness(m_totalThickness_gcm2/(std::fabs(std::cos(theta))));

		return catima::reverse_integrate_energyloss(m_projectile, m_material);

		//return EnergyLoss::GetReverseEnergyLoss(m_params);
	}

	/*Calculates the reverse energy loss for traveling some fraction through the target*/
	double Target::GetReverseEnergyLossFractionalDepth(int zp, int ap, double finalEnergy, double theta, double percent_depth)
	{
		m_params.ZP = zp;
		m_params.massP = MassLookup::GetInstance().FindMass(zp, ap)*EnergyLoss::mev2u;

		if(theta == M_PI/2.)
			return finalEnergy;
		else if (theta > M_PI/2.)
			theta = M_PI-theta;

		m_params.energy = finalEnergy;
		m_params.thickness = m_totalThickness*percent_depth/(std::fabs(std::cos(theta)));

		m_projectile.A = MassLookup::GetInstance().FindMassU(zp, ap);
		m_projectile.Z = zp;
		m_projectile.Q = zp;
		m_projectile.T = finalEnergy/m_projectile.A;
		m_material.thickness(m_totalThickness_gcm2*percent_depth/(std::fabs(std::cos(theta))));

		return catima::reverse_integrate_energyloss(m_projectile, m_material);

		//return EnergyLoss::GetReverseEnergyLoss(m_params);
	}

}
