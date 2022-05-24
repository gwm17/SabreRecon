#include "CutHandler.h"
#include <iostream>

namespace SabreRecon {

	CutHandler::CutHandler() :
		m_eventPtr(nullptr), m_isValid(false), m_eventInit(false), m_cutInit(false)
	{
	}

	CutHandler::CutHandler(const std::vector<ReconCut>& cuts, CalEvent* event) :
		m_eventPtr(nullptr), m_isValid(false), m_eventInit(false), m_cutInit(false)
	{
		InitCuts(cuts);
		InitEvent(event);
	}

	CutHandler::~CutHandler()
	{
		for(auto& cut : m_cuts)
		{
			if(cut.file_ptr->IsOpen())
				cut.file_ptr->Close();
		}
	}

	void CutHandler::InitCuts(const std::vector<ReconCut>& cuts)
	{
		m_cuts = cuts;

		for(auto& cut : m_cuts)
		{
			cut.file_ptr = TFile::Open(cut.filename.c_str(), "READ");
			if(!cut.file_ptr->IsOpen())
			{
				std::cerr<<"Unable to open cutfile "<<cut.filename<<std::endl;
				m_cutInit = false;
				m_isValid = false;
				return;
			}

			cut.cut_ptr = (TCutG*) cut.file_ptr->Get("CUTG");
			if(cut.cut_ptr == nullptr)
			{
				std::cerr<<"Unable to get CUTG from cutfile "<<cut.filename<<std::endl;
				m_cutInit = false;
				m_isValid = false;
				return;
			}

			cut.cut_ptr->SetName(cut.cutname.c_str());
		}

		m_cutInit = true;
		if(m_eventInit)
			m_isValid = true;
	}

	void CutHandler::InitEvent(CalEvent* event)
	{
		m_eventPtr = event;

		if(m_eventPtr == nullptr)
		{
			std::cerr<<"Invalid event pointer given to CutHandler::InitEvent!"<<std::endl;
			m_eventInit = false;
			m_isValid = false;
			return;
		}

		m_varMap["xavg"] = &m_eventPtr->xavg;
		m_varMap["x1"] = &m_eventPtr->x1;
		m_varMap["x2"] = &m_eventPtr->x2;
		m_varMap["scintE"] = &m_eventPtr->scintE;
		m_varMap["cathodeE"] = &m_eventPtr->cathodeE;
		m_varMap["anodeFrontE"] = &m_eventPtr->anodeFrontE;
		m_varMap["anodeBackE"] = &m_eventPtr->anodeBackE;
		m_varMap["theta"] = &m_eventPtr->theta;
		m_varMap["scintT"] = &m_eventPtr->scintT;

		m_eventInit = true;
		if(m_cutInit)
			m_isValid = true;
	}

	bool CutHandler::IsInside()
	{
		for(auto& cut : m_cuts)
		{
			auto iterX = m_varMap.find(cut.xparam);
			auto iterY = m_varMap.find(cut.yparam);
			if(iterX == m_varMap.end() || iterY == m_varMap.end())
			{
				std::cerr<<"Bad variables at CutHandler::IsInside for cut "<<cut.cutname<<" with x: "<<cut.xparam<<" y: "<<cut.yparam<<std::endl;
				return false;
			}
			else if(!cut.cut_ptr->IsInside(*(iterX->second), *(iterY->second)))
				return false;
		}

		return true;
	}
}