#ifndef CUT_HANDLER_H
#define CUT_HANDLER_H

#include <vector>
#include <string>
#include <unordered_map>
#include "TFile.h"
#include "TCutG.h"
#include "CalDict/DataStructs.h"

namespace SabreRecon {

	struct ReconCut
	{
		ReconCut() {}
		ReconCut(const std::string& file, const std::string& cut, const std::string& x, const std::string& y) :
			xparam(x), yparam(y), filename(file), cutname(cut)
		{
		}

		TFile* file_ptr = nullptr;
		TCutG* cut_ptr = nullptr;
		std::string xparam = "";
		std::string yparam = "";
		std::string filename = "";
		std::string cutname = "";
	};

	class CutHandler
	{
	public:
		CutHandler();
		CutHandler(const std::vector<ReconCut>& cuts, CalEvent* event);
		~CutHandler();

		void InitCuts(const std::vector<ReconCut>& cuts);
		void InitEvent(CalEvent* event);
		inline const bool IsValid() const { return m_isValid; }

		bool IsInside();

	private:
		std::vector<ReconCut> m_cuts;
		CalEvent* m_eventPtr;
		bool m_isValid;

		bool m_eventInit;
		bool m_cutInit;

		std::unordered_map<std::string, double*> m_varMap;
	};
}

#endif