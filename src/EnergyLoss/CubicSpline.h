/*
CubicSpline.h 
Class for generating cubic splines of data in tables or held in memory. Cubic splines are a form of interpolation,
somewhat more advanced than linear interpolation, but not so complicated that it significantly slows down a calculation.
For more information see Wikipedia or Num. Rec. in C.

Gordon M. May 2021
*/
#ifndef CUBICSPLINE_H
#define CUBICSPLINE_H

#include <vector>
#include <string>

namespace PunchTable {

	//struct holding final spline info
	struct Spline
	{
		double y1=0, y2=0;
		double x1=0, x2=0;
		double k1=0, k2=0;
	};

	class CubicSpline
	{
	public:
		CubicSpline();
		CubicSpline(const std::string& filename);
		CubicSpline(const std::vector<double>& x, const std::vector<double>& y);
		~CubicSpline();
		void ReadFile(const std::string& filename);
		inline void ReadData(const std::vector<double>& x, const std::vector<double>& y)
		{ 
			m_dataX = x; 
			m_dataY = y; 
			m_validFlag=true;
			MakeSplines();
		}
		bool IsValid() { return m_validFlag; }
		double Evaluate(double x);
		double EvaluateROOT(double* x, double* p); //for plotting as ROOT function

	private:
		void MakeSplines();

		std::vector<double> m_dataX;
		std::vector<double> m_dataY;
		std::vector<Spline> m_splines;

		bool m_validFlag;
	};

}
#endif