#pragma once

//ofc include fit engines wouldn't be good otherwise
#include"fit_engines.h"

//stuff specific to ceres
//macro for windows machines
#if defined(_WIN32)||defined(WIN32)
#define GLOG_NO_ABBREVIATED_SEVERITIES
#endif

#include<ceres/ceres.h>
//algebra library
#include<Eigen/Dense>

class CeresFitEngine :public FitEngine<double>
{
private:
	//interface implemention stuff
private:
	ceres::Solver::Options options;
	ceres::Solver::Summary summary;

public:
	CeresFitEngine();
	~CeresFitEngine()
	{}

	CGaussParams<double> fitVortCoreCGauss(const std::vector<CPair>& pts, const std::array<double, 2>& initPos, const double& initR);

	EGaussParams<double> fitVortCoreEGauss(const std::vector<CPair>& pts, const std::array<double, 2>& initPos,
		const std::array<double, 2>& initR, const double& initAng);
};
