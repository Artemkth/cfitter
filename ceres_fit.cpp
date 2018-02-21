//core interface definitions
#include"ceres_fit.h" 
//for core algebraic functions
#include<cmath>

typedef std::pair<std::array<double, 3>, std::vector<double>> CPair;

//inline function to return rotation matrix(for fit with elliptic core)
inline Eigen::Matrix<double, 2, 2> rotMatrix(const double& phi) {
	Eigen::Matrix<double, 2, 2> retMatrix;
	retMatrix <<	cos(phi), -sin(phi),
					sin(phi), cos(phi);
	return retMatrix;
}

#ifndef M_PI
constexpr double M_PI = 3.14159265358979323846;
#endif

//functors for ceres fitting
//radial
class RadialGaussFtor : public ceres::SizedCostFunction<1, 2, 1, 1> {
private:
	CPair VFPoint;
public:
	RadialGaussFtor(const CPair& pnt) :VFPoint(pnt)
	{}
	bool Evaluate(double const* const* parameters,
		double* residuals,
		double** jacobians) const
	{
		//symbols for block parameters
		const double& x = VFPoint.first[0];
		const double& y = VFPoint.first[1];

		//first row of parameters is to be supposed coordinate of the peak
		const double& x0 = parameters[0][0];
		const double& y0 = parameters[0][1];

		//second row is only radius
		const double& r = parameters[1][0];
		//only failure scenario is one of the radii being lower or equal to zero
		if (r <= 0.)
			return false;

		//third row is "amplitude"
		const double& am = parameters[2][0];

		//useful thing to precalculate
		double expFactor = exp(-((x - x0) * (x - x0) + (y - y0) * (y - y0)) / (2. * r * r));

		//manual promises that residuals is never NULL
		//and there should be only one :p
		residuals[0] = VFPoint.second[2] - am * expFactor ;

		//if jacobian array is nullptr there is no jacobian terms to compute
		if (jacobians == nullptr)
			return true;

		//first row is derivatives due to assumed position
		if (jacobians[0] != nullptr)
		{
			jacobians[0][0] = - (x - x0) * expFactor / (2. * r * r);
			jacobians[0][1] = - (y - y0) * expFactor / (2. * r * r);
		}

		//radius row derivatives
		if (jacobians[1] != nullptr)
			jacobians[1][0] = -((x - x0) * (x - x0) + (y - y0) * (y - y0)) * expFactor / (2. * r * r * r);

		//amplitude derivative
		if (jacobians[2] != nullptr)
			jacobians[2][0] = - expFactor;

		//nothing here should fail, right?
		return true;
	}
};

//functors for ceres fitting
//elliptical
class ElipticalGaussFtor : public ceres::SizedCostFunction<1, 2, 2, 1, 1> {
private:
	CPair VFPoint;
public:
	ElipticalGaussFtor(const CPair& pnt) :VFPoint(pnt)
	{}
	bool Evaluate(double const* const* parameters,
		double* residuals,
		double** jacobians) const
	{
		//first row of parameters is to be supposed coordinate of the peak
		const double& x0 = parameters[0][0];
		const double& y0 = parameters[0][1];
		//relative position vector
		Eigen::Matrix<double, 2, 1> dispVect;
		dispVect << (VFPoint.first[0] - x0), (VFPoint.first[1] - y0);

		//second row is radii of minor and major axis
		const double& rx = parameters[1][0];
		const double& ry = parameters[1][1];
		//only failure scenario is one of the radii being exact zero(avoiding div 0)
		if (rx <= 0.||ry <= 0.)
			return false;
		//initialize diagonalized dispersion matrix
		Eigen::DiagonalMatrix<double, 2> rDiag;
		rDiag.diagonal() <<	1. / (2. * rx * rx), 1. / (2. * ry * ry);

		//third row is rotation angle
		const double& phi = parameters[2][0];
		auto rmat = rotMatrix(phi);
		//rotated dispersion matrix
		auto DispMat = rmat * rDiag * rmat.transpose();
		//exp parameter
		Eigen::Matrix<double, 1, 1> eparam = dispVect.transpose() * DispMat * dispVect;
		//useful thing to precalculate
		double expFactor = exp(-eparam(0, 0));

		//fourth row is "amplitude"
		const double& am = parameters[3][0];


		//manual promises that residuals is never NULL
		//and there should be only one :p
		residuals[0] = VFPoint.second[2] - am * expFactor;

		//if jacobian array is nullptr there is no jacobian terms to compute
		if (jacobians == nullptr)
			return true;

		//first row is derivatives due to assumed position
		if (jacobians[0] != nullptr)
		{
			auto grad = 2. * expFactor * DispMat * dispVect;
			jacobians[0][0] = - grad(0, 0);
			jacobians[0][1] = - grad(1, 0);
		}

		//radius row derivatives
		if (jacobians[1] != nullptr)
		{
			Eigen::DiagonalMatrix<double, 2> tmp;
			tmp.diagonal() << 1. / (rx * rx * rx), 0.;
			Eigen::Matrix<double,1,1> res = expFactor * dispVect.transpose() * rmat * tmp * rmat.transpose() * dispVect;
			jacobians[1][0] = - res(0, 0);
			tmp.diagonal() << 0., 1. / (ry * ry * ry);
			res = expFactor * dispVect.transpose() * rmat * tmp * rmat.transpose() * dispVect;
			jacobians[1][1] = - res(0, 0);
		}

		//angle row derivate
		if (jacobians[2] != nullptr)
		{
			auto drotMat = rotMatrix(phi + M_PI / 2.0);
			auto res = expFactor * dispVect.transpose() * ( drotMat * rDiag * rmat.transpose() + rmat * rDiag * drotMat.transpose() ) * dispVect;
			jacobians[2][0] = - res(0, 0);
		}

		//amplitude derivative
		if (jacobians[3] != nullptr)
			jacobians[3][0] = - expFactor;

		//nothing here should fail, right?
		return true;
	}
};

CeresFitEngine::CeresFitEngine()
{
	options.max_num_iterations = 100;
	options.linear_solver_type = ceres::DENSE_QR;
	options.minimizer_type = ceres::TRUST_REGION;
}

CGaussParams<double> CeresFitEngine::fitVortCoreCGauss(const std::vector<CPair>& pts, const std::array<double, 2>& initPos, const double& initR)
{
	CGaussParams<double> result;
	result.pos = initPos;
	result.r = initR;

	double am = 1.;

	ceres::Problem problem;

	for (const auto& x : pts)
		problem.AddResidualBlock(
			new RadialGaussFtor(x),
			nullptr,
			result.pos.data(),
			&result.r,
			&am
		);

	ceres::Solve(options, &problem, &summary);

	return result;
}

EGaussParams<double> CeresFitEngine::fitVortCoreEGauss(const std::vector<CPair>& pts, const std::array<double, 2>& initPos,
	const std::array<double, 2>& initR, const double& initAng)
{
	EGaussParams<double> result;
	result.pos = initPos;
	result.prAxis = initR;
	result.angle = initAng;

	double am = 1.;

	ceres::Problem problem;

	for (const auto& x : pts)
		problem.AddResidualBlock(
			new ElipticalGaussFtor(x),
			nullptr,
			result.pos.data(),
			result.prAxis.data(),
			&result.angle,
			&am
		);

	ceres::Solve(options, &problem, &summary);

	return result;
}
