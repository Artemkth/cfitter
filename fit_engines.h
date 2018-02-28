#pragma once
#include<vector>
#include<array>
#include<ostream>


//container for circular gauss fit params
template<typename T>
struct CGaussParams {
	std::array<T, 2> pos;
	T r;
};

template<typename T>
std::ostream& operator<<(std::ostream& os, const CGaussParams<T>& par)
{
	os<<par.pos[0]<<'\t'<<par.pos[1]<<'\t'<<par.r;
	return os;
}

//container for elyptical gauss fit params
template<typename T>
struct EGaussParams {
	std::array<T, 2> pos;
	std::array<T, 2> prAxis;
	T angle;
};

template<typename T>
std::ostream& operator<<(std::ostream& os, const EGaussParams<T>& par)
{
	os<<par.pos[0]<<'\t'<<par.pos[1]<<'\t'<<par.prAxis[0]<<'\t'<<par.prAxis[1]<<'\t'<<par.angle;
	return os;
}

template<typename T>
class FitEngine {
public:
	typedef  std::pair<std::array<T, 3>, std::vector<T>> CPair;

	//these should create and establish context needed by fitter
	FitEngine() = default;
	~FitEngine() = default;

	//prohibit copy construction
	FitEngine(const FitEngine&) = delete;
	FitEngine& operator=(const FitEngine&) = delete;

	

	//common interfaces to suit magnetization fit needs
	//fit vortex core with circular 2D gauss distribution
	virtual CGaussParams<T> fitVortCoreCGauss(const std::vector<CPair>& pts, const std::array<T,2>& initPos, const T& initR) = 0;

	//fit vortex core to ellyptical 2D gauss
	virtual EGaussParams<T> fitVortCoreEGauss(const std::vector<CPair>& pts, const std::array<T, 2>& initPos,
		const std::array<T,2>& initR, const T& initAng) = 0;
};
