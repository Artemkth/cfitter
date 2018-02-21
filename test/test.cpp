#include<iostream>
#include"OVFReader.h"
#include"ceres_fit.h"
#include<cstdlib>
#include"mathematica_fit.h"
#include<cmath>

bool peakPred(const std::array<double,3>& x, const std::vector<double>& y)
{
	if (y[2] > .25)
		return true;
	return false;
}

template<typename T>
void PrintVF(const std::vector<std::pair<std::array<T, 3>, std::vector<T>>>& data)
{
	std::array<std::string, 3> xlabel = { "x","y","z" };
	std::array<std::string, 3> mlabel = { "mx","my","mz" };
	for (const auto& x : data)
	{
		//Print x
		for (int i = 0; i < 3; i++)
			std::cout << xlabel[i] << "=" << x.first[i] << (i == 2 ? '\n' : '\t');
		//Print m
		for (int i = 0; i < 3; i++)
			std::cout << mlabel[i] << "=" << x.second[i] << (i == 2 ? '\t' : '\t');
		std::cout << "m^2: " << x.second[0] * x.second[0] + x.second[1] * x.second[1] + x.second[2] * x.second[2] << std::endl;
	}
}

int main(int argc, char* argv[])
{
	std::cout << "Launched from: " << argv[0] << std::endl;
	auto data = readOVF("../test_files/test_large.ovf");
	std::cout << "read successfully huh? read " << data.size() << " blocks" << std::endl;
	std::cout << "first block with description: \"" << data[0].descr << "\"" << std::endl;
	std::cout << " Trying to pick reasonable ammount of points just for lulz" << std::endl;
	auto vfield = data[0].cast<double>();
	auto exPoint = vfield.getPoint(0, 64, 0);
	PrintVF<double>({ exPoint });

	auto slice = vfield.selectInLayer(peakPred, 0);
	std::cout << "Took a slice of data points in 1st layer with mz > .25" << std::endl;
	PrintVF(slice);
	//try to find maximum
	auto ctr = slice[0];
	for (const auto& x : slice)
		if (fabs(ctr.second[2]) < fabs(x.second[2]))
			ctr = x;
	std::cout << "Highest point with magnitude mz=" << ctr.second[2] << ", found at x=" << ctr.first[0] << ", y=" << ctr.first[1] << std::endl;
	std::array<double, 2> initPos = { ctr.first[0], ctr.first[1] };

	std::cout << "Trying to start a mathlink" << std::endl;
	MathematicaFitter<double> fit;
	std::cout << "Connection to math kernel established" << std::endl;

	std::cout << "Trying to fit the core shape for radial shape" << std::endl;
	auto cfitRes = fit.fitVortCoreCGauss(slice, initPos, 5.e-9);
	std::cout << "Results are: x0=" << cfitRes.pos[0] << ", y0=" << cfitRes.pos[1] << ", r=" << cfitRes.r << std::endl;

	std::cout << "Trying to fit the core shape for elliptical shape" << std::endl;
	auto cfitRes2 = fit.fitVortCoreEGauss(slice, initPos, { 5.e-9 ,5.e-9 }, 0);
	std::cout << "Results are: x0=" << cfitRes2.pos[0] << ", y0=" << cfitRes2.pos[1] << 
		", rx=" << cfitRes2.prAxis[0] << ", ry=" << cfitRes2.prAxis[1] << ", phi=" << cfitRes2.angle << std::endl;

	std::cout << "Trying ceres backend for lulz" << std::endl;
	CeresFitEngine cfit;

	initPos[0] *= 1.e9;
	initPos[1] *= 1.e9;

	for (auto& x : slice)
		for (auto& y : x.first)
			y *= 1.e9;

	std::cout << "Trying to fit the core shape for radial shape" << std::endl;
	auto cfitCeres = cfit.fitVortCoreCGauss(slice, initPos, 5);
	std::cout << "Results are: x0=" << cfitCeres.pos[0] << ", y0=" << cfitCeres.pos[1] << ", r=" << cfitCeres.r << std::endl;

	std::cout << "Trying to fit the core shape for elliptical shape" << std::endl;
	auto cfitCeres2 = cfit.fitVortCoreEGauss(slice, initPos, { 5 ,5 }, 0);
	std::cout << "Results are: x0=" << cfitCeres2.pos[0] << ", y0=" << cfitCeres2.pos[1] <<
		", rx=" << cfitCeres2.prAxis[0] << ", ry=" << cfitCeres2.prAxis[1] << ", phi=" << cfitCeres2.angle << std::endl;

	return 0;
}
