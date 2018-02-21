#pragma once
#include"fit_engines.h"
#include<string.h>
#include<stdexcept>

//includes for mathematica
#include<wstp.h>

template<typename T>
class MathematicaFitter : FitEngine<T>
{
private:
	//mathematica env variables
	WSENV env;
	WSLINK link;
	int error;

	const std::string CGaussFitModel = "Times[am,Power[E,Times[Rational[1,2],Power[r,-2],Plus[Times[-1,Power[Plus[x,Times[-1,x0]],2]],Times[-1,Power[Plus[y,Times[-1,y0]],2]]]]]]";
	const std::string EGaussFitModel = "Times[am,Power[E,Times[Rational[1,4],Power[a,-2],Power[b,-2],Plus[Times[-1,Plus[Power[a,2],Power[b,2]],Plus[Power[Plus[x,Times[-1,x0]],2],Power[Plus[y,Times[-1,y0]],2]]],Times[Plus[a,Times[-1,b]],Plus[a,b],Plus[Times[Plus[x,Times[-1,x0],y,Times[-1,y0]],Plus[x,Times[-1,x0],Times[-1,y],y0],Cos[Times[2,fi]]],Times[2,Plus[x,Times[-1,x0]],Plus[y,Times[-1,y0]],Sin[Times[2,fi]]]]]]]]]";

	inline void PutReal(const T& x)
	{
		if (std::is_same<T, double>::value)
			WSPutReal64(link, x);
		if (std::is_same<T, float>::value)
			WSPutReal32(link, x);
	}
public:
	MathematicaFitter()
	{
		env = WSInitialize((WSEnvironmentParameter)0);
		if (env == (WSENV)0)
			throw std::runtime_error("Could not set up mathematica enviroment");

		link = WSOpenString(env, "-linkmode launch -linkname \'math -mathlink\'", &error);
		if (link == (WSLINK)0 || error != WSEOK)
			throw std::runtime_error("Could not open math link");
	}

	~MathematicaFitter()
	{
		WSPutFunction(link, "Exit", 0);

		WSClose(link);
		WSDeinitialize(env);
	}

	EGaussParams<T> fitVortCoreEGauss(const std::vector<std::pair<std::array<T,3>,std::vector<T>>>& pts, const std::array<T, 2>& initPos,
		const std::array<T, 2>& initR, const T& initAng)
	{
		struct EGaussParams<T> result;
		//sending the package to evaluate for mathematica
		//EvaluatePacket[
		//	ToExpression["{x0, y0 ...}"]/.NonlinearModelFit[{},Evaluate[ToExpression["model"],{{x0, 9e-9},...},{x,y}]]@@{"BestFitParameters"}
		//]
		WSPutFunction(link, "EvaluatePacket", 1);
		{
			WSPutFunction(link, "ReplaceAll", 2);
			WSPutFunction(link, "ToExpression", 1);
			{
				WSPutString(link, "{x0,y0,a,b,fi,am}");
			}
			WSPutFunction(link, "Apply", 2);
			{
				WSPutFunction(link, "NonlinearModelFit", 4);
				{
					//First parameter is list of tripplets
					WSPutFunction(link, "List", pts.size());
					for (auto it = pts.begin(); it != pts.end(); ++it)
					{
						WSPutFunction(link, "List", 3);
						//put x and y
						for (size_t i = 0; i < 2; i++)
							PutReal(it->first[i]);
						//put mz
						PutReal(it->second[2]);
					}
					//Second parameter is model for fitting
					WSPutFunction(link, "Evaluate", 1);
					{
						WSPutFunction(link, "ToExpression", 1);
						{
							WSPutString(link, EGaussFitModel.c_str());
						}
					}

					//copied from full form in mathematica
					//Third parameter is model parameter descriptions
					WSPutFunction(link, "List", 6);
					{
						WSPutFunction(link, "List", 2);
						{
							WSPutSymbol(link, "x0");
							PutReal(initPos[0]);
						}
						WSPutFunction(link, "List", 2);
						{
							WSPutSymbol(link, "y0");
							PutReal(initPos[1]);
						}
						WSPutFunction(link, "List", 2);
						{
							WSPutSymbol(link, "a");
							PutReal(initR[0]);
						}
						WSPutFunction(link, "List", 2);
						{
							WSPutSymbol(link, "b");
							PutReal(initR[1]);
						}
						WSPutFunction(link, "List", 2);
						{
							WSPutSymbol(link, "fi");
							PutReal(initAng);
						}
						WSPutFunction(link, "List", 2);
						{
							WSPutSymbol(link, "am");
							WSPutReal64(link, 1.);
						}
					}

					//Fourth argument is list of model free variables
					WSPutFunction(link, "List", 2);
					{
						WSPutSymbol(link, "x"); WSPutSymbol(link, "y");
					}
				}
				WSPutFunction(link, "List", 1);
				{
					WSPutString(link, "BestFitParameters");
				}
			}
		}
		WSEndPacket(link);
		WSFlush(link);

		//receive an answer
		int pkt;

		while ((pkt = WSNextPacket(link), pkt) && pkt != RETURNPKT) {
			WSNewPacket(link);
			if (WSError(link))
				throw std::runtime_error("Error reading a response: " );
		}

		int len;
		WSTestHead(link, "List", &len);
		if (len != 6)
			throw std::runtime_error("Got unexpected response");

		double* retArray = new double[len];
		for (long i = 0; i < len; i++)
			WSGetReal64(link, retArray + i);

		//positions
		result.pos[0] = retArray[0];
		result.pos[1] = retArray[1];
		//core radius
		result.prAxis[0] = retArray[2];
		result.prAxis[1] = retArray[3];
		//angle
		result.angle = retArray[4];

		free(reinterpret_cast<void*>(retArray));

		return result;
	}

	CGaussParams<T> fitVortCoreCGauss(const std::vector<std::pair<std::array<T, 3>, std::vector<T>>>& pts, const std::array<T, 2>& initPos, const T& initR)
	{
		struct CGaussParams<T> result;
		//sending the package to evaluate for mathematica
		//EvaluatePacket[
		//	ToExpression["{x0, y0 ...}"]/.NonlinearModelFit[{},Evaluate[ToExpression["model"],{{x0, 9e-9},...},{x,y}]]@@{"BestFitParameters"}
		//]
		WSPutFunction(link, "EvaluatePacket", 1);
		{
			WSPutFunction(link, "ReplaceAll", 2);
			WSPutFunction(link, "ToExpression", 1);
			{
				WSPutString(link, "{x0,y0,r,am}");
			}
			WSPutFunction(link, "Apply", 2);
			{
				WSPutFunction(link, "NonlinearModelFit", 4);
				{
					//First parameter is list of tripplets
					WSPutFunction(link, "List", pts.size());
					for (auto it = pts.begin(); it != pts.end(); ++it)
					{
						WSPutFunction(link, "List", 3);
						//put x and y
						for (size_t i = 0; i < 2; i++)
							PutReal(it->first[i]);
						//put mz
						PutReal(it->second[2]);
					}
					//Second parameter is model for fitting
					WSPutFunction(link, "Evaluate", 1);
					{
						WSPutFunction(link, "ToExpression", 1);
						{
							WSPutString(link, CGaussFitModel.c_str());
						}
					}

					//copied from full form in mathematica
					//Third parameter is model parameter descriptions
					WSPutFunction(link, "List", 4);
					{
						WSPutFunction(link, "List", 2);
						{
							WSPutSymbol(link, "x0");
							PutReal(initPos[0]);
						}
						WSPutFunction(link, "List", 2);
						{
							WSPutSymbol(link, "y0");
							PutReal(initPos[1]);
						}
						WSPutFunction(link, "List", 2);
						{
							WSPutSymbol(link, "r");
							PutReal(initR);
						}
						WSPutFunction(link, "List", 2);
						{
							WSPutSymbol(link, "am");
							WSPutReal64(link, 1.);
						}
					}

					//Fourth argument is list of model free variables
					WSPutFunction(link, "List", 2);
					{
						WSPutSymbol(link, "x"); WSPutSymbol(link, "y");
					}
				}
				WSPutFunction(link, "List", 1);
				{
					WSPutString(link, "BestFitParameters");
				}
			}
		}
		WSEndPacket(link);
		WSFlush(link);

		//receive an answer
		int pkt;

		while ((pkt = WSNextPacket(link), pkt) && pkt != RETURNPKT) {
			WSNewPacket(link);
			if (WSError(link))
				throw std::runtime_error("Error reading a response: ");
		}

		int len;
		WSTestHead(link, "List", &len);
		if (len != 4)
			throw std::runtime_error("Got unexpected response");

		double* retArray = new double[len];
		for (long i = 0; i < len; i++)
			WSGetReal64(link, retArray + i);

		//positions
		result.pos[0] = retArray[0];
		result.pos[1] = retArray[1];
		//core radius
		result.r = retArray[2];

		free(reinterpret_cast<void*>(retArray));

		return result;
	}
};
