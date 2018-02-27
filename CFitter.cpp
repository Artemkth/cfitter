//standard libraries
#include<iostream>
#include<fstream>
#include<thread>
#include<vector>
#include<string>
#include<regex>
#include<map>
#include<queue>
#include<condition_variable>

//importer library
#include"OVFReader.h"

//fitting backends
#include"ceres_fit.h"
#include"mathematica_fit.h"

//logging backend
#include<spdlog/spdlog.h>

//program options reader
#include<boost/program_options.hpp>
//boost filesystem
#include<boost/filesystem.hpp>

//math headers
#include<cmath>

struct Layer {
	size_t n{ 0 };
	std::array<double, 3> max{ 0., 0., 0. };
	std::vector<
		std::pair<
				std::array<double, 3>, std::vector<double>
			>
		> pnts{};
};


//in-plane deviation norm
inline double pl_norm(const std::array<double, 3>& x, const std::array<double, 3>& y)
{
	std::array<double, 2> d{ x[0] - y[0],x[1] - y[1] };
	return sqrt(d[0] * d[0] + d[1] * d[1]);
}

template<size_t n>
inline std::array<double,n> convBase(const std::array<double, n>& ref, double base)
{
	std::array<double, n> res;
	for (auto i = 0; i < n; i++)
		res[i] = ref[i] / base;
	return res;
}

std::vector<std::vector<Layer>> getData(const std::string& fName, const std::vector<unsigned int>& layers, double thresh, double rad, double ubase)
{
	auto rLog = spdlog::get("reader");
	rLog->debug("Trying to read a file: {}, cutoff height is {}, cutoff radius is {}", fName, thresh, rad);
	
	std::vector<OVFSegment> data{};

	try {
		data = readOVF(fName);
	}
	catch (std::exception& e)
	{
		rLog->error("ERROR reading \"{}\", exception occured: \"{}\"", fName, e.what());
		return {};
	}

	rLog->debug("{}: read {} segment/s", fName, data.size());

	std::vector<std::vector<Layer>> result{  };

	size_t segCnt = 0;
	for (auto& seg : data)
	{
		auto field = seg.cast<double>();

		std::vector<Layer> segJobs{};
		for (auto& layer : layers)
		{
			std::vector<std::pair<std::array<double, 3>, std::vector<double>>> slice{};
			try {
				slice = field.selectInLayer([&thresh](const std::array<double, 3>& x, const std::vector<double>& y) {return fabs(y[2]) > thresh ? true : false; }, layer);
			}
			catch (std::exception& e)
			{
				rLog->error("ERROR reading \"{}\" segment {} layer {}, exception occured: \"{}\"", fName, segCnt, layer, e.what());
				return {};
			}

			auto ctr = slice[0];
			for (const auto& x : slice)
				if (fabs(ctr.second[2]) < fabs(x.second[2]))
					ctr = x;

			struct Layer ret;
			ret.n = layer;
			ret.max = convBase(ctr.first, ubase);
			for (const auto& x : slice)
				if (rad < 0 || pl_norm(x.first, ctr.first) <= rad * ubase)
					ret.pnts.push_back({ convBase(x.first, ubase), x.second });

			rLog->debug("{}: Selected the points in layer {} of segment {}, highest value is {} at {} {}, {} points in total.", fName, layer, segCnt, ctr.second[2], ctr.first[0] / ubase, ctr.first[1] / ubase, ret.pnts.size());
			segJobs.push_back(ret);
		}
		segCnt++;
		result.push_back(segJobs);
	}
	return result;
}

//template for return result
template<typename T>
using ResList =  std::vector<
	std::pair<
		std::string,
		std::vector<std::map<size_t, T>>
	>
>;

//worker class, template params T -- proper return structure, U -- class name of fitter
template<typename T, typename U>
ResList<T> worker(const std::vector<std::string>& files, const std::vector<unsigned int>& layers, bool cont, const size_t tcount,
	const double thresh, const double radius, const double ubase)
{
	//check if the function is called for correct 
	static_assert(std::is_same<T, CGaussParams<double>>::value || std::is_same<T, EGaussParams<double>>::value, "called a worker() for an unsupported fit function class");
	static_assert(std::is_base_of<FitEngine<double>, U>, "provided wrong incompatible fitter backend");

	//get the logger
	auto conLog = spdlog::get("main");
	
	//make map with all the files
	ResList<T> results{};
	for (const auto& file : files)
		results.push_back({ file, {} });
	//and its mutexes
	std::mutex resultMutex;
	std::condition_variable resultCV;

	//make queue for importing files
	std::queue<std::string> importJobs{};
	std::vector<std::string> importsDone{};
	for (const auto& file : files)
		importJobs.push(file);
	//and its mutex
	std::mutex impJobMutex;

	//ticket structure
	struct Ticket {
		//file number, nothing to do with Linking Park songs :p
		std::string fileName{ "" };
		//supposed maximum
		std::array<double, 3> max{ 0.,0.,0. };
		//segment number
		size_t segN{ 0 };
		//layer number
		size_t layerN{ 0 };
		//data points
		std::vector<std::pair<std::array<double, 3>, std::vector<double>>> pts{};
	};

	//queue for tickets
	std::vector<Ticket> tickets{};
	std::mutex ticketMutex;
	std::condition_variable ticketCV;
	std::atomic_bool moreTicketsExpected{ true };

	//checking if dependency exists
	auto existDep = [&](const Ticket& t) {
		std::string prev{ "" };
		//first find which file is previous
		for(auto it = files.begin(); it!= files.end(); ++it)
			if (it == t.fileName)
			{
				if (it == files.begin())
					return true;

				prev = *std::prev(it);
			}
		//find result column
		auto it = results.begin();
		for (; it != results.end(); ++it)
			if (it->first == prev)
				break;

		if (it->second.empty())
			return false;

		return it->second.find(t.layerN).second;
	};

	//definition of worker thread
	auto workerThread = [&]() {
		//attach and initialize a fitter to a thread
		U fitEngine;

		while (true) 
		{
			//ticket for current thread
			Ticket curT;
			std::string curFile;

			//learn the current state and take decision
			//don't take ownership of import lock before ticket lock because tickets can be accessed seprately
			ticketMutex.lock();
			impJobMutex.lock();
			//state of both queues
			bool impJobsAvailable = !importJobs.empty();
			size_t ticketN = tickets.size();

			//status flags for the thread, look for cubo maps in notebook for verification
			bool waitForTicket{ !impJobsAvailable && ticketN == 0 && moreTicketsExpected };
			bool isFitting{ ticketN >= tcount || waitForTicket || (!impJobsAvailable&&ticketN != 0) };
			bool isImporting{ impJobsAvailable && ticketN < tcount };

			if (isFitting && !waitForTicket && !cont)
			{
				curT = std::move(tickets.front());
				tickets.erase(tickets.begin());
			}
			if (isImporting)
			{
				curFile = std::move(importJobs.front());
				inportJobs.pop();
			}

			//don't forget to unlock everything
			ticketMutex.unlock();
			impJobMutex.unlock();

			if (isImporting)
			{
				//import code to be implemented here
				//tickets to be merged
				std::queue<Ticket> fileTickets{};

				//data
				auto rawData = getData(curFile, layers, thresh, rad, ubase);
				for (size_t seg = 0; seg < rawData.size(); seg++)
					for (auto& x : rawData[seg])
					{
						Ticket t;
						t.fileName = curFile;
						t.segN = seg;
						t.layerN = x.n;
						t.max = std::move(x.max);
						t.pts = std::move(x.pnts);

						//push new tickets into buffer
						fileTickets.push(std::move(t));
					}

				ticketMutex.lock();
				while (!fileTickets.empty())
				{
					tickets.push_back(std::move(fileTickets.front()));
					fileTickets.pop();
				}

				//check if there should be any more jobs
				impJobMutex.lock();
				importsDone.push_back(curFile);
				moreTicketsExpected = importsDone.size() != files.size();
				impJobMutex.unlock();

				ticketCV.notify_all();
				ticketMutex.unlock();
				
			}
			else if (isFitting)
			{
				//fitting code here

				//handling for getting ticket in case it is needed(either there was none or there is continuity condition)
				//first wait for tickets to appear if there is none
				if (waitForTicket) {
					//unique lock should lock upon creation
					std::unique_lock<std::mutex> ticketLock(ticketMutex);
					//if before this point new tickets appeared, wait for tickets
					if(tickets.empty())
						ticketCV.wait(ticketLock);
				}
				//getting ticket
				bool gotTicket{ !cont && waitForTicket };
				while (!gotTicket) {
					ticketMutex.lock();
					if (cont)
						resultMutex.lock();
					//break if there is no reason for waiting
					if (!moreTicketsExpected&&tickets.empty())
					{
						if (cont)
							resultMutex.unlock();
						ticketMutex.unlock();
						break;
					}

					for (auto it = tickets.begin(); it != tickets.end(); ++it)
						if (!cont || existDep(*it))
						{
							curT = *it;
							gotTicket = true;
							break;
						}

					if (cont)
						resultMutex.unlock();
					ticketMutex.unlock();
				}

				//if no ticket received anyway continue onto next cycle
				if (!gotTicket)
					continue;

				//import code after this point
			}
			else
			{
				break;
			}
		}
	};



	return results;
}


int main(int argc, char* argv[])
{
	//initialize main spdlog logger
	auto conLog = spdlog::stdout_color_st("main");
	auto rLog = spdlog::stdout_color_mt("reader");

	//thread count constant for future paralelism
	unsigned int threadCnt{ 0 };
	//is sequence of files to be treated as continuous?
	bool isCont{ false };
	//if log needs to be verbous?
	bool isVerb{ false };
	//list of input files
	std::vector<std::string> fileNames{};
	//list of layer to be fit
	std::vector<unsigned int> layers{};
	//second cut parameter
	double cutR{ -1. };
	//unit base
	double unit_base{ 1.e-9 };
	//cut thresh
	double thresh{ .1 };

	//boost program options initialization
	boost::program_options::options_description desc("Supported options");
	desc.add_options()
		("help,h", "display this message")
		("input-file", boost::program_options::value<std::vector<std::string>>(&fileNames)->required(), "input OVF file names")
		("layers,l", boost::program_options::value<std::vector<unsigned int>>(&layers)->required()->multitoken(), "list of layer number to be fit")
		("backend,b", boost::program_options::value<std::string>()->default_value("ceres"), "Backend to use for fitting core-position. Available options are: ceres/mathematica")
		("type", boost::program_options::value<std::string>()->default_value("elliptical", "Fit function type: elliptical/circular"))
		//due to includes from wstp.h fucking up everything
#if defined(_WIN32)||defined(WIN32)
		("threads,t", boost::program_options::value<unsigned int>(&threadCnt)->default_value(min(4, std::thread::hardware_concurrency())))
#elif
		("threads", boost::program_options::value<unsigned int>(&threadCnt)->default_value(std::min(4, std::thread::hardware_concurrency())))
#endif
		("output,o", boost::program_options::value<std::string>()->default_value(""), "output text file name, by default table is spit into stdout")
		("continuous,c", boost::program_options::bool_switch(&isCont)->default_value(false), "treat files as continuous series to improve fit speed")
		("verbose,v", boost::program_options::bool_switch(&isVerb)->default_value(false), "enable debug output")
		("sec-cut-r", boost::program_options::value<double>(&cutR)->default_value(-1.0), "secondary cut radius")
		("unit-base", boost::program_options::value<double>(&unit_base)->default_value(1.e-9), "unit base")
		("core-thresh", boost::program_options::value<double>(&thresh)->default_value(.3), "threshold for points to be considered core")
		;

	//set special status to input file
	boost::program_options::positional_options_description posDesc;
	posDesc.add("input-file", -1);

	//parse commandline
	boost::program_options::variables_map varMap;
	boost::program_options::store(
		boost::program_options::command_line_parser(argc, argv).options(desc).positional(posDesc).run(),
		varMap
	);

	if (varMap.count("help"))
	{
		//if help is required do it and quit
		std::cout << desc << '\n';
		return 0;
	}

	//TODO: optional exception handling here 
	//check if required arguments are provided
	try {
		boost::program_options::notify(varMap);
	}
	catch (boost::program_options::error& e) {
		conLog->error(e.what());
		return -1;
	}

	if (isVerb)
	{
		conLog->set_level(spdlog::level::debug);
		rLog->set_level(spdlog::level::debug);
	}

	//internal variables
	enum class Backend {
		Mathematica,
		Ceres,
		Unknown
	};

	Backend bend{ Backend::Unknown };

	enum class FitFunction {
		EllipticalGauss,
		CircularGauss,
		Unknown
	};

	FitFunction ffunc{ FitFunction::Unknown };

	//set few more internal switches
	{
		std::string beString = varMap["backend"].as<std::string>();
		if (beString == "ceres")
			bend = Backend::Ceres;
		else if (beString == "mathematica")
			bend = Backend::Mathematica;
		else
		{
			conLog->critical("Unknown backend type: {}", beString);
			return -1;
		}
	}
	{
		std::string tString = varMap["type"].as<std::string>();
		if (tString == "elliptical")
			ffunc = FitFunction::EllipticalGauss;
		else if (tString == "circular")
			ffunc = FitFunction::CircularGauss;
		else
		{
			conLog->critical("Unknown fit function type: {}", tString);
			return -1;
		}
	}

	conLog->debug("Launch parameters parsed:");
	//TODO implement exploding wildcards for windows
	for(const auto& file: fileNames)
		if (!boost::filesystem::exists(file) || !boost::filesystem::is_regular_file(file))
		{
			conLog->error("File \"{}\" does not exist, aborting.", file);
			return -1;
		}

	//block for logging file names
	{
		std::string cumulNames{ "" };
		for (const auto& file : fileNames)
		{
			cumulNames += file;
			cumulNames += " ";
		}

		conLog->debug("File names to open: " + cumulNames);
	}

	//block for logging layers
	{
		std::string cumulNames{ "" };
		for (const auto& layer : layers)
		{
			cumulNames += std::to_string(layer);
			cumulNames += " ";
		}

		conLog->debug("Layers to fit: " + cumulNames);
	}

	//check other fit options
	conLog->info("Starting fit using {} fitter threads of {} using {} vortex core model", threadCnt, bend == Backend::Ceres ? "Ceres" : "Mathematica", 
		ffunc == FitFunction::CircularGauss ? "circular" : "elliptical");


	//all the parallel stuff begins here
	for (const auto& file : fileNames)
		getData(file, layers, thresh, cutR, unit_base);


	return 0;
}
