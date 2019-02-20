//standard libraries
#include<iostream>
#include<iomanip>
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
#include "spdlog/sinks/stdout_color_sinks.h"

//program options reader
#include<boost/program_options.hpp>
//boost filesystem
#include<boost/filesystem.hpp>

//math headers
#include<cmath>

struct Layer {
	size_t n{ 0 };
	std::array<double, 2> max{ 0., 0. };
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
			ret.max = convBase(std::array<double,2>{ ctr.first[0], ctr.first[1] }, ubase);
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
using ResList =  std::vector<std::pair<
		std::string,
		std::vector<
			std::map<size_t, T>		//vector of segments(size_t indexed maps of 
		>
>>;

//fitter helper, hack arround strong typization
template<typename T>
T fitHelper(FitEngine<double>& engine,const std::vector<std::pair<std::array<double, 3>, std::vector<double>>>& pts, const std::array<double, 2> &maxPos, bool cont, const T* const prev, double ubase)
{
	//static_assert(std::is_base_of<FitEngine<double>, U>::value, "Engine should be derived from FitEngine class");
	bool useCont = cont&&prev!=nullptr;
	//FitEngine<double>& engine{eng};
	void* result{nullptr};
	CGaussParams<double> resCVal;
	if(std::is_same<T, CGaussParams<double>>::value){
		resCVal = engine.fitVortCoreCGauss(pts, useCont? (reinterpret_cast<const CGaussParams<double>* const>(prev)->pos):maxPos,
										   useCont? (reinterpret_cast<const CGaussParams<double>* const>(prev)->r): (5.e-9/ubase));
		
		result = reinterpret_cast<void*>(&resCVal);
	}
	EGaussParams<double> resEVal;
	if(std::is_same<T, EGaussParams<double>>::value){
		std::array<double, 2> prAxisBUp{5.e-9/ubase, 5.e-9/ubase};
		resEVal = engine.fitVortCoreEGauss(pts, useCont? (reinterpret_cast<const EGaussParams<double>* const>(prev)->pos):maxPos,
										   useCont? (reinterpret_cast<const EGaussParams<double>* const>(prev)->prAxis):prAxisBUp,
										   useCont? (reinterpret_cast<const EGaussParams<double>* const>(prev)->angle): 0);
		result = reinterpret_cast<void*>(&resEVal);
	}
	return *reinterpret_cast<T*>(result);
}

//worker class, template params T -- proper return structure, U -- class name of fitter
template<typename T, typename U>
ResList<T> worker(const std::vector<std::string>& files, const std::vector<unsigned int>& layers, bool cont, const size_t tcount,
	const double thresh, const double radius, const double ubase, const bool verbose)
{
	//check if the function is called for correct 
	static_assert(std::is_same<T, CGaussParams<double>>::value || std::is_same<T, EGaussParams<double>>::value, "called a worker() for an unsupported fit function class");
	static_assert(std::is_base_of<FitEngine<double>, U>::value, "provided wrong incompatible fitter backend");

	//get the logger
	auto conLog = spdlog::get("main");
	
	//datatypes used
	//output result
	struct Result{
		//flag to set when the result is done
		bool ready {false};
		//flag to set if previous file is not found
		bool orphaned{ false };
		
		//indexing stuff
		size_t fileNumb	{ 0 };
		
		//result storage
		T fit_res;
		
		//reference to required segment
		Result* req {nullptr};
	};
	
	//ticket structure
	struct Ticket {
		//file index number
		size_t fileNumb{ 0 };
		//supposed maximum
		std::array<double, 2> max{ 0., 0. };
		//data read from importer
		std::vector<std::pair<std::array<double, 3>, std::vector<double>>> pts{};
		//result drop off point
		Result *dropRes{nullptr};
	};
	
	//and result list
	std::vector<
		std::pair<
			bool,
			std::vector<
				std::map<
					size_t,
					Result
				>
			>
		>
	> results;
	
	//initialize it for array of 
	for(size_t i = 0; i < files.size(); i++)
		results.push_back({false, {}});
	//and its mutexes
	std::mutex resultMutex;
	std::condition_variable resultCV;

	//make queue for importing files
	std::queue<std::pair<size_t, std::string>> importJobs{};
	for (size_t i = 0; i < files.size(); i++)
		importJobs.push({i, files[i]});
	//and its mutex
	std::mutex impJobMutex;

	//queue for tickets
	std::vector<Ticket> tickets{};
	std::mutex ticketMutex;
	std::condition_variable ticketCV;
	std::atomic_bool moreTicketsExpected{ true };

	//checking if dependency exists
	auto existDep = [&](const Ticket& t) {
		if (t.fileNumb == 0)
			return true;
		if(!results[t.fileNumb - 1].first)
			return false;
		return t.dropRes->orphaned || (t.dropRes->req != nullptr && t.dropRes->req->ready);
	};
	
	//imported counter
	size_t impCount { 0 };

	//definition of worker thread
	auto workerThread = [&](const size_t threadN) {
		//attach and initialize a fitter to a thread
		U fitEngine;
		//and open a logger for it
		auto log = spdlog::stdout_color_st("worker#"+std::to_string(threadN));
		if(verbose)
			log->set_level(spdlog::level::debug);

		while (true) 
		{
			//ticket for current thread
			Ticket curT;
			std::pair<size_t, std::string> curFile;

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
				importJobs.pop();
			}

			//don't forget to unlock everything
			ticketMutex.unlock();
			impJobMutex.unlock();

			if (isImporting)
			{
				log->debug("decided to import a file #{}, name is: \"{}\"", curFile.first, curFile.second);
				//import code to be implemented here
				//tickets to be merged
				std::queue<Ticket> fileTickets{};

				//data
				//TODO implement exception handling
				auto rawData = getData(curFile.second, layers, thresh, radius, ubase);
				log->debug("Reading #{}: \"{}\" done, {} segments found.", curFile.first, curFile.second, rawData.size());
				
				//lock in result mutex
				resultMutex.lock();
				//reserve space for results
				auto& curResult = results[curFile.first];
				//set read status to true
				curResult.first = true;
				//fill the spaces
				auto& curSegs = curResult.second;
				size_t sCnt = 0;
				for (auto& seg: rawData)
				{
					curSegs.push_back({});
					auto& lastSeg = curSegs.back();
					//reserve places
					
					//temporary container for tickets
					std::queue<Ticket> tmpTickets{};
					
					for(auto& layer: seg)
					{
						//check if result has 0 points in slice
						if(layer.pnts.size() == 0)
						{
							log->error("No vortex core found in segment {} layer {} of file: \"{}\"",
								sCnt, layer.n, curFile.second
							);
							continue;
						}
						
						//fill out ticket and result
						Ticket t;
						t.fileNumb = curFile.first;
						t.max = layer.max;
						t.pts = std::move(layer.pnts);
						//result
						Result r;
						r.ready = false;
						r.fileNumb = t.fileNumb;
						
						//push it into result map
						lastSeg[layer.n] = r;
						Result* thisLayer = &lastSeg[layer.n];
						t.dropRes = thisLayer;
						fileTickets.push(std::move(t));
						
						//if the whole deal is required to be continuous crosslink
						if(cont)
						{
							//first try linking this file to previous(only if next file is read and 
							//this file isn't first)
							if(t.fileNumb != 0&&results[curFile.first - 1].first){
								if (results[curFile.first - 1].second.size() <= sCnt)
								{
									log->warn("Current segment {} of \"{}\", didn't exist in previous file", sCnt, curFile.second);
									thisLayer->orphaned = true;
								}
								else{
									auto& prevFileMap = results[curFile.first - 1].second[sCnt];
									auto prev = prevFileMap.find(layer.n);
									if(prev != prevFileMap.end())
										thisLayer->req = &prev->second;
									else
									{
										log->warn("Current layer {} in segment {} of \"{}\", didn't exist in previous file",
											layer.n, sCnt, curFile.second);
										thisLayer->orphaned = true;
									}
								}
							}
							
							//then try linking next file, if needed
							if(t.fileNumb != files.size() - 1 && results[curFile.first + 1].first){
								if(results[curFile.first + 1].second.size() <= sCnt)
									log->warn("Current segment {} of \"{}\", doesn't exist in next file", sCnt, curFile.second);
								else{
									auto& nextFileMap = results[curFile.first + 1].second[sCnt];
									auto next = nextFileMap.find(layer.n);
									if(next != nextFileMap.end() && next->second.req == nullptr)
										next->second.req =thisLayer;
									else
										log->warn("Current layer {} in segment {} of \"{}\", doesn't exist in next file",
												  layer.n, sCnt, curFile.second);
								}
							}
						}
						tmpTickets.push(std::move(t));
					}
					sCnt++;
				}
				impCount++;
				if(impCount == files.size())
					moreTicketsExpected = false;
				resultMutex.unlock();

				ticketMutex.lock();
				while (!fileTickets.empty())
				{
					tickets.push_back(std::move(fileTickets.front()));
					fileTickets.pop();
				}
				ticketCV.notify_all();
				ticketMutex.unlock();
			}
			else if (isFitting)
			{
				//fitting code here
				log->debug("Fitting this round");

				//handling for getting ticket in case it is needed(either there was none or there is continuity condition)
				//first wait for tickets to appear if there is none
				if (waitForTicket) {
					log->debug("Waiting for a file to fit to appear");
					//unique lock should lock upon creation
					std::unique_lock<std::mutex> ticketLock(ticketMutex);
					//if before this point new tickets appeared, wait for tickets
					if(tickets.empty())
						ticketCV.wait(ticketLock);
				}
				//getting ticket
				log->debug("Trying to aquire a ticket");
				bool gotTicket{ !cont && !waitForTicket };
				while (!gotTicket) {
					ticketMutex.lock();
					if (cont)
						resultMutex.lock();
					//break if there is no reason for waiting
					if (!moreTicketsExpected&&tickets.empty())
					{
						log->debug("No more tickets expected, braking");
						if (cont)
							resultMutex.unlock();
						ticketMutex.unlock();
						break;
					}

					for (auto it = tickets.begin(); it != tickets.end(); ++it)
						if (!cont || existDep(*it))
						{
							curT = *it;
							tickets.erase(it);
							gotTicket = true;
							break;
						}

					if (cont)
						resultMutex.unlock();
					ticketMutex.unlock();
					
					//TODO check!!
					//if this is cont run wait for something to change
					if(!gotTicket){
						ticketMutex.lock();
						if(!moreTicketsExpected&&tickets.empty())
						{
							ticketMutex.unlock();
							break;
						}
						ticketMutex.unlock();
						std::unique_lock<std::mutex> resLock(resultMutex);
						resultCV.wait(resLock);
					}
				}

				//if no ticket received anyway continue onto next cycle
				if (!gotTicket)
					continue;

				//import code after this point
				
				//calling custom helper
				curT.dropRes->fit_res = fitHelper<T>(fitEngine, curT.pts, curT.max, cont && curT.fileNumb != 0, &curT.dropRes->req->fit_res, ubase);
				log->debug("fit succeeded");
				curT.dropRes->ready = true;
				resultCV.notify_all();
			}
			else
			{
				log->debug("terminating!");
				break;
			}
		}//main while(true) loop
	};
	
	
	std::vector<std::thread> workerArray{};
	for(size_t i = 0; i<tcount; i++)
		workerArray.push_back(std::thread(workerThread, i));
	
	for(auto& thr: workerArray)
		thr.join();
	
	//make map with all the files
	ResList<T> retArray{};
	for(size_t i = 0; i < files.size(); i++)
	{
		auto& file = results[i].second;
		
		retArray.push_back({files[i], {}});
		auto& retFile = retArray.back().second;
		for(auto& seg: file)
		{
			retFile.push_back({});
			auto& cSeg = retFile.back();
			for(auto& layer: seg)
				if(layer.second.ready)
					cSeg[layer.first] = layer.second.fit_res;
		}
	}

	return retArray;
}

//templated printer
template<typename T>
void pout(const ResList<T>& res, std::ostream& str)
{
	for(const auto& file: res)
		for(size_t sCnt = 0; sCnt < file.second.size(); sCnt++){
			str<<'\t'<<"File: \""<<file.first<<"\", segment "<<sCnt<<"\n";
			for(const auto& layer: file.second[sCnt]){
				str<<layer.first<<'\t'<<layer.second<<'\n';
			}
		}
}

int main(int argc, char* argv[])
{
	//initialize main spdlog logger
	auto conLog = spdlog::stdout_color_mt("main");
	auto rLog = spdlog::stdout_color_mt("reader");

	//thread count constant for future paralelism
	size_t threadCnt{ 0 };
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
	//output file name
	std::string ofile{ "" };
	
	//boost program options initialization
	boost::program_options::options_description desc("Supported options");
	desc.add_options()
		("help,h", "display this message")
		("input-file", boost::program_options::value<std::vector<std::string>>(&fileNames)->required(), "input OVF file names")
		("layers,l", boost::program_options::value<std::vector<unsigned int>>(&layers)->required()->multitoken(), "list of layer number to be fit")
		("backend,b", boost::program_options::value<std::string>()->default_value("ceres"), "Backend to use for fitting core-position. Available options are: ceres/mathematica")
		("type", boost::program_options::value<std::string>()->default_value("elliptical", "Fit function type: elliptical/circular"))
		//due to includes from wstp.h fucking up everything
		//thanks MICROSOFT
#if defined(_WIN32)||defined(WIN32)
		("threads,t", boost::program_options::value<size_t>(&threadCnt)->default_value(min(4, std::thread::hardware_concurrency())))
#else
		//'u' in the literal is required because of how std::min is defined
		("threads,t", boost::program_options::value<size_t>(&threadCnt)->default_value(std::min(4u, std::thread::hardware_concurrency())))
#endif
		("output,o", boost::program_options::value<std::string>(&ofile)->default_value(""), "output text file name, by default table is spit into stdout")
		("continuous,c", boost::program_options::bool_switch(&isCont)->default_value(false), "treat files as continuous series to improve fit speed, fit results from previous file is used as initial conditions for next fit")
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
	try {
		boost::program_options::store(
			boost::program_options::command_line_parser(argc, argv).options(desc).positional(posDesc).run(),
			varMap
		);
	} catch(boost::program_options::error& e){
		conLog->error(e.what());
		return -1;
	}

	if (varMap.count("help"))
	{
		//if help is required do it and quit
		std::cout << desc << '\n';
		return 0;
	}
	
	//validate layer list
	std::sort(layers.begin(), layers.end());
	if(std::unique(layers.begin(), layers.end()) != layers.end())
	{
		conLog->error("layer list contains duplicates!");
		return -1;
	}

	//TODO: optional exception handling here 
	//check if required arguments are provided
	try {
		boost::program_options::notify(varMap);
	}
	catch (boost::program_options::error& e) {
		conLog->error("Wrong command-line usage!\n{}", e.what());
		std::cout << desc << '\n';
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
	//has to be done with win api, so it can SUCC
	for(const auto& file: fileNames)
		if (!boost::filesystem::is_regular_file(file))
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

	std::ofstream fout;
	//try to open the output file
	if (ofile != "")
	{
		//open file for writing with modifier to discard old file content before starting a write
		fout.open(ofile, std::fstream::out|std::fstream::trunc);
		if (!fout.good())
		{
			conLog->error("Could not open the file \"{}\" for writing", ofile);
			return -1;
		}
	}

	//output stream of the whole deal
	std::ostream& outs{ (ofile != "") ? fout : std::cout };
	if(ofile != "")outs << "//Backend: " << (bend == Backend::Ceres ? "Ceres" : "Mathematica") << "; Model: " << (ffunc == FitFunction::CircularGauss ? "circular gauss core" : "elliptical gauss core") << '\n';
		
	//set precision
	outs << std::setprecision(10);

	if (ffunc == FitFunction::CircularGauss) {
		pout<CGaussParams<double>>(
			(bend == Backend::Mathematica)?
			worker<CGaussParams<double>, MathematicaFitter<double>>(fileNames, layers, isCont, threadCnt, thresh, cutR, unit_base, isVerb):
			worker<CGaussParams<double>, CeresFitEngine>(fileNames, layers, isCont, threadCnt, thresh, cutR, unit_base, isVerb),
			outs
			);
	}
	else {
		pout<EGaussParams<double>>(
			(bend == Backend::Mathematica) ?
			worker<EGaussParams<double>, MathematicaFitter<double>>(fileNames, layers, isCont, threadCnt, thresh, cutR, unit_base, isVerb) :
			worker<EGaussParams<double>, CeresFitEngine>(fileNames, layers, isCont, threadCnt, thresh, cutR, unit_base, isVerb),
			outs
			);
	}
	return 0;
}
