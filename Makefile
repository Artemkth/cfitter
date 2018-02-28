 
CXX = g++
CC = gcc
EIGEN_INCLUDES = -I/usr/include/eigen3
BOOST_LIBS = -lboost_system -lboost_program_options -lboost_filesystem
BUILD_DIR = ./build
CXX_FLAGS = -O2

test: ceres_fit.o OVFReader.o test.o
	cd build;$(CXX) -luuid -lWSTP64i4 -lceres $^ -o test
	
cfitter: ceres_fit.o OVFReader.o CFitter.o
	cd build;$(CXX) -luuid -lWSTP64i4 -lceres -lpthread $(BOOST_LIBS) $^ -o cfitter

ceres_fit.o:
	$(CXX) ceres_fit.cpp $(EIGEN_INCLUDES) $(CXX_FLAGS) -c -o $(BUILD_DIR)/ceres_fit.o

OVFReader.o:
	$(CXX) OVFReader.cpp $(CXX_FLAGS) -c -o $(BUILD_DIR)/OVFReader.o

test.o:
	$(CXX) ./test/test.cpp -I./ $(EIGEN_INCLUDES) $(CXX_FLAGS) -c -o $(BUILD_DIR)/test.o
	
CFitter.o:
	$(CXX) CFitter.cpp $(EIGEN_INCLUDES) $(CXX_FLAGS) -c -o $(BUILD_DIR)/CFitter.o

clean:
	rm -f build/*
