 
CXX = g++
CC = gcc
EIGEN_INCLUDES = -I/usr/include/eigen3
BUILD_DIR = ./build
CXX_FLAGS = -O2

test: ceres_fit.o OVFReader.o test.o
	cd build;$(CXX) -luuid -lWSTP64i4 -lceres $^ -o test

ceres_fit.o:
	$(CXX) ceres_fit.cpp $(EIGEN_INCLUDES) $(CXX_FLAGS) -c -o $(BUILD_DIR)/ceres_fit.o

OVFReader.o:
	$(CXX) OVFReader.cpp $(CXX_FLAGS) -c -o $(BUILD_DIR)/OVFReader.o

test.o:
	$(CXX) ./test/test.cpp -I./ $(EIGEN_INCLUDES) $(CXX_FLAGS) -c -o $(BUILD_DIR)/test.o

clean:
	rm -f build/*
