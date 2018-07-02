GIT_COMMIT_INFO			= $(shell git log -1 --format='%ci  commit: %h')

CXX				= g++
ICPC				= icpc

CXXFLAGS			= -Wall -std=c++0x -O3 -march=native -D NDEBUG
CXXFLAGS_INTEL			= -Wall -std=c++0x -O3 -march=native -D NDEBUG -static
CXXFLAGS_GIT_COMMIT_INFO	= -D GIT_COMMIT_INFO="\"$(GIT_COMMIT_INFO)\""


all:	cbust

cbust:	*.cpp *.hpp
	$(CXX) $(CXXFLAGS) $(CXXFLAGS_GIT_COMMIT_INFO) -o cbust *.cpp

cbust_intel:	*.cpp *.hpp
	$(ICPC) $(CXXFLAGS_INTEL) $(CXXFLAGS_GIT_COMMIT_INFO) -o cbust_intel *.cpp

clean:
	rm -f cbust cbust_intel
