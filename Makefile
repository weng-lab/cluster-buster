GIT_COMMIT_INFO			= $(shell git log -1 --format='%ci  commit: %h')

CXX				= g++
ICPC				= icpc

CXXFLAGS			= -Wall -std=c++0x -O3 -march=x86-64 -D NDEBUG
CXXFLAGS_INTEL			= -Wall -std=c++0x -O3 -march=native -D NDEBUG -static
CXXFLAGS_GIT_COMMIT_INFO	= -D GIT_COMMIT_INFO="\"$(GIT_COMMIT_INFO)\""

# Optionally compile with AMD Math Library (LibM) for a slightly faster binary:
#   https://developer.amd.com/amd-aocl/amd-math-library-libm/
#
#     mkdir ../aocl-libm
#     cd ../aocl-libm
#
#     AMD_LIB_VERSION=3.1.0
#
#     tar xzf aocl-libm-linux-aocc-${AMD_LIB_VERSION}.tar.gz
#     mv amd-libm amd-libm-aocc
#
#     tar xzf aocl-libm-linux-gcc-${AMD_LIB_VERSION}.tar.gz
#     mv amd-libm amd-libm-gcc
AMD_LIBM_AOCC			= ../aocl-libm/amd-libm-aocc
AMD_LIBM_GCC			= ../aocl-libm/amd-libm-gcc


all:	cbust

cbust:	*.cpp *.hpp
	$(CXX) $(CXXFLAGS) $(CXXFLAGS_GIT_COMMIT_INFO) -o cbust *.cpp

cbust_static:	*.cpp *.hpp
	$(CXX) $(CXXFLAGS) $(CXXFLAGS_GIT_COMMIT_INFO) -o cbust_static *.cpp -static

cbust_amd_libm_aocc:	*.cpp *.hpp
	$(CXX) $(CXXFLAGS) $(CXXFLAGS_GIT_COMMIT_INFO) -I $(AMD_LIBM_AOCC)/include -o cbust_amd_libm_aocc *.cpp -L $(AMD_LIBM_AOCC)/lib -static -lalm -lm

cbust_amd_libm_gcc:	*.cpp *.hpp
	$(CXX) $(CXXFLAGS) $(CXXFLAGS_GIT_COMMIT_INFO) -I $(AMD_LIBM_GCC)/include -o cbust_amd_libm_gcc *.cpp -L $(AMD_LIBM_GCC)/lib -static -lalm -lm

cbust_intel:	*.cpp *.hpp
	$(ICPC) $(CXXFLAGS_INTEL) $(CXXFLAGS_GIT_COMMIT_INFO) -o cbust_intel *.cpp

clean:
	rm -f cbust cbust_static cbust_amd_libm_aocc cbust_amd_libm_gcc cbust_intel
