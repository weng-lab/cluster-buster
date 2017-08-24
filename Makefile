GIT_COMMIT_INFO = $(shell git log -1 --format='%ci  commit: %h')

CXX		= g++
CXXFLAGS	= -Wall -O3 -march=native -D NDEBUG -D GIT_COMMIT_INFO="\"$(GIT_COMMIT_INFO)\""


all :	cbust

cbust :	*.cpp *.hpp
	$(CXX) $(CXXFLAGS) -o cbust *.cpp

clean :
	rm -f cbust
