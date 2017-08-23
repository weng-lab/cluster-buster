CXX		= g++
CXXFLAGS	= -Wall -O3 -DNDEBUG


all :	cbust

cbust :	*.cpp *.hpp
	$(CXX) $(CXXFLAGS) -o cbust *.cpp

clean :
	rm -f cbust
