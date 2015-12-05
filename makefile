CXX = g++
CXXFLAGS = -Wall

ugh: flow.cc
	$(CXX) $(CXXFLAGS) -o $@ flow.cc

clean:
	rm -f ugh flow.cc

