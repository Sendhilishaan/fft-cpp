CXX = g++
CXXFLAGS = -std=c++11 -O3 -Wall

all: fft bluestein

fft: fft.cpp
	$(CXX) $(CXXFLAGS) -o fft fft.cpp

bluestein: bluestein.cpp
	$(CXX) $(CXXFLAGS) -o bluestein bluestein.cpp

clean:
	rm -f fft bluestein *.o