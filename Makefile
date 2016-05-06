all:
	g++ -Wall -O3 -o MaxwellBloch MaxwellBloch.cc -I/usr/local/boost_1_58_0 -std=c++11 -fopenmp
