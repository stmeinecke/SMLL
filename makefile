all:
	g++ -std=c++11 SML_Laser.cpp -o SML_Laser -Ofast -I $$PWD/incl/ -lfftw3 -lm
