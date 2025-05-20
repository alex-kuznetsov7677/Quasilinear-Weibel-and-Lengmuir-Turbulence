#pragma once
#include<fstream>
#include <vector>
#include<sstream> 
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <complex>
#include <omp.h>
#include <time.h>
#include <stdio.h>
#include "mpi.h"
#include <algorithm> 
#include <random>

using namespace std;

extern const int NUM_NODES;       // Количество узлов
extern const int NUM_THREADS;     // Количество потоков на каждый процесс


//start level of modes
extern double start_level_F_modes;
extern double start_level_TS_modes;
//numerous of modes
extern int Ngarmonik_TS_y;
extern int Ngarmonik_TS_x;
extern int Ngarmonik_TS;

extern int Ngarmonik_F_x;
extern int Ngarmonik_F_y;
extern int Ngarmonik_F;

// parameters of the distribution function
extern double V_stream;
extern double Beta_bkg;
extern double Beta_stream;
extern double ratio;
extern double A;

//grid parameters in the velocity space
extern double BetaperpsqrtA;
extern int setkaBB;
extern int setkaBBkvadr;
extern int setkaBBminus;
extern double dBBX;
extern double dBBY;
extern double dBBY2;
extern double dBBYkvadr;
extern double dBBXkvadr;
extern double BBXmax;
extern double BBYmax;
extern double dBBX2;


extern double Tmax;				
extern double dt;



extern complex<double> I; 

// useful constants
extern double C1;
extern double C2;
extern double C3;
extern double C4;
