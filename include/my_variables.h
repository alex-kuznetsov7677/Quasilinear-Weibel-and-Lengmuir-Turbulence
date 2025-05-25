#pragma once

#include <complex>


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
extern double Vstepx;
extern double Vstepy;
extern double Vstepy2;
extern double Vstepykvadr;
extern double Vstepxkvadr;
extern double Vminx;
extern double Vminy;
extern double Vstepx2;




extern double Tmax;				
extern double dt;



extern std::complex<double> I; 

// useful constants
extern double C1;
extern double C2;
extern double C3;
extern double C4;
