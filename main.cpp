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

#include "my_variables.h" 
#include "Preparation.h" 
#include "Fields_changing.h" 
#include "Distribution_function_changing.h" 
#include "Integrals.h" 
using namespace std;

//
const int NUM_NODES = 20;       // Количество узлов
const int NUM_THREADS = 12;     // Количество потоков на каждый процесс

//start level of modes
double start_level_F_modes = 0.0045;
double start_level_TS_modes = 0.00003;

//numerous of modes
int Ngarmonik_TS_y = 60;
int Ngarmonik_TS_x = 10;
int Ngarmonik_TS = Ngarmonik_TS_y * Ngarmonik_TS_x;

int Ngarmonik_F_x = 30;
int Ngarmonik_F_y = 10;
int Ngarmonik_F = Ngarmonik_F_x * Ngarmonik_F_y;

// parameters of the distribution function
double V_stream = 0.1 * 4.5 / sqrt(2);
double Beta_bkg = 0.1;
double Beta_stream = 0.1;
double ratio = 0.01;
double A = (ratio * (pow(V_stream, 2) + pow(Beta_stream, 2) / 2.) + (pow(V_stream * ratio, 2) + pow(Beta_bkg, 2) / 2.)) / (pow(Beta_bkg, 2) / 2.) - 1;

//grid parameters in the velocity space
double BetaperpsqrtA = Beta_bkg * sqrt(1 + A);
int setkaBB = 300;
int setkaBBkvadr = setkaBB * setkaBB;
int setkaBBminus = setkaBB - 1;
double dBBX = 8 * Beta_bkg / setkaBB;
double dBBY = 1.2 * dBBX;
double dBBY2 = dBBY * 2;
double dBBYkvadr = dBBY * dBBY;
double dBBXkvadr = dBBX * dBBX;
double BBXmax = Beta_bkg * 4 - (dBBX / 2.);
double BBYmax = BBXmax;
double dBBX2 = dBBX * 2;

double Tmax = 4000;
double dt = 0.2;

complex<double> I(0., 1.);

double C1 = 2. * 2. / ((1. + A) * pow(Beta_bkg, 3.) * pow(2. * acos(-1.0), 3. / 2.));
double C2 = Beta_bkg * (1. + A) * Beta_bkg / (8. * acos(-1.0));
double C3 = Beta_bkg * sqrt(1. + A) / (2 * pow(2., 3. / 2.) * sqrt(acos(-1.0)));
double C4 = 2. * sqrt(2.) * acos(-1.0) / (sqrt(acos(-1.0)) * Beta_bkg * sqrt(1. + A));



int main(int argc, char** argv)
{
	double endtime;

	double starttime;


	vector <double> Kvector_F;
	vector <double> Kvector_TS;
	Creating_grid_wavevectors(Kvector_F, Kvector_TS);


	vector <complex<double> > BEXTvector;
	vector <complex<double> > IEyEXT;
	vector <complex<double> > IExEXT;

	vector <complex<double> > parallel_mas1_F;
	vector <complex<double> > parallel_mas1_TS;
	vector <complex<double> > parallel_mas2;
	vector <complex<double> > parallel_mas3;

	vector <complex<double> > B;
	vector <complex<double> > IEX_TS;
	vector <complex<double> > IEY_TS;
	vector <complex<double> > B_TS;
	vector <complex<double> > B1;


	vector <complex<double> > b_F;	vector <complex<double> > b1_F; 	vector <complex<double> > blast_F;	vector <complex<double> > b1last_F;
	vector <complex<double> > IEx_F; 	vector <complex<double> > IEx1_F;  vector <complex<double> > IExlast_F; 	vector <complex<double> > IEx1last_F; 	
	vector <complex<double> > IEy_F; 	vector <complex<double> > IEy1_F;  vector <complex<double> > IEylast_F;		vector <complex<double> > IEy1last_F; 
	vector <complex<double> > IEx_TS;	vector <complex<double> > IEx1_TS;	vector <complex<double> > IExlast_TS;	vector <complex<double> > IEx1last_TS;
	vector <complex<double> > IEy_TS;	vector <complex<double> > IEy1_TS;	vector <complex<double> > IEylast_TS;	vector <complex<double> > IEy1last_TS;
	vector <complex<double> > b_TS;	vector <complex<double> > b1_TS;	vector <complex<double> > blast_TS;	vector <complex<double> > b1last_TS;

	vector <complex<double> > fk_F; 	vector <complex<double> > fkcel_F;
	vector <complex<double> > fk_TS;	vector <complex<double> > fkcel_TS;


	vector <complex<double> > f0k;	vector <complex<double> > f0k2;	vector <complex<double> > f0kcel;	vector <complex<double> > f0k2cel;


	ofstream CMFt100000("DESCRIPTIONRe[B]__365.txt");
	CMFt100000 << "ratio=0.01; V_stream=0.1*4.5/sqrt(2); T_bkg=0.1; Beta_stream=0.1;dt=0.2(cout100*0.2),Beta_perp=0.1,30x10garmoniki.start_without_pause_ravnyeNU(0.12 * pow(2., i * 0.06))(K=-0.06+ 0.012 * i))фазанорм. начальная каппа=infty ФР, setka == 300, BBYmax = BBXmax;double BBXmax = Beta_perp * 4;dBBY=1.2*dBBX" << '\n';
	CMFt100000 << '\n' << '\n' << "1D_PUCHOK_setka30x10=Kyvector.push_back(2.5 + 0.05* i)Kxvector_puchok.push_back(-0.5.+0.05+0.1*i)";
	CMFt100000 << '\n' << "(bez bz(ky=0)not impact on two-stream modes)NU=0";
	CMFt100000 << '\n' << "без уждвоенной моды";
	CMFt100000 << '\n' << "1p.version"; cout << '\n' << "1p.version";
	ofstream CMFt2000("abs[B]___365MATLAB.txt");
	ofstream CMFt2003("abs[Bts]___365MATLAB.txt");
	ofstream CMFt2001("abs[Ex]___365MATLAB.txt");
	ofstream CMFt2002("abs[Ey]___365MATLAB.txt");
	ofstream CMFt1301("Wy0__365.txt");
	ofstream CMFt1303("Wx0__365.txt");
	ofstream CMFt1302("A__365.txt");
	ofstream CMFt30("sr_b___365.txt");
	ofstream CMFt31("sr_bts___365.txt");
	ofstream CMFt13045("NU_eff__365.txt");
	ofstream CMFt32("Re[Time]__365.txt");
	ofstream CMFt1("sr_ey___365.txt");
	ofstream CMFt2("sr_ex___365.txt");
	ofstream CMFt_108("ReFUNC0_365.txt");
	ofstream CMFt_109("ImFUNC111_30000_157____631_365.txt");
	ofstream CMFt_110("ImFUNC111_30000_47___631_365.txt");
	ofstream CMFt_161("ImFUNC111_160000_47____631_365.txt");



	//double trimax = 0;
	double max_with_stream = 0;

	for (double i = 0; i < setkaBB; i++) {
		double BBnomer = i * setkaBB;
		double BBXrl = -BBXmax + i * dBBX;
		for (double j = 0; j < setkaBB; j++) {
			double BBYrl = -BBYmax + j * dBBY;
			max_with_stream = 1 / (acos(-1.0) * Beta_bkg * Beta_bkg) * exp(-pow(BBXrl / Beta_bkg, 2.) - pow((BBYrl + V_stream * ratio) / Beta_bkg, 2.)) + ratio / (acos(-1.0) * Beta_stream * Beta_stream) * exp(-pow((BBYrl - V_stream) / Beta_stream, 2.) - pow(BBXrl / Beta_stream, 2.));
			f0kcel.push_back(max_with_stream);
			f0k2cel.push_back(max_with_stream);
			f0k.push_back(max_with_stream);
			f0k2.push_back(max_with_stream);

		}
	}
	double VVxzero = real(ENERGY_X(f0k2));
	double VVyzero = real(ENERGY_Y(f0k2));
	double Konc1 = real(DENSITY(f0k2, 0));


	for (double i = 0; i < setkaBB; i++) {
		double BBnomer = i * setkaBB;
		for (double j = 0; j < setkaBB; j++) {
			f0kcel[BBnomer + j] = f0kcel[BBnomer + j] / Konc1;
			f0k2cel[BBnomer + j] = f0k2cel[BBnomer + j] / Konc1;
			f0k[BBnomer + j] = f0k[BBnomer + j] / Konc1;
			f0k2[BBnomer + j] = f0k2[BBnomer + j] / Konc1;
		}
	}

	CMFt100000 << VVyzero << '\n';
	CMFt100000 << VVxzero << '\n';

	for (int i = 0; i < Ngarmonik_F; i++) {
		B.push_back((0., 0.));
		B1.push_back((0., 0.));
	}
	for (int i = 0; i < Ngarmonik_TS; i++) {
		IEX_TS.push_back((0., 0.));
		IEY_TS.push_back((0., 0.));
		B_TS.push_back((0., 0.));

	}

	int rank, size;
	MPI_Init(&argc, &argv);
	starttime = MPI_Wtime();
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (Ngarmonik_F % size != 0)
	{
		printf("This application is meant to be run with Ngarmonik/2 MPI processes.\n");
		MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
	}
	double B_MOD; double DBprod_MOD; double Kprod; double NU_effective; double sqrt_VVyzero;
	for (int i = 0; i < Ngarmonik_F / size; i++) {
		b_F.push_back((0., 0.));		b1_F.push_back((0., 0.)); 		blast_F.push_back((0., 0.));		b1last_F.push_back((0., 0.));
		IEx_F.push_back((0., 0.));		IEx1_F.push_back((0., 0.));		IExlast_F.push_back((0., 0.));		IEx1last_F.push_back((0., 0.));
		IEy_F.push_back((0., 0.));		IEy1_F.push_back((0., 0.));		IEylast_F.push_back((0., 0.));		IEy1last_F.push_back((0., 0.));
	}

	int start_ampl = -2;
	int end_ampl = 2;

	int start_phase = -2;
	int end_phase = 2;

	for (int p2 = 0; p2 < Ngarmonik_TS / size; p2++) {

		int preal = p2 + rank * Ngarmonik_TS / size;
		int p3 = preal % Ngarmonik_TS_y;
		int p1 = (preal - p3) / Ngarmonik_TS_y;

		int x = rand() % (end_ampl - start_ampl + 1) + start_ampl;
		int y = rand() % (end_phase - start_phase + 1) + start_phase;
		complex< double > z(x, y);
		complex< double > start_condition = start_level_TS_modes * z / Kvector_TS[Ngarmonik_TS_x + p3] * Kvector_TS[p1];
		IEx_TS.push_back(start_condition);
		IEx1_TS.push_back(start_condition);
		IExlast_TS.push_back(start_condition);
		IEx1last_TS.push_back(start_condition);

		start_condition = 0.00003 * z;
		IEy_TS.push_back(start_condition);
		IEy1_TS.push_back(start_condition);
		IEylast_TS.push_back(start_condition);
		IEy1last_TS.push_back(start_condition);

		b_TS.push_back(0);
		b1_TS.push_back(0);
		blast_TS.push_back(0);
		b1last_TS.push_back(0);

	}


	for (int i = 0; i < Ngarmonik_F / size; i++) {
		complex< double > z(exp(0));
		BEXTvector.push_back(start_level_F_modes * z);
	}

	for (int i = 0; i < setkaBBkvadr * Ngarmonik_F / size; i++) {
		parallel_mas1_F.push_back((0., 0.));
	}
	for (int i = 0; i < setkaBBkvadr * Ngarmonik_TS / size; i++) {
		parallel_mas1_TS.push_back((0., 0.));
	}

	for (int i = 0; i < setkaBBkvadr * Ngarmonik_F / size; i++) {
		fk_F.push_back((0., 0.));
		fkcel_F.push_back((0., 0.));

	}
	for (int i = 0; i < setkaBBkvadr * Ngarmonik_TS / size; i++) {
		fk_TS.push_back((0., 0.));
		fkcel_TS.push_back((0., 0.));

	}

	for (int i = 0; i < setkaBBkvadr; i++) {
		parallel_mas2.push_back((0., 0.));
		parallel_mas3.push_back((0., 0.));
	}
	for (int i = 2; i < Tmax / dt; i++) {

		if (i == 100) {
			start_level_F_modes = 0;
		}

		f0kcel = f0k2cel;
		fill(parallel_mas2.begin(), parallel_mas2.end(), 0);
		PERTURBATION_OF_UNIFORM_DISTRIBUTION(rank, size, Ngarmonik_F, Kvector_F, IEy1last_F, IEx1last_F, b1last_F, fk_F, parallel_mas2);
		MPI_Reduce(parallel_mas2.data(), parallel_mas3.data(), setkaBBkvadr, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);
		if (rank == 0) {
			transform(f0kcel.begin(), f0kcel.end(), parallel_mas3.begin(), f0k2cel.begin(), std::plus<complex<double> >());
		}

		fill(parallel_mas2.begin(), parallel_mas2.end(), 0);
		PERTURBATION_OF_UNIFORM_DISTRIBUTION(rank, size, Ngarmonik_TS, Kvector_TS, IEy1last_TS, IEx1last_TS, b1last_TS, fk_TS, parallel_mas2);
		MPI_Reduce(parallel_mas2.data(), parallel_mas3.data(), setkaBBkvadr, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);
		if (rank == 0) {
			transform(f0k2cel.begin(), f0k2cel.end(), parallel_mas3.begin(), f0k2cel.begin(), std::plus<complex<double> >());
		}


		MPI_Bcast(f0k2cel.data(), setkaBBkvadr, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
		PERTURBATION_OF_MODES_DISTRIBUTION(rank, size, Ngarmonik_F, Kvector_F, IEylast_F, IExlast_F, blast_F, fkcel_F, fk_F, f0k2cel,  parallel_mas1_F);
		fk_F = parallel_mas1_F;
		PERTURBATION_OF_MODES_DISTRIBUTION(rank, size, Ngarmonik_TS, Kvector_TS, IEylast_TS, IExlast_TS, blast_TS, fkcel_TS, fk_TS, f0k2cel, parallel_mas1_TS);
		fk_TS = parallel_mas1_TS;





		FILAMENTATION_FIELDS(rank, size, Ngarmonik_F, Kvector_F, IEylast_F, IExlast_F, blast_F, IEy_F, IEx_F, b_F, IEy1_F, IEx1_F, b1_F, IEy1last_F, IEx1last_F, b1last_F, fk_F, fkcel_F);
		MPI_Gather(b_F.data(), Ngarmonik_F / size, MPI_DOUBLE_COMPLEX, B.data(), Ngarmonik_F / size, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
		MPI_Gather(b1_F.data(), Ngarmonik_F / size, MPI_DOUBLE_COMPLEX, B1.data(), Ngarmonik_F / size, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);


		TWO_STREAM_FIELDS(rank, size, Ngarmonik_TS, Kvector_TS, IEylast_TS, IExlast_TS, blast_TS, IEy_TS, IEx_TS, b_TS, IEy1_TS, IEx1_TS, b1_TS, IEy1last_TS, IEx1last_TS, b1last_TS, fk_TS, fkcel_TS);
		MPI_Gather(IEx_TS.data(), Ngarmonik_TS / size, MPI_DOUBLE_COMPLEX, IEX_TS.data(), Ngarmonik_TS / size, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
		MPI_Gather(IEy_TS.data(), Ngarmonik_TS / size, MPI_DOUBLE_COMPLEX, IEY_TS.data(), Ngarmonik_TS / size, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
		MPI_Gather(b_TS.data(), Ngarmonik_TS / size, MPI_DOUBLE_COMPLEX, B_TS.data(), Ngarmonik_TS / size, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

		int pk = i;
		if (rank == 0) {

			if (pk % 100 == 0) {
				cout << pk * dt << '\t';
				double Bmodkvadr = 0.;
				for (int it = 0; it < Ngarmonik_F; it++) {

					Bmodkvadr = Bmodkvadr + pow(abs(B[it]), 2);
					
				}
				double IExkvadr = 0.;  double IEykvadr = 0.;  double Bkvadr = 0.;
				for (double j = 0; j < Ngarmonik_TS; j++) {
					IExkvadr = IExkvadr + pow(abs(IEX_TS[j]), 2);
					IEykvadr = IEykvadr + pow(abs(IEY_TS[j]), 2);
					Bkvadr = Bkvadr + pow(abs(B_TS[j]), 2);

				}

				double Bmod = sqrt(Bmodkvadr);
				double Bmodts = sqrt(Bkvadr);
				double VVx0 = real(ENERGY_X(f0k2cel));
				double VVy0 = real(ENERGY_Y(f0k2cel));
				double Areal = (VVy0) / (VVx0)-1;
				double IExmod = sqrt(IExkvadr);
				double IEymod = sqrt(IEykvadr);
				cout << real(DENSITY(f0k2cel, 0)) << '\t' << Areal << '\t' << Bmod << '\t' << IExmod << '\t' << IEymod << '\t' << Bmodts << '\t' << NU_effective << '\t' << Kprod << '\t' << '\n';

				for (double p = 0; p < Ngarmonik_F; p++) {
					CMFt2000 << abs(B[p]) << " ";
				}
				for (double p = 0; p < Ngarmonik_TS; p++) {
					CMFt2002 << abs(IEY_TS[p]) << " ";
					CMFt2001 << abs(IEX_TS[p]) << " ";
					CMFt2003 << abs(B_TS[p]) << " ";
				}

				CMFt2000 << ";"; CMFt2001 << ";"; CMFt2002 << ";";
				CMFt30 << Bmod << ",";
				CMFt31 << Bmodts << ",";
				CMFt1 << IEymod << ",";
				CMFt2 << IExmod << ",";
				CMFt1302 << Areal << ",";
				CMFt1301 << VVy0 << ",";
				CMFt1303 << VVx0 << ",";
				CMFt13045 << NU_effective << ",";
			}
		}

		if (rank == 0) {
			if (pk % 400 == 0) {
				for (double j = 0; j < setkaBB; j++) {
					double BBXnomer = j * setkaBB;
					for (double k = 0; k < setkaBB; k++) {
						double BBXnomerk = BBXnomer + k;
						CMFt_108 << real(f0k[BBXnomerk]) << " ";
					}

				}
				CMFt_108 << ";";

			}
		}


		f0k = f0k2;
		fill(parallel_mas2.begin(), parallel_mas2.end(), 0);
		PERTURBATION_OF_UNIFORM_DISTRIBUTION(rank, size, Ngarmonik_F, Kvector_F, IEylast_F, IExlast_F, blast_F, fkcel_F, parallel_mas2);
		MPI_Reduce(parallel_mas2.data(), parallel_mas3.data(), setkaBBkvadr, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);
		if (rank == 0) {
			transform(f0k.begin(), f0k.end(), parallel_mas3.begin(), f0k2.begin(), std::plus<complex<double> >());
		}
		fill(parallel_mas2.begin(), parallel_mas2.end(), 0);
		PERTURBATION_OF_UNIFORM_DISTRIBUTION(rank, size, Ngarmonik_TS, Kvector_TS, IEylast_TS, IExlast_TS, blast_TS, fkcel_TS, parallel_mas2);
		MPI_Reduce(parallel_mas2.data(), parallel_mas3.data(), setkaBBkvadr, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);
		if (rank == 0) {
			transform(f0k2.begin(), f0k2.end(), parallel_mas3.begin(), f0k2.begin(), std::plus<complex<double> >());
		}

		
		MPI_Bcast(f0k2.data(), setkaBBkvadr, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

		PERTURBATION_OF_MODES_DISTRIBUTION(rank, size, Ngarmonik_F, Kvector_F, IEy1_F, IEx1_F, b1_F, fkcel_F, fk_F, f0k2, parallel_mas1_F);
		fkcel_F = parallel_mas1_F;
		PERTURBATION_OF_MODES_DISTRIBUTION(rank, size, Ngarmonik_TS, Kvector_TS, IEy1_TS, IEx1_TS, b1_TS, fk_TS, fkcel_TS, f0k2,parallel_mas1_TS);
		fkcel_TS = parallel_mas1_TS;


		b1last_F = b1_F;
		IEy1last_F = IEy1_F;
		blast_F = b_F;
		IEylast_F = IEy_F;
		IExlast_F = IEx_F;
		IEx1last_F = IEx1_F;

		IEy1last_TS = IEy1_TS;
		IEylast_TS = IEy_TS;
		IExlast_TS = IEx_TS;
		IEx1last_TS = IEx1_TS;
		blast_TS = b_TS;
		b1last_TS = b1_TS;

	}

	CMFt100000.close();
	CMFt2000.close();
	CMFt30.close();
	CMFt31.close();
	CMFt1301.close();
	CMFt1303.close();
	CMFt1302.close();
	CMFt13045.close();
	CMFt2001.close();
	CMFt2002.close();
	CMFt2003.close();

	CMFt1.close();
	CMFt2.close();
	CMFt_108.close();

	if (rank == 0) {
		endtime = MPI_Wtime();
		cout << "final" << endtime - starttime << '\n';
		CMFt32 << endtime - starttime << '\n';
	}
	CMFt32.close();
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Finalize();



}