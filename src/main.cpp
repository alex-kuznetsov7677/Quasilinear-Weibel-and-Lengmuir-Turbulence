#include <cmath>
#include <complex>
#include <vector>
#include <cstdlib>
#include <random>
#include <iostream>
#include <algorithm> 

#include <mpi.h>
// 4. Локальные заголовки проекта
#include "my_variables.h" 
#include "Preparation.h" 
#include "Fields_changing.h" 
#include "Distribution_function_changing.h" 
#include "Integrals.h" 
#include "ConfigReader.h"
#include "OutputWriter.h"

using namespace std;


auto config = ReadConfig("config.txt");

//start level of modes
double start_level_F_modes = config["start_level_F_modes"];
double start_level_TS_modes =config["start_level_TS_modes"];



//numerous of modes
int Ngarmonik_TS_x = static_cast<int>(config["Ngarmonik_TS_x"]);
int Ngarmonik_TS_y = static_cast<int>(config["Ngarmonik_TS_y"]);
int Ngarmonik_TS = Ngarmonik_TS_x * Ngarmonik_TS_y;


int Ngarmonik_F_x = static_cast<int>(config["Ngarmonik_F_x"]);
int Ngarmonik_F_y = static_cast<int>(config["Ngarmonik_F_y"]);
int Ngarmonik_F = Ngarmonik_F_x * Ngarmonik_F_y;

// parameters of the distribution function

double V_stream = config["V_stream"];
double Beta_bkg = config["Beta_bkg"];
double Beta_stream = config["Beta_stream"];
double ratio = config["ratio"];
double A = (ratio * (pow(V_stream, 2) + pow(Beta_stream, 2) / 2.) + (pow(V_stream * ratio, 2) + pow(Beta_bkg, 2) / 2.)) / (pow(Beta_bkg, 2) / 2.) - 1;

//grid parameters in the velocity space
double BetaperpsqrtA = Beta_bkg * sqrt(1 + A);
int setkaBB = static_cast<int>(config["setkaBB"]);
int setkaBBkvadr = setkaBB * setkaBB;
int setkaBBminus = setkaBB - 1;


double Vstepx =  config["ratio_Vstepx_multiple_setkaBB_to_Beta_bkg"] * Beta_bkg / setkaBB;
double Vstepy = config["ratio_Vstepy_to_Vstepx"] * Vstepx;
double Vstepy2 = Vstepy * 2;
double Vstepykvadr = Vstepy * Vstepy;
double Vstepxkvadr = Vstepx * Vstepx;
double Vminx = Beta_bkg * config["ratio_Vminx_to_Beta_bkg"] + (Vstepx / 2.);
double Vminy = Beta_bkg * config["ratio_Vminy_to_Beta_bkg"] + (Vstepy / 2.);
double Vstepx2 = Vstepx * 2;

double Tmax = config["Tmax"];
double dt = config["dt"];

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


	//double trimax = 0;
	double max_with_stream = 0;

	for (double i = 0; i < setkaBB; i++) {
		double BBnomer = i * setkaBB;
		double BBXrl = Vminx + i * Vstepx;
		for (double j = 0; j < setkaBB; j++) {
			double BBYrl = Vminy + j * Vstepy;
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
	OutputFiles out;
	if (rank==0){
		out.log_number_MPI_processes(size);
		InitializeOutputFiles(out,config);

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
	if (rank == 0) {
		starttime = MPI_Wtime();
	}



	for (int pk = 2; pk < Tmax / dt; pk++) {

		if (pk == 100) {
			start_level_F_modes = 0;
		}


		f0kcel = f0k2cel;

		fill(parallel_mas2.begin(), parallel_mas2.end(), complex<double>(0., 0.));

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
		TWO_STREAM_FIELDS(rank, size, Ngarmonik_TS, Kvector_TS, IEylast_TS, IExlast_TS, blast_TS, IEy_TS, IEx_TS, b_TS, IEy1_TS, IEx1_TS, b1_TS, IEy1last_TS, IEx1last_TS, b1last_TS, fk_TS, fkcel_TS);

		if (pk % 10 == 0) {
			MPI_Gather(b_F.data(), Ngarmonik_F / size, MPI_DOUBLE_COMPLEX, B.data(), Ngarmonik_F / size, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
			MPI_Gather(b1_F.data(), Ngarmonik_F / size, MPI_DOUBLE_COMPLEX, B1.data(), Ngarmonik_F / size, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);



			MPI_Gather(IEx_TS.data(), Ngarmonik_TS / size, MPI_DOUBLE_COMPLEX, IEX_TS.data(), Ngarmonik_TS / size, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
			MPI_Gather(IEy_TS.data(), Ngarmonik_TS / size, MPI_DOUBLE_COMPLEX, IEY_TS.data(), Ngarmonik_TS / size, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
			MPI_Gather(b_TS.data(), Ngarmonik_TS / size, MPI_DOUBLE_COMPLEX, B_TS.data(), Ngarmonik_TS / size, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
		}

		
		if (rank == 0 && pk % 10 == 0) {
			double time = pk*dt;
            		out.log_fields(time, B, IEX_TS, IEY_TS, B_TS, Ngarmonik_F, Ngarmonik_TS,f0k2cel, NU_effective, Kprod);
		}
		if (rank == 0 && pk % 400 == 0) {
            		out.log_distribution(setkaBB, f0k);
		}



		f0k = f0k2;
		fill(parallel_mas2.begin(), parallel_mas2.end(), complex<double>(0., 0.));
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

	void close_all();
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0) {
		endtime = MPI_Wtime();
		double time_ex=endtime - starttime;
		out.log_execution_time(time_ex);
	}
	MPI_Finalize();



}