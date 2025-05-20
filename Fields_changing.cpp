#include <vector>
#include <cmath>
#include <iostream>

#include "my_variables.h"


complex <double> CURRENT_DENSITY_Y(std::vector<complex<double> >& fka, int p) {

	complex<double> k(0, 0);
	for (int i = 0; i < setkaBB; i++) {
		int BBnomer = p * setkaBBkvadr + i * setkaBB;
		double BBreal = -BBXmax + i * dBBX;
		for (int j = 0; j < setkaBB; j++) {
			double BBYrl = -BBYmax + j * dBBY;
			k = k + C4 * fka[BBnomer + j] * BBYrl * dBBY * dBBX;
		}
	}
	return k;
}
complex <double> CURRENT_DENSITY_X(std::vector<complex<double> >& fka, int p3) {

	complex<double> k(0, 0);
	for (int i = 0; i < setkaBB; i++) {
		int BBnomer = p3 * setkaBBkvadr + i * setkaBB;
		double BBreal = -BBXmax + i * dBBX;
		for (int j = 0; j < setkaBB; j++) {

			k = k + C4 * fka[BBnomer + j] * BBreal * dBBY * dBBX;
		}
	}
	return k;
}

void FILAMENTATION_FIELDS(int& oblast, int& size, int& Ngarmonik,std::vector<double>& Kvector, std::vector<complex<double> >& IEylast, std::vector<complex<double> >& IExlast, std::vector<complex<double> >& blast,
	std::vector<complex<double> >& IEy, std::vector<complex<double> >& IEx, std::vector<complex<double> >& b, std::vector<complex<double> >& IEy1, std::vector<complex<double> >& IEx1, std::vector<complex<double> >& b1,
	std::vector<complex<double> >& IEy1last, std::vector<complex<double> >& IEx1last, std::vector<complex<double> >& b1last, std::vector<complex<double> >& fk, std::vector<complex<double> >& fkcel) {

	for (int p2 = 0; p2 < Ngarmonik_F / size; p2++) {

		int preal = p2 + oblast * Ngarmonik_F / size;
		int p3 = preal % Ngarmonik_F_y;
		int p1 = (preal - p3) / Ngarmonik_F_y;
		int Ng = Ngarmonik_F_x;

		b1[p2] = b1last[p2] + (Kvector[Ng + p3] * IExlast[p2] - Kvector[p1] * IEylast[p2]) * dt;
		complex <double> Integ = CURRENT_DENSITY_X(fk, p2);
		IEx[p2] = IExlast[p2] + dt * (-Kvector[Ng + p3] * b1[p2] - I * (Integ));
		Integ = CURRENT_DENSITY_Y(fk, p2);
		IEy[p2] = IEylast[p2] + dt * (Kvector[p1] * b1[p2] - I * (Integ));


		Integ = CURRENT_DENSITY_Y(fkcel, p2);
		IEy1[p2] = IEy1last[p2] + dt * (Kvector[p1] * blast[p2] - I * (Integ));
		Integ = CURRENT_DENSITY_X(fkcel, p2);
		IEx1[p2] = IEx1last[p2] + dt * (-Kvector[Ng + p3] * blast[p2] - I * (Integ));
		b[p2] = blast[p2] + (Kvector[Ng + p3] * IEx1[p2] - Kvector[p1] * IEy1[p2]) * dt;

	}
}

void TWO_STREAM_FIELDS(int& oblast, int& size, int& Ngarmonik_puc,  std::vector<double>& Kvector_puc, std::vector<complex<double> >& IEylast_puc, std::vector<complex<double> >& IExlast_puc, std::vector<complex<double> >& blast_puc,
	std::vector<complex<double> >& IEy_puc, std::vector<complex<double> >& IEx_puc, std::vector<complex<double> >& b_puc, std::vector<complex<double> >& IEy1_puc, std::vector<complex<double> >& IEx1_puc, std::vector<complex<double> >& b1_puc,
	std::vector<complex<double> >& IEy1last_puc, std::vector<complex<double> >& IEx1last_puc, std::vector<complex<double> >& b1last_puc, std::vector<complex<double> >& fk_puc, std::vector<complex<double> >& fkcel_puc) {

	for (int p2 = 0; p2 < Ngarmonik_TS / size; p2++) {

		int preal = p2 + oblast * Ngarmonik_TS / size;
		int p3 = preal % Ngarmonik_TS_y;
		int p1 = (preal - p3) / Ngarmonik_TS_y;
		int Ng = Ngarmonik_TS_x;

		b1_puc[p2] = b1last_puc[p2] + (Kvector_puc[Ng + p3] * IExlast_puc[p2] - Kvector_puc[p1] * IEylast_puc[p2]) * dt;
		complex <double> Integ = CURRENT_DENSITY_Y(fk_puc, p2);
		IEy_puc[p2] = IEylast_puc[p2] + dt * (Kvector_puc[p1] * b1_puc[p2] - I * (Integ));
		Integ = CURRENT_DENSITY_X(fk_puc, p2);
		IEx_puc[p2] = IExlast_puc[p2] + dt * (-Kvector_puc[Ng + p3] * b1_puc[p2] - I * (Integ));

		IEy1_puc[p2] = (IEy_puc[p2] + IEylast_puc[p2]) / 2.;
		IEx1_puc[p2] = (IEx_puc[p2] + IExlast_puc[p2]) / 2.;
		b_puc[p2] = blast_puc[p2] + (Kvector_puc[Ng + p3] * IEx1_puc[p2] - Kvector_puc[p1] * IEy1_puc[p2]) * dt;
	}
}

