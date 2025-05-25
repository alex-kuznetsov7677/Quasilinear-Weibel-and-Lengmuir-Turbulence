#include "my_variables.h"
#include "Distribution_function_changing.h"

#include <vector>
#include <cmath>
#include <complex> 

using std::complex;
using std::vector;

complex<double>DifY(int& N, std::vector<complex<double> >& fka) {
	complex<double> k;

	if (N - setkaBB != 0 && (N - setkaBBminus) % setkaBB != 0) {
		return k = (fka[N + 1] - fka[N - 1]) / Vstepy2;
	}
	else {
		return k = 0;
	}
}
complex<double>DifX(int& N, std::vector<complex<double> >& fka) {
	complex<double> k;

	int fk = (N - N % setkaBB);
	if (fk % setkaBBkvadr != 0 && (fk - setkaBBminus * setkaBB) % setkaBBkvadr != 0) {
		return k = (fka[N + setkaBB] - fka[N - setkaBB]) / Vstepx2;
	}
	else {
		return k = 0;
	}

}
void PERTURBATION_OF_UNIFORM_DISTRIBUTION(int& oblast, int& size, int& Ngarmonik, std::vector<double>& Kvector, std::vector<complex<double> >& IEy, std::vector<complex<double> >& IEx, std::vector<complex<double> >& b1, std::vector<complex<double> >& fk,
	std::vector<complex<double> >& parallel_mas2) {
	double KOEF;
	if (Ngarmonik == Ngarmonik_TS) {
		KOEF = dt * C3/2.;
	}
	else {
		KOEF = dt * C3;
	}

	for (int p2 = 0; p2 < Ngarmonik / size; p2++) {
		for (int j = 0; j < setkaBB; j++) {
			int BBXnomer = j * setkaBB;
			double BBXrl = Vminx + j * Vstepx;

			for (int k = 0; k < setkaBB; k++) {
				int BBXnomerkF0 = BBXnomer + k;
				int BBXnomerk = p2 * setkaBBkvadr + BBXnomer + k;

				double BBYrl = Vminy + k * Vstepy;
				complex<double> dfx = DifX(BBXnomerk, fk);
				complex<double> dfy = DifY(BBXnomerk, fk);
				complex<double> dfa = (BBYrl * dfx - BBXrl * dfy);

				complex<double> yh1 = -I * IEx[p2] * conj(dfx) - I * IEy[p2] * conj(dfy) + b1[p2] * conj(dfa);
				parallel_mas2[BBXnomerkF0] = parallel_mas2[BBXnomerkF0] - KOEF * (yh1 + conj(yh1));
			}

		}
	}
}






void PERTURBATION_OF_MODES_DISTRIBUTION(int& oblast, int& size, int& Ngarmonik,  std::vector<double>& Kvector, std::vector<complex<double> >& IEy, std::vector<complex<double> >& IEx, std::vector<complex<double> >& b, std::vector<complex<double> >& fkcel, std::vector<complex<double> >& fk,
	std::vector<complex<double> >& f0k2cel, std::vector<complex<double> >& parallel_mas1) {

	double KOEF1 = 2 * dt * C3;


	for (int p2 = 0; p2 < Ngarmonik / size; p2++) {
		int preal = p2 + oblast * Ngarmonik / size;

		double B_EXT_start;
		double Ngarmonikx;
		int p1; int p3;
		if (Ngarmonik == Ngarmonik_TS) {
			B_EXT_start = 0;
			Ngarmonikx = Ngarmonik_TS_x;
			p3 = preal % Ngarmonik_TS_y;
			p1 = (preal - p3) / Ngarmonik_TS_y;
		}
		else {
			B_EXT_start = start_level_F_modes;
			Ngarmonikx = Ngarmonik_F_x;
			p3 = preal % Ngarmonik_F_y;
			p1 = (preal - p3) / Ngarmonik_F_y;
		}
		complex<double> BBB =b[p2]+ B_EXT_start;

		for (int j = 0; j < setkaBB; j++) {

			int BBXnomer = j * setkaBB;
			double BBXrl = Vminx + j * Vstepx; 

			for (double k = 0; k < setkaBB; k++) {
				int BBXnomerk = p2 * setkaBBkvadr + BBXnomer + k;
				int BBXnomerkreal = preal * setkaBBkvadr + BBXnomer + k;
				int BBXnomerkF0 = BBXnomer + k;
				double BBYrl = Vminy + k * Vstepy;

				complex<double> df0x = DifX(BBXnomerkF0, f0k2cel);
				complex<double> df0y = DifY(BBXnomerkF0, f0k2cel);
				complex<double> df0a = (BBYrl * df0x - BBXrl * df0y);

				parallel_mas1[BBXnomerk] = fk[BBXnomerk] - dt * I * (BBXrl * Kvector[p1] * fkcel[BBXnomerk] + BBYrl * Kvector[Ngarmonikx + p3] * fkcel[BBXnomerk]) + KOEF1 * (I * IEx[p2] * df0x + I * IEy[p2] * df0y - BBB * df0a);	//основной цикл


			}


		}
	}
}
