#include <vector>
#include <cmath>
#include <iostream>

#include "my_variables.h"


std::complex <double> DENSITY(const std::vector<std::complex<double> >& fka, int p3) {


	complex<double> k(0, 0);

	for (double i = 0; i < setkaBB; i++) {
		double BBnomer = p3 * setkaBBkvadr + i * setkaBB;
		for (double j = 0; j < setkaBB; j++) {
			k = k + fka[BBnomer + j] * Vstepy * Vstepx;
		}
	}
	return k;
}
std::complex <double> ENERGY_Y(const std::vector<std::complex<double> >& fka) {

	complex<double> k(0, 0);
	for (int i = 0; i < setkaBB; i++) {
		double BBnomer = i * setkaBB;
		for (int j = 0; j < setkaBB; j++) {
			double BBYrl = Vminy + j * Vstepy;
			k = k + fka[BBnomer + j] * BBYrl * BBYrl * Vstepy * Vstepx;
		}
	}
	return k;
}

std::complex <double> ENERGY_X(const std::vector<std::complex<double> >& fka) {


	complex<double> k(0, 0);
	for (int i = 0; i < setkaBB; i++) {
		int BBnomer = i * setkaBB;
		double BBXrl = Vminx + i * Vstepx;
		for (int j = 0; j < setkaBB; j++) {
			k = k + fka[BBnomer + j] * BBXrl * BBXrl * Vstepy * Vstepx;
		}
	}
	return k;
}
