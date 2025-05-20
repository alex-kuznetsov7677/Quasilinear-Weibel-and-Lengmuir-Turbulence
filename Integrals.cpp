#include <vector>
#include <cmath>
#include <iostream>

#include "my_variables.h"


complex <double> DENSITY(std::vector<complex<double> >& fka, int p3) {


	complex<double> k(0, 0);

	for (double i = 0; i < setkaBB; i++) {
		double BBnomer = p3 * setkaBBkvadr + i * setkaBB;
		double BBreal = -BBXmax + i * dBBX;
		for (double j = 0; j < setkaBB; j++) {

			k = k + fka[BBnomer + j] * dBBY * dBBX;
		}
	}
	return k;
}
complex <double> ENERGY_Y(std::vector<complex<double> >& fka) {

	complex<double> k(0, 0);
	for (int i = 0; i < setkaBB; i++) {
		double BBnomer = i * setkaBB;
		for (int j = 0; j < setkaBB; j++) {
			double BBYrl = -BBYmax + j * dBBY;
			k = k + fka[BBnomer + j] * BBYrl * BBYrl * dBBY * dBBX;
		}
	}
	return k;
}
complex <double> ENERGY_X(std::vector<complex<double> >& fka) {


	complex<double> k(0, 0);
	for (int i = 0; i < setkaBB; i++) {
		int BBnomer = i * setkaBB;
		double BBreal = -BBXmax + i * dBBX;
		for (int j = 0; j < setkaBB; j++) {
			k = k + fka[BBnomer + j] * BBreal * BBreal * dBBY * dBBX;
		}
	}
	return k;
}
