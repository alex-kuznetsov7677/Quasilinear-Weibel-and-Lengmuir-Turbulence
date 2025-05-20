#include <vector>
#include <cmath>
#include <iostream>

#include "my_variables.h"

void Creating_grid_wavevectors(std::vector<double>& Kvector_F, std::vector<double>& Kvector_TS) {
	for (double i = 0; i < Ngarmonik_F_x; i++) {
		Kvector_F.push_back(0.12 * pow(2., i * 0.06));
	}
	for (double i = 0; i < Ngarmonik_F_y; i++) {
		Kvector_F.push_back(-0.06 + 0.012 * i);
	}

	for (double i = 0; i < Ngarmonik_TS_x; i++) {
		Kvector_TS.push_back(-0.5 + 0.05 + 0.1 * i);
	}
	for (double i = 0; i < Ngarmonik_TS_y; i++) {
		Kvector_TS.push_back(2.5 + 0.05 * i);
	}
}
