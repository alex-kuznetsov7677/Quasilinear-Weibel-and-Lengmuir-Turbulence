#include <vector>
#include <cmath>


#include "my_variables.h"
#include "ConfigReader.h"

void Creating_grid_wavevectors(std::vector<double>& Kvector_F, std::vector<double>& Kvector_TS) {

	auto config = ReadConfig("config.txt");
	for (double i = 0; i < Ngarmonik_F_x; i++) {
		Kvector_F.push_back(config["kmin_x_F"]+i*config["kstep_x_F"]);
	}
	for (double i = 0; i < Ngarmonik_F_y; i++) {
		Kvector_F.push_back(config["kmin_y_F"]+i*config["kstep_y_F"]);
	}

	for (double i = 0; i < Ngarmonik_TS_x; i++) {
		Kvector_TS.push_back(config["kmin_x_TS"]+i*config["kstep_x_TS"]);
	}
	for (double i = 0; i < Ngarmonik_TS_y; i++) {
		Kvector_TS.push_back(config["kmin_y_TS"]+i*config["kstep_y_TS"]);
	}
}
