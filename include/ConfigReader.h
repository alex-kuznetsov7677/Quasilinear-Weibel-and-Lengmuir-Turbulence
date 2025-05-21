#ifndef CONFIG_READER_H
#define CONFIG_READER_H

#include <map>
#include <string>
#include <fstream>
#include <sstream>

std::map<std::string, double> ReadConfig(const std::string& filename);

#endif
