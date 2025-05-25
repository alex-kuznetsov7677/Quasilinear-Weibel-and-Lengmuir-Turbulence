#pragma once

// ћинимальный набор заголовков (только то, что реально используетс€)
#include <fstream>   // дл€ std::ofstream
#include <map>       // дл€ std::map
#include <string>    // дл€ std::string
#include <vector>    // дл€ std::vector
#include <complex>   // дл€ std::complex
struct OutputFiles {
    std::ofstream CMFt100000;
    std::ofstream CMFt2000, CMFt2003, CMFt2001, CMFt2002;
    std::ofstream CMFt1301, CMFt1303, CMFt1302;
    std::ofstream CMFt30, CMFt31, CMFt13045;
    std::ofstream CMFt32, CMFt1, CMFt2;
    std::ofstream CMFt_108, CMFt_109, CMFt_110, CMFt_161;
    void log_number_MPI_processes(int number);
    void log_execution_time(double time);
    void close_all();
    void log_fields(double time,
                    const std::vector<std::complex<double>>& B,
                    const std::vector<std::complex<double>>& IEX_TS,
                    const std::vector<std::complex<double>>& IEY_TS,
                    const std::vector<std::complex<double>>& B_TS,
                    int Ngarmonik_F,
                    int Ngarmonik_TS,
                    const std::vector<std::complex<double>>& f0k2cel,
                    double NU_effective,
                    double Kprod);

    void log_distribution(int setkaBB,  const std::vector<std::complex<double>>& f0k);
};


void InitializeOutputFiles(OutputFiles &files, const std::map<std::string, double>& config);

