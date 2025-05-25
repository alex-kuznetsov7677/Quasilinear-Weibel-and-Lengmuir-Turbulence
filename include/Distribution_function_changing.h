#pragma once
#include <complex>
#include <vector>

std::complex<double> DifY(int& N, std::vector<std::complex<double>>& fka);
std::complex<double> DifX(int& N, std::vector<std::complex<double>>& fka);

void PERTURBATION_OF_UNIFORM_DISTRIBUTION_COMBINED(
    int& oblast, int& size,
    int Ngarmonik_F, std::vector<double>& Kvector_F,
    std::vector<std::complex<double>>& IEy_F, std::vector<std::complex<double>>& IEx_F,
    std::vector<std::complex<double>>& b1_F, std::vector<std::complex<double>>& fk_F,  
    std::vector<std::complex<double>>& parallel_mas2,
    int Ngarmonik_TS, std::vector<double>& Kvector_TS,
    std::vector<std::complex<double>>& IEy_TS, std::vector<std::complex<double>>& IEx_TS,
    std::vector<std::complex<double>>& b1_TS, std::vector<std::complex<double>>& fk_TS
);

void PERTURBATION_OF_UNIFORM_DISTRIBUTION(
    int& oblast,
    int& size,
    int& Ngarmonik,
    std::vector<double>& Kvector,
    std::vector<std::complex<double>>& IEy,
    std::vector<std::complex<double>>& IEx,
    std::vector<std::complex<double>>& b1,
    std::vector<std::complex<double>>& fk,
    std::vector<std::complex<double>>& parallel_mas2
);

void PERTURBATION_OF_MODES_DISTRIBUTION(
    int& oblast,
    int& size,
    int& Ngarmonik,
    std::vector<double>& Kvector,
    std::vector<std::complex<double>>& IEy,
    std::vector<std::complex<double>>& IEx,
    std::vector<std::complex<double>>& b,
    std::vector<std::complex<double>>& fkcel,
    std::vector<std::complex<double>>& fk,
    std::vector<std::complex<double>>& f0k2cel,
    std::vector<std::complex<double>>& parallel_mas1
);