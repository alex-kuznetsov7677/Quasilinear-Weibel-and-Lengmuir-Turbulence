#pragma once
#include <complex>
#include <vector>

std::complex<double> DifY(int& N, std::vector<std::complex<double>>& fka);
std::complex<double> DifX(int& N, std::vector<std::complex<double>>& fka);

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