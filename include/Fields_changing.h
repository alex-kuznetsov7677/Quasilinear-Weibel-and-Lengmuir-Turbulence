#pragma once 
#include <vector>
#include <cmath>
#include <complex>

std::complex<double> CURRENT_DENSITY_Y(std::vector<std::complex<double>>& fka, int p);
std::complex<double> CURRENT_DENSITY_X(std::vector<std::complex<double>>& fka, int p3);

void FILAMENTATION_FIELDS(int& oblast, int& size, int& Ngarmonik, 
    std::vector<double>& Kvector, 
    std::vector<std::complex<double>>& IEylast, std::vector<std::complex<double>>& IExlast, std::vector<std::complex<double>>& blast,
    std::vector<std::complex<double>>& IEy, std::vector<std::complex<double>>& IEx, std::vector<std::complex<double>>& b, 
    std::vector<std::complex<double>>& IEy1, std::vector<std::complex<double>>& IEx1, std::vector<std::complex<double>>& b1,
    std::vector<std::complex<double>>& IEy1last, std::vector<std::complex<double>>& IEx1last, std::vector<std::complex<double>>& b1last, 
    std::vector<std::complex<double>>& fk, std::vector<std::complex<double>>& fkcel);

void TWO_STREAM_FIELDS(int& oblast, int& size, int& Ngarmonik_puc, 
    std::vector<double>& Kvector_puc, 
    std::vector<std::complex<double>>& IEylast_puc, std::vector<std::complex<double>>& IExlast_puc, std::vector<std::complex<double>>& blast_puc,
    std::vector<std::complex<double>>& IEy_puc, std::vector<std::complex<double>>& IEx_puc, std::vector<std::complex<double>>& b_puc, 
    std::vector<std::complex<double>>& IEy1_puc, std::vector<std::complex<double>>& IEx1_puc, std::vector<std::complex<double>>& b1_puc,
    std::vector<std::complex<double>>& IEy1last_puc, std::vector<std::complex<double>>& IEx1last_puc, std::vector<std::complex<double>>& b1last_puc, 
    std::vector<std::complex<double>>& fk_puc, std::vector<std::complex<double>>& fkcel_puc);