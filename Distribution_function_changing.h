#pragma once
complex<double>DifY(int& N, std::vector<complex<double> >& fka);
complex<double>DifX(int& N, std::vector<complex<double> >& fka);
void PERTURBATION_OF_UNIFORM_DISTRIBUTION(int& oblast, int& size, int& Ngarmonik, std::vector<double>& Kvector, std::vector<complex<double> >& IEy, std::vector<complex<double> >& IEx, std::vector<complex<double> >& b1, std::vector<complex<double> >& fk, std::vector<complex<double> >& parallel_mas2);
void PERTURBATION_OF_MODES_DISTRIBUTION(int& oblast, int& size, int& Ngarmonik, std::vector<double>& Kvector, std::vector<complex<double> >& IEy, std::vector<complex<double> >& IEx, std::vector<complex<double> >& b, std::vector<complex<double> >& fkcel, std::vector<complex<double> >& fk, std::vector<complex<double> >& f0k2cel,  std::vector<complex<double> >& parallel_mas1);

