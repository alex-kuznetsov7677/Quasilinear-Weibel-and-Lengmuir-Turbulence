#pragma once
complex <double> CURRENT_DENSITY_Y(std::vector<complex<double> >& fka, int p);
complex <double> CURRENT_DENSITY_X(std::vector<complex<double> >& fka, int p3);
void FILAMENTATION_FIELDS(int& oblast, int& size, int& Ngarmonik, std::vector<double>& Kvector, std::vector<complex<double> >& IEylast, std::vector<complex<double> >& IExlast, std::vector<complex<double> >& blast,
	std::vector<complex<double> >& IEy, std::vector<complex<double> >& IEx, std::vector<complex<double> >& b, std::vector<complex<double> >& IEy1, std::vector<complex<double> >& IEx1, std::vector<complex<double> >& b1,
	std::vector<complex<double> >& IEy1last, std::vector<complex<double> >& IEx1last, std::vector<complex<double> >& b1last, std::vector<complex<double> >& fk, std::vector<complex<double> >& fkcel);
void TWO_STREAM_FIELDS(int& oblast, int& size, int& Ngarmonik_puc,  std::vector<double>& Kvector_puc, std::vector<complex<double> >& IEylast_puc, std::vector<complex<double> >& IExlast_puc, std::vector<complex<double> >& blast_puc,
	std::vector<complex<double> >& IEy_puc, std::vector<complex<double> >& IEx_puc, std::vector<complex<double> >& b_puc, std::vector<complex<double> >& IEy1_puc, std::vector<complex<double> >& IEx1_puc, std::vector<complex<double> >& b1_puc,
	std::vector<complex<double> >& IEy1last_puc, std::vector<complex<double> >& IEx1last_puc, std::vector<complex<double> >& b1last_puc, std::vector<complex<double> >& fk_puc, std::vector<complex<double> >& fkcel_puc);