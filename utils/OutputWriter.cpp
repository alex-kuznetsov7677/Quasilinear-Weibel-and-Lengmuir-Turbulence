// 1. Главный заголовок (должен быть первым)
#include "OutputWriter.h"

// 2. Стандартные библиотеки (только то, что реально используется)
#include <cmath>       // для std::sqrt, std::norm, std::abs
#include <complex>     // для std::complex, std::real

// 3. Локальные заголовки
#include "Integrals.h"  // для ENERGY_X, ENERGY_Y, DENSITY

void InitializeOutputFiles(OutputFiles &files, const std::map<std::string, double>& config) {
    files.CMFt100000.open("DESCRIPTIONRe[B]__366.txt");
    
    files.CMFt100000 << "=== CONFIGURATION PARAMETERS ===\n";
  	for (const auto& pair : config) {
        files.CMFt100000 << pair.first << " = " << pair.second << "\n";
    }
	
    files.CMFt2000.open("abs[B]___366MATLAB.txt");
    files.CMFt2003.open("abs[Bts]___366MATLAB.txt");
    files.CMFt2001.open("abs[Ex]___366MATLAB.txt");
    files.CMFt2002.open("abs[Ey]___366MATLAB.txt");
    files.CMFt1301.open("Wy0__366.txt");
    files.CMFt1303.open("Wx0__366.txt");
    files.CMFt1302.open("A__366.txt");
    files.CMFt30.open("sr_b___366.txt");
    files.CMFt31.open("sr_bts___366.txt");
    files.CMFt13045.open("NU_eff__366.txt");
    files.CMFt32.open("Re[Time]__366.txt");
    files.CMFt1.open("sr_ey___366.txt");
    files.CMFt2.open("sr_ex___366.txt");
    files.CMFt_108.open("ReFUNC0_366.txt");
    files.CMFt_109.open("output366.txt");
    files.CMFt_110.open("ImFUNC111_30000_47___631_366.txt");
    files.CMFt_161.open("ImFUNC111_160000_47____631_366.txt");
}


void OutputFiles::log_fields(
    double time,
    const std::vector<std::complex<double>>& B,
    const std::vector<std::complex<double>>& IEX_TS,
    const std::vector<std::complex<double>>& IEY_TS,
    const std::vector<std::complex<double>>& B_TS,
    int Ngarmonik_F,
    int Ngarmonik_TS,
    const std::vector<std::complex<double>>& f0k2cel,
    double NU_effective,
    double Kprod
) {
  
    CMFt_109 << time << '\t';

    double Bmodkvadr = 0.;
    for (int it = 0; it < Ngarmonik_F; ++it) {
        Bmodkvadr += std::norm(B[it]);
    }

    double IExkvadr = 0., IEykvadr = 0., Bkvadr = 0.;
    for (int j = 0; j < Ngarmonik_TS; ++j) {
        IExkvadr += std::norm(IEX_TS[j]);
        IEykvadr += std::norm(IEY_TS[j]);
        Bkvadr += std::norm(B_TS[j]);
    }

    double Bmod = std::sqrt(Bmodkvadr);
    double Bmodts = std::sqrt(Bkvadr);
    double VVx0 = real(ENERGY_X(f0k2cel));
    double VVy0 = real(ENERGY_Y(f0k2cel));
    double Areal = (VVy0 / VVx0) - 1;
    double IExmod = std::sqrt(IExkvadr);
    double IEymod = std::sqrt(IEykvadr);

    CMFt_109 << real(DENSITY(f0k2cel, 0)) << '\t' << Areal << '\t' << Bmod << '\t'
              << IExmod << '\t' << IEymod << '\t' << Bmodts << '\t'
              << NU_effective << '\t' << Kprod << '\n';

    for (int p = 0; p < Ngarmonik_F; ++p) {
        CMFt2000 << std::abs(B[p]) << " ";
    }
    for (int p = 0; p < Ngarmonik_TS; ++p) {
        CMFt2002 << std::abs(IEY_TS[p]) << " ";
        CMFt2001 << std::abs(IEX_TS[p]) << " ";
        CMFt2003 << std::abs(B_TS[p]) << " ";
    }

    CMFt2000 << ";"; CMFt2001 << ";"; CMFt2002 << ";";

    CMFt30 << Bmod << ",";
    CMFt31 << Bmodts << ",";
    CMFt1 << IEymod << ",";
    CMFt2 << IExmod << ",";
    CMFt1302 << Areal << ",";
    CMFt1301 << VVy0 << ",";
    CMFt1303 << VVx0 << ",";
    CMFt13045 << NU_effective << ",";
}

void OutputFiles::log_distribution(
    int setkaBB,
    const std::vector<std::complex<double>>& f0k
) {
    for (int j = 0; j < setkaBB; ++j) {
        int BBXnomer = j * setkaBB;
        for (int k = 0; k < setkaBB; ++k) {
            int BBXnomerk = BBXnomer + k;
            CMFt_108 << real(f0k[BBXnomerk]) << " ";
        }
    }
    CMFt_108 << ";";
}
void OutputFiles::close_all() {
    CMFt100000.close();
    CMFt2000.close();
    CMFt2003.close();
    CMFt2001.close();
    CMFt2002.close();
    CMFt1301.close();
    CMFt1303.close();
    CMFt1302.close();
    CMFt30.close();
    CMFt31.close();
    CMFt13045.close();
    CMFt32.close();
    CMFt1.close();
    CMFt2.close();
    CMFt_108.close();
    CMFt_109.close();
    CMFt_110.close();
    CMFt_161.close();
}