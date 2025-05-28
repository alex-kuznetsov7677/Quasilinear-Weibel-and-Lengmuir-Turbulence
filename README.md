# Quasilinear Interaction of Langmuir and Weibel Turbulence

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)
[![arXiv](https://img.shields.io/badge/arXiv-2504.18859-b31b1b.svg)](https://arxiv.org/abs/2504.18859)

Numerical implementation of the quasilinear model for coupled evolution of Langmuir (two-stream) and Weibel (filamentation) turbulence in beam-plasma systems from the article: 
**A.A. Kuznetsov, Vl.V. Kocharovsky**  
*"Quasilinear interaction between Langmuir and Weibel turbulence in a beam-plasma system"*  
https://arxiv.org/abs/2504.18859

## 📌 Key Features
- Solves coupled Vlasov-Maxwell equations in 2D velocity space
- Models simultaneous evolution of:
  - Langmuir (two-stream) turbulence
  - Weibel (filamentation) turbulence
- Tracks spectral dynamics and field energy transfer
- Implements Stormer-Verlet (Leapfrog) numerical scheme
- Parallelized for efficient computation

## Compilation and Execution
  To compile the program, run:
  ```bash
  make

  After compilation, run the program with:
  ```bash
  mpirun -np <number_of_processes> -maxtime <time_of_execution_in_minutes> build/exe9_c

## Configuration
Simulation parameters are set in config.txt. Key parameters include:

- start_level_F_modes - initial perturbation level for each filamentation mode
- start_level_TS_modes - initial perturbation level for each two-stream mode
- Ngarmonik_F_x, Ngarmonik_F_y - number of filamentation modes in x and y
- Ngarmonik_TS_x, Ngarmonik_TS_y - number of two-stream modes in x and y
- V_stream - directional beam velocity
- Beta_bkg, Beta_stream - thermal velocities of background and beam plasma
- ratio - relative beam density
- setkaBB - velocity space grid size
- Tmax - maximum simulation time (in plasma frequency)
- dt - time step (in plasma frequency)
## License

Apache License 2.0 - See [LICENSE](LICENSE) for details.

## ❓ Support

For questions and bug reports:

- [Open an issue](https://github.com/alex-kuznetsov7677/Quasilinear-Weibel-and-Lengmuir-Turbulence)
- Email: [kuznetsov.alexey@ipfran.ru](mailto:kuznetsov.alexey@ipfran.ru)
