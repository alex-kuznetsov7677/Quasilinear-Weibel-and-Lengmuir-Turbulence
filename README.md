# Quasilinear Interaction of Langmuir and Weibel Turbulence

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)
[![arXiv](https://img.shields.io/badge/arXiv-2504.18859-b31b1b.svg)](https://arxiv.org/abs/2504.18859)

Numerical implementation of the quasilinear model for coupled evolution of Langmuir (two-stream) and Weibel (filamentation) turbulence in beam-plasma systems from the article: 
**A.A. Kuznetsov, Vl.V. Kocharovsky**  
*"Quasilinear interaction between Langmuir and Weibel turbulence in a beam-plasma system"*  
https://arxiv.org/abs/2504.18859

## 📌 Key Features
- Solves the Vlasov-Maxwell equations in 2D velocity-space and wavevector-space grids.
- Simulates quasilinear interaction between Langmuir (two-stream) and Weibel (filamentation) instabilities in a collisionless, unmagnetized plasma-beam system.
- Tracks spectral dynamics and field energy transfer
- Implements Stormer-Verlet (Leapfrog) numerical scheme
- MPI-parallelized for efficient computation of mode dynamics across distributed systems (e.g., KIAM K100 supercomputer).

## Compilation and Execution

This code was developed on the [KIAM K100 supercomputer (Moscow, Russia)](https://www.kiam.ru/MVS/).

  To compile the program, run:
  ```bash
  make
  ```
  After compilation, run the program with:
  ```bash
  mpirun -np <number_of_processes> -maxtime <time_of_execution_in_minutes> build/exe9_c
  ```
## Configuration
Simulation parameters are set in config.txt. Key parameters include:

- start_level_F_modes  &mdash;  initial perturbation level for each filamentation mode
- start_level_TS_modes  &mdash;  initial perturbation level for each two-stream mode
  
- Ngarmonik_F_x, Ngarmonik_F_y  &mdash;  number of filamentation modes in x and y
- Ngarmonik_TS_x, Ngarmonik_TS_y  &mdash;  number of two-stream modes in x and y
  
- V_stream  &mdash;  directional beam velocity
- Beta_bkg, Beta_stream  &mdash;  thermal velocities of background and beam plasma
- ratio  &mdash;  relative beam density

- setkaBB  &mdash;  velocity space grid size
- ratio_Vminx_to_Beta_bkg  &mdash;  X-min boundary - Lower limit of velocity grid in x-direction (in units of β<sub>bkg</sub>)
- ratio_Vminy_to_Beta_bkg  &mdash;  Y-min boundary - Lower limit of velocity grid in y-direction (in units of β<sub>bkg</sub>)
- ratio_Vstepx_multiple_setkaBB_to_Beta_bkg  &mdash;  X-step scaling - Determines velocity grid spacing in x-direction: Δv<sub>x</sub> = (β<sub>bkg</sub> × 8)/setkaBB
- ratio_Vstepy_to_Vstepx  &mdash;  Y-step ratio - Velocity grid spacing in y-direction: Δv<sub>y</sub> = 1.2 × Δv<sub>x</sub>

- Tmax  maximum simulation time (in plasma frequency)
- dt  time step (in plasma frequency)

- kmin_x_TS  &mdash;  Minimum k<sub>x</sub> - Lower bound of wavevector range in x-direction
- kstep_x_TS  &mdash;  k<sub>x</sub> step - Spacing between wavevectors in x-direction
- kmin_y_TS  &mdash;  Minimum k<sub>y</sub> - Lower bound of wavevector range in y-direction
- kstep_y_TS  &mdash;  k<sub>y</sub> step - Spacing between wavevectors in y-direction

- kmin_x_F  &mdash;  Minimum k<sub>x</sub> - Lower bound of wavevector range in x-direction
- kstep_x_F  &mdash;  k<sub>x</sub> step - Spacing between wavevectors in x-direction
- kmin_y_F  &mdash;  Minimum k<sub>y</sub> - Lower bound of wavevector range in y-direction
- kstep_y_F  &mdash;  k<sub>y</sub> step - Spacing between wavevectors in y-direction

## Requirements

- C++ compiler (e.g., g++)
- MPI (e.g., OpenMPI)
- Standard math libraries

## License

Apache License 2.0 - See [LICENSE](LICENSE) for details.

## ❓ Support

For questions and bug reports:

- [Open an issue](https://github.com/alex-kuznetsov7677/Quasilinear-Weibel-and-Lengmuir-Turbulence)
- Email:
  - [kuznetsov.alexey@ipfran.ru](mailto:kuznetsov.alexey@ipfran.ru)
  - [a.kuznetsov7677@gmail.com](mailto:a.kuznetsov7677@gmail.com)
