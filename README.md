# SplineAdaptiveFilter-Verilog
Synthesizable RTL code for the Hammerstein-type spline adaptive filter. Code accompanying the following IEEE TVLSI paper: https://doi.org/10.1109/TVLSI.2025.3621655

This repo consists of the RTL and test bench files of the following designs

1. HSAF - The original HSAF algorithm implementation without delayed adaptation
2. DHSAF-I - Retimed version of HSAF with a critical path delay of one multiplier and 4 adders
3. DHSAF-II - Retimed version of HSAF with a critical path delay of one multiplier.
4. DHSAF-I-CG - Retimed DHSAF-I with data and clock gating
5. DHSAF-II-CG - Retimed DHSAF-II with data and clock gating

The test bench and its inputs are available in the simulation folder. Multiple test inputs are available in DHSAF-I folder.

Instructions to run simulation
------------------------------

* Run the test bench using the simulator of choice. (Tested with Cadence Incisive 15.2)
* The outputs are written to the outputs folder.
* Plot the error using `plot\_rtl\_mse.m` MATLAB file. Ensure file paths are updated correctly.
