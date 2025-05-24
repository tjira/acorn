# TODO List

- [ ] The resampling within unSMASH algorithm should probably not resample all the spheres if resampling during jump.
- [x] Implement the TDC using the logarithm of U.
- [ ] Implement the function that generates all combinations using GSL instead of my shitty function.
- [x] Make the OpenBLAS library compile with arbitrary number of cores, not just the number of cores of the current machine.
- [x] Fix the bug when loading integrals from disk, the basis still needs to be specified.
- [x] Make the Eigen tensor contraction parallel.
- [ ] Velocity and acceleration in lambda TSH should be in the direction of motion for higher dimensional systems.
- [x] Check if i really need to store both J_AO and J_AO_A for the HF method.
- [x] The coefficient matrix from GHF is not necessarily in a,b,a,b order, which might mess up following calcualtions.
- [ ] Better integrator for the Bohmian trajectories so that the trajectories would not lag behind.
- [ ] Bohmian dynamics should have an option to output only the mean trajectories.
