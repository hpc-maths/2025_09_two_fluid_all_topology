This repository contains a numerical solver for the all-topology model presented in ``An all-topology two-fluid model for two-phase flows derived through Hamilton's Stationary Action Principle'' (https://hal.science/hal-05249139).
The implementation is carried out in the framework of [samurai](https://github.com/hpc-maths/samurai).
The spatial discretization is based on the finite volume method using a non-conservative Rusanov flux.

In order to build the executable, please first install ```samurai``` folliwing the instructions at https://github.com/hpc-maths/samurai. Then, move into the directory with the test case of interest run

```bash
source configure.sh
```

Finally, to run the program, execute

```bash
./build/all_topology_test
```

The parameters declared in ```input.json``` will be used to solve the corresponding Riemann problem. More specifically, ```1D_RIEMANN``` simulates a 1D Riemann problem, whereas ```2D_RIEMANN``` simulates a 2D Riemann problem. The solver is dimension-independent. The main difference regards only the initial and boundary conditions, which are obviously test-case dependent. This is the reason of the two different ```All_topology_solver.hpp```, i.e. to handle initial conditions (also through ```containers.hpp```) and boundary conditions (see ```user_bc.hpp``` for the 2D Riemann problem). Moreover, the spatial dimension is a template parameter, hence to be known at compile-time (it is fixed in ```.cpp``` file). If you take into account these constraints, you can easily design you own test case to simulate the all-topology model.

If one is interested in the results only, you can find them into the subfolders ```RESULTS```. The results are in ```.h5``` format. For 2D problems, post-processing and visualization can be easily done in [paraview](https://www.paraview.org/), where for the 1D test cases, please use the provided routine ```plot_results.py```. You need ```python``` with ```matplotlib``` and ```h5py```. Then simply run, e.g.,

```bash
python plot_results.py RESULTS/Rusanov_Flux_order1_ite_0 --field conserved_0
```

to visualize the volume fraction of phase 1. The order of the parameters can be find in the struct ```EquationData```in ```include/flux_base.hpp```. For the sake of completeness, we mentnion that it is: volume fraction phase 1, partial mass phase 1, momentum phase 1, specific total energy of phase 1, and then the same order for the fields associated to phase 2.\
The parameter ```field``` can be filled with the conserved variables (from ```conserved_0``` to ```conserved_6```) or with primtive variables, such as phasic densities, pressures, velocities, speeds of sound, specific entropies, temperatures etc... Please refer to the file ```All_topology_solver.hpp``` and in particular to the function ```create_fields``` for an overview of the different options. As a matter of example, you can execute

```bash
python plot_results.py RESULTS/Rusanov_Flux_order1_ite_0 --field rho1
```

to visualize the density of phase 1.

In order to further facilitate the experiences for the visualization of 1D Riemann problems, one can execute the jupyter notebook ```post_process_1D_Riemann
.ipnyb``` where the routine ```plot_1D_Riemann_results``` provides visualization tools. It is sufficient to pass in input the base filename of the output, the file with the time saving (both given by the code), and the name of the fields for which post-processing is deisred (up to 2 simultaneously)

For the reference Riemann problem already available in the repository, one can play around with data at

https://hpc-maths.github.io/2025_09_two_fluid_all_topology/
