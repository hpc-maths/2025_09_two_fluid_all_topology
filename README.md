This repository contains a numerical solver for the novel model presented in ``An all-topology two-fluid model for two-phase flows derived through Hamilton's Stationary Action Principle'' (https://hal.science/hal-05249139).
The implementation is carried out in the framework of [samurai](https://github.com/hpc-maths/samurai).
The spatial discretization is based on the finite volume method using a non-conservative Rusanov flux.

In order to build the executable, please first install ```samurai``` following the instructions at https://github.com/hpc-maths/samurai. Please, include ```nlohmann```in your ```conda```environment

```bash
conda install nlohmann_json
```

To reproduce the figures and videos shown on the website, a ```paraview``` installation is required as well as ```ffmpeg``` library

```bash
conda install ffmpeg
```

We provide a conda environment file (```conda/environment_samurai.yml```) to facilitate the setup process. You can create and activate the conda environment by running the following commands in your terminal:

```bash
conda env create -f conda/environment_samurai.yml
conda activate samurai-0-27
```

Then, move into the directory with the test case of interest and run

```bash
source configure.sh
```

Finally, to run the program, execute

```bash
./All_topology_test
```

The parameters declared in ```input.json``` will be used to solve the corresponding Riemann problem. More specifically, ```1D_RIEMANN``` simulates a 1D Riemann problem, whereas ```2D_RIEMANN``` simulates a 2D Riemann problem. The solver is dimension-independent. The main difference regards only the initial and boundary conditions, which are obviously test-case dependent. This is the reason of the two different ```All_topology_solver.hpp```, i.e. to handle initial conditions (also through ```containers.hpp```) and boundary conditions (see ```user_bc.hpp``` for the 2D Riemann problem). Moreover, the spatial dimension is a template parameter, hence to be known at compile-time (it is fixed in ```.cpp``` file). If you take into account these constraints, you can easily design you own test case to simulate the all-topology model.

If one is interested in the results only, you can find them into the subfolders ```RESULTS``` for the 1D test case. The results are in ```.h5``` format. For 2D problems, due to its size, the data is not provided. However, compiling and executing the code yields ```.xdmf``` and ```.h5``` files. Executing ```export_data.py``` in the subfolder ```RESULTS``` with ```pvpython``` (native ```paraview``` python executable) will generate the video

```bash
pvpython export_data.py
````
Please, notice that ```paraview``` has to be in your ```PATH``` so as to avoid the error

```bash
-bash: paraview: command not found
```

In order to further facilitate the experiences for the visualization for 1D problems, one can execute the jupyter notebook ```post_process_1D_Riemann
.ipnyb``` where the routine ```plot_1D_Riemann_results``` provides visualization tools. It is sufficient to pass in input the base filename of the output, the file with the time saving (both given by the code), and the name of the fields for which post-processing is desired (up to 2 simultaneously)

For the reference Riemann problems already available in the repository, one can play around with data at

https://hpc-maths.github.io/2025_09_two_fluid_all_topology/
