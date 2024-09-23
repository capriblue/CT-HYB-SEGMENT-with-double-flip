#### CT-HYB-SEGMENT
`CT-HYB-SEGMENT` is the hybridization expansion continuous-time quantum Monte Carlo code, implemented and described in:
H. Hafermann, P. Werner, E. Gull, CPC 184, 1280 (2013).
This version of the code is its adaptation to *ALPSCore* library. 

#### Dependencies
1. ALPSCore (http://alpscore.org)
2. LAPACK
3. NFFT3 (https://www-user.tu-chemnitz.de/~potts/nfft/)

#### Usage 
```
alps_cthyb parameters
```
See `alps_cthyb --help` for the list of parameters

Python bindings are temporarily disabled

#### Installation

Replace the variables with `${}` to the proper values for your system.

```ShellSession
$ git clone https://github.com/ALPSCore/CT-HYB-SEGMENT.git
$ cd CT-HYB-SEGMENT
$ CTHYB_INSTALL_DIR=${where-to-install}
$ mkdir build
$ cd build
$ NFFT3_DIR=${NFFT3-install-dir} CC=${C-compiler} CXX=${C++-compiler} cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DALPSCore_DIR=${ALPSCore_INSTALL_DIR}/share/ALPSCore \
    -DCMAKE_INSTALL_PREFIX=$CTHYB_INSTALL_DIR
$ make
$ make install
```

#### MY MODIFICATION for this fork. 

I added the directory structure to a original CT-HYB-SEGMENT repo and modified the cmake files to be able to compile this sources and silence cmake warning when compiling. But I did not change any critical codes. If you find it difficult to read original code, check my repository. 

I diveided codes into 5 parts: `evaluation`, `model`, `monte_carlo_definition`, `state`, `utility`.

- `evaluation`: After Monte Carlo sampling, this codes accumulate whole data from each machines and do some statistic calculation. I summarized these function into this folder. you can add new features in `evaluation/evaluator` directory. 
- `model`: This part deals with the model of the CT-QMC. As you know, CT-HYB solves the impurity problem which is mainly composed of 2 parts: hybridization function and original lattice parameters such as U/t, mu etc. The code in this part transforms data you input the runner. 
- `monte_carlo_definiton`: This application is based on the ALPSCore library providing the Monte Carlo calculation templates. The codes in this folder charactaraize a variety of setting for Monte Carlo calculation, like the initialize and finalize definitons, the progress of the calculation, concrete way of calculation which is accutally  deveided to the funcitons and defined other folders.
- `state` : This folder has the state definitions. Changing the state (configuration; those who familier with the correct expressions) of the model, we can obtain the natural distribution and thermo dynamical values. In CT-QMC, the state is correspoinding to the segments. but for the faster calculation, the implementation of CTQMC also have the Matrix which is reduced from the segments. `state/hyblocal` contains the segment definition and utility for changing the collection of segments. `state/hybmatrix` contais matrix definiton and utility for changing the states.
- `utility` : The files I could not classify are in this folder. 