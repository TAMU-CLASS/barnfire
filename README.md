# barnfire #
## Cross-Section Generation Framework <br> Texas A&M University: <br> Center for LArge-Scale Scientific Simulations


This repository houses several scripts to generate nuclear cross sections used in computational neutron transport. Beginning with ENDF/B files, the scripts allow a user to generate consistent MG (multigroup), FEDS (see below), and continuous-energy cross sections. Currently, [NJOY](http://t2.lanl.gov/nis/codes/njoy99/) is used to do most of the heavy-lifting. MG/FEDS cross sections are generated in PDT (see below) format. Continuous-energy cross sections are generated in the ACE format. 

At the present time, capabilities are being coalesced and enhanced. Future capabilities will include generation of cross section covariance matrices in MG/FEDS format and ability to generate custom thermal treatments, such as bound-thermal (S(alpha,beta)) treatments with custom phonon spectra.

### FEDS ###

FEDS, or Finite-Element-with-Discontiguous-Support, is a new energy discretization method for particle transport. It is essentially MG with discontiguous energy groups. FEDS shares similarities to MB (multiband) methods, but has several distinct characteristics. A library of high-energy-resolution spectra "snapshots" is used to generate the FEDS energy mesh. FEDS cross sections can be used in existing MG transport codes without modification, provided those codes can handle upscattering.

References for FEDS:

* A. T. Till. *Finite Elements with Discontiguous Support for Energy Discretization in Particle Transport*. Ph.D. thesis, Texas A&M University, College Station, TX (2015).

### PDT ###

PDT is a massively-parallel discrete-ordinates research transport code under development at Texas A&M University. For more information, see

* [Efficient Massively-Parallel Transport Sweeps](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.642.622&rep=rep1&type=pdf)
* [Provably Optimal Parallel Transport Sweeps on Regular Grids](https://e-reports-ext.llnl.gov/pdf/718952.pdf)
* [Validation of Full-Domain Massively Parallel Transport Sweep Algorithms](https://e-reports-ext.llnl.gov/pdf/777348.pdf)


### License ###

barnfire is distributed under the MIT license



### Installation ###

There are three major components in the repository.

* [src](src), which is a folder of Python scripts used to generate the NJOY input decks and to create the final MG/FEDS cross sections
* Cross section / data files, which are either downloaded automatically from [LANL](t2.lanl.gov), or are in the repository [here](dat/thermal_endf)
* Dependencies, whose sources are again not in the repository, but which are necessary to run the Python scripts; links may be found in the following section

### Dependencies ###

Using the Python codes requires you have the following installed on your system:

* [parallel](http://www.gnu.org/software/parallel/)
* [wget](http://www.gnu.org/software/wget/)
* Python 2 (2.7.3 or later)
* Numpy (1.8.1 or later)
* Scipy (0.13.0 or later)
* [sci-kit learn](http://scikit-learn.org/) (0.16.1 or later)
* [uncertainties](https://pypi.python.org/pypi/uncertainties/) (2.4.6 or later)
* [nuclide-data](https://github.com/attom/nuclide-data)
* [iapws](https://github.com/jjgomera/iapws/) (steam tables, used in [materials\_materials.py](src/materials_materials.py))
* [NJOY99.364](http://t2.lanl.gov/nis/codes/njoy99/) (used with [materials\_njoy.py](src/materials_njoy.py), [write\_njoy.py](src/write_njoy.py), [Readgroupr.py](src/Readgroupr.py))

### Examples ###

An example of how to run and produce FEDS XS is given in [examples/create\_xs.sh](examples/create_xs.sh). Each input has several phases, to wit, specification of PDT materials and creation of NJOY input decks, determination of FEDS energy mesh, running of NJOY, and condensation of MG cross sections (XS) from NJOY into FEDS XS. RunProb calls combinations of [materials](src/materials.py), [indicators](src/indicators.py), and [indicators\_clustering](src/indicators_clustering.py), which call [Readgroupr](src/Readgroupr.py) and [write\_njoy](src/write_njoy.py).


### Python scripts ###

[This](src) directory hosts several Python codes and a few shell scripts. Of note:

* [RunProb.sh](src/RunProb.sh) -- List of several example problem inputs. 
* [directories.py](src/directories.py) -- specifies directory locations where input data, problems, output data, and output figures may be found. You may need to __export__ the locations of your scratch directory (`SCRATCH_BARN`), ENDF/B data (`ENDF`), and NJOY executable (`NJOY`).
* [materials.py](src/materials.py) -- Caller for a set of scripts with names `materials_*.py` that handle material specification (`--njoy` option) and cross section generation (`--bonderanko` option). Material specification mixes nuclides to generate a materials and creates NJOY input scripts for each unique nuclide in the problem. Cross section generation uses material specifications and interpolates MG GENDF files produced by NJOY into FEDS XS in PDT format. To get help, run `./materials.py -h`
* [indicators.py](src/indicators.py) -- Generates spectra that are used
  either as inputs to the minimization problem, which produces the FEDS energy
  mesh, or as the basis functions, which are used in XS condensation. Spectra
  come in three flavors: the macroscopic total cross section, and solutions to
  the infinite-medium, fixed-source, slowing-down equation, with or without
  escape cross section. To get help, run `./indicators.py -h`
* [indicators\_clustering.py](src/indicators_clustering.py) -- Solves the minimization problem using clustering algorithms from machine learning. Can produce discontiguous (FEDS) or contiguous (MG) energy meshes. To get help, run `./indicators_clustering.h`
* [Readgroupr.py](src/Readgroupr.py) -- Set of functions for reading
  nuclear data files, including ENDF, PENDF, and GENDF. Can read GENDF files,
  interpolate in temperature and background cross section, condense from
  subelements to (dis)contiguous elements, and print in PDT format. Called by
  several `materials_*.py` modules. To get help, run `./Readgroupr.py -h`
* [write\_njoy.py](src/write_njoy.py) -- Creates NJOY input files for one nuclide.
  NJOY processing is done in two steps: ENDF to PENDF and PENDF to GENDF.
  Called by [materials\_njoy.py](src/materials_njoy.py)


### How do I get set up? ###

* First, download nuclide-data and put it in `dat`:

```
#!bash

cd dat
git clone https://github.com/attom/nuclide-data.git
cd ..
```

* Second, download and make NJOY2016:

```
#!bash

git clone https://github.com/njoy/NJOY2016.git
mkdir NJOY2016-build
cd NJOY2016-build
cmake -DCMAKE_BUILD_TYPE=Release ../NJOY2016
make
export NJOY=`pwd`
cd ..
```

* Then, be sure all the dependencies are met. IAPWS can be installed with `sudo pip install iapws` If you install it yourself, be sure that Python can find it. Specifically, be sure [materials\_materials.py](src/materials_materials.py) can find it.

* Change your `.bashrc`, `.cshrc`, etc. file to export the variables `SCRATCH_BARN`, `ENDF`, and `NJOY`

  * `SCRATCH_BARN` will house your output cross sections and other temporary files
  * `ENDF` will be where the raw ENDF/B and thermal ENDF/B files are housed
  * `NJOY` points to the directory with `njoy` (If using NJOY2012, rename binary to `njoy` and use v82 or later)

* Peruse [directories.py](src/directories.py) and be sure it is pointing to the proper locations on your system

### How do run these scripts on "cluster"? ###

cluster, the local compute cluster for the Nuclear Engineering Department at Texas A&M University, has been set up with all required dependencies. You need to do the following to run on cluster:

* Compile NJOY99 (the one on cluster is currently bad) and note the path to your NJOY executable, `xnjoy`. Call this location `[NJOY_PATH]`

* Log onto cluster, and put the following near the beginning of your `.cshrc` or `.bashrc` file:

```
module load gcc-4.9.1 lapack-3.5.0-gcc-4.9.1 python-2.7.8 gnu-parallel
```

For `.cshrc` files, use:

```
setenv SCRATCH_BARN "/scratch/[username]/barnfire"
setenv ENDF "/scratch/[username]/endf"
setenv NJOY "[NJOY_PATH]
```

For `.bashrc` files, use:

```
export SCRATCH_BARN=/scratch/[username]/barnfire
export ENDF=/scratch/[username]/endf
export NJOY=[NJOY_PATH]
```

Of course, you should replace `[username]` with your username and `[NJOY_PATH]` with the actual path to your NJOY executable.


### Who do I talk to? ###

* [Yunhuang Zhang](mailto:greatfrog@tamu.edu)
* [Andrew Till](mailto:attom@tamu.edu)

