# incompressible 𝜇(I) rheology for OpenFOAM

this is a viscosity module for OpenFOAM with a focus on granular flow.
Currently I work on implementing:

* a standard 𝜇(I) rheology for incompressible solvers (e.g. interFoam)
* a well-posed 𝜇(I) rheology for incompressible solvers
* a compressible 𝜇(I) rheology (one day...)

as a scientific basis for this work I follow the following papers:

* [Barker, T., & Gray, J. (2017). Partial regularisation of the incompressible 𝜇(I)-rheology for granular flow. Journal of Fluid Mechanics, 828, 5-32.](https://doi.org/10.1017/jfm.2017.428)
* [Barker, T., Rauter, M., Maguire, E., Johnson, C., & Gray, J. (2021). Coupling rheology and segregation in granular flows. Journal of Fluid Mechanics, 909, A22.](https://doi.org/10.1017/jfm.2020.973)

## Installation

after having installed OpenFOAM on your system, clone this git repository to a directory on your system:

```
git clone https://github.com/alexjarosch/OpenFOAM-muI.git
```
now change into the folder just created and compile the code

```
cd OpenFOAM-muI
wmake libso
```
If you want to clean the build, use
```
wclean
```

## Usage

The viscosity library will be placed into `$FOAM_USER_LIBBIN` and is called `libviscositymuI.so`.
To make use of the library with a solver place the following code at the bottom of `system/controlDict`:
```
libs
(
    "libviscositymuI.so"
);
```

Parameters of the regularized 𝜇(I) rheology are defined in `constant/transportProperties`:
```
    transportModel muIreg;
    muIregCoeffs
    {
        mus      0.342;
        mud      0.557;
        muInf    0.05;
        alphaReg 1.9;
        I0       0.069;
        dg       5.0e-4;
        rhog     2500;
        nuMax    1e3;
        nuMin    1e-5;
        pMin     10.0;
    }
    rho          2500;
```

## Examples

in the `run` directory you can find two tutorial cases, which are modified version of the OpenFOAM interFoam tutorial case "damBreak". You can run a regularized and a non-regularized version of the examples.

## Which Implementation

Two versions of a 𝜇(I) rheology are implemented:
* muI follows the classical, unregularized 𝜇(I) rheology, implemented in accordance with equation 2.21 in Barker & Gray 2017.
* muIreg is the regularized version of a 𝜇(I) rheology, implemented accordance equation 6.3 in Barker & Gray 2017.


## Disclaimer

this code is under development and not yet benchmarked against known granular flow examples. **Use at your own risk!**
