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

after having installed OpenFOAM on your system, go into this folder and
```
wmake libso
```
to clean the build, use
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

Parameters of the 𝜇(I) rheology are defined in `constant/transportProperties`:
```
transportModel muI;
muICoeffs
{
    mus      0.3;
    mud      0.5;
    I0       0.05;
    dg       1.0e-4;
    rhog     400;
    nuMax    1e3;
    nuMin    1e-5;
}

```

## Disclaimer

this code is under development and not yet benchmarked against known granular flow examples. **Use at your own risk!**
