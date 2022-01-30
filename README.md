# Incompressible 𝜇(I) rheology for OpenFOAM

this is a viscosity module for OpenFOAM with a focus on granular flow.
Currently I work on implementing:

* a standard 𝜇(I) rheology for incompressible solvers (e.g. interFoam)
* a regularized 𝜇(I) rheology for incompressible solvers

as a scientific basis for this work I follow the following papers:

* [Barker, T., & Gray, J. (2017). Partial regularisation of the incompressible 𝜇(I)-rheology for granular flow. Journal of Fluid Mechanics, 828, 5-32.](https://doi.org/10.1017/jfm.2017.428)
* [Barker, T., Rauter, M., Maguire, E., Johnson, C., & Gray, J. (2021). Coupling rheology and segregation in granular flows. Journal of Fluid Mechanics, 909, A22.](https://doi.org/10.1017/jfm.2020.973)


[![Release](https://img.shields.io/badge/release-1.0.1-blue.svg)](https://github.com/alexjarosch/OpenFOAM-muI)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4964189.svg)](https://doi.org/10.5281/zenodo.4964189)
[![OpenFOAM 8](https://img.shields.io/badge/OpenFOAM-8-brightgreen.svg)](https://openfoam.org/)
[![OpenFOAM 9](https://img.shields.io/badge/OpenFOAM-9-brightgreen.svg)](https://openfoam.org/)

## Citing the code

If you use this OpenFOAM-muI module in your work/research/simulations, please cite the code in your publications by using the dedicated DOI number for a given code version: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4964189.svg)](https://doi.org/10.5281/zenodo.4964189)

The DOI number displayed here always refers to the latest release of OpenFOAM-muI. If you require a DOI for an older release, please refer to the Zenodo homepage of this code (by clicking the DOI badge above) and find a list of all DOI's for all code releases.

## Installation

Current development on the main branch is for OpenFOAM version 9. After having installed [OpenFOAM](https://openfoam.org/download/) version 9 on your system, clone this git repository to a directory on your system:

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

For an explanation what the parameters mean and recommended values for a given granular medium, please consult the published 𝜇(I) rheology literature.

## Examples

In the `tutorials` directory you can find two tutorial cases, which are modified version of the OpenFOAM interFoam tutorial case ["damBreak"](https://github.com/OpenFOAM/OpenFOAM-9/tree/master/tutorials/multiphase/interFoam/laminar/damBreak/damBreak) and are inspired by the benchmark of [Balmforth, N. J., & Kerswell, R. R. (2005). Granular collapse in two dimensions. Journal of Fluid Mechanics, 538, 399-428.](https://doi.org/10.1017/S0022112005005537). You can run a regularized and a non-regularized version of the examples.

**Note:** To reproduce the actual benchmark data of Balmforth and Kerswell, a more advanced configuration of the example cases is required.

## Which Implementation

Two versions of a 𝜇(I) rheology are implemented:
* muI follows the classical, unregularized 𝜇(I) rheology, implemented in accordance with equation 2.21 in Barker & Gray 2017.
* muIreg is the regularized version of a 𝜇(I) rheology, implemented accordance equation 6.3 in Barker & Gray 2017.

## Regularized 𝜇(I) Python utility

In the 'utils' directory there is a Python script which you can use to estimate the I threshold for the 𝜇(I) regularization. It is the same algorithm the muIreg viscosity module uses.

## Acknowledgements

[OpenFOAM](https://github.com/OpenFOAM/) is a free, open-source software for computational fluid dynamics (CFD). The development is primarily done by [CFD Direct](https://cfd.direct/) on behalf of the [OpenFOAM](https://openfoam.org/) Foundation.

## Disclaimer

This code is under development and not yet fully benchmarked against known granular flow examples. **Use at your own risk!**
