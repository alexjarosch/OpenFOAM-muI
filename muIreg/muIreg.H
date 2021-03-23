/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::viscosityModels::muI

Description
    Regularized mu(I) granular rheology. Following
    Barker, T., Rauter, M., Maguire, E., Johnson, C., & Gray, J. (2021).
    Coupling rheology and segregation in granular flows.
    Journal of Fluid Mechanics, 909, A22. doi:10.1017/jfm.2020.973

Author
    Alexander Jarosch research@alexj.at

SourceFiles
    muI.C

\*---------------------------------------------------------------------------*/

#ifndef muIreg_H
#define muIreg_H

#include "viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class muIreg Declaration
\*---------------------------------------------------------------------------*/

class muIreg
:
    public viscosityModel
{
    // Private Data

        dictionary muIregCoeffs_;

        dimensionedScalar mus_;
        dimensionedScalar mud_;
        dimensionedScalar I0_;
        dimensionedScalar I1_;
        dimensionedScalar dg_;
        dimensionedScalar rhog_;
        dimensionedScalar nuMax_;
        dimensionedScalar nuMin_;

        volScalarField nu_;
        volScalarField muI_;
        volScalarField I_;


    // Private Member Functions

        //- Calculate and return the laminar viscosity
        tmp<volScalarField> calcNu() const;
        //- Calculate and return the mu(I) rheology number
        tmp<volScalarField> calcMuI() const;
        //- Calculate and return the inertial I number
        tmp<volScalarField> calcI() const;


public:

    //- Runtime type information
    TypeName("muIreg");


    // Constructors

        //- Construct from components
        muIreg
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~muIreg()
    {}


    // Member Functions

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
            return nu_;
        }

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }

        //- Correct the laminar viscosity
        virtual void correct()
        {
            nu_ = calcNu();
            muI_ = calcMuI();
            I_ = calcI();
        }

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& viscosityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
