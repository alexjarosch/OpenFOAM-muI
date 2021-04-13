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
    Foam::viscosityModels::muJreg

Description
    Standard mu(I) granular rheology. Following
    Barker, T., Rauter, M., Maguire, E., Johnson, C., & Gray, J. (2021).
    Coupling rheology and segregation in granular flows.
    Journal of Fluid Mechanics, 909, A22. doi:10.1017/jfm.2020.973

Author
    Alexander Jarosch research@alexj.at

SourceFiles
    muJreg.C

\*---------------------------------------------------------------------------*/

#ifndef muJreg_H
#define muJreg_H

#include "viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class muJreg Declaration
\*---------------------------------------------------------------------------*/

class muJreg
:
    public viscosityModel
{
    // Private Data

        dictionary muJregCoeffs_;

        dimensionedScalar mus_;
        dimensionedScalar mud_;
        dimensionedScalar muInf_;
        dimensionedScalar I0_;
        dimensionedScalar dg_;
        dimensionedScalar rhog_;
        dimensionedScalar nuMax_;
        dimensionedScalar nuMin_;
        dimensionedScalar pMin_;

        volScalarField nu_;
        mutable volScalarField mu_;
        mutable volScalarField I_;
        mutable volScalarField peff_;
        mutable volScalarField normD_;


    // Private Member Functions

        //- Calculate and return the laminar viscosity
        tmp<volScalarField> calcNu() const;
        //- Calculate and return the mu(I) rheology number
        tmp<volScalarField> calcMu() const;
        //- Calculate and return the inertial I number
        tmp<volScalarField> calcI() const;
        //- Calculate and return the effective pressure
        tmp<volScalarField> calcPeff() const;
        //- Calculate and return the effective pressure
        tmp<volScalarField> calcNormD() const;


public:

    //- Runtime type information
    TypeName("muJreg");


    // Constructors

        //- Construct from components
        muJreg
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~muJreg()
    {}


    // Member Functions

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
            return nu_;
        }

        //- Return the mu(J) friction
        virtual tmp<volScalarField> mu() const
        {
            return mu_;
        }

        //- Return the inertial number I
        virtual tmp<volScalarField> I() const
        {
            return I_;
        }

        //- Return the effective pressure
        virtual tmp<volScalarField> peff() const
        {
            return peff_;
        }

        //- Return the norm of the deviatoric strain rate
        virtual tmp<volScalarField> normD() const
        {
            return normD_;
        }

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }

        //- Correct the laminar viscosity
        virtual void correct()
        {
            normD_ = calcNormD();
            peff_ = calcPeff();
            I_ = calcI();
            mu_ = calcMu();
            nu_ = calcNu();
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
