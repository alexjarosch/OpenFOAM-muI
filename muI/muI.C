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

Author
    Alexander Jarosch research@alexj.at

\*---------------------------------------------------------------------------*/

#include "muI.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(muI, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        muI,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::muI::calcNu() const
{
    Info<< "Calculate mu(I) based on pressure\n" << endl;
    const volScalarField& ptot = U_.mesh().lookupObject<volScalarField>("p_rgh");
    tmp<volScalarField> sr(max(strainRate(), dimensionedScalar ("vSmall", dimless/dimTime, VSMALL)));
    return
    (
        max(
            min(
                (calcMuI()*ptot)/(2.0*rhog_*sr()), nuMax_
            ), nuMin_
        )
    );
}

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::muI::calcMuI() const
{
    return
    (
        (mus_*I0_ + mud_*calcI())/
        (I0_ + calcI())
    );
}

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::muI::calcI() const
{
    Info<< "Calculate I based on pressure\n" << endl;
    const volScalarField& pp = U_.mesh().lookupObject<volScalarField>("p_rgh");
    tmp<volScalarField> sr(max(strainRate(), dimensionedScalar ("vSmall", dimless/dimTime, VSMALL)));
    tmp<volScalarField> sr2(mag(symm(fvc::grad(U_))));
    return
    (
        2.0*dg_*sr()*pow((max(pp, dimensionedScalar ("pmin", dimPressure, 1.0)))/rhog_, -0.5)
    );
    // return  tmp<volScalarField>
    // (
    //     new volScalarField
    //     (
    //         IOobject
    //         (
    //             "II0",
    //             U_.time().timeName(),
    //             U_.db(),
    //             IOobject::NO_READ,
    //             IOobject::NO_WRITE,
    //             false
    //         ),
    //        U_.mesh(),
    //        dimensionedScalar("II0", dimless, 0.5)
    //    )
    // );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::muI::muI
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    muICoeffs_(viscosityProperties.optionalSubDict(typeName + "Coeffs")),
    mus_("mus", dimless, muICoeffs_),
    mud_("mud", dimless, muICoeffs_),
    I0_("I0", dimless, muICoeffs_),
    dg_("dg", dimLength, muICoeffs_),
    rhog_("rhog", dimDensity, muICoeffs_),
    nuMax_("nuMax", dimViscosity, muICoeffs_),
    nuMin_("nuMin", dimViscosity, muICoeffs_),
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcNu()
    ),
    muI_
    (
        IOobject
        (
            "muI",
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcMuI()
    ),
    I_
    (
        IOobject
        (
            "I",
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcI()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::muI::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    muICoeffs_ = viscosityProperties.optionalSubDict(typeName + "Coeffs");

    muICoeffs_.lookup("mus") >> mus_;
    muICoeffs_.lookup("mud") >> mud_;
    muICoeffs_.lookup("I0") >> I0_;
    muICoeffs_.lookup("dg") >> dg_;
    muICoeffs_.lookup("rhog") >> rhog_;
    muICoeffs_.lookup("nuMax") >> nuMax_;
    muICoeffs_.lookup("nuMin") >> nuMin_;

    return true;
}


// ************************************************************************* //
