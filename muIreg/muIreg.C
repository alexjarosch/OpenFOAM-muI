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

#include "muIreg.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(muIreg, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        muIreg,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::muIreg::calcNu() const
{
    const objectRegistry& db = U_.db();
    if (db.foundObject<volScalarField>("p")) {
        Info<< "Calculate mu(I) based on pressure" << endl;
        const volScalarField& ptot = U_.mesh().lookupObject<volScalarField>("p");
        tmp<volScalarField> sr(max(strainRate(), dimensionedScalar ("vSmall", dimless/dimTime, VSMALL)));
        return
        (
            max(
                min(
                    (calcMuI()*max(ptot, dimensionedScalar ("pmin", dimPressure, 1.0)))/(2.0*rhog_*sr()), nuMax_
                ), nuMin_
            )
        );
    } else{
        Info<< "Pressure not found for mu(I), return zero" << endl;
        return  tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "nuis0",
                    U_.time().timeName(),
                    U_.db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
               U_.mesh(),
               dimensionedScalar("nuis0", dimViscosity, 0.0)
           )
        );
    }

}

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::muIreg::calcMuI() const
{
    return
    (
        (mus_*I0_ + mud_*calcI() + I1_*pow(calcI(), 2))/
        (I0_ + calcI())
    );
}

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::muIreg::calcI() const
{
    const objectRegistry& db = U_.db();
    if (db.foundObject<volScalarField>("p")) {
        // Info<< "Calculate I based on pressure" << endl;
        const volScalarField& pp = U_.mesh().lookupObject<volScalarField>("p");
        tmp<volScalarField> sr(max(strainRate(), dimensionedScalar ("vSmall", dimless/dimTime, VSMALL)));
        tmp<volScalarField> sr2(mag(symm(fvc::grad(U_))));
        return
        (
            2.0*dg_*sr()*pow((max(pp, dimensionedScalar ("pmin", dimPressure, 1.0)))/rhog_, -0.5)
        );
    } else {
        Info<< "Pressure not found for I, return zero" << endl;
        return  tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "II0",
                    U_.time().timeName(),
                    U_.db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
               U_.mesh(),
               dimensionedScalar("II0", dimless, 0.0)
           )
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::muIreg::muIreg
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    muIregCoeffs_(viscosityProperties.optionalSubDict(typeName + "Coeffs")),
    mus_("mus", dimless, muIregCoeffs_),
    mud_("mud", dimless, muIregCoeffs_),
    I0_("I0", dimless, muIregCoeffs_),
    I1_("I1", dimless, muIregCoeffs_),
    dg_("dg", dimLength, muIregCoeffs_),
    rhog_("rhog", dimDensity, muIregCoeffs_),
    nuMax_("nuMax", dimViscosity, muIregCoeffs_),
    nuMin_("nuMin", dimViscosity, muIregCoeffs_),
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

bool Foam::viscosityModels::muIreg::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    muIregCoeffs_ = viscosityProperties.optionalSubDict(typeName + "Coeffs");

    muIregCoeffs_.lookup("mus") >> mus_;
    muIregCoeffs_.lookup("mud") >> mud_;
    muIregCoeffs_.lookup("I0") >> I0_;
    muIregCoeffs_.lookup("I1") >> I1_;
    muIregCoeffs_.lookup("dg") >> dg_;
    muIregCoeffs_.lookup("rhog") >> rhog_;
    muIregCoeffs_.lookup("nuMax") >> nuMax_;
    muIregCoeffs_.lookup("nuMin") >> nuMin_;

    return true;
}


// ************************************************************************* //
