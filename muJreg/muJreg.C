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

#include "muJreg.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(muJreg, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        muJreg,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::muJreg::calcNu() const
{
    const objectRegistry& db = U_.db();
    if (db.foundObject<volScalarField>("p")) {
        Info<< "Calculate mu(J) based on pressure" << endl;
        tmp<volScalarField> normDlim(normD_+dimensionedScalar ("vSmall", dimless/dimTime, VSMALL));
        return
        (
            max(
                min(
                    (mu_*peff_)/(2.0*rhog_*normDlim()), nuMax_
                ), nuMin_
            )
        );
    } else{
        Info<< "Pressure not found for mu(J), return zero" << endl;
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
Foam::viscosityModels::muJreg::calcMu() const
{
    const objectRegistry& db = U_.db();
    if (db.foundObject<volScalarField>("p")) {
        return
        (
            (mus_*I0_ + mud_*I_ + muInf_*pow(I_, 2))/
            (I0_ + I_)
        );
    } else {
        Info<< "Pressure not found for mu(J), return mus" << endl;
        return mus_*calcI();
    }
}

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::muJreg::calcI() const
{
    const objectRegistry& db = U_.db();
    if (db.foundObject<volScalarField>("p")) {
        // Info<< "Calculate I based on pressure" << endl;
        return
        (
            2.0*dg_*normD_*pow((peff_ + dimensionedScalar ("psmall", dimPressure, SMALL))/rhog_, -0.5)
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

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::muJreg::calcPeff() const
{
    const objectRegistry& db = U_.db();
    if (db.foundObject<volScalarField>("p")) {
        // Info<< "Calculate I based on pressure" << endl;
        const volScalarField& ptot = U_.mesh().lookupObject<volScalarField>("p");
        return max(ptot, pMin_);
    } else {
        Info<< "Effective pressure not calculated, return zero" << endl;
        return  tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "peff0",
                    U_.time().timeName(),
                    U_.db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
               U_.mesh(),
               dimensionedScalar("peff0", dimPressure, 0.0)
           )
        );
    }
}

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::muJreg::calcNormD() const
{
    // note this is different than the classical OpenFOAM strainRate
    return mag(symm(fvc::grad(U_)))/sqrt(2.0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::muJreg::muJreg
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    muJregCoeffs_(viscosityProperties.optionalSubDict(typeName + "Coeffs")),
    mus_("mus", dimless, muJregCoeffs_),
    mud_("mud", dimless, muJregCoeffs_),
    muInf_("muInf", dimless, muJregCoeffs_),
    I0_("I0", dimless, muJregCoeffs_),
    dg_("dg", dimLength, muJregCoeffs_),
    rhog_("rhog", dimDensity, muJregCoeffs_),
    nuMax_("nuMax", dimViscosity, muJregCoeffs_),
    nuMin_("nuMin", dimViscosity, muJregCoeffs_),
    pMin_("pMin", dimPressure, muJregCoeffs_),
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
    mu_
    (
        IOobject
        (
            "mu",
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcMu()
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
    ),
    peff_
    (
        IOobject
        (
            "peff",
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcPeff()
    ),
    normD_
    (
        IOobject
        (
            "normD",
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcNormD()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::muJreg::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    muJregCoeffs_ = viscosityProperties.optionalSubDict(typeName + "Coeffs");

    muJregCoeffs_.lookup("mus") >> mus_;
    muJregCoeffs_.lookup("mud") >> mud_;
    muJregCoeffs_.lookup("muInf") >> muInf_;
    muJregCoeffs_.lookup("I0") >> I0_;
    muJregCoeffs_.lookup("dg") >> dg_;
    muJregCoeffs_.lookup("rhog") >> rhog_;
    muJregCoeffs_.lookup("nuMax") >> nuMax_;
    muJregCoeffs_.lookup("nuMin") >> nuMin_;
    muJregCoeffs_.lookup("pMin") >> pMin_;

    return true;
}


// ************************************************************************* //
