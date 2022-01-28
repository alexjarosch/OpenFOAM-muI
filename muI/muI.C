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
    const objectRegistry& db = U_.db();
    if (db.foundObject<volScalarField>("p")) {
        Info<< "Calculate mu(I) based on pressure" << endl;
        return
        (
            max(
                min(
                    (mu_*peff_)/(2.0*rhog_*normD_), nuMax_
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
Foam::viscosityModels::muI::calcMu() const
{
    const objectRegistry& db = U_.db();
    if (db.foundObject<volScalarField>("p")) {
        return
        (
            (mus_*I0_ + mud_*I_)/
            (I0_ + I_)
        );
    } else {
        Info<< "Pressure not found for mu(I), return mus" << endl;
        return mus_*calcI();
    }
}

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::muI::calcI() const
{
    const objectRegistry& db = U_.db();
    if (db.foundObject<volScalarField>("p")) {
        // Info<< "Calculate I based on pressure" << endl;
        return
        (
            2.0*dg_*normD_*pow(peff_/rhog_, -0.5)
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
Foam::viscosityModels::muI::calcPeff() const
{
    const objectRegistry& db = U_.db();
    const Time& runTime= db.time();
    if (db.foundObject<volScalarField>("p") && runTime.timeIndex() > 1) {
        // Info<< "Calculate I based on pressure" << endl;
        const volScalarField& ptot = U_.mesh().lookupObject<volScalarField>("p");
        const volScalarField& gh = U_.mesh().lookupObject<volScalarField>("gh");
        if (rmHydAirP_) {
            Info<< "Hydrostatic pressure of air phase removed" << endl;
            return max(ptot - rhoAir_*gh, pMin_);
        } else {
            return max(ptot, pMin_);
        }

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
               dimensionedScalar("peff0", dimPressure, SMALL)
           )
        );
    }
}

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::muI::calcNormD() const
{
    // note this is different than the classical OpenFOAM strainRate
    return max(mag(symm(fvc::grad(U_)))/sqrt(2.0), dimensionedScalar ("vSmall", dimless/dimTime, VSMALL));
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
    pMin_("pMin", dimPressure, muICoeffs_),
    rhoAir_("rhoAir", dimDensity, muICoeffs_),
    rmHydAirP_(muICoeffs_.lookupOrDefault<Switch>("rmHydAirP", false)),
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
    muICoeffs_.lookup("pMin") >> pMin_;
    muICoeffs_.lookup("rmHydAirP_") >> rmHydAirP_;
    muICoeffs_.lookup("rhoAir") >> rhoAir_;
    return true;
}


// ************************************************************************* //
