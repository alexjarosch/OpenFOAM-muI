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

#include "muIregp.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(muIregp, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        muIregp,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::muIregp::calcNu() const
{
    const objectRegistry& db = U_.db();
    const volScalarField& alphag = U_.mesh().lookupObject<volScalarField>("alpha.granul");

    if (db.foundObject<volScalarField>("p")) {
        Info<< "Calculate reg mu(I) based on pressure" << endl;
        return
        (
            max(
                min(
                    (mu_*peff_)/(2.0*rhog_*normD_*max(alphag, alphaSmall_)), nuMax_
                ), nuMin_
            )
        );
    } else{
        Info<< "Pressure not found for reg mu(I), return zero" << endl;
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
Foam::viscosityModels::muIregp::calcMu() const
{
    const objectRegistry& db = U_.db();
    dimensionedScalar IVsmall("IVsmall", dimless, VSMALL);

    if (db.foundObject<volScalarField>("p")) {
        return
        (
            // assume that I and A_m are positive
            pos(IN1_ - I_)*sqrt(alphaReg_/max((log(A_m_) - log(I_ + IVsmall)), IVsmall))
            +
            pos(I_ - IN1_)*(mus_*I0_ + mud_*I_ + muInf_*pow(I_, 2))/(I0_ + I_)
        );
    } else {
        Info<< "Pressure not found for reg mu(I), return mus" << endl;
        return mus_*calcI();
    }
}

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::muIregp::calcI() const
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
Foam::viscosityModels::muIregp::calcPeff() const
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
Foam::viscosityModels::muIregp::calcNormD() const
{
    // note this is different than the classical OpenFOAM strainRate
    return max(mag(symm(fvc::grad(U_)))/sqrt(2.0), dimensionedScalar ("vSmall", dimless/dimTime, VSMALL));
}

void Foam::viscosityModels::muIregp::initRegParameter()
{
    // I use Z instead of the lowercase zeta in Barker & Gray 2017
    scalar tanZs = tan((atan(mus_) + atan(mud_))/2.).value();
    // lowest I value, close to zero
    scalar Ilower = 1e-9;
    // upper I value
    scalar Iupper = (I0_*(tanZs - mus_)/(mud_ - tanZs)).value();
    // I for now, I call all values reg here
    scalar Ireg = (Ilower + Iupper)/2.;

    // solve the equation in a while loop
    while(Iupper - Ilower > 1e-9)
    {
        Ireg = (Ilower + Iupper)/2.;

        scalar muIreg = (mus_ + Ireg*(mud_ - mus_)/(Ireg + I0_)).value();
        // muPrime in Barker & Gray 2017, derivative of mu with respect to I
        scalar muPrime = ((I0_*(mud_ - mus_))/(Ireg + I0_)/(Ireg + I0_)).value();
        // Inupnu is the fraction used in the C eq 3.9 in Barker & Gray
        scalar Inupnu = muPrime*Ireg/muIreg;
        // C equation 3.9 in Barker & Gray 2017
        scalar C = 4*Inupnu*Inupnu - 4*Inupnu + muIreg*muIreg*(1 - Inupnu/2.)*(1 - Inupnu/2.);
        if (C < 0.)
            Iupper = Ireg;
        else
            Ilower = Ireg;
    }
    // the estimated I becomes I^N_1
    IN1_.value() = Ireg;
    // A minus according to eq 6.4
    A_m_ = IN1_*exp(alphaReg_*pow(I0_ + IN1_, 2)/pow(mus_*I0_ + mud_*IN1_ + muInf_*IN1_*IN1_, 2));

    // Info << "IN1 = " << IN1_.value() << endl
    //      << "A_minus = " << A_m_.value() << endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::muIregp::muIregp
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    muIregpCoeffs_(viscosityProperties.optionalSubDict(typeName + "Coeffs")),
    mus_("mus", dimless, muIregpCoeffs_),
    mud_("mud", dimless, muIregpCoeffs_),
    muInf_("muInf", dimless, muIregpCoeffs_),
    I0_("I0", dimless, muIregpCoeffs_),
    dg_("dg", dimLength, muIregpCoeffs_),
    rhog_("rhog", dimDensity, muIregpCoeffs_),
    nuMax_("nuMax", dimViscosity, muIregpCoeffs_),
    nuMin_("nuMin", dimViscosity, muIregpCoeffs_),
    pMin_("pMin", dimPressure, muIregpCoeffs_),
    alphaReg_("alphaReg", dimless, muIregpCoeffs_),
    IN1_("IN1", dimless, 0.),
    A_m_("A_m", dimless, 0.),
    rhoAir_("rhoAir", dimDensity, muIregpCoeffs_),
    alphaSmall_("alphaSmall", dimless, muIregpCoeffs_),
    rmHydAirP_(muIregpCoeffs_.lookupOrDefault<Switch>("rmHydAirP", false)),
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
            IOobject::NO_WRITE
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
            IOobject::NO_WRITE
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
            IOobject::NO_WRITE
        ),
        calcNormD()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::muIregp::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    muIregpCoeffs_ = viscosityProperties.optionalSubDict(typeName + "Coeffs");

    muIregpCoeffs_.lookup("mus") >> mus_;
    muIregpCoeffs_.lookup("mud") >> mud_;
    muIregpCoeffs_.lookup("muInf") >> muInf_;
    muIregpCoeffs_.lookup("I0") >> I0_;
    muIregpCoeffs_.lookup("dg") >> dg_;
    muIregpCoeffs_.lookup("rhog") >> rhog_;
    muIregpCoeffs_.lookup("nuMax") >> nuMax_;
    muIregpCoeffs_.lookup("nuMin") >> nuMin_;
    muIregpCoeffs_.lookup("pMin") >> pMin_;
    muIregpCoeffs_.lookup("alphaReg") >> alphaReg_;
    muIregpCoeffs_.lookup("rmHydAirP_") >> rmHydAirP_;
    muIregpCoeffs_.lookup("rhoAir") >> rhoAir_;
    muIregpCoeffs_.lookup("alphaSmall") >> alphaSmall_;

    initRegParameter();

    return true;
}


// ************************************************************************* //
