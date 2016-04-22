/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "turbulentPotentialJoTomega.H"
#include "addToRunTimeSelectionTable.H"
#include "backwardsCompatibilityWallFunctions.H"
#include "components.H"
#include "fvCFD.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(turbulentPotentialJoTomega, 0);
addToRunTimeSelectionTable(RASModel, turbulentPotentialJoTomega, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


tmp<volScalarField> turbulentPotentialJoTomega::Ts() const
{
	if(tslimiter_ == "true")
	{
        return max(1.0/omega_, minTS());
	}
	
    return (k_/epsilon_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentPotentialJoTomega::turbulentPotentialJoTomega
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel
)
:
    RASModel(typeName, U, phi, lamTransportModel),


    cEp1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cEp1",
            coeffDict_,
            1.45
        )
    ),
    cEp2con_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cEp2con",
            coeffDict_,
            1.83
        )
    ),
    cEp3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cEp3",
            coeffDict_,
            0.15
        )
    ),
    cD1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cD1",
       	    coeffDict_,
            0.5
        )
    ),
    cD2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cD2",
       	    coeffDict_,
            0.88
        )
    ),
    cVv1_
    (
     	dimensioned<scalar>::lookupOrAddToDict
        (
            "cVv1",
            coeffDict_,
            2.0
        )
    ),
    cTv1_
    (
     	dimensioned<scalar>::lookupOrAddToDict
        (
            "cTv1",
            coeffDict_,
            0.0
        )
    ),
    cP1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cP1",
            coeffDict_,
            2.0
        )
    ),
    cP2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cP2",
            coeffDict_,
            0.6
        )
    ),
    cP3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cP3",
            coeffDict_,
            0.12
        )
    ),
    cP4_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cP4",
            coeffDict_,
            0.85714
        )
    ),
    cPphi_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cPphi",
            coeffDict_,
            2.0
        )
    ),
    cMu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cMu",
            coeffDict_,
            0.21
        )
    ),

    cT_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cT",
            coeffDict_,
            0.0033
        )
    ),

    cPr_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cPr",
            coeffDict_,
            1.0
        )
    ),

    cEhm_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cEhm",
            coeffDict_,
            10.0
        )
    ),
	
    cEhR_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cEhR",
            coeffDict_,
            1.0
        )
    ),

	gT1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gT1",
            coeffDict_,
            0.0
        )
    ),
	
	gT2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gT2",
            coeffDict_,
            0.0
        )
    ),
	
	cNF_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cNF",
            coeffDict_,
            10.0
        )
    ),
    
    	cPw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cPw",
            coeffDict_,
            18.0
        )
    ),
    
    sigmaKInit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaKInit",
            coeffDict_,
            1.0
        )
    ),

    sigmaOmInit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaOmInit",
            coeffDict_,
            0.833
        )
    ),
    
    sigmaEpsInit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEpsInit",
            coeffDict_,
            0.833
        )
    ),

    sigmaEpsVisc_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEpsVisc",
            coeffDict_,
            1.0
        )
    ),

    sigmaPhiInit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaPhiInit",
            coeffDict_,
            0.33
        )
    ),

    sigmaPsiInit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaPsiInit",
            coeffDict_,
            1.0
        )
    ),

    psiNuFrac_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "psiNuFrac",
            coeffDict_,
            1.0
        )
    ),


   solveK_
   (
       coeffDict_.lookup("solveK")
   ),

   solveOmega_
   (
       coeffDict_.lookup("solveOmega")
   ),

   solveEps_
   (
       coeffDict_.lookup("solveEps")
   ),

   solvePsi_
   (
       coeffDict_.lookup("solvePsi")
   ),

   solvePhi_
   (
       coeffDict_.lookup("solvePhi")
   ),

   solveNut_
   (
       coeffDict_.lookup("solveNut")
   ),

   eqnSigmaK_
   (
       coeffDict_.lookup("eqnSigmaK")
   ),

   eqnSigmaEps_
   (
       coeffDict_.lookup("eqnSigmaEps")
   ),

   eqnSigmaPhi_
   (
       coeffDict_.lookup("eqnSigmaPhi")
   ),

   eqnSigmaPsi_
   (
       coeffDict_.lookup("eqnSigmaPsi")
   ),

   eqncEp2_
   (
       coeffDict_.lookup("eqncEp2")
   ),

   eqnEpsHat_
   (
       coeffDict_.lookup("eqnEpsHat")
   ),
   
   timeScaleEps_
   (
       coeffDict_.lookup("timeScaleEps")
   ),
   prodType_
   (
       coeffDict_.lookup("prodType")
   ),
   debugWrite_
   (
       coeffDict_.lookup("debugWrite")
   ),
   tslimiter_
   (
       coeffDict_.lookup("tslimiter")
   ),
   psiProd_
   (
       coeffDict_.lookup("psiProd")
   ),
    y_(mesh_),


    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    gradk_
    (
        IOobject
        (
            "gradk",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::grad(k_)
    ),
    omega_
    (
        IOobject
        (
            "omega",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        omega_*k_
    ),
    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    nutNorm_
    (
        IOobject
        (
            "nutNorm",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (nut_/max(nut_))
    ),
    tpphi_
    (
        IOobject
        (
            "tpphi",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    tpphisqrt_
    (
        IOobject
        (
            "tpphi",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt(tpphi_)
    ),
    vorticity_
    (
        IOobject
        (
            "vorticity",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::curl(U_)
    ),
	phis_
    (
        IOobject
        (
            "phis",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    vorticityTmp_
    (
        IOobject
        (
            "vorticityTmp",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::curl(U_)
    ),
    ivorticity_
    (
        IOobject
        (
            "ivorticity",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("iv", dimensionSet(0,0,1,0,0,0,0), vector(1.0,1.0,1.0))
    ),
    tppsi_
    (
        IOobject
        (
            "tppsi",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    uGrad_
    (
        IOobject
        (
            "uGrad",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::grad(U_)
    ),
	epsHat_
    (
        IOobject
        (
            "epsHat",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("epsHat", dimensionSet(0,0,-1,0,0,0,0), 1.0)
    ),
    kol_
    (
        IOobject
        (
            "kol",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("kol", dimensionSet(0,0,-1,0,0,0,0), 0.0)
    ),
    kSafe_
    (
        IOobject
        (
            "kSafe",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (k_)
    ),
    kSqrt_
    (
        IOobject
        (
            "kSqrt",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    gradkSqrt_
    (
        IOobject
        (
            "gradkSqrt",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::grad(kSqrt_)
    ),
    nutSafe_
    (
        IOobject
        (
            "nutSafe",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (nut_)
    ),
    epsilonSafe_
    (
        IOobject
        (
            "epsilonSafe",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (epsilon_)
    ),
    sigmaK_
    (
        IOobject
        (
            "sigmaK",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        sigmaKInit_
    ),   
    sigmaOm_
    (
        IOobject
        (
            "sigmaOm",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        sigmaOmInit_
    ),
    sigmaEps_
    (
        IOobject
        (
            "sigmaEps",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        sigmaEpsInit_
    ),
    sigmaPhi_
    (
        IOobject
        (
            "sigmaPhi",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        sigmaPhiInit_
    ),
    sigmaPsi_
    (
        IOobject
        (
            "sigmaPsi",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        sigmaPsiInit_
    ),
    cEp2_
    (
        IOobject
        (
            "cEp2",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (cEp2con_ - 0.16*exp(-0.1*sqr(k_)/(nu()*epsilon_)))
    ),
    tpProd_
    (
        IOobject
        (
            "tpProd",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        cPr_*(2*nut_*magSqr(symm(fvc::grad(U_)))/k_)
    ),
    cP1eqn_
    (
        IOobject
        (
            "cP1eqn",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        2.0*(0.5+0.5*((tpProd_*k_)/epsilon_))
    ),
    dimRat_
    (
        IOobject
        (
            "dimRat",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (psiReal() & psiReal())/(k_*phiReal())
    ),
    gradTpphi_
    (
        IOobject
        (
            "gradTpphi",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::grad(tpphi_)
    ),
    gradTppsi_
    (
        IOobject
        (
            "gradTppsi",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::grad(tppsi_)
    ),
    tpProdSqr_
    (
        IOobject
        (
            "tpProdSqr",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sqr(tppsi_ & vorticity_)
    ),
    tpProd3d_
    (
        IOobject
        (
            "tpProd3d",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (mag(psiReal() ^ vorticity_))
    )
{

    Info<< "Made it past constructors " << endl;

    // Calculate eddy viscosity
    if(solveNut_ == "true")
    {
		if(timeScaleEps_ == "epsilon" || timeScaleEps_ != "epsHat")
		{
            nut_ = cMu_*tpphi_*Ts();
        }
               
    }

    kSafe_ = max(k_, dimensionedScalar("minK", k_.dimensions(), 1.0e-15));

    Info<< "solveK is: " <<solveK_ <<endl;
    Info<< "solveEps is: " <<solveEps_ <<endl;
    Info<< "solvePhi is: " <<solvePhi_ <<endl;
    Info<< "solvePsi is: " <<solvePsi_ <<endl;

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Not used but necessary for RAS Model
tmp<volSymmTensorField> turbulentPotentialJoTomega::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            (
			  (2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)
			)
        )
    );
}

// Not used but necessary for RAS Model
tmp<volSymmTensorField> turbulentPotentialJoTomega::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}

// Term that is directly added to the momentum equation
tmp<fvVectorMatrix> turbulentPotentialJoTomega::divDevReff(volVectorField& U) const
{
    return
    (
       fvc::grad(tpphi_)
     + fvc::curl(tppsi_)
     + fvc::laplacian(nut_, U, "laplacian(nuEff,U)")
     - fvm::laplacian(nuEff(), U)
    );
}



bool turbulentPotentialJoTomega::read()
{
    if (RASModel::read())
    {
        cEp1_.readIfPresent(coeffDict());
        cEp2con_.readIfPresent(coeffDict());
        cEp3_.readIfPresent(coeffDict());
        cP1_.readIfPresent(coeffDict());
        cP2_.readIfPresent(coeffDict());
        cP3_.readIfPresent(coeffDict());
        cMu_.readIfPresent(coeffDict());
		cPphi_.readIfPresent(coeffDict());
		cEhm_.readIfPresent(coeffDict());
		cEhR_.readIfPresent(coeffDict());
		cPr_.readIfPresent(coeffDict());
		cD1_.readIfPresent(coeffDict());
		cD2_.readIfPresent(coeffDict());
        cVv1_.readIfPresent(coeffDict());
        cTv1_.readIfPresent(coeffDict());
		cT_.readIfPresent(coeffDict());
		sigmaKInit_.readIfPresent(coeffDict());
        sigmaEpsInit_.readIfPresent(coeffDict());
        sigmaEpsVisc_.readIfPresent(coeffDict());
        sigmaPhiInit_.readIfPresent(coeffDict());
		sigmaPsiInit_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void turbulentPotentialJoTomega::correct()
{

    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    if (mesh_.changing())
    {
        y_.correct();
        bound(k_, dimensionedScalar("minK", k_.dimensions(), 1.0e-15));
        bound(epsilon_, dimensionedScalar("minEps", epsilon_.dimensions(), 1.0e-15));
		bound(tpphi_,dimensionedScalar("minTpphi", tpphi_.dimensions(), 1.0e-15));
    }


    // Set the time scale using epsilon 
    volScalarField T("TimeScale",Ts());
		
    // Bound and Display log output if production is negative
    bound(T, dimensionedScalar("minT", T.dimensions(), 1.0e-15));
        	
    // Vorticity
    vorticity_ = fvc::curl(U_);
	uGrad_ = fvc::grad(U_);
    Info<< "Made it past vorticity" << endl;
	
	// Production Terms
	volScalarField S2 = magSqr(dev(symm(uGrad_)));
    volScalarField G("RASModel::G", nut_*2*S2);
    Info<< "Made it production" << endl;

	// Sigma equations 
    if(eqnSigmaK_ == "true")
    {
	    sigmaK_ = 0.67 + 0.33*(tpProd_/epsHat_);
	}

    if(eqnSigmaEps_ == "true")
    {
	    sigmaEps_ = 0.33 + 0.5*(tpProd_/epsHat_);
    }

    if(eqnSigmaPhi_ == "true")
    {
	    sigmaPhi_ = 0.21 + 0.12*(tpProd_/epsHat_);
	}

    if(eqnSigmaPsi_ == "true")
    {
	    sigmaEps_ = 0.33 + 0.4*(tpProd_/epsHat_);
    }

    epsilonSafe_ = max(epsilon_, dimensionedScalar("minEps", epsilon_.dimensions(), 1.0e-15));

    if(eqncEp2_ == "true")
    {
        cEp2_ = cEp2con_ - 0.16*exp(-0.25*sqr(k_)/(nu()*epsilonSafe_));
    }
    else
    {
        cEp2_ =  cEp2con_;
    }

    cP1eqn_ = 2.0*(0.5+0.5*((tpProd_*k_)/epsilonSafe_));


    //Dissipation equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(omega_)
      + fvm::div(phi_, omega_)
      + fvm::SuSp(-fvc::div(phi_), omega_)
      - fvm::laplacian(DomegaEff(), omega_)
     ==
       (cEp1_-1.0)*G/(k_*Ts())
     + fvm::Sp(-1.0*(cEp2_-1.0)*omega_,omega_)
    );

    if(solveOmega_ == "true")
    {
    omegaEqn().relax();
    solve(omegaEqn);
    bound(omega_,dimensionedScalar("minOmega", epsilon_.dimensions(), 1.0e-15));
    }


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (

        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      + fvm::SuSp(-fvc::div(phi_), k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G
      - fvm::Sp(omega_,k_)
    );


    if(solveK_ == "true")
    {
    kEqn().relax();
    solve(kEqn);
    bound(k_,dimensionedScalar("minK", k_.dimensions(), 1.0e-15));
    }

    kSqrt_ = sqrt(k_);
    bound(kSqrt_,dimensionedScalar("minKsqrt", kSqrt_.dimensions(), 3.16e-8));
    kSqrt_.correctBoundaryConditions();

    gradk_ = fvc::grad(k_);
    gradkSqrt_ = fvc::grad(kSqrt_);
    tpphisqrt_ = sqrt(tpphi_);
    volVectorField gradPhiSqrt("gradPhiSqrt", fvc::grad(tpphisqrt_));
       
    volScalarField GdK("GdK", G/k_);
    
    volVectorField gradPhiOverK("gradPhiOverK",fvc::grad(PhiOverK()));

    // Phi equation
    tmp<fvScalarMatrix> tpphiEqn
    (
        fvm::ddt(tpphi_)
      + fvm::div(phi_, tpphi_)
      + fvm::SuSp(-fvc::div(phi_), tpphi_)
      - fvm::laplacian(DphiEff(), tpphi_)
      ==
        //3.0*cP1_*nutFrac()*epsHat_*(1.0 - Alpha())*(2.0*k_/3.0 - tpphi_) 
        (1.0 + cP1_)*2*nutFrac()*omega_*(2*Alpha() - 1.0)*tpphi_		
      + fvm::Sp(-1.0*(cP2_+cP4_)*GdK,tpphi_) 
      + (cP2_+cP4_)*Alpha()*((tppsi_ & tppsi_)/((nut_*k_)*(1.0+cPw_/reTau())))*tpphi_
      + cP2_*GdK*tpphi_
      + fvm::Sp(-2.0*Alpha()*omega_,tpphi_)
      + gT1_*fvm::Sp(-1.0*nut_*(gradk_ & gradPhiOverK)/tpphi_,tpphi_) 	  
      + gT2_*fvm::Sp(-2.0*nu()*(gradPhiSqrt & gradPhiSqrt)/tpphi_,tpphi_)
      + cT_*(1.5*tpphi_-k_)*GdK*sqrt((nut_/nu()))
    );

    if(solvePhi_ == "true")
    {
    tpphiEqn().relax();
    solve(tpphiEqn);
    bound(tpphi_,dimensionedScalar("minTpphi", tpphi_.dimensions(), 1.0e-15));
    }

    gradTpphi_ = fvc::grad(tpphi_);

	volTensorField gradPsiOverK("gradPsiOverK",fvc::grad(PsiOverK()));
	
	
	// Choice of production term for psi equation
	volScalarField psProd("psProd",GdK);
	
	if(psiProd_ == "psi")
	{
		volScalarField psProd("psProd",(tppsi_ & vorticity_)/k_);
	}


    // Psi Equation
    tmp<fvVectorMatrix> tppsiEqn
    (
        fvm::ddt(tppsi_)
      + fvm::div(phi_, tppsi_)
      + fvm::Sp(-fvc::div(phi_), tppsi_)
      - fvm::laplacian(DpsiEff(), tppsi_)

      ==
        (1.0 - cP2_)*tpphi_*vorticity_
      - (2*Alpha() - cP2_)*psProd*tppsi_ 
      + 0.3*(2*Alpha() - 1.0)*tpphi_*vorticity_
      - gT1_*2.0*nut_*(gradk_ & gradPsiOverK)     
	  + fvm::Sp(-1.0*(cP1_*nutFrac()*(1.0-Alpha()))*omega_,tppsi_)   
      + fvm::Sp(-1.0*Alpha()*omega_,tppsi_)
      + fvm::Sp(-1.0*omega_/((1+cP3_*(nut_/nu()))),tppsi_)
      + gT2_*fvm::Sp(-2.0*nu()*(gradkSqrt_ & gradPhiSqrt)/(sqrt(k_*tpphi_)),tppsi_)
      + cT_*sqrt((nut_/nu()))*vorticity_*k_
    );

    if(solvePsi_ == "true")
    {
    tppsiEqn().relax();
    solve(tppsiEqn);
    }

    gradTppsi_ = fvc::grad(tppsi_);
    
    
    volScalarField psiZ(tppsi_.component(2));
	
	
	epsilon_ = omega_*k_;
	
	    // Calculate eddy viscosity
    if(solveNut_ == "true")
    {
        if(timeScaleEps_ == "epsilon" || timeScaleEps_ != "epsHat")
		{
            nut_ = cMu_*tpphi_*Ts();
            nut_.correctBoundaryConditions();
        }
        
    }
		
	
	Info<< "Maximum nut: " << gMax(nut_) << " Maximum K: " << gMax(k_) << " Maximum Epsilon: " << gMax(epsilon_) <<endl;
    Info<< "Maximum Phi: " << gMax(tpphi_) << " Maximum Psi_z: " << gMax(psiZ) << " Maximum Production: " << gMax(G) <<endl;
    Info<< "Maximum Omega: " << gMax(omega_) <<endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
