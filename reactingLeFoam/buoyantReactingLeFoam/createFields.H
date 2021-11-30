#include "createRDeltaT.H"

Info<< "Reading thermophysical properties\n" << endl;
autoPtr<fluidReactionThermo> pThermo(fluidReactionThermo::New(mesh));
fluidReactionThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

// Added code
const word inertSpecie(thermo.properties().lookup("inertSpecie"));
if (!composition.species().found(inertSpecie))
{
    FatalIOErrorIn(args.executable().c_str(), thermo.properties())
        << "Inert specie " << inertSpecie << " not found in available species "
        << composition.species() << exit(FatalIOError);
}
// ----------

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh
    ),
    thermo.rho()
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


volScalarField& p = thermo.p();

#include "compressibleCreatePhi.H"

mesh.setFluxRequired(p.name());


Info << "Creating turbulence model.\n" << nl;
autoPtr<compressible::momentumTransportModel> turbulence
(
    compressible::momentumTransportModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

Info<< "Creating thermophysical transport model\n" << endl;
autoPtr<fluidReactionThermophysicalTransportModel>
thermophysicalTransport
(
    fluidReactionThermophysicalTransportModel::New
    (
        turbulence(),
        thermo
    )
);

Info<< "Creating reaction model\n" << endl;
autoPtr<combustionModel> reaction(combustionModel::New(thermo, turbulence()));

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"
#include "readpRef.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

pressureReference pressureReference
(
    p,
    p_rgh,
    pimple.dict(),
    thermo.incompressible()
);

mesh.setFluxRequired(p_rgh.name());

hydrostaticInitialisation
(
    p_rgh,
    p,
    rho,
    U,
    gh,
    ghf,
    pRef,
    thermo,
    pimple.dict()
);

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

dimensionedScalar initialMass = fvc::domainIntegrate(rho);


multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());


// Added code

volVectorField heatFluxCorr
(
    IOobject
    (
        "heatFluxCorr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("heatFluxCorr", thermo.he().dimensions()/dimLength, Foam::vector(0,0,0))
);

volVectorField gradY
(
    IOobject
    (
        "gradY",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("gradY", dimless/dimLength, Foam::vector(0,0,0))
);


dimensionedScalar Hs("Hs",thermo.he().dimensions(),0.0);
label Ysize = Y.size();

scalarList LewisNumber(Ysize,1.0);

Info<< "\nReading Lewis numbers\n" << endl;
IOdictionary chemistryProperties
(
    IOobject
    (
        "chemistryProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);


dictionary LewisNumberDict(chemistryProperties.subDict("LewisNumbers"));


for (label i=0; i<Ysize; i++)
{
    if (LewisNumberDict.found(Y[i].member()))
    {
        LewisNumber[i] = readScalar(LewisNumberDict.lookup(Y[i].member()));
    }
    Info<<"Lewis number of specie "<<Y[i].name()<<" is: "<<LewisNumber[i]<<endl;
}

// -----------


#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"
