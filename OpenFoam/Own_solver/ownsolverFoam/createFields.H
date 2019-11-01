Info<< "Reading field T\n" << endl;

volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);


Info<< "Reading diffusivity DT\n" << endl;

dimensionedScalar DT
(
    transportProperties.lookup("DT")
);


Info<< "Reading density rho\n" << endl;

dimensionedScalar rho
(
    transportProperties.lookup("rho")
);

Info<< "Reading specific heat c\n" << endl;

dimensionedScalar c
(
    transportProperties.lookup("c")
);

Info<< "Reading power P\n" << endl;

dimensionedScalar P
(
    transportProperties.lookup("P")
);

Info<< "Reading volume affected by the source V\n" << endl;

dimensionedScalar V
(
    transportProperties.lookup("V")
);

#include "createFvOptions.H"
