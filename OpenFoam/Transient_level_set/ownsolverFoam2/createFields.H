Info << "Reading field G\n" << endl;

volScalarField G
(
    IOobject
    (
	"G",
	runTime.timeName(),
	mesh,
	IOobject::MUST_READ,
	IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField G0
(
    IOobject
    (
	"G0",
	runTime.timeName(),
	mesh,
	IOobject::MUST_READ,
	IOobject::AUTO_WRITE
    ),
    G
);
