surfaceVectorField gradGf = fvc::interpolate(fvc::grad(G));
surfaceScalarField w = sign(fvc::interpolate(G))*(mesh.Sf() & gradGf)/(mag(gradGf) + SMALL);

solve
(
fvm::ddt(G) + fvc::div(w,G) + fvc::SuSp(-fvc::div(w), G) - sign(G0)
);
