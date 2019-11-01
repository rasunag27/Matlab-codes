volScalarField x = mesh.C().component(vector::X);
volScalarField y = mesh.C().component(vector::Y);
volScalarField r = sqrt(x*x + y*y);
scalar R = 1;
forAll(r,icell)
{
if (r[icell] > R) G[icell] = 1;
else G[icell] = -1;
}
G0=G;
