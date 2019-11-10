singlePhaseTransportModel laminarTransport(U, phi);

// Thermal expansion coefficient [1/K]
dimensionedScalar beta
(
    "beta",
    dimless/dimTemperature,
    laminarTransport
);

// Reference temperature [K]
dimensionedScalar TRef("TRef", dimTemperature, laminarTransport);

// Laminar Prandtl number
dimensionedScalar Pr("Pr", dimless, laminarTransport);

// Turbulent Prandtl number
dimensionedScalar Prt("Prt", dimless, laminarTransport);

// Density
dimensionedScalar rho("rho", dimDensity, laminarTransport);

// Specific Heat Capacity
dimensionedScalar c("c", (dimLength*dimLength)/(dimTime*dimTime)/dimTemperature, laminarTransport);

// Thermal Conductivity
dimensionedScalar K("K", dimMass*dimLength/(dimTime*dimTime*dimTime)/dimTemperature, laminarTransport);
