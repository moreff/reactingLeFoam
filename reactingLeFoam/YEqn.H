tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.schemes().div("div(phi,Yi_h)")
    )
);

reaction->correct();
volScalarField Yt(0.0*Y[0]);

forAll(Y, i)
{
    if (i != inertIndex && composition.active(i))
    {
        volScalarField& Yi = Y[i];

        fvScalarMatrix YiEqn
        (
    	fvm::ddt(rho, Yi)
          + mvConvection->fvmDiv(phi, Yi)
          - fvm::laplacian(thermophysicalTransport->alphaEff()/LewisNumber[i], Yi)
         ==
    	reaction->R(Yi)
          + fvModels.source(rho, Yi)
        );

        YiEqn.relax();

        fvConstraints.constrain(YiEqn);

        YiEqn.solve("Yi");

        fvConstraints.constrain(Yi);

        Yi.max(0.0);
        Yt += Yi;
    }
}

Y[inertIndex] = scalar(1) - Yt;
Y[inertIndex].max(0.0);

// composition.normalise();
