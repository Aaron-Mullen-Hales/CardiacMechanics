/*---------------------------------------------------------------------------*\
  Case-local Aróstica benchmark time-dependent endocardial pressure
  boundary condition.
\*---------------------------------------------------------------------------*/

#include "arosticaActivationPressureTractionFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "IStringStream.H"

#include <cmath>

namespace Foam
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

dictionary arosticaActivationPressureTractionFvPatchVectorField::withDefaults
(
    const dictionary& input
)
{
    dictionary result(input);
    if (!result.found("traction"))
    {
        IStringStream stream("traction uniform (0 0 0);");
        entry::New(result, stream);
    }
    if
    (
        !result.found("pressure")
     && !result.found("pressureSeries")
     && !result.found("pressureField")
    )
    {
        IStringStream stream("pressure uniform 0;");
        entry::New(result, stream);
    }
    return result;
}


void arosticaActivationPressureTractionFvPatchVectorField::
validateParameters() const
{
    const dimensionSet rateDims(dimless/dimTime);

    if
    (
        alphaMin_.dimensions() != rateDims || !std::isfinite(alphaMin_.value())
    )
    {
        FatalErrorInFunction
            << "alphaMin must have dimensions of inverse time and be "
            << "finite; got " << alphaMin_ << exit(FatalError);
    }

    if
    (
        alphaMax_.dimensions() != rateDims || !std::isfinite(alphaMax_.value())
    )
    {
        FatalErrorInFunction
            << "alphaMax must have dimensions of inverse time and be "
            << "finite; got " << alphaMax_ << exit(FatalError);
    }

    if
    (
        alphaPre_.dimensions() != rateDims || !std::isfinite(alphaPre_.value())
    )
    {
        FatalErrorInFunction
            << "alphaPre must have dimensions of inverse time and be "
            << "finite; got " << alphaPre_ << exit(FatalError);
    }

    if
    (
        alphaMid_.dimensions() != rateDims || !std::isfinite(alphaMid_.value())
    )
    {
        FatalErrorInFunction
            << "alphaMid must have dimensions of inverse time and be "
            << "finite; got " << alphaMid_ << exit(FatalError);
    }

    if
    (
        sigmaPre_.dimensions() != dimPressure
     || !std::isfinite(sigmaPre_.value())
     || sigmaPre_.value() < 0.0
    )
    {
        FatalErrorInFunction
            << "sigmaPre must have dimensions of pressure and be finite "
            << "and non-negative; got " << sigmaPre_ << exit(FatalError);
    }

    if
    (
        sigmaMid_.dimensions() != dimPressure
     || !std::isfinite(sigmaMid_.value())
     || sigmaMid_.value() < 0.0
    )
    {
        FatalErrorInFunction
            << "sigmaMid must have dimensions of pressure and be finite "
            << "and non-negative; got " << sigmaMid_ << exit(FatalError);
    }

    if
    (
        tSysPre_.dimensions() != dimTime
     || !std::isfinite(tSysPre_.value())
     || tSysPre_.value() < 0.0
    )
    {
        FatalErrorInFunction
            << "tSysPre must have dimensions of time and be finite and "
            << "non-negative; got " << tSysPre_ << exit(FatalError);
    }

    if
    (
        tDiasPre_.dimensions() != dimTime
     || !std::isfinite(tDiasPre_.value())
     || tDiasPre_.value() <= tSysPre_.value()
    )
    {
        FatalErrorInFunction
            << "tDiasPre must have dimensions of time, be finite and be "
            << "greater than tSysPre; got tDiasPre = " << tDiasPre_
            << ", tSysPre = " << tSysPre_ << exit(FatalError);
    }

    if
    (
        gamma_.dimensions() != dimTime
     || !std::isfinite(gamma_.value())
     || gamma_.value() <= 0.0
    )
    {
        FatalErrorInFunction
            << "gamma must have dimensions of time and be finite and "
            << "positive; got " << gamma_ << exit(FatalError);
    }
}


scalar arosticaActivationPressureTractionFvPatchVectorField::sPlus
(
    const scalar deltaT
) const
{
    return 0.5*(1.0 + std::tanh(deltaT/gamma_.value()));
}


scalar arosticaActivationPressureTractionFvPatchVectorField::sMinus
(
    const scalar deltaT
) const
{
    return 0.5*(1.0 - std::tanh(deltaT/gamma_.value()));
}


scalar arosticaActivationPressureTractionFvPatchVectorField::gPre
(
    const scalar t
) const
{
    // Eq. (8): gPre(t) = S-(t - tDiasPre)
    return sMinus(t - tDiasPre_.value());
}


scalar arosticaActivationPressureTractionFvPatchVectorField::
activationFunction(const scalar t) const
{
    // Eq. (8): b(t) = aPre(t) + alphaPre gPre(t) + alphaMid
    const scalar fPre = sPlus(t - tSysPre_.value())*sMinus(t - tDiasPre_.value());
    const scalar aPre = alphaMax_.value()*fPre + alphaMin_.value()*(1.0 - fPre);

    return aPre + alphaPre_.value()*gPre(t) + alphaMid_.value();
}


scalar arosticaActivationPressureTractionFvPatchVectorField::pDerivative
(
    const scalar t,
    const scalar p
) const
{
    // Eq. (7): dp/dt = -|b(t)| p + sigmaMid |b(t)|_+ + sigmaPre |gPre(t)|_+
    const scalar b = activationFunction(t);
    const scalar g = gPre(t);

    return -mag(b)*p + sigmaMid_.value()*max(b, 0.0) + sigmaPre_.value()*max(g, 0.0);
}


scalar arosticaActivationPressureTractionFvPatchVectorField::integrateP
(
    const scalar tEnd
) const
{
    if (tEnd <= 0.0)
    {
        return 0.0;
    }

    // p(t) is spatially uniform and depends only on time, so it is
    // integrated here with a fixed sub-step fine enough to resolve the
    // gamma_ transition time scale, independent of the outer solver time
    // step
    const scalar subStep = min(gamma_.value()/20.0, 1.0e-4);

    const label nSteps = max(label(tEnd/subStep) + 1, 1);

    if (nSteps > 20000000)
    {
        FatalErrorInFunction
            << "Excessive number of sub-steps (" << nSteps << ") required "
            << "to integrate the Aróstica activation pressure to t = "
            << tEnd << "; check the gamma parameter" << exit(FatalError);
    }

    const scalar h = tEnd/nSteps;

    scalar t = 0.0;
    scalar p = 0.0;

    for (label stepI = 0; stepI < nSteps; ++stepI)
    {
        const scalar k1 = pDerivative(t, p);
        const scalar k2 = pDerivative(t + 0.5*h, p + 0.5*h*k1);
        const scalar k3 = pDerivative(t + 0.5*h, p + 0.5*h*k2);
        const scalar k4 = pDerivative(t + h, p + h*k3);

        p += (h/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4);
        t += h;

        if (!std::isfinite(p))
        {
            FatalErrorInFunction
                << "Non-finite pressure p = " << p << " encountered while "
                << "integrating the Aróstica activation pressure model at "
                << "t = " << t << exit(FatalError);
        }
    }

    return p;
}


scalar arosticaActivationPressureTractionFvPatchVectorField::
currentPressure() const
{
    const label timeIndex = db().time().timeIndex();

    if (timeIndex != cachedTimeIndex_)
    {
        p_ = integrateP(db().time().value());
        cachedTimeIndex_ = timeIndex;
    }

    return p_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

arosticaActivationPressureTractionFvPatchVectorField::
arosticaActivationPressureTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(p, iF),
    alphaMin_("alphaMin", dimless/dimTime, 0.0),
    alphaMax_("alphaMax", dimless/dimTime, 0.0),
    alphaPre_("alphaPre", dimless/dimTime, 0.0),
    alphaMid_("alphaMid", dimless/dimTime, 0.0),
    sigmaPre_("sigmaPre", dimPressure, 0.0),
    sigmaMid_("sigmaMid", dimPressure, 0.0),
    tSysPre_("tSysPre", dimTime, 0.0),
    tDiasPre_("tDiasPre", dimTime, 1.0),
    gamma_("gamma", dimTime, 1.0),
    cachedTimeIndex_(-1),
    p_(0.0)
{}


arosticaActivationPressureTractionFvPatchVectorField::
arosticaActivationPressureTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    solidTractionFvPatchVectorField(p, iF, withDefaults(dict)),
    // Read with the (name, dimensions, dict) constructor rather than
    // dict.lookup(). lookup() hands the value stream to
    // dimensionedScalar(Istream&), which expects the name to still be in the
    // stream and therefore rejects the dimensioned form
    //
    //     alphaMin  [0 0 -1 0 0 0 0]  -30;
    //
    // that write() itself emits -- so a case could be written but never read
    // back, which blocked every restart. This form accepts both the plain
    // and the dimensioned spelling and validates the dimensions on the way
    // in. The values obtained are identical; nothing in the residual changes.
    alphaMin_("alphaMin", dimless/dimTime, dict),
    alphaMax_("alphaMax", dimless/dimTime, dict),
    alphaPre_("alphaPre", dimless/dimTime, dict),
    alphaMid_("alphaMid", dimless/dimTime, dict),
    sigmaPre_("sigmaPre", dimPressure, dict),
    sigmaMid_("sigmaMid", dimPressure, dict),
    tSysPre_("tSysPre", dimTime, dict),
    tDiasPre_("tDiasPre", dimTime, dict),
    gamma_("gamma", dimTime, dict),
    cachedTimeIndex_(-1),
    p_(0.0)
{
    validateParameters();

    Info<< "Creating " << type() << " on " << patch().name() << nl
        << "    alphaMin = " << alphaMin_ << nl
        << "    alphaMax = " << alphaMax_ << nl
        << "    alphaPre = " << alphaPre_ << nl
        << "    alphaMid = " << alphaMid_ << nl
        << "    sigmaPre = " << sigmaPre_ << nl
        << "    sigmaMid = " << sigmaMid_ << nl
        << "    tSysPre = " << tSysPre_ << nl
        << "    tDiasPre = " << tDiasPre_ << nl
        << "    gamma = " << gamma_ << endl;
}


arosticaActivationPressureTractionFvPatchVectorField::
arosticaActivationPressureTractionFvPatchVectorField
(
    const arosticaActivationPressureTractionFvPatchVectorField& pvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidTractionFvPatchVectorField(pvf, p, iF, mapper),
    alphaMin_(pvf.alphaMin_),
    alphaMax_(pvf.alphaMax_),
    alphaPre_(pvf.alphaPre_),
    alphaMid_(pvf.alphaMid_),
    sigmaPre_(pvf.sigmaPre_),
    sigmaMid_(pvf.sigmaMid_),
    tSysPre_(pvf.tSysPre_),
    tDiasPre_(pvf.tDiasPre_),
    gamma_(pvf.gamma_),
    cachedTimeIndex_(-1),
    p_(0.0)
{}


#ifndef OPENFOAM_ORG
arosticaActivationPressureTractionFvPatchVectorField::
arosticaActivationPressureTractionFvPatchVectorField
(
    const arosticaActivationPressureTractionFvPatchVectorField& pvf
)
:
    solidTractionFvPatchVectorField(pvf),
    alphaMin_(pvf.alphaMin_),
    alphaMax_(pvf.alphaMax_),
    alphaPre_(pvf.alphaPre_),
    alphaMid_(pvf.alphaMid_),
    sigmaPre_(pvf.sigmaPre_),
    sigmaMid_(pvf.sigmaMid_),
    tSysPre_(pvf.tSysPre_),
    tDiasPre_(pvf.tDiasPre_),
    gamma_(pvf.gamma_),
    cachedTimeIndex_(pvf.cachedTimeIndex_),
    p_(pvf.p_)
{}
#endif


arosticaActivationPressureTractionFvPatchVectorField::
arosticaActivationPressureTractionFvPatchVectorField
(
    const arosticaActivationPressureTractionFvPatchVectorField& pvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(pvf, iF),
    alphaMin_(pvf.alphaMin_),
    alphaMax_(pvf.alphaMax_),
    alphaPre_(pvf.alphaPre_),
    alphaMid_(pvf.alphaMid_),
    sigmaPre_(pvf.sigmaPre_),
    sigmaMid_(pvf.sigmaMid_),
    tSysPre_(pvf.tSysPre_),
    tDiasPre_(pvf.tDiasPre_),
    gamma_(pvf.gamma_),
    cachedTimeIndex_(pvf.cachedTimeIndex_),
    p_(pvf.p_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void arosticaActivationPressureTractionFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    traction() = vector::zero;
    pressure() = currentPressure();

    solidTractionFvPatchVectorField::updateCoeffs();
}


void arosticaActivationPressureTractionFvPatchVectorField::write
(
    Ostream& os
) const
{
    solidTractionFvPatchVectorField::write(os);

    alphaMin_.writeEntry("alphaMin", os);
    alphaMax_.writeEntry("alphaMax", os);
    alphaPre_.writeEntry("alphaPre", os);
    alphaMid_.writeEntry("alphaMid", os);
    sigmaPre_.writeEntry("sigmaPre", os);
    sigmaMid_.writeEntry("sigmaMid", os);
    tSysPre_.writeEntry("tSysPre", os);
    tDiasPre_.writeEntry("tDiasPre", os);
    gamma_.writeEntry("gamma", os);
}


makePatchTypeField
(
    fvPatchVectorField,
    arosticaActivationPressureTractionFvPatchVectorField
);

} // End namespace Foam
