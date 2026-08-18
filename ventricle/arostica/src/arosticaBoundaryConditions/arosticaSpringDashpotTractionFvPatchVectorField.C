/*---------------------------------------------------------------------------*\
  Case-local Aróstica benchmark spring/dashpot traction conditions.
\*---------------------------------------------------------------------------*/

#include "arosticaSpringDashpotTractionFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "IStringStream.H"

#include <cmath>

namespace Foam
{

dictionary arosticaSpringDashpotTractionFvPatchVectorField::withDefaults
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
    if (!result.found("pressure"))
    {
        IStringStream stream("pressure uniform 0;");
        entry::New(result, stream);
    }
    return result;
}


dimensionedScalar
arosticaSpringDashpotTractionFvPatchVectorField::readCoefficient
(
    const dictionary& dict,
    const word& dimensionedName,
    const word& legacyName,
    const dimensionSet& expectedDimensions
)
{
    dimensionedScalar result(dimensionedName, expectedDimensions, 0.0);
    if (dict.found(dimensionedName))
    {
        // Same defect as the activation-pressure BC: dict.lookup() feeds the
        // value stream to dimensionedScalar(Istream&), which expects the name
        // still to be present and so cannot read back the dimensioned form
        // that write() emits. The (name, dimensions, dict) constructor reads
        // both spellings and the value obtained is identical.
        const dimensionedScalar input(dimensionedName, expectedDimensions, dict);
        if (input.dimensions() != expectedDimensions)
        {
            FatalIOErrorInFunction(dict)
                << dimensionedName << " has dimensions " << input.dimensions()
                << "; expected " << expectedDimensions << exit(FatalIOError);
        }
        result = dimensionedScalar
        (
            dimensionedName, expectedDimensions, input.value()
        );
    }
    else if (dict.found(legacyName))
    {
        result.value() = readScalar(dict.lookup(legacyName));
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Missing " << dimensionedName << " (legacy alias "
            << legacyName << ')' << exit(FatalIOError);
    }
    if (!std::isfinite(result.value()) || result.value() < 0.0)
    {
        FatalIOErrorInFunction(dict)
            << dimensionedName << " must be finite and non-negative"
            << exit(FatalIOError);
    }
    return result;
}


void arosticaSpringDashpotTractionFvPatchVectorField::requireReferenceArea
(
    const dictionary& dict
)
{
    if (!dict.lookupOrDefault<Switch>("useUndeformedArea", false))
    {
        FatalIOErrorInFunction(dict)
            << "Aróstica support traction requires explicit "
            << "useUndeformedArea true" << exit(FatalIOError);
    }
}


arosticaSpringDashpotTractionFvPatchVectorField::
arosticaSpringDashpotTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(p, iF),
    springCoefficient_("springCoefficient", dimPressure/dimLength, 0.0),
    dashpotCoefficient_
    (
        "dashpotCoefficient", dimPressure*dimTime/dimLength, 0.0
    ),
    writeDiagnostics_(false),
    displacementMode_("ownerCell")
{}


arosticaSpringDashpotTractionFvPatchVectorField::
arosticaSpringDashpotTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    solidTractionFvPatchVectorField(p, iF, withDefaults(dict)),
    springCoefficient_
    (
        readCoefficient
        (
            dict, "springCoefficient", "alpha", dimPressure/dimLength
        )
    ),
    dashpotCoefficient_
    (
        readCoefficient
        (
            dict,
            "dashpotCoefficient",
            "beta",
            dimPressure*dimTime/dimLength
        )
    ),
    writeDiagnostics_(dict.lookupOrDefault<Switch>("writeDiagnostics", false)),
    displacementMode_
    (
        dict.lookupOrDefault<word>("supportDisplacementMode", "ownerCell")
    )
{
    requireReferenceArea(dict);
    Info<< "Creating " << type() << " on " << patch().name() << nl
        << "    springCoefficient = " << springCoefficient_ << nl
        << "    dashpotCoefficient = " << dashpotCoefficient_ << nl
        << "    useUndeformedArea = true" << endl;
}


arosticaSpringDashpotTractionFvPatchVectorField::
arosticaSpringDashpotTractionFvPatchVectorField
(
    const arosticaSpringDashpotTractionFvPatchVectorField& pvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidTractionFvPatchVectorField(pvf, p, iF, mapper),
    springCoefficient_(pvf.springCoefficient_),
    dashpotCoefficient_(pvf.dashpotCoefficient_),
    writeDiagnostics_(pvf.writeDiagnostics_),
    displacementMode_(pvf.displacementMode_)
{}

#ifndef OPENFOAM_ORG
arosticaSpringDashpotTractionFvPatchVectorField::
arosticaSpringDashpotTractionFvPatchVectorField
(
    const arosticaSpringDashpotTractionFvPatchVectorField& pvf
)
:
    solidTractionFvPatchVectorField(pvf),
    springCoefficient_(pvf.springCoefficient_),
    dashpotCoefficient_(pvf.dashpotCoefficient_),
    writeDiagnostics_(pvf.writeDiagnostics_),
    displacementMode_(pvf.displacementMode_)
{}
#endif

arosticaSpringDashpotTractionFvPatchVectorField::
arosticaSpringDashpotTractionFvPatchVectorField
(
    const arosticaSpringDashpotTractionFvPatchVectorField& pvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(pvf, iF),
    springCoefficient_(pvf.springCoefficient_),
    dashpotCoefficient_(pvf.dashpotCoefficient_),
    writeDiagnostics_(pvf.writeDiagnostics_),
    displacementMode_(pvf.displacementMode_)
{}


void arosticaSpringDashpotTractionFvPatchVectorField::updateCoeffs()
{
    if (updated()) return;

    const volVectorField& D = db().lookupObject<volVectorField>("D");
    const label patchI = patch().index();

    // The fixed-gradient patch value may belong to the preceding PETSc
    // trial. Use the adjacent cell value consistently for the spring and
    // dashpot, including all old-time levels in the ddt scheme.
    Field<vector>::operator=(patchInternalField());

    const tmp<volVectorField> tDdot(fvc::ddt(D));

    vectorField displacement(patchInternalField());
    vectorField velocity
    (
        tDdot().boundaryField()[patchI].patchInternalField()
    );

    if (displacementMode_ == "cellReconstructedFace")
    {
        // Reconstruct the SURFACE displacement the benchmark's Robin
        // condition actually refers to, from CELL data only:
        //
        //     D_f = D_c + (Cf - Cc) & grad(D)_c
        //
        // grad() here is the case's leastSquaresS4f scheme, which sets
        // useBoundaryFaceValues = false on every patch, so the gradient is
        // built from owner and internal-neighbour cell values and never
        // inherits this patch's own fixedGradient value.
        //
        // The velocity is reconstructed with the SAME operator applied to the
        // cell-based ddt field rather than by differencing reconstructed
        // levels. The two are algebraically identical: reconstruction is
        // linear in the field, the BDF2 coefficients are constants and the
        // mesh is fixed, so
        //
        //     rec(ddt D) = ddt(rec D).
        //
        // Every quantity below is therefore cell-derived at every time level;
        // no old-time boundary values enter, which is what the Change-1 fix
        // required.
        const vectorField delta(patch().delta());

        const volTensorField gradD(fvc::grad(D));
        const tensorField gradDc(gradD.boundaryField()[patchI].patchInternalField());
        displacement += (delta & gradDc);

        const volTensorField gradDdot(fvc::grad(tDdot()));
        const tensorField gradDdotc
        (
            gradDdot.boundaryField()[patchI].patchInternalField()
        );
        velocity += (delta & gradDdotc);
    }
    else if (displacementMode_ != "ownerCell")
    {
        FatalErrorInFunction
            << "Unknown supportDisplacementMode " << displacementMode_
            << "; expected ownerCell or cellReconstructedFace"
            << exit(FatalError);
    }

    vectorField spring(patch().size(), vector::zero);
    vectorField dashpot(patch().size(), vector::zero);
    calculateContributions(spring, dashpot, displacement, velocity);
    traction() = spring + dashpot;
    pressure() = scalar(0.0);

    if (writeDiagnostics_)
    {
        const vectorField displacementOld
        (
            D.oldTime().boundaryField()[patchI].patchInternalField()
        );
        const vectorField displacementOldOld
        (
            D.oldTime().oldTime().boundaryField()[patchI]
                .patchInternalField()
        );
        const vectorField mixedMeasureVelocity
        (
            tDdot().boundaryField()[patchI]
        );
        const scalarField dashpotPower(dashpot & velocity);

        Info<< "Aróstica support diagnostics, patch " << patch().name()
            << ": max|Dcell|=" << max(mag(displacement))
            << ", max|DcellOld|=" << max(mag(displacementOld))
            << ", max|DcellOldOld|=" << max(mag(displacementOldOld))
            << ", max|DdotCell|=" << max(mag(velocity))
            << ", max|DdotMixed-DdotCell|="
            << max(mag(mixedMeasureVelocity - velocity))
            << ", max|dashpotTraction|=" << max(mag(dashpot))
            << ", dashpotPowerRange=[" << min(dashpotPower)
            << ',' << max(dashpotPower) << ']'
            << ", max|traction|=" << max(mag(traction())) << endl;
    }
    solidTractionFvPatchVectorField::updateCoeffs();
}


void arosticaSpringDashpotTractionFvPatchVectorField::write(Ostream& os) const
{
    solidTractionFvPatchVectorField::write(os);
    dimensionedScalar springOut
    (
        "springCoefficientValue",
        springCoefficient_.dimensions(),
        springCoefficient_.value()
    );
    dimensionedScalar dashpotOut
    (
        "dashpotCoefficientValue",
        dashpotCoefficient_.dimensions(),
        dashpotCoefficient_.value()
    );
    os.writeEntry("supportDisplacementMode", displacementMode_);
    springOut.writeEntry("springCoefficient", os);
    dashpotOut.writeEntry("dashpotCoefficient", os);
    os.writeKeyword("writeDiagnostics")
        << writeDiagnostics_ << token::END_STATEMENT << nl;
}


arosticaNormalSpringDashpotTractionFvPatchVectorField::
arosticaNormalSpringDashpotTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    arosticaSpringDashpotTractionFvPatchVectorField(p, iF)
{}

arosticaNormalSpringDashpotTractionFvPatchVectorField::
arosticaNormalSpringDashpotTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    arosticaSpringDashpotTractionFvPatchVectorField(p, iF, dict)
{
    validateTangentialTraction(dict);
}

arosticaNormalSpringDashpotTractionFvPatchVectorField::
arosticaNormalSpringDashpotTractionFvPatchVectorField
(
    const arosticaNormalSpringDashpotTractionFvPatchVectorField& pvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    arosticaSpringDashpotTractionFvPatchVectorField(pvf, p, iF, mapper)
{}

#ifndef OPENFOAM_ORG
arosticaNormalSpringDashpotTractionFvPatchVectorField::
arosticaNormalSpringDashpotTractionFvPatchVectorField
(
    const arosticaNormalSpringDashpotTractionFvPatchVectorField& pvf
)
:
    arosticaSpringDashpotTractionFvPatchVectorField(pvf)
{}
#endif

arosticaNormalSpringDashpotTractionFvPatchVectorField::
arosticaNormalSpringDashpotTractionFvPatchVectorField
(
    const arosticaNormalSpringDashpotTractionFvPatchVectorField& pvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    arosticaSpringDashpotTractionFvPatchVectorField(pvf, iF)
{}

void arosticaNormalSpringDashpotTractionFvPatchVectorField::
validateTangentialTraction(const dictionary& dict) const
{
    if (!dict.found("tangentialTraction")) return;
    const vectorField tangential("tangentialTraction", dict, patch().size());
    if (max(mag(tangential)) > SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "Aróstica epicardial tangentialTraction must be zero"
            << exit(FatalIOError);
    }
}

void arosticaNormalSpringDashpotTractionFvPatchVectorField::
calculateContributions
(
    vectorField& spring,
    vectorField& dashpot,
    const vectorField& displacement,
    const vectorField& velocity
) const
{
    const vectorField N(patch().nf());
    spring = -springCoefficient_.value()*(displacement & N)*N;
    dashpot = -dashpotCoefficient_.value()*(velocity & N)*N;
}


arosticaVectorSpringDashpotTractionFvPatchVectorField::
arosticaVectorSpringDashpotTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    arosticaSpringDashpotTractionFvPatchVectorField(p, iF)
{}

arosticaVectorSpringDashpotTractionFvPatchVectorField::
arosticaVectorSpringDashpotTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    arosticaSpringDashpotTractionFvPatchVectorField(p, iF, dict)
{}

arosticaVectorSpringDashpotTractionFvPatchVectorField::
arosticaVectorSpringDashpotTractionFvPatchVectorField
(
    const arosticaVectorSpringDashpotTractionFvPatchVectorField& pvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    arosticaSpringDashpotTractionFvPatchVectorField(pvf, p, iF, mapper)
{}

#ifndef OPENFOAM_ORG
arosticaVectorSpringDashpotTractionFvPatchVectorField::
arosticaVectorSpringDashpotTractionFvPatchVectorField
(
    const arosticaVectorSpringDashpotTractionFvPatchVectorField& pvf
)
:
    arosticaSpringDashpotTractionFvPatchVectorField(pvf)
{}
#endif

arosticaVectorSpringDashpotTractionFvPatchVectorField::
arosticaVectorSpringDashpotTractionFvPatchVectorField
(
    const arosticaVectorSpringDashpotTractionFvPatchVectorField& pvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    arosticaSpringDashpotTractionFvPatchVectorField(pvf, iF)
{}

void arosticaVectorSpringDashpotTractionFvPatchVectorField::
calculateContributions
(
    vectorField& spring,
    vectorField& dashpot,
    const vectorField& displacement,
    const vectorField& velocity
) const
{
    spring = -springCoefficient_.value()*displacement;
    dashpot = -dashpotCoefficient_.value()*velocity;
}

Foam::scalar
Foam::arosticaSpringDashpotTractionFvPatchVectorField::ddtCoefficient() const
{
    const volVectorField& D = db().lookupObject<volVectorField>("D");
    const fvMesh& m = patch().boundaryMesh().mesh();
    const scalar dt = m.time().deltaTValue();

#ifdef OPENFOAM_NOT_EXTEND
    const word ddtName(m.ddtScheme("ddt(" + D.name() + ')'));
#else
    const word ddtName(m.schemesDict().ddtScheme("ddt(" + D.name() + ')'));
#endif

    if (ddtName == "backward")
    {
        // backwardDdtScheme falls back to Euler on the first step, when the
        // old and old-old time indices coincide
        const scalar dt0 =
        (
            D.oldTime().timeIndex() == D.oldTime().oldTime().timeIndex()
          ? GREAT
          : m.time().deltaT0Value()
        );

        return (1.0 + dt/(dt + dt0))/dt;
    }

    // Euler, and a safe conservative default for anything else
    return 1.0/dt;
}


Foam::scalar
Foam::arosticaSpringDashpotTractionFvPatchVectorField::effectiveStiffness() const
{
    return
        springCoefficient_.value()
      + dashpotCoefficient_.value()*ddtCoefficient();
}


bool Foam::arosticaNormalSpringDashpotTractionFvPatchVectorField::
supportTangent(scalarField& kEff, bool& normalOnly) const
{
    kEff.setSize(patch().size(), effectiveStiffness());
    normalOnly = true;
    return true;
}


bool Foam::arosticaVectorSpringDashpotTractionFvPatchVectorField::
supportTangent(scalarField& kEff, bool& normalOnly) const
{
    kEff.setSize(patch().size(), effectiveStiffness());
    normalOnly = false;
    return true;
}


makePatchTypeField
(
    fvPatchVectorField,
    arosticaNormalSpringDashpotTractionFvPatchVectorField
);
makePatchTypeField
(
    fvPatchVectorField,
    arosticaVectorSpringDashpotTractionFvPatchVectorField
);

} // End namespace Foam
