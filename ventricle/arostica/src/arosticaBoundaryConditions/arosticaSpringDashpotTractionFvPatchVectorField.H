/*---------------------------------------------------------------------------*\
  Case-local Aróstica benchmark spring/dashpot traction conditions.
\*---------------------------------------------------------------------------*/

#ifndef arosticaSpringDashpotTractionFvPatchVectorField_H
#define arosticaSpringDashpotTractionFvPatchVectorField_H

#include "solidTractionFvPatchVectorField.H"
#include "volFields.H"

namespace Foam
{

class arosticaSpringDashpotTractionFvPatchVectorField
:
    public solidTractionFvPatchVectorField
{
protected:
    const dimensionedScalar springCoefficient_;
    const dimensionedScalar dashpotCoefficient_;
    const Switch writeDiagnostics_;

    //- Where the Robin support evaluates the displacement it acts on:
    //
    //    ownerCell             D_c              (first order in h)
    //    cellReconstructedFace D_c + delta & grad(D)_c
    //
    //  The benchmark's condition is a SURFACE condition in u|_Gamma, so the
    //  face value is the consistent one. With alpha_epi = 1e8 Pa/m the O(h)
    //  difference between the two reaches ~28 kPa at t=0.26, against a peak
    //  chamber pressure of 16.1 kPa, so the distinction is not academic.
    //
    //  The reconstruction uses the production leastSquaresS4f gradient, which
    //  hardcodes useBoundaryFaceValues = false on every patch, so it is built
    //  from owner and internal-neighbour CELL values only and never inherits
    //  this patch's own fixedGradient value. That is what keeps it free of the
    //  traction/impK closure involved in the original Change-1 defect.
    const word displacementMode_;

    static dictionary withDefaults(const dictionary&);
    static dimensionedScalar readCoefficient
    (
        const dictionary&,
        const word& dimensionedName,
        const word& legacyName,
        const dimensionSet&
    );
    static void requireReferenceArea(const dictionary&);

        //- Current-time coefficient of the ddt scheme the residual uses for
        //  the dashpot velocity. Replicates backwardDdtScheme's own logic,
        //  including its first-step fallback to Euler
        scalar ddtCoefficient() const;

        //- Effective support stiffness per unit reference area,
        //  kEff = alpha + beta*ddtCoefficient()
        scalar effectiveStiffness() const;

    virtual void calculateContributions
    (
        vectorField& spring,
        vectorField& dashpot,
        const vectorField& displacement,
        const vectorField& velocity
    ) const = 0;

public:
    arosticaSpringDashpotTractionFvPatchVectorField
    (
        const fvPatch&,
        const DimensionedField<vector, volMesh>&
    );
    arosticaSpringDashpotTractionFvPatchVectorField
    (
        const fvPatch&,
        const DimensionedField<vector, volMesh>&,
        const dictionary&
    );
    arosticaSpringDashpotTractionFvPatchVectorField
    (
        const arosticaSpringDashpotTractionFvPatchVectorField&,
        const fvPatch&,
        const DimensionedField<vector, volMesh>&,
        const fvPatchFieldMapper&
    );
#ifndef OPENFOAM_ORG
    arosticaSpringDashpotTractionFvPatchVectorField
    (
        const arosticaSpringDashpotTractionFvPatchVectorField&
    );
#endif
    arosticaSpringDashpotTractionFvPatchVectorField
    (
        const arosticaSpringDashpotTractionFvPatchVectorField&,
        const DimensionedField<vector, volMesh>&
    );

    virtual void updateCoeffs();
    virtual void write(Ostream&) const;
};


class arosticaNormalSpringDashpotTractionFvPatchVectorField
:
    public arosticaSpringDashpotTractionFvPatchVectorField
{
    virtual void calculateContributions
    (
        vectorField&,
        vectorField&,
        const vectorField&,
        const vectorField&
    ) const;
    void validateTangentialTraction(const dictionary&) const;

public:
    TypeName("arosticaNormalSpringDashpotTraction");

    //- Support tangent: normal-only, kEff*(N outer N)
    virtual bool supportTangent(scalarField& kEff, bool& normalOnly) const;

    arosticaNormalSpringDashpotTractionFvPatchVectorField
    (
        const fvPatch&,
        const DimensionedField<vector, volMesh>&
    );
    arosticaNormalSpringDashpotTractionFvPatchVectorField
    (
        const fvPatch&,
        const DimensionedField<vector, volMesh>&,
        const dictionary&
    );
    arosticaNormalSpringDashpotTractionFvPatchVectorField
    (
        const arosticaNormalSpringDashpotTractionFvPatchVectorField&,
        const fvPatch&,
        const DimensionedField<vector, volMesh>&,
        const fvPatchFieldMapper&
    );
#ifndef OPENFOAM_ORG
    arosticaNormalSpringDashpotTractionFvPatchVectorField
    (
        const arosticaNormalSpringDashpotTractionFvPatchVectorField&
    );
#endif
    arosticaNormalSpringDashpotTractionFvPatchVectorField
    (
        const arosticaNormalSpringDashpotTractionFvPatchVectorField&,
        const DimensionedField<vector, volMesh>&
    );

    virtual tmp<fvPatchVectorField> clone() const
    {
        return tmp<fvPatchVectorField>
        (
            new arosticaNormalSpringDashpotTractionFvPatchVectorField(*this)
        );
    }
    virtual tmp<fvPatchVectorField> clone
    (
        const DimensionedField<vector, volMesh>& iF
    ) const
    {
        return tmp<fvPatchVectorField>
        (
            new arosticaNormalSpringDashpotTractionFvPatchVectorField(*this, iF)
        );
    }
};


class arosticaVectorSpringDashpotTractionFvPatchVectorField
:
    public arosticaSpringDashpotTractionFvPatchVectorField
{
    virtual void calculateContributions
    (
        vectorField&,
        vectorField&,
        const vectorField&,
        const vectorField&
    ) const;

public:
    TypeName("arosticaVectorSpringDashpotTraction");

    //- Support tangent: full vector, kEff*I
    virtual bool supportTangent(scalarField& kEff, bool& normalOnly) const;

    arosticaVectorSpringDashpotTractionFvPatchVectorField
    (
        const fvPatch&,
        const DimensionedField<vector, volMesh>&
    );
    arosticaVectorSpringDashpotTractionFvPatchVectorField
    (
        const fvPatch&,
        const DimensionedField<vector, volMesh>&,
        const dictionary&
    );
    arosticaVectorSpringDashpotTractionFvPatchVectorField
    (
        const arosticaVectorSpringDashpotTractionFvPatchVectorField&,
        const fvPatch&,
        const DimensionedField<vector, volMesh>&,
        const fvPatchFieldMapper&
    );
#ifndef OPENFOAM_ORG
    arosticaVectorSpringDashpotTractionFvPatchVectorField
    (
        const arosticaVectorSpringDashpotTractionFvPatchVectorField&
    );
#endif
    arosticaVectorSpringDashpotTractionFvPatchVectorField
    (
        const arosticaVectorSpringDashpotTractionFvPatchVectorField&,
        const DimensionedField<vector, volMesh>&
    );

    virtual tmp<fvPatchVectorField> clone() const
    {
        return tmp<fvPatchVectorField>
        (
            new arosticaVectorSpringDashpotTractionFvPatchVectorField(*this)
        );
    }
    virtual tmp<fvPatchVectorField> clone
    (
        const DimensionedField<vector, volMesh>& iF
    ) const
    {
        return tmp<fvPatchVectorField>
        (
            new arosticaVectorSpringDashpotTractionFvPatchVectorField(*this, iF)
        );
    }
};

} // End namespace Foam

#endif

