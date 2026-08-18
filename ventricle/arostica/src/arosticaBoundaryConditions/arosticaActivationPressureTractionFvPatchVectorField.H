/*---------------------------------------------------------------------------*\
  Case-local Aróstica benchmark time-dependent endocardial pressure
  boundary condition.

  Description
      Solves the Aróstica et al. pressure evolution equation (Eq. (7)-(8) in
      "A software benchmark for cardiac elastodynamics", Comput. Methods
      Appl. Mech. Engrg. 435 (2025) 117485):

          dp/dt = -|b(t)| p(t) + sigmaMid |b(t)|_+ + sigmaPre |gPre(t)|_+

          b(t)    = aPre(t) + alphaPre gPre(t) + alphaMid
          aPre(t) = alphaMax fPre(t) + alphaMin (1 - fPre(t))
          fPre(t) = S+(t - tSysPre) S-(t - tDiasPre)
          gPre(t) = S-(t - tDiasPre)
          S+-(x)  = 0.5 (1 +- tanh(x/gamma))

      p(t) is spatially uniform (a pure function of time), so it is
      integrated once per accepted time step with a fixed-substep RK4
      scheme fine enough to resolve the gamma time scale, independent of
      the outer solver time step, then applied as a uniform endocardial
      pressure via the inherited solidTraction pressure mechanism. This
      replaces reading a pre-tabulated pressure curve.
\*---------------------------------------------------------------------------*/

#ifndef arosticaActivationPressureTractionFvPatchVectorField_H
#define arosticaActivationPressureTractionFvPatchVectorField_H

#include "solidTractionFvPatchVectorField.H"
#include "volFields.H"

namespace Foam
{

class arosticaActivationPressureTractionFvPatchVectorField
:
    public solidTractionFvPatchVectorField
{
    // Private Data

        //- Minimum activation rate, alphaMin in Eq. (8)
        const dimensionedScalar alphaMin_;

        //- Maximum activation rate, alphaMax in Eq. (8)
        const dimensionedScalar alphaMax_;

        //- Pre-ejection activation rate, alphaPre in Eq. (8)
        const dimensionedScalar alphaPre_;

        //- Mid-cycle activation rate, alphaMid in Eq. (8)
        const dimensionedScalar alphaMid_;

        //- Pre-ejection driving pressure, sigmaPre in Eq. (7)
        const dimensionedScalar sigmaPre_;

        //- Mid-cycle driving pressure, sigmaMid in Eq. (7)
        const dimensionedScalar sigmaMid_;

        //- Systole onset time, tSysPre in Eq. (8)
        const dimensionedScalar tSysPre_;

        //- Diastole onset time, tDiasPre in Eq. (8)
        const dimensionedScalar tDiasPre_;

        //- Activation function transition sharpness, gamma in Eq. (8)
        const dimensionedScalar gamma_;

        //- Mesh time index at which p_ was last integrated
        mutable label cachedTimeIndex_;

        //- Pressure p(t) at the cached time index
        mutable scalar p_;


    // Private Member Functions

        //- Fill in default traction/pressure entries expected by the base
        //  solidTraction class; the pressure value is overwritten every
        //  time step by the ODE solution
        static dictionary withDefaults(const dictionary&);

        //- Validate the pressure model parameters
        void validateParameters() const;

        //- S+-(deltaT) = 0.5*(1 +- tanh(deltaT/gamma_))
        scalar sPlus(const scalar deltaT) const;
        scalar sMinus(const scalar deltaT) const;

        //- gPre(t), Eq. (8)
        scalar gPre(const scalar t) const;

        //- Activation function b(t), Eq. (8)
        scalar activationFunction(const scalar t) const;

        //- Time derivative of p, Eq. (7)
        scalar pDerivative(const scalar t, const scalar p) const;

        //- Integrate p(t) from 0 to tEnd using fixed-substep RK4
        scalar integrateP(const scalar tEnd) const;

        //- Return p(t) at the current time, integrating and caching it
        //  once per accepted time step
        scalar currentPressure() const;


public:

    //- Runtime type information
    TypeName("arosticaActivationPressureTraction");


    // Constructors

        arosticaActivationPressureTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        arosticaActivationPressureTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        arosticaActivationPressureTractionFvPatchVectorField
        (
            const arosticaActivationPressureTractionFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

#ifndef OPENFOAM_ORG
        arosticaActivationPressureTractionFvPatchVectorField
        (
            const arosticaActivationPressureTractionFvPatchVectorField&
        );
#endif

        arosticaActivationPressureTractionFvPatchVectorField
        (
            const arosticaActivationPressureTractionFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new arosticaActivationPressureTractionFvPatchVectorField(*this)
            );
        }

        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new arosticaActivationPressureTractionFvPatchVectorField
                (
                    *this, iF
                )
            );
        }


    // Member Functions

        virtual void updateCoeffs();

        virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif
