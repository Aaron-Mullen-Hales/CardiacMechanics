# Aróstica support boundary conditions

This case-local library provides the benchmark's two nominal-traction support
laws:

- `arosticaNormalSpringDashpotTraction`: epicardial normal-only support;
- `arosticaVectorSpringDashpotTraction`: basal full-vector support.

The implementation is stateless with respect to nonlinear trials. It evaluates
the current displacement and the configured OpenFOAM `ddt(D)` and inserts the
result through `solidTraction`.

Every patch dictionary must explicitly contain `useUndeformedArea true`. This
is essential with `s4f-cardiac-clean`: its `solidTraction` accessor is
non-virtual, so reference-area integration is selected from the base-class
dictionary value, not from a derived override.

The implementation is adapted from the previously tested archived Aróstica
conditions in this workspace; it is kept local because those runtime classes
are absent from the clean solids4foam build.

