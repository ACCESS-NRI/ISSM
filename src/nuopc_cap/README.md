# ISSM NUOPC Cap

This directory contains a first, minimal `NUOPC` cap for running ISSM as an
in-process coupled component.

## Scope

This first version is intentionally small:

- geometry: `ESMF_Mesh` built from the native ISSM unstructured mesh
- import: one nodal field, `floatingIceMeltRate`
- exports: `iceThickness`, `iceSurface`, `iceMask`
- runtime: initialize from a prepared ISSM case directory containing
  `<model_name>.bin` and `<model_name>.toolkits`

## Expected Workflow

1. Prepare an ISSM case directory with the usual ISSM input files.
2. Ensure the model is marshalled with `TransientSolution`.
3. Ensure `md.transient.isoceancoupling = 0`.
   This cap owns the coupling exchange and does not use the legacy
   `issm_ocean` / `OceanExchangeDatax` runtime path.
4. Configure the NUOPC component with the attributes:
   - `case_dir`
   - `model_name`
   - `solution_name` (optional, defaults to `TransientSolution`)
   - `write_restart` (optional, `true` or `false`)

## Limitations

- only transient coupling is supported
- only 2-D triangular meshes are supported
- only nodal import/export fields are supported
- the imported melt field is written directly into
  `BasalforcingsFloatingiceMeltingRateEnum`, so the caller must supply values
  using ISSM's expected sign and units
- PETSc ownership still follows the ISSM runtime pattern, so this first version
  assumes ISSM owns PETSc initialization/finalization inside the coupled job
