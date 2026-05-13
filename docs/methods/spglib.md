# spglib k-point IBZ symmetry path

Tracker for the spglib-based irreducible-Brillouin-zone reduction and
density-matrix reconstruction in `kpoint_density_spglib.F`. This page is the
single source of truth for what is wired up, what is gated by which guard,
and what is queued next. Keep it in sync with the code; do not let the
`Status` table drift.

## Status

| Path                                  | Trigger                                                       | Works at | Notes                                                                  |
|---------------------------------------|---------------------------------------------------------------|---------:|------------------------------------------------------------------------|
| spglib DM reconstruction              | `kpoint%use_spglib_dm` and `kpoint%scoord` associated         | any rank | `D~ P D~^dagger` runs through `cp_cfm_gemm` (PZGEMM) on cp_cfm tiles cloned from `fmwork(1)%matrix_struct`; replicated dense buffers feed `cp_cfm_set_submatrix` and capture the gathered result |
| spglib MO rotation                    | `use_spglib` plus `my_kpgrp` plus `kp%mos` plus `(.NOT. do_ext)` | any rank | `cp_fm_get_submatrix` replicates IBZ MO FM tile on every rank, then the per-MO ZGEMM runs redundantly |
| spglib IBZ generation (`cryssym.F`)   | `SYMMETRY ON` in `&KPOINTS`                                   | any rank | `spg_get_ir_reciprocal_mesh`; populates `csym%rt/vt/f0/nrtot/ibrot` |
| legacy `symtrans` fallback            | dispatch with neither `use_spglib` nor `use_mo_rot`           |    rank 1 | aborts at `kpoint_methods.F:2130` for distributed matrices; do not rely on it |

## Regtest coverage

Two dirs in `tests/QS/`, both tagged `spglib`:

- `regtest-kp-spglib-sym/` -- `mpiranks==1`, exercises the spglib MO path
  for the SYMMETRY ON inputs and the DM path as fallback.
- `regtest-kp-spglib-sym-mpi/` -- `mpiranks>=2`, exercises the multi-rank
  DM path. Both dirs share the same 26 sym/nosym pairs (one stable
  Materials Project structure per crystallographic point group, excluding
  622, m-3, 432). Tolerance 1e-8 Ha; references taken from the
  no-symmetry full-grid runs.

Latest local run on `cp2k.psmp`, gfortran 14.2, this branch:

```
mpiranks=1   52 / 52 OK   (8 min)   regtest-kp-spglib-sym
mpiranks=4   52 / 52 OK   (3 min)   regtest-kp-spglib-sym-mpi  (MO path)
mpiranks=12  52 / 52 OK   (3 min)   regtest-kp-spglib-sym-mpi  (MO path)
```

The DM path itself was separately verified by forcing `use_mo_rot=.FALSE.`
in `kpoint_methods.F` and re-running the same suites:

```
mpiranks=1   52 / 52 OK   (8 min)   DM path via cp_cfm_gemm
mpiranks=4   52 / 52 OK   (3 min)   DM path via cp_cfm_gemm
mpiranks=12  52 / 52 OK   (3 min)   DM path via cp_cfm_gemm
```

## Implementation pointers

- `src/kpoint_density_spglib.F`
  - `spglib_dm_workspace_init` -- pulls P(k_IBZ) into a dense complex AO x AO buffer (`pmat_kibz`) replicated on every rank, then distributes it into a cp_cfm tile (`cfm_p`) cloned from the caller's `fm_struct_template` (passed down as `fmwork(1)%matrix_struct`).
  - `dbcsr_to_dense_complex` -- iterates each rank's locally-owned blocks, then `group%sum(pmat)` so every rank holds the full matrix; still used by the workspace init for the initial replicated copy.
  - `spglib_dm_rotate_one_star` -- one star member: `P(k_star) = D~ P(k_IBZ) D~^dagger`. The atom permutation `f0` and the non-symmorphic phase `exp{2 pi i (k_rot . tau_i - k_IBZ . tau_{f0(i)})}` are absorbed into a single dense `D~`. The product chain runs via two `cp_cfm_gemm` calls (PZGEMM) on `cfm_d / cfm_p / cfm_w / cfm_r`; the final result is gathered to the replicated `pmat_kstar` for Hermitization and the DBCSR write-back.
  - `dense_complex_to_dbcsr` -- writes only the destination DBCSR's locally-owned blocks; missing dense entries are silently dropped.
  - `spglib_mo_*` family -- mirrors the FHI-aims MO-side path; multi-rank capable (cp_fm_get_submatrix replicates the IBZ MO tile on every rank, then the per-MO ZGEMM runs redundantly).
- `src/cryssym.F` -- spglib dispatch for the IBZ generation; populates `csym%rt/vt/f0/nrtot/ibrot`.
- `src/kpoint_methods.F`
  - DM/MO dispatch around `kpoint_density_transform`, in the `apply_symmetry` branch (search `use_spglib`).
  - Legacy `symtrans` (`kpoint_methods.F:2024`) is reachable only when the spglib path is opted out; its multi-rank `cp_abort` lives at `:2130`.
- `src/kpoint_types.F` -- adds `scoord`, `symop_rt`, `symop_vt`, and the `use_spglib_dm` toggle to `kpoint_type`.

## Future plan

In rough order:

1. **Remove the replicated dense buffers in the DM path.** The compute is
   already distributed (PZGEMM), but `pmat_kibz`, `pmat_kstar`, `delta_full`
   are still O(N_basis^2) per rank. Building `cfm_p` directly via
   `copy_dbcsr_to_fm` + `cp_fm_to_cfm` (mirroring `dbcsr_to_cfm` in
   `rpa_gw.F`) and writing back via `cp_cfm_to_fm` + `copy_fm_to_dbcsr`
   removes the gather/replicate steps; only the `delta_full` build remains
   replicated (cheap relative to P). After that, the DM path is fully
   distributed in both memory and compute.
2. **Distribute the MO path too.** The MO rotation still replicates the
   (nao, nmo) coefficient matrix on every rank via `cp_fm_get_submatrix`.
   The same cp_cfm-based product chain can be reused.
3. **Retire the legacy `symtrans` multi-rank abort.** Once the spglib path
   covers every code site that currently dispatches to `symtrans`, replace
   the abort at `kpoint_methods.F:2130` with a hard fallback to the
   spglib path or remove the branch.
4. **Cover the missing point groups.** Regtests skip 622, m-3, 432 because
   the discovery pass on Materials Project did not return a single small
   stable representative for those groups. Re-run `discover.py` with a
   larger candidate window if needed.
5. **Coverage extensions.** Aux-fit (`for_aux_fit=.TRUE.`) and the
   energy-weighted DM via `pmat_ext` (`do_ext=.TRUE.`) currently use the
   DM path; spot-check that the regtests touch both branches, otherwise
   add a small dedicated case.
6. **Performance baseline.** Take wall-time numbers vs the dense-replicated
   variant and the legacy `symtrans` path on a representative cell (e.g.
   the 80-atom NASICON used in the battery campaign) so the trade-off is
   documented.

## DBCSR detour (not in tree)

An earlier attempt rewrote the DM rotation entirely on top of DBCSR
(`dbcsr_multiply` chain on a sparse `delta` and a desymmetrized P).
Single-rank runs were correct but the multi-rank path produced ~0.3 Ha
energy drift that resisted every variant tried (ownership-filtered vs.
fully-reserved delta, set+filter vs. release+create resets, with and
without explicit Hermitization, `transpose_distribution=.FALSE.` on
`dbcsr_transposed`, src- vs dst-iterating projection). Root cause is
suspected to be a SUMMA-level interaction with the extremely sparse
`delta` operand. The branch was discarded; the cp_cfm/PZGEMM path here
is what landed.

## Source bench

Original prototype lives outside the CP2K tree at
`/workspace/cp2k_spglib_integration/sym_regtest/` (29 of the 32
crystallographic point groups, sym/nosym pairs, results.csv). The
two regtest dirs in `tests/QS/` are derived from that bench; keep the
prototype around for any future point-group additions or for re-deriving
references after a change in DFT defaults.
