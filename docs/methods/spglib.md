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

## Outlook

Where we are: IBZ generation, DM reconstruction via `cp_cfm_gemm` (PZGEMM),
and the MO-side rotation are all in tree on master and pass 52 / 52 regtest
cases at mpiranks 1, 4, 12. Compute is distributed. What remains splits
into two categories: closing the memory-distribution gap (items 1 and 2,
which together still touch capability), and polish (items 3 through 6).

Effort and value are rated relative to the integration that landed up to
commit `9a8f26587b` -- roughly two weeks of design, the DBCSR detour, the
ScaLAPACK rewrite, and 156 regtest runs.

| #   | Item                                                                 | Scope     | Risk   | Value  | When        |
|-----|----------------------------------------------------------------------|-----------|--------|--------|-------------|
| 1   | Drop replicated dense buffers in DM path (`copy_dbcsr_to_fm` + `cp_fm_to_cfm`, mirror `dbcsr_to_cfm` in `rpa_gw.F:6448`; symmetric write-back) | ~10%      | low    | high   | do first    |
| 2   | Distribute the MO path the same way (currently replicates the (nao, nmo) coefficient tile via `cp_fm_get_submatrix`) | ~10%      | low    | medium | after #1    |
| 3   | Retire the `symtrans` multi-rank abort at `kpoint_methods.F:2130` once every dispatch site is covered | ~2%       | medium | medium | after #1/#2 |
| 4   | Cover missing point groups 622, m-3, 432 (re-run `discover.py` with a larger Materials Project window) | ~3%       | low    | low    | nice-to-have |
| 5   | Aux-fit (`for_aux_fit=.TRUE.`) + energy-weighted DM (`pmat_ext`, `do_ext=.TRUE.`) regtest spot-check | ~3%       | low    | low-med| nice-to-have |
| 6   | Performance baseline vs replicated and legacy `symtrans` on the 80-atom NASICON cell | ~5%       | low    | medium | after #1+#2 |

Roughly 30% of the original integration effort remains, almost all of the
capability gain in items 1 and 2. After those, the spglib path is on par
with the rest of the k-point code for both compute and memory; the rest
is documentation-grade.

### Items in detail

1. **Remove the replicated dense buffers in the DM path.** The compute is
   already distributed (PZGEMM), but `pmat_kibz`, `pmat_kstar`, `delta_full`
   are still O(N_basis^2) per rank. Build `cfm_p` directly via
   `copy_dbcsr_to_fm` + `cp_fm_to_cfm` (mirroring `dbcsr_to_cfm` in
   `rpa_gw.F`) and write back via `cp_cfm_to_fm` + `copy_fm_to_dbcsr`.
   Only `delta_full` stays replicated (cheap relative to P). After that,
   the DM path is fully distributed in both memory and compute.
2. **Distribute the MO path too.** The MO rotation still replicates the
   (nao, nmo) coefficient matrix on every rank via `cp_fm_get_submatrix`.
   The same cp_cfm-based product chain can be reused, modulo the real
   FM type. MO is the production default for the SYMMETRY ON cases, so
   this is the bigger memory win in practice once #1 proves the pattern.
3. **Retire the legacy `symtrans` multi-rank abort.** Once the spglib
   path covers every code site that currently dispatches to `symtrans`,
   replace the abort at `kpoint_methods.F:2130` with a hard fallback to
   the spglib path or remove the branch outright. Verification is the
   bulk of the work, not the edit.
4. **Cover the missing point groups.** Regtests skip 622, m-3, 432
   because the discovery pass on Materials Project did not return a
   single small stable representative for those groups. Re-run
   `discover.py` with a larger candidate window. Completeness, not
   correctness; the existing 29 groups already exercise every code path.
5. **Coverage extensions.** Aux-fit (`for_aux_fit=.TRUE.`) and the
   energy-weighted DM via `pmat_ext` (`do_ext=.TRUE.`) currently use the
   DM path; confirm regtests touch both branches, otherwise add a small
   dedicated case. Protects against future regressions on branches that
   aren't currently asserted.
6. **Performance baseline.** Take wall-time numbers vs the dense
   replicated variant and the legacy `symtrans` path on a representative
   cell (e.g. the 80-atom NASICON used in the battery campaign). Run
   this *after* #1 and #2 so the numbers reflect the final state, not
   an intermediate one.

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
