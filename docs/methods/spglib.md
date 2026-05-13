# spglib k-point IBZ symmetry path

Tracker for the spglib-based irreducible-Brillouin-zone reduction and
density-matrix reconstruction in `kpoint_density_spglib.F`. This page is the
single source of truth for what is wired up, what is gated by which guard,
and what is queued next. Keep it in sync with the code; do not let the
`Status` table drift.

## Status

| Path                                  | Trigger                                                       | Works at | Notes                                                                  |
|---------------------------------------|---------------------------------------------------------------|---------:|------------------------------------------------------------------------|
| spglib DM reconstruction              | `kpoint%use_spglib_dm` and `kpoint%scoord` associated         | any rank | replicate-rotate-scatter via `group%sum(pmat)` inside `dbcsr_to_dense_complex` |
| spglib MO rotation                    | `use_spglib` plus single-rank DBCSR plus `(.NOT. do_ext)`     |    rank 1 | needs flat IBZ MO FM tile; rank > 1 falls through to DM path           |
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
mpiranks=1   52 / 52 OK   ( 5 min)
mpiranks=4   52 / 52 OK   ( 2 min)
```

## Implementation pointers

- `src/kpoint_density_spglib.F`
  - `spglib_dm_workspace_init` -- pulls P(k_IBZ) into a dense complex AO x AO buffer.
  - `dbcsr_to_dense_complex` -- iterates each rank's locally-owned blocks, then `group%sum(pmat)` so every rank holds the full matrix.
  - `spglib_dm_rotate_one_star` -- one star member: `P(k_star) = D~ P(k_IBZ) D~^dagger`. The atom permutation `f0` and the non-symmorphic phase `exp{2 pi i (k_rot . tau_i - k_IBZ . tau_{f0(i)})}` are absorbed into a single dense `D~`.
  - `dense_complex_to_dbcsr` -- writes only the destination DBCSR's locally-owned blocks; missing dense entries are silently dropped.
  - `spglib_mo_*` family -- mirrors the FHI-aims MO-side path; gated single-rank for now, primarily kept for dump-and-diff debugging.
- `src/cryssym.F` -- spglib dispatch for the IBZ generation; populates `csym%rt/vt/f0/nrtot/ibrot`.
- `src/kpoint_methods.F`
  - DM/MO dispatch around `kpoint_density_transform`, in the `apply_symmetry` branch (search `use_spglib`).
  - Legacy `symtrans` (`kpoint_methods.F:2024`) is reachable only when the spglib path is opted out; its multi-rank `cp_abort` lives at `:2130`.
- `src/kpoint_types.F` -- adds `scoord`, `symop_rt`, `symop_vt`, and the `use_spglib_dm` toggle to `kpoint_type`.

## Future plan

In rough order:

1. **Multi-rank MO path.** Generalise `spglib_mo_workspace_init` to accept a distributed FM (gather to dense via `cp_fm_get_submatrix` or a transpose-reduce; run the per-MO rotation locally; redistribute on write). Lift the `(spglib_numnodes == 1)` gate from `use_mo_rot` once that lands.
2. **Distributed `D~ P D~^dagger`.** The current DM path replicates the full dense P on every rank (memory O(N_basis^2)). Move to a ScaLAPACK / `cp_fm`-distributed product when memory becomes the bottleneck (large primitive cells, dense basis sets).
3. **Retire the legacy `symtrans` multi-rank abort.** Once the DM path covers every code site that currently dispatches to `symtrans`, replace the abort at `kpoint_methods.F:2130` with a hard fallback to the spglib path or remove the branch.
4. **Cover the missing point groups.** Regtests skip 622, m-3, 432 because the discovery pass on Materials Project did not return a single small stable representative for those groups. Re-run `discover.py` with a larger candidate window if needed.
5. **Coverage extensions.** Aux-fit (`for_aux_fit=.TRUE.`) and the energy-weighted DM via `pmat_ext` (`do_ext=.TRUE.`) currently use the DM path; spot-check that the regtests touch both branches, otherwise add a small dedicated case.
6. **Performance baseline.** Once the distributed product lands, take wall-time numbers vs the legacy path on a representative cell (e.g. the 80-atom NASICON used in the battery campaign) so the trade-off is documented.

## Source bench

Original prototype lives outside the CP2K tree at
`/workspace/cp2k_spglib_integration/sym_regtest/` (29 of the 32
crystallographic point groups, sym/nosym pairs, results.csv). The
two regtest dirs in `tests/QS/` are derived from that bench; keep the
prototype around for any future point-group additions or for re-deriving
references after a change in DFT defaults.
