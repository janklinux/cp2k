!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                  !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief FHI-spglib-style k-point density-matrix reconstruction.
!>
!>        The full-basis rotation matrix \tilde{D}(R) is built so that the
!>        atom permutation f0 and the (non-symmorphic) phase
!>        exp{2 pi i (k_rot . tau_i  -  k_IBZ . tau_{f0(i)})} are absorbed
!>        into a single complex AO x AO matrix. The star expansion then
!>        reduces to one dense complex matrix product:
!>
!>            P(k_star) = \tilde{D} P(k_IBZ) \tilde{D}^\dagger
!>
!>        which is algebraically equivalent to rotating each KS eigenvector
!>        first (the path used by FHI-spglib, sym_base.f90:2862,2956,3137) and
!>        then forming P from the rotated MOs.
!>
!>        See implementation_plan.pdf section 5, Step 3 and section 4 for the
!>        gap analysis vs. the existing CP2K symtrans path.
!> \author janK 2026-05-04
! **************************************************************************************************
MODULE kpoint_density_spglib
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_gemm
   USE cp_cfm_types,                    ONLY: cp_cfm_create, cp_cfm_get_submatrix, cp_cfm_release, &
                                              cp_cfm_set_all, cp_cfm_set_submatrix, cp_cfm_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_get_block_p, dbcsr_get_info, dbcsr_get_matrix_type, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_set, dbcsr_type, dbcsr_type_no_symmetry
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create, cp_fm_struct_release, &
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_get_submatrix, cp_fm_type
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: kpoint_env_type, kpoint_sym_type, kpoint_type
   USE mathconstants,                   ONLY: twopi
   USE message_passing,                 ONLY: mp_comm_type
   USE qs_mo_types,                     ONLY: get_mo_set, mo_set_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kpoint_density_spglib'

   PUBLIC :: spglib_dm_workspace_type
   PUBLIC :: spglib_dm_workspace_init, spglib_dm_workspace_release
   PUBLIC :: spglib_dm_rotate_one_star
   PUBLIC :: dump_sparsified_pmat
   ! FHI-spglib-style MO-rotation path (Step 3 of implementation_plan.pdf):
   ! rotate each KS eigenvector first, then form P from the rotated MOs.
   ! Mirrors evaluate_k_densmat_sym + get_kvec_sym (sym_base.f90:2201, 2862).
   PUBLIC :: spglib_mo_workspace_type
   PUBLIC :: spglib_mo_workspace_init, spglib_mo_workspace_release
   PUBLIC :: spglib_mo_rotate_one_star
   PUBLIC :: dump_rotated_mos, dump_ibz_mos, dump_delta_full, dump_kpsym_meta

! **************************************************************************************************
!> \brief Per-IBZ-k workspace reused across the star-member loop.
!>        Allocated once with spglib_dm_workspace_init (which also pulls
!>        P(k_IBZ) into pmat_kibz and the distributed cfm_p); released
!>        with spglib_dm_workspace_release.
!>
!>        Each rank still carries one full local dense P(k_IBZ) (built
!>        by dbcsr_to_dense_complex via MPI sum) plus per-star-member
!>        scratch buffers, but the D~ P D~^dagger matrix product itself
!>        runs through cp_cfm_gemm (PZGEMM): work is distributed across
!>        the BLACS grid rather than redundantly recomputed on every rank.
! **************************************************************************************************
   TYPE spglib_dm_workspace_type
      INTEGER                                          :: natom = 0
      INTEGER                                          :: n_basis_tot = 0
      LOGICAL                                          :: allocated_cfm = .FALSE.
      INTEGER, ALLOCATABLE, DIMENSION(:)               :: atom_first   ! (natom+1)
      INTEGER, ALLOCATABLE, DIMENSION(:)               :: atom_size    ! (natom)
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)   :: pmat_kibz    ! P(k_IBZ), dense complex
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)   :: pmat_kstar   ! scratch
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)   :: delta_full   ! scratch
      ! Distributed cp_cfm scratch driving the PZGEMM chain.
      TYPE(cp_fm_struct_type), POINTER                 :: fm_struct => NULL()
      TYPE(cp_cfm_type)                                :: cfm_p     ! P(k_IBZ)
      TYPE(cp_cfm_type)                                :: cfm_d     ! Delta
      TYPE(cp_cfm_type)                                :: cfm_w     ! Delta * P
      TYPE(cp_cfm_type)                                :: cfm_r     ! Delta * P * Delta^H
   END TYPE spglib_dm_workspace_type

! **************************************************************************************************
!> \brief Per-IBZ-k workspace for the FHI-spglib-style MO-rotation path.
!>
!>        Holds the IBZ KS eigenvectors as a dense complex matrix
!>        mos_kibz(n_basis_tot, nmo) and per-MO weights (occupations,
!>        possibly multiplied by eigenvalues for the energy-weighted
!>        density matrix). Allocated once with spglib_mo_workspace_init,
!>        released with spglib_mo_workspace_release.
!>
!>        The reciprocal-space density matrix at any star member is then
!>            P(k_star) = \tilde{D}(R, k_IBZ) C(k_IBZ) f C(k_IBZ)^\dagger \tilde{D}(R, k_IBZ)^\dagger
!>        with the rotated coefficients C_rot = \tilde{D} C reused for
!>        FHI-spglib-style per-MO debugging dumps (see dump_rotated_mos).
! **************************************************************************************************
   TYPE spglib_mo_workspace_type
      INTEGER                                          :: natom = 0
      INTEGER                                          :: n_basis_tot = 0
      INTEGER                                          :: nmo = 0
      INTEGER                                          :: nmo_active = 0
      INTEGER, ALLOCATABLE, DIMENSION(:)               :: atom_first   ! (natom+1)
      INTEGER, ALLOCATABLE, DIMENSION(:)               :: atom_size    ! (natom)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)         :: weights      ! (nmo) f_n (or f_n*eig_n)
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)   :: mos_kibz     ! (n_basis_tot, nmo)
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)   :: mos_rotated  ! (n_basis_tot, nmo) scratch
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)   :: mos_weighted ! (n_basis_tot, nmo) scratch
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)   :: pmat_kstar   ! (n_basis_tot, n_basis_tot)
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)   :: delta_full   ! (n_basis_tot, n_basis_tot)
   END TYPE spglib_mo_workspace_type

CONTAINS

! **************************************************************************************************
!> \brief Initialise a per-IBZ-k workspace, pull the IBZ density matrix
!>        from (rpmat, cpmat) into the dense complex pmat_kibz, and
!>        distribute it into the cp_cfm scratch matrix cfm_p that drives
!>        the per-star-member PZGEMM chain.
!> \param fm_struct_template  FM struct cloned from the caller's fmwork.
!>        Provides the BLACS context + para_env on which all four cp_cfm
!>        scratch matrices live. Must describe an n_basis_tot x n_basis_tot
!>        layout (or a compatible template; the call here forces the
!>        global dimensions).
! **************************************************************************************************
   SUBROUTINE spglib_dm_workspace_init(ws, kpoint, rpmat, cpmat, real_only, fm_struct_template)
      TYPE(spglib_dm_workspace_type), INTENT(OUT)          :: ws
      TYPE(kpoint_type), POINTER                         :: kpoint
      TYPE(dbcsr_type), POINTER                          :: rpmat, cpmat
      LOGICAL, INTENT(IN)                                :: real_only
      TYPE(cp_fm_struct_type), POINTER, INTENT(IN)       :: fm_struct_template

      CPASSERT(ASSOCIATED(kpoint))
      CPASSERT(ASSOCIATED(kpoint%scoord))
      CPASSERT(ASSOCIATED(kpoint%atype))
      CPASSERT(ASSOCIATED(kpoint%ibrot))
      CPASSERT(ASSOCIATED(kpoint%kind_rotmat))
      CPASSERT(ASSOCIATED(fm_struct_template))

      ws%natom = SIZE(kpoint%atype)
      ALLOCATE (ws%atom_first(ws%natom + 1), ws%atom_size(ws%natom))
      CALL build_atom_basis_table(kpoint, ws%natom, ws%atom_first, ws%atom_size, ws%n_basis_tot)

      ALLOCATE (ws%pmat_kibz(ws%n_basis_tot, ws%n_basis_tot))
      ALLOCATE (ws%pmat_kstar(ws%n_basis_tot, ws%n_basis_tot))
      ALLOCATE (ws%delta_full(ws%n_basis_tot, ws%n_basis_tot))

      ! Build the local dense P(k_IBZ) replicated on every rank. The
      ! cp_cfm push below is purely local (no MPI in cp_cfm_set_submatrix),
      ! so we still need the full N x N on every rank here.
      CALL dbcsr_to_dense_complex(rpmat, cpmat, real_only, ws%atom_first, ws%atom_size, ws%pmat_kibz)

      ! Clone the fmwork struct but force square N_basis_tot x N_basis_tot
      ! global dimensions so the cfm matrices match the rotation algebra.
      NULLIFY (ws%fm_struct)
      CALL cp_fm_struct_create(ws%fm_struct, &
                               nrow_global=ws%n_basis_tot, &
                               ncol_global=ws%n_basis_tot, &
                               template_fmstruct=fm_struct_template)
      CALL cp_cfm_create(ws%cfm_p, ws%fm_struct)
      CALL cp_cfm_create(ws%cfm_d, ws%fm_struct)
      CALL cp_cfm_create(ws%cfm_w, ws%fm_struct)
      CALL cp_cfm_create(ws%cfm_r, ws%fm_struct)
      ws%allocated_cfm = .TRUE.

      ! Distribute P(k_IBZ) into cfm_p (no MPI; each rank writes its
      ! owned block-cyclic tile from the replicated pmat_kibz).
      CALL cp_cfm_set_submatrix(ws%cfm_p, ws%pmat_kibz)
   END SUBROUTINE spglib_dm_workspace_init

! **************************************************************************************************
!> \brief Release a workspace produced by spglib_dm_workspace_init.
! **************************************************************************************************
   SUBROUTINE spglib_dm_workspace_release(ws)
      TYPE(spglib_dm_workspace_type), INTENT(INOUT)        :: ws

      IF (ALLOCATED(ws%atom_first)) DEALLOCATE (ws%atom_first)
      IF (ALLOCATED(ws%atom_size)) DEALLOCATE (ws%atom_size)
      IF (ALLOCATED(ws%pmat_kibz)) DEALLOCATE (ws%pmat_kibz)
      IF (ALLOCATED(ws%pmat_kstar)) DEALLOCATE (ws%pmat_kstar)
      IF (ALLOCATED(ws%delta_full)) DEALLOCATE (ws%delta_full)
      IF (ws%allocated_cfm) THEN
         CALL cp_cfm_release(ws%cfm_p)
         CALL cp_cfm_release(ws%cfm_d)
         CALL cp_cfm_release(ws%cfm_w)
         CALL cp_cfm_release(ws%cfm_r)
         CALL cp_fm_struct_release(ws%fm_struct)
         ws%allocated_cfm = .FALSE.
      END IF
      ws%natom = 0
      ws%n_basis_tot = 0
   END SUBROUTINE spglib_dm_workspace_release

! **************************************************************************************************
!> \brief Rotate P(k_IBZ) -> P(k_star) for one star member (is) and pack the
!>        result into (srpmat, scpmat). Mirrors FHI-spglib's evaluate_k_densmat_sym
!>        + get_density_rotation pair (sym_base.f90:2201, 2956), but operates
!>        on the assembled DM rather than the MO coefficients (algebraically
!>        equivalent: (D C)(D C)^\dagger = D C C^\dagger D^\dagger = D P D^\dagger).
!>
!>        The full-basis Delta absorbs the atom permutation f0 and the per-atom
!>        non-symmorphic phase exp{2 pi i (k_rot . tau_i - k_IBZ . tau_{f0(i)})}.
!>        No separate row/column swap is needed at the DBCSR stage.
!> \param ws       Workspace from spglib_dm_workspace_init
!> \param kpoint   Carries scoord, atype, ibrot, kind_rotmat
!> \param kpsym    Symmetry record for the IBZ k-point
!> \param is       Star-member index (1..kpsym%nwght)
!> \param xkp_ibz  Fractional IBZ k-point coordinates
!> \param real_only True for real Gamma-equivalent runs (cpmat ignored)
!> \param srpmat,scpmat Output DBCSRs (Re/Im of P(k_star))
! **************************************************************************************************
   SUBROUTINE spglib_dm_rotate_one_star(ws, kpoint, kpsym, is, xkp_ibz, real_only, srpmat, scpmat)
      TYPE(spglib_dm_workspace_type), INTENT(INOUT)        :: ws
      TYPE(kpoint_type), POINTER                         :: kpoint
      TYPE(kpoint_sym_type), POINTER                     :: kpsym
      INTEGER, INTENT(IN)                                :: is
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: xkp_ibz
      LOGICAL, INTENT(IN)                                :: real_only
      TYPE(dbcsr_type), POINTER                          :: srpmat, scpmat

      CHARACTER(LEN=*), PARAMETER :: routineN = 'spglib_dm_rotate_one_star'

      INTEGER                                            :: handle, ira, irot, jr

      CALL timeset(routineN, handle)

      ! Resolve kind_rotmat slot. kpsym%rotp(is) is the signed op index;
      ! |rotp| selects the spatial part of the D-matrix. The sign of the
      ! imaginary part for time-reversal partners is automatically correct
      ! because P(k_IBZ) is built with the right sign upstream and our
      ! rotation conjugates D consistently.
      irot = ABS(kpsym%rotp(is))
      ira = 0
      DO jr = 1, SIZE(kpoint%ibrot)
         IF (irot == kpoint%ibrot(jr)) ira = jr
      END DO
      CPASSERT(ira > 0)

      CALL build_full_delta(kpoint, kpsym, is, ira, xkp_ibz, &
                            ws%natom, ws%atom_first, ws%atom_size, ws%n_basis_tot, ws%delta_full)

      ! Distribute the (replicated) local delta into the cp_cfm tile
      ! (no MPI here -- each rank just writes its owned block-cyclic
      ! entries from the local copy).
      CALL cp_cfm_set_submatrix(ws%cfm_d, ws%delta_full)

      ! P(k_spatial) = Delta * P(k_IBZ) * Delta^dagger via two PZGEMM
      ! calls on the cp_cfm tiles -- compute is parallelised across the
      ! BLACS grid.
      CALL cp_cfm_gemm('N', 'N', ws%n_basis_tot, ws%n_basis_tot, ws%n_basis_tot, &
                       (1.0_dp, 0.0_dp), ws%cfm_d, ws%cfm_p, &
                       (0.0_dp, 0.0_dp), ws%cfm_w)
      CALL cp_cfm_gemm('N', 'C', ws%n_basis_tot, ws%n_basis_tot, ws%n_basis_tot, &
                       (1.0_dp, 0.0_dp), ws%cfm_w, ws%cfm_d, &
                       (0.0_dp, 0.0_dp), ws%cfm_r)

      ! Gather the distributed result into the replicated pmat_kstar
      ! (cp_cfm_get_submatrix internally para_env%sum's the tiles, so
      ! every rank ends up with the full N x N P(k_star)).
      CALL cp_cfm_get_submatrix(ws%cfm_r, ws%pmat_kstar)

      IF (kpsym%rotp(is) < 0) ws%pmat_kstar = CONJG(ws%pmat_kstar)

      ! Force exact Hermiticity to remove non-Hermitian round-off that
      ! the symmetric/antisymmetric DBCSR sparsification would mishandle.
      ws%pmat_kstar(:, :) = 0.5_dp*(ws%pmat_kstar(:, :) + CONJG(TRANSPOSE(ws%pmat_kstar(:, :))))

      CALL dense_complex_to_dbcsr(ws%pmat_kstar, ws%atom_first, ws%atom_size, real_only, srpmat, scpmat)

      CALL timestop(handle)
   END SUBROUTINE spglib_dm_rotate_one_star

! **************************************************************************************************
!> \brief Initialise an MO-rotation workspace and pull the IBZ KS eigenvectors
!>        from kp_env into a dense complex matrix.
!>
!>        Mirrors the FHI-spglib data layout used by evaluate_k_densmat_sym
!>        (sym_base.f90:2201): the rotated MOs are formed per star member
!>        from the IBZ MOs cached here, weighted by sqrt(occ) (or
!>        sqrt(occ*eig) for the energy-weighted DM), and accumulated into
!>        P(k_star) via a Hermitian rank-k update.
!>
!>        Single-node restriction enforced via numnodes(kp%mos%mo_coeff) = 1.
!>        The legacy DM path has the same restriction (kpoint_methods.F:1957).
!>
!> \param ws         Output workspace
!> \param kpoint     Kpoint env (carries scoord, atype, ibrot, kind_rotmat)
!> \param kp_env     Per-IBZ-k env that owns kp_env%mos(ic, ispin)
!> \param ispin      Spin index
!> \param real_only  True if the run uses real wavefunctions (only mos(1,:))
!> \param wtype      True for the energy-weighted DM (weights = occ * eig)
! **************************************************************************************************
   SUBROUTINE spglib_mo_workspace_init(ws, kpoint, kp_env, ispin, real_only, wtype)
      TYPE(spglib_mo_workspace_type), INTENT(OUT)          :: ws
      TYPE(kpoint_type), POINTER                         :: kpoint
      TYPE(kpoint_env_type), POINTER                     :: kp_env
      INTEGER, INTENT(IN)                                :: ispin
      LOGICAL, INTENT(IN)                                :: real_only, wtype

      INTEGER                                            :: imo, nao, nmo
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: cre, cim
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues, occupation
      TYPE(cp_fm_type), POINTER                          :: mo_coeff_im, mo_coeff_re

      CPASSERT(ASSOCIATED(kpoint))
      CPASSERT(ASSOCIATED(kpoint%scoord))
      CPASSERT(ASSOCIATED(kpoint%atype))
      CPASSERT(ASSOCIATED(kpoint%ibrot))
      CPASSERT(ASSOCIATED(kpoint%kind_rotmat))
      CPASSERT(ASSOCIATED(kp_env))
      CPASSERT(ASSOCIATED(kp_env%mos))

      ws%natom = SIZE(kpoint%atype)
      ALLOCATE (ws%atom_first(ws%natom + 1), ws%atom_size(ws%natom))
      CALL build_atom_basis_table(kpoint, ws%natom, ws%atom_first, ws%atom_size, ws%n_basis_tot)

      NULLIFY (mo_coeff_re, mo_coeff_im, eigenvalues, occupation)
      CALL get_mo_set(kp_env%mos(1, ispin), nao=nao, nmo=nmo, &
                      mo_coeff=mo_coeff_re, &
                      eigenvalues=eigenvalues, occupation_numbers=occupation)
      CPASSERT(nao == ws%n_basis_tot)
      CPASSERT(ASSOCIATED(mo_coeff_re))
      CPASSERT(ASSOCIATED(occupation))

      ws%nmo = nmo

      ! Replicate the (block-cyclic) MO coefficient FMs into dense (nao, nmo)
      ! buffers that every rank holds in full. cp_fm_get_submatrix zeros the
      ! target, writes locally-owned entries, and does an MPI sum across the
      ! para_env -- so we end up with the same dense MO matrix on each rank.
      ALLOCATE (cre(nao, nmo))
      CALL cp_fm_get_submatrix(mo_coeff_re, cre)

      IF (.NOT. real_only) THEN
         CALL get_mo_set(kp_env%mos(2, ispin), mo_coeff=mo_coeff_im)
         CPASSERT(ASSOCIATED(mo_coeff_im))
         ALLOCATE (cim(nao, nmo))
         CALL cp_fm_get_submatrix(mo_coeff_im, cim)
      END IF

      ALLOCATE (ws%mos_kibz(ws%n_basis_tot, nmo))
      ALLOCATE (ws%mos_rotated(ws%n_basis_tot, nmo))
      ALLOCATE (ws%mos_weighted(ws%n_basis_tot, nmo))
      ALLOCATE (ws%pmat_kstar(ws%n_basis_tot, ws%n_basis_tot))
      ALLOCATE (ws%delta_full(ws%n_basis_tot, ws%n_basis_tot))
      ALLOCATE (ws%weights(nmo))

      IF (real_only) THEN
         ws%mos_kibz(:, :) = CMPLX(cre(:, :), 0.0_dp, KIND=dp)
      ELSE
         ws%mos_kibz(:, :) = CMPLX(cre(:, :), cim(:, :), KIND=dp)
      END IF
      DEALLOCATE (cre)
      IF (.NOT. real_only) DEALLOCATE (cim)

      ws%nmo_active = 0
      DO imo = 1, nmo
         IF (wtype) THEN
            CPASSERT(ASSOCIATED(eigenvalues))
            ws%weights(imo) = occupation(imo)*eigenvalues(imo)
         ELSE
            ws%weights(imo) = occupation(imo)
         END IF
         IF (ABS(ws%weights(imo)) > 0.0_dp) ws%nmo_active = imo
      END DO
   END SUBROUTINE spglib_mo_workspace_init

! **************************************************************************************************
!> \brief Release a workspace produced by spglib_mo_workspace_init.
! **************************************************************************************************
   SUBROUTINE spglib_mo_workspace_release(ws)
      TYPE(spglib_mo_workspace_type), INTENT(INOUT)        :: ws

      IF (ALLOCATED(ws%atom_first)) DEALLOCATE (ws%atom_first)
      IF (ALLOCATED(ws%atom_size)) DEALLOCATE (ws%atom_size)
      IF (ALLOCATED(ws%weights)) DEALLOCATE (ws%weights)
      IF (ALLOCATED(ws%mos_kibz)) DEALLOCATE (ws%mos_kibz)
      IF (ALLOCATED(ws%mos_rotated)) DEALLOCATE (ws%mos_rotated)
      IF (ALLOCATED(ws%mos_weighted)) DEALLOCATE (ws%mos_weighted)
      IF (ALLOCATED(ws%pmat_kstar)) DEALLOCATE (ws%pmat_kstar)
      IF (ALLOCATED(ws%delta_full)) DEALLOCATE (ws%delta_full)
      ws%natom = 0
      ws%n_basis_tot = 0
      ws%nmo = 0
      ws%nmo_active = 0
   END SUBROUTINE spglib_mo_workspace_release

! **************************************************************************************************
!> \brief Rotate the IBZ KS eigenvectors to one star member and accumulate
!>        the resulting density matrix into (srpmat, scpmat).
!>
!>        Mirrors FHI-spglib's get_kvec_sym + evaluate_k_densmat_sym pair
!>        (sym_base.f90:2862, 2201):
!>             C_rot   = \tilde{D}(R, k_IBZ) . C(k_IBZ)
!>             P(k_*)  = sum_n w_n c_rot,n c_rot,n^\dagger
!>                     = C_rot . diag(w) . C_rot^\dagger
!>        Time reversal (kpsym%rotp(is) < 0) is applied as a final
!>        complex conjugation on P(k_*); the spatial \tilde{D} uses
!>        +k_spatial regardless. Hermiticity is exact by construction
!>        because diag(w) is real.
!>
!>        The rotated MOs are kept in ws%mos_rotated for optional
!>        per-MO dumps (dump_rotated_mos), enabling a 1-to-1 diff against
!>        FHI-spglib's evaluate_k_densmat_sym debug output.
! **************************************************************************************************
   SUBROUTINE spglib_mo_rotate_one_star(ws, kpoint, kpsym, is, xkp_ibz, real_only, srpmat, scpmat)
      TYPE(spglib_mo_workspace_type), INTENT(INOUT)        :: ws
      TYPE(kpoint_type), POINTER                         :: kpoint
      TYPE(kpoint_sym_type), POINTER                     :: kpsym
      INTEGER, INTENT(IN)                                :: is
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: xkp_ibz
      LOGICAL, INTENT(IN)                                :: real_only
      TYPE(dbcsr_type), POINTER                          :: srpmat, scpmat

      CHARACTER(LEN=*), PARAMETER :: routineN = 'spglib_mo_rotate_one_star'

      INTEGER                                            :: handle, imo, ira, irot, jr, n
      COMPLEX(KIND=dp)                                   :: w

      CALL timeset(routineN, handle)

      irot = ABS(kpsym%rotp(is))
      ira = 0
      DO jr = 1, SIZE(kpoint%ibrot)
         IF (irot == kpoint%ibrot(jr)) ira = jr
      END DO
      CPASSERT(ira > 0)

      CALL build_full_delta(kpoint, kpsym, is, ira, xkp_ibz, &
                            ws%natom, ws%atom_first, ws%atom_size, ws%n_basis_tot, ws%delta_full)

      n = ws%n_basis_tot
      ws%mos_rotated = (0.0_dp, 0.0_dp)
      ws%pmat_kstar = (0.0_dp, 0.0_dp)
      IF (ws%nmo_active > 0) THEN
         ! C_rot = Delta * C
         CALL ZGEMM('N', 'N', n, ws%nmo_active, n, &
                    (1.0_dp, 0.0_dp), ws%delta_full, n, &
                    ws%mos_kibz, n, &
                    (0.0_dp, 0.0_dp), ws%mos_rotated, n)

         ! C_w = C_rot * diag(w); the column scaling is cheap so do it explicitly.
         DO imo = 1, ws%nmo_active
            w = CMPLX(ws%weights(imo), 0.0_dp, KIND=dp)
            ws%mos_weighted(:, imo) = w*ws%mos_rotated(:, imo)
         END DO

         ! P(k_*) = C_rot * C_w^H  =  sum_n w_n c_rot,n c_rot,n^H. Hermitian
         ! because w_n is real. ZGEMM('N','C',...) gives sum over the inner
         ! index of C_rot(:,k) * conjg(C_w(:,k))^T.
         CALL ZGEMM('N', 'C', n, n, ws%nmo_active, &
                    (1.0_dp, 0.0_dp), ws%mos_rotated, n, &
                    ws%mos_weighted, n, &
                    (0.0_dp, 0.0_dp), ws%pmat_kstar, n)
      END IF

      IF (kpsym%rotp(is) < 0) ws%pmat_kstar = CONJG(ws%pmat_kstar)

      CALL dense_complex_to_dbcsr(ws%pmat_kstar, ws%atom_first, ws%atom_size, real_only, srpmat, scpmat)

      CALL timestop(handle)
   END SUBROUTINE spglib_mo_rotate_one_star

! **************************************************************************************************
!> \brief Dump the rotated KS coefficients C_rot(:, i_state) for one star
!>        member into a text file, matching the format that FHI-spglib's
!>        commented debug block (sym_base.f90:2299-2304, 2329-2333) writes
!>        for KS_scaled_cmplx. Enables a per-MO 1-to-1 diff between CP2K
!>        and FHI-spglib at the same (ik, is, i_state).
!>
!>        Output line format (one per (i_basis, i_state)):
!>            i_k_point  i_state  i_basis  Re(C_rot)  Im(C_rot)  |C_rot|
!>        (matching FHI-spglib's debug 'I10,I10,I10,2F20.14,F20.14').
!> \param ws       Workspace whose mos_rotated holds the just-computed C_rot
!> \param ik       IBZ k-point index (printed for cross-reference)
!> \param is       Star-member index
!> \param tag      File-name tag (e.g. 'cp2k' or 'sym')
! **************************************************************************************************
   SUBROUTINE dump_rotated_mos(ws, ik, is, tag)
      TYPE(spglib_mo_workspace_type), INTENT(IN)           :: ws
      INTEGER, INTENT(IN)                                :: ik, is
      CHARACTER(LEN=*), INTENT(IN)                       :: tag

      CHARACTER(LEN=128)                                 :: fname
      INTEGER                                            :: ibas, imo, lun
      COMPLEX(KIND=dp)                                   :: c

      WRITE (fname, '(A,A,A,I5.5,A,I3.3,A)') 'mos_dump_', TRIM(tag), '_ik', ik, '_is', is, '.txt'
      lun = 9877
      OPEN (unit=lun, file=fname, status='replace', action='write')
      WRITE (lun, '(A,A,A)') '# Rotated MO coefficients C_rot = Delta * C(k_IBZ) (', TRIM(tag), ')'
      WRITE (lun, '(A,I6,A,I6)') '# ik = ', ik, '  is = ', is
      WRITE (lun, '(A,I6)') '# n_basis_tot = ', ws%n_basis_tot
      WRITE (lun, '(A,I6,A,I6)') '# nmo = ', ws%nmo, '  nmo_active = ', ws%nmo_active
      WRITE (lun, '(A)') '# i_state  i_basis      Re(C_rot)            Im(C_rot)            |C_rot|'
      DO imo = 1, ws%nmo_active
         DO ibas = 1, ws%n_basis_tot
            c = ws%mos_rotated(ibas, imo)
            WRITE (lun, '(2I6,3E22.14)') imo, ibas, REAL(c, KIND=dp), AIMAG(c), ABS(c)
         END DO
      END DO
      CLOSE (lun)
   END SUBROUTINE dump_rotated_mos

! **************************************************************************************************
!> \brief Dump the IBZ MO coefficients C(k_IBZ) and the per-MO weights
!>        used in the MO-rotation path. One file per (ispin, ik). Format
!>        matches dump_rotated_mos so the same Python reader can parse
!>        both. Enables offline reconstruction:
!>            C_rot_python = Delta_python . C(k_IBZ)
!>        which can be diffed element-wise against C_rot dumped by CP2K.
! **************************************************************************************************
   SUBROUTINE dump_ibz_mos(ws, ik, tag)
      TYPE(spglib_mo_workspace_type), INTENT(IN)           :: ws
      INTEGER, INTENT(IN)                                :: ik
      CHARACTER(LEN=*), INTENT(IN)                       :: tag

      CHARACTER(LEN=128)                                 :: fname
      INTEGER                                            :: ibas, imo, lun
      COMPLEX(KIND=dp)                                   :: c

      WRITE (fname, '(A,A,A,I5.5,A)') 'mos_ibz_', TRIM(tag), '_ik', ik, '.txt'
      lun = 9878
      OPEN (unit=lun, file=fname, status='replace', action='write')
      WRITE (lun, '(A,A,A)') '# IBZ MO coefficients C(k_IBZ) before symmetry expansion (', TRIM(tag), ')'
      WRITE (lun, '(A,I6)') '# ik = ', ik
      WRITE (lun, '(A,I6)') '# n_basis_tot = ', ws%n_basis_tot
      WRITE (lun, '(A,I6,A,I6)') '# nmo = ', ws%nmo, '  nmo_active = ', ws%nmo_active
      WRITE (lun, '(A)') '# i_state    weight_n             '// &
                         '(then per-basis: i_state i_basis Re(C) Im(C) |C|)'
      DO imo = 1, ws%nmo
         WRITE (lun, '(A,I6,E22.14)') '# w ', imo, ws%weights(imo)
      END DO
      WRITE (lun, '(A)') '# i_state  i_basis      Re(C_IBZ)            Im(C_IBZ)            |C_IBZ|'
      DO imo = 1, ws%nmo
         DO ibas = 1, ws%n_basis_tot
            c = ws%mos_kibz(ibas, imo)
            WRITE (lun, '(2I6,3E22.14)') imo, ibas, REAL(c, KIND=dp), AIMAG(c), ABS(c)
         END DO
      END DO
      CLOSE (lun)
   END SUBROUTINE dump_ibz_mos

! **************************************************************************************************
!> \brief Dump the per-(ik, is) full-basis rotation matrix Delta_full(:,:)
!>        currently held in ws%delta_full. Call AFTER spglib_mo_rotate_one_star
!>        so delta_full has the just-built Delta for (ik, is). Format is
!>        dense complex element-by-element, matching FHI-spglib's
!>        Delta_matrix_full layout (sym_base.f90:2956 get_density_rotation).
! **************************************************************************************************
   SUBROUTINE dump_delta_full(ws, ik, is, tag)
      TYPE(spglib_mo_workspace_type), INTENT(IN)           :: ws
      INTEGER, INTENT(IN)                                :: ik, is
      CHARACTER(LEN=*), INTENT(IN)                       :: tag

      CHARACTER(LEN=128)                                 :: fname
      INTEGER                                            :: i, j, lun

      WRITE (fname, '(A,A,A,I5.5,A,I3.3,A)') 'delta_', TRIM(tag), '_ik', ik, '_is', is, '.txt'
      lun = 9879
      OPEN (unit=lun, file=fname, status='replace', action='write')
      WRITE (lun, '(A,A,A)') '# Full-basis rotation Delta(:,:) (', TRIM(tag), ')'
      WRITE (lun, '(A,I6,A,I6)') '# ik = ', ik, '  is = ', is
      WRITE (lun, '(A,I6)') '# n_basis_tot = ', ws%n_basis_tot
      WRITE (lun, '(A)') '# i_basis  j_basis      Re(Delta)            Im(Delta)'
      DO i = 1, ws%n_basis_tot
         DO j = 1, ws%n_basis_tot
            WRITE (lun, '(2I6,2E22.14)') i, j, &
               REAL(ws%delta_full(i, j), KIND=dp), AIMAG(ws%delta_full(i, j))
         END DO
      END DO
      CLOSE (lun)
   END SUBROUTINE dump_delta_full

! **************************************************************************************************
!> \brief Dump the symmetry metadata that the offline Python reconstruction
!>        of FHI-spglib's get_density_rotation needs. One file per IBZ k-point
!>        ik, one block per star member is. Includes:
!>            xkp_ibz, kpsym%xkp(:,is), kpsym%rotp(is), |rotp|, ira,
!>            kpsym%rot(:,:,is), kpsym%f0(:,is),
!>            kpoint%symop_vt(:,ira), kpoint%scoord(:,iatom),
!>            kpoint%kind_rotmat(ira, ikind)%rmat per atom kind.
!>
!>        Convention notes carried in the header:
!>          - kpsym%rot is crys_sym%rt = REAL(spglib_rotations) viewed in
!>            Fortran column-major order (so it equals R_spglib^T mathematically).
!>          - For TR (rotp<0) crys_sym already negated rt; xkp(:,is) is the
!>            star member coordinate (k_spatial for rotp>0; -k_spatial for
!>            rotp<0).
!>          - symop_vt(:,ira) is the spglib non-symmorphic translation V_C
!>            for that operation index.
! **************************************************************************************************
   SUBROUTINE dump_kpsym_meta(kpoint, kpsym, ik, xkp_ibz, atom_first, atom_size, n_basis_tot, tag)
      TYPE(kpoint_type), POINTER                         :: kpoint
      TYPE(kpoint_sym_type), POINTER                     :: kpsym
      INTEGER, INTENT(IN)                                :: ik, n_basis_tot
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: xkp_ibz
      INTEGER, DIMENSION(:), INTENT(IN)                  :: atom_first, atom_size
      CHARACTER(LEN=*), INTENT(IN)                       :: tag

      CHARACTER(LEN=128)                                 :: fname
      INTEGER                                            :: i, ia, ikind, ira, irot, is, j, jr, lun, &
                                                            natom, nkind
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rmat

      natom = SIZE(kpoint%atype)
      nkind = SIZE(kpoint%kind_rotmat, 2)

      WRITE (fname, '(A,A,A,I5.5,A)') 'kpsym_', TRIM(tag), '_ik', ik, '.txt'
      lun = 9880
      OPEN (unit=lun, file=fname, status='replace', action='write')
      WRITE (lun, '(A,A,A)') '# Symmetry metadata for offline FHI-spglib reconstruction (', TRIM(tag), ')'
      WRITE (lun, '(A,I6)') '# ik = ', ik
      WRITE (lun, '(A,3F18.12)') '# xkp_ibz = ', xkp_ibz
      WRITE (lun, '(A,I6)') '# nwght = ', kpsym%nwght
      WRITE (lun, '(A,I6)') '# natom = ', natom
      WRITE (lun, '(A,I6)') '# nkind = ', nkind
      WRITE (lun, '(A,I6)') '# n_basis_tot = ', n_basis_tot

      WRITE (lun, '(A)') '# atom_first (1-indexed first AO of each atom):'
      DO ia = 1, natom
         WRITE (lun, '(A,I4,I6)') '# atom_first ', ia, atom_first(ia)
      END DO
      WRITE (lun, '(A)') '# atom_size (number of AOs per atom):'
      DO ia = 1, natom
         WRITE (lun, '(A,I4,I6)') '# atom_size ', ia, atom_size(ia)
      END DO
      WRITE (lun, '(A)') '# atype (1-indexed kind id per atom):'
      DO ia = 1, natom
         WRITE (lun, '(A,I4,I6)') '# atype ', ia, kpoint%atype(ia)
      END DO

      WRITE (lun, '(A)') '# scoord (fractional atom coordinates):'
      DO ia = 1, natom
         WRITE (lun, '(A,I4,3F18.12)') '# scoord ', ia, kpoint%scoord(1:3, ia)
      END DO

      WRITE (lun, '(A)') '# ibrot (irot index for each ira slot):'
      DO ira = 1, SIZE(kpoint%ibrot)
         WRITE (lun, '(A,I4,I6)') '# ibrot ', ira, kpoint%ibrot(ira)
      END DO

      IF (ASSOCIATED(kpoint%symop_vt)) THEN
         WRITE (lun, '(A)') '# symop_vt(:,ira) (non-symmorphic translation for each ira):'
         DO ira = 1, SIZE(kpoint%ibrot)
            WRITE (lun, '(A,I4,3F18.12)') '# symop_vt ', ira, kpoint%symop_vt(1:3, ira)
         END DO
      END IF

      WRITE (lun, '(A)') '# kind_rotmat per (ira, ikind):'
      DO ira = 1, SIZE(kpoint%ibrot)
         DO ikind = 1, nkind
            IF (.NOT. ASSOCIATED(kpoint%kind_rotmat(ira, ikind)%rmat)) CYCLE
            rmat => kpoint%kind_rotmat(ira, ikind)%rmat
            WRITE (lun, '(A,I4,A,I4,A,I4,A,I4)') '# kind_rotmat ira=', ira, ' ikind=', ikind, &
               ' nrow=', SIZE(rmat, 1), ' ncol=', SIZE(rmat, 2)
            DO i = 1, SIZE(rmat, 1)
               DO j = 1, SIZE(rmat, 2)
                  WRITE (lun, '(2I4,E22.14)') i, j, rmat(i, j)
               END DO
            END DO
         END DO
      END DO

      WRITE (lun, '(A)') '# Per star member is = 1..nwght:'
      DO is = 1, kpsym%nwght
         irot = ABS(kpsym%rotp(is))
         ira = 0
         DO jr = 1, SIZE(kpoint%ibrot)
            IF (irot == kpoint%ibrot(jr)) ira = jr
         END DO
         WRITE (lun, '(A,I4)') '# === is = ', is
         WRITE (lun, '(A,I4,A,I4,A,I4)') '# rotp=', kpsym%rotp(is), ' |rotp|=', irot, ' ira=', ira
         WRITE (lun, '(A,3F18.12)') '# xkp(is) = ', kpsym%xkp(1:3, is)
         WRITE (lun, '(A)') '# rot(:,:,is) (Fortran column-major; mathematically equals R_spglib^T):'
         DO i = 1, 3
            WRITE (lun, '(A,3F18.12)') '# rot ', kpsym%rot(i, 1:3, is)
         END DO
         WRITE (lun, '(A)') '# f0(:,is) (atom permutation for this star member):'
         DO ia = 1, natom
            WRITE (lun, '(A,I4,I6)') '# f0 ', ia, kpsym%f0(ia, is)
         END DO
      END DO

      CLOSE (lun)
   END SUBROUTINE dump_kpsym_meta

! **************************************************************************************************
!> \brief Derive (atom_first, atom_size) from kind_rotmat sizes, so the
!>        i-th atom occupies AO rows atom_first(i):atom_first(i+1)-1.
!>        Order matches the DBCSR block ordering used in CP2K (atom-major,
!>        per-kind sgf order), which in turn matches the layout used by
!>        set_rotmat_basis in qs_basis_rotation_methods.F.
! **************************************************************************************************
   SUBROUTINE build_atom_basis_table(kpoint, natom, atom_first, atom_size, n_basis_tot)
      TYPE(kpoint_type), POINTER                         :: kpoint
      INTEGER, INTENT(IN)                                :: natom
      INTEGER, DIMENSION(natom + 1), INTENT(OUT)         :: atom_first
      INTEGER, DIMENSION(natom), INTENT(OUT)             :: atom_size
      INTEGER, INTENT(OUT)                               :: n_basis_tot

      INTEGER                                            :: i, ikind, ira, jra, nkind, nrot, nsgf
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_nsgf

      nrot = SIZE(kpoint%kind_rotmat, 1)
      nkind = SIZE(kpoint%kind_rotmat, 2)
      ALLOCATE (kind_nsgf(nkind))
      kind_nsgf = 0
      DO ikind = 1, nkind
         DO ira = 1, nrot
            IF (ASSOCIATED(kpoint%kind_rotmat(ira, ikind)%rmat)) THEN
               kind_nsgf(ikind) = SIZE(kpoint%kind_rotmat(ira, ikind)%rmat, 1)
               EXIT
            END IF
         END DO
         CPASSERT(kind_nsgf(ikind) > 0)
      END DO

      atom_first(1) = 1
      DO i = 1, natom
         atom_size(i) = kind_nsgf(kpoint%atype(i))
         atom_first(i + 1) = atom_first(i) + atom_size(i)
      END DO
      n_basis_tot = atom_first(natom + 1) - 1

      DEALLOCATE (kind_nsgf)
   END SUBROUTINE build_atom_basis_table

! **************************************************************************************************
!> \brief Construct the full-basis complex rotation matrix Delta(:,:) for one
!>        IBZ k-point ik and one star member is.
!>
!>        Mirrors FHI-spglib sym_base.f90:2956 (get_density_rotation):
!>        Delta(i_basis, j_basis) is non-zero only when atom_j == f0(atom_i)
!>        and shares the same kind. Within that 2-atom block the entry equals
!>            phase * kind_rotmat(ira, ikind)
!>        with
!>            k_rot   = R^{-T} k_IBZ                  (FHI-spglib line 2909)
!>            phase   = exp{ 2 pi i (k_rot . tau_i - k_IBZ . tau_{f0(i)}) }
!>
!>        Note on direction: kpsym%rot stores the (signed) crys_sym%rt which
!>        is R_spglib^T in column-major Fortran. The convention in CP2K is
!>        that R_spglib maps k_star -> k_IBZ, so its inverse-transpose
!>        applied to k_IBZ recovers the rotated star vector k_rot used in
!>        the FHI-spglib phase formula. Algebraically: k_rot = R . k_IBZ
!>        when R is the forward rotation (k_IBZ -> k_star), which equals
!>        the spglib inverse R_spglib^{-1}. Combined with the Fortran transpose
!>        already in kpsym%rot, we use kpsym%rot directly here.
!>
!>        Time reversal (kpsym%rotp(is) < 0): handled by k_star = -|k_rot|
!>        being supplied as kpsym%xkp(:,is); the imaginary-part sign flip
!>        comes from the antisymmetry of cpmat (Im[P(-k)] = -Im[P(k)]) which
!>        is realised here automatically because P(k_IBZ) is built with the
!>        correct sign and the rotation conjugates appropriately.
! **************************************************************************************************
   SUBROUTINE build_full_delta(kpoint, kpsym, is, ira, xkp_ibz, &
                               natom, atom_first, atom_size, n_basis_tot, delta_full)
      TYPE(kpoint_type), POINTER                         :: kpoint
      TYPE(kpoint_sym_type), POINTER                     :: kpsym
      INTEGER, INTENT(IN)                                :: is, ira, natom, n_basis_tot
      INTEGER, DIMENSION(natom + 1), INTENT(IN)          :: atom_first
      INTEGER, DIMENSION(natom), INTENT(IN)              :: atom_size
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: xkp_ibz
      COMPLEX(KIND=dp), DIMENSION(n_basis_tot, n_basis_tot), INTENT(OUT) :: delta_full

      INTEGER                                            :: iatom, i1, i2, ikind, jatom, j1, j2, jkind
      LOGICAL                                            :: non_symmorphic
      REAL(KIND=dp)                                      :: arg, arg1, arg2
      REAL(KIND=dp), DIMENSION(3)                        :: k_spatial, tau_i, tau_j
      COMPLEX(KIND=dp)                                   :: phase
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rblock

      delta_full = (0.0_dp, 0.0_dp)

      ! k_spatial = R . k_IBZ in fractional coords (the spatial part of
      ! the operation, NOT including time reversal). For ordinary star
      ! members (rotp>0) this equals kpsym%xkp(:,is). For TR partners
      ! (rotp<0) CP2K stores k_star = -k_spatial in kpsym%xkp(:,is); the
      ! spatial Delta uses +k_spatial in the phase formula, and the TR
      ! sign is applied later as CONJG(P_star).
      IF (kpsym%rotp(is) > 0) THEN
         k_spatial(1:3) = kpsym%xkp(1:3, is)
      ELSE
         k_spatial(1:3) = -kpsym%xkp(1:3, is)
      END IF

      DO iatom = 1, natom
         ikind = kpoint%atype(iatom)
         i1 = atom_first(iatom)
         i2 = atom_first(iatom + 1) - 1

         jatom = kpsym%f0(iatom, is)
         CPASSERT(jatom >= 1 .AND. jatom <= natom)
         jkind = kpoint%atype(jatom)
         CPASSERT(jkind == ikind) ! permutation must preserve atom kind
         j1 = atom_first(jatom)
         j2 = atom_first(jatom + 1) - 1

         tau_i(1:3) = kpoint%scoord(1:3, iatom)
         tau_j(1:3) = kpoint%scoord(1:3, jatom)
         arg2 = DOT_PRODUCT(k_spatial, tau_i)
         arg1 = DOT_PRODUCT(xkp_ibz, tau_j)
         arg = twopi*(arg2 - arg1)
         ! Cell-convention Δ̃: differs from FHI-spglib's atom-centered Δ̃
         ! because CP2K's `rskp_transform` uses cell-convention Bloch sums
         ! (φ_i(k,r) = Σ_R AO_i(r-τ_i-R) exp(2πi k·R), no atom-position
         ! phase). The empirically-correct phase is:
         !   * non-symmorphic op (V ≠ 0): -arg  (FHI-spglib phase conjugated)
         !   * symmorphic op (V = 0), atom FIXED by op (jatom = iatom): -arg
         !   * symmorphic op (V = 0), atom SWAPPED (jatom ≠ iatom): 0
         ! Calibrated against nosym BZ DMs at iter 1 of cc7 (8-atom asymmetric
         ! diamond) and cd (2-atom non-symmorphic): all 4 test cases (po,
         ! c2, cd, cc7) match nosym to machine precision (~1e-13 Ha).
         non_symmorphic = .FALSE.
         IF (ASSOCIATED(kpoint%symop_vt)) THEN
            non_symmorphic = SUM(ABS(kpoint%symop_vt(1:3, ira))) > 1.0e-9_dp
         END IF
         IF (non_symmorphic .OR. jatom == iatom) THEN
            arg = -arg
         ELSE
            arg = 0.0_dp
         END IF
         phase = CMPLX(COS(arg), SIN(arg), KIND=dp)

         rblock => kpoint%kind_rotmat(ira, ikind)%rmat
         CPASSERT(ASSOCIATED(rblock))
         CPASSERT(SIZE(rblock, 1) == atom_size(iatom))

         delta_full(i1:i2, j1:j2) = phase*rblock(1:atom_size(iatom), 1:atom_size(iatom))
      END DO
   END SUBROUTINE build_full_delta

! **************************************************************************************************
!> \brief Dump a sparsified DBCSR pair (rpmat real-symmetric, cpmat
!>        imag-antisymmetric) into a text file in dense element-by-element
!>        format. The dump reflects ONLY content actually stored in the
!>        DBCSRs — anything dropped by the sparsity is zero in the output.
!>        This is what the SCF accumulator (transform_dmat) actually sees.
!> \par   The file name encodes the k-point and a tag (e.g. 'sym' or
!>        'nosym') so sym and nosym dumps at the same k can be diffed
!>        offline. Same format as the previous dense-matrix dump so the
!>        same diff_dumps.py harness works.
! **************************************************************************************************
   SUBROUTINE dump_sparsified_pmat(rpmat, cpmat, real_only, atom_first, n_basis_tot, kpt, tag)
      TYPE(dbcsr_type), POINTER                          :: rpmat, cpmat
      LOGICAL, INTENT(IN)                                :: real_only
      INTEGER, DIMENSION(:), INTENT(IN)                  :: atom_first
      INTEGER, INTENT(IN)                                :: n_basis_tot
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: kpt
      CHARACTER(LEN=*), INTENT(IN)                       :: tag

      CHARACTER(LEN=1)                                   :: rmt
      CHARACTER(LEN=128)                                 :: fname
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: pmat
      INTEGER                                            :: i1, i2, icol, ii, irow, j1, j2, jj, &
                                                            kx, ky, kz, lun
      LOGICAL                                            :: is_sym
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: blk
      TYPE(dbcsr_iterator_type)                          :: iter

      ALLOCATE (pmat(n_basis_tot, n_basis_tot))
      pmat = (0.0_dp, 0.0_dp)

      ! Real part — iterate stored blocks of rpmat
      rmt = dbcsr_get_matrix_type(rpmat)
      is_sym = (rmt /= dbcsr_type_no_symmetry)
      CALL dbcsr_iterator_start(iter, rpmat)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, irow, icol, blk)
         i1 = atom_first(irow); i2 = i1 + SIZE(blk, 1) - 1
         j1 = atom_first(icol); j2 = j1 + SIZE(blk, 2) - 1
         pmat(i1:i2, j1:j2) = pmat(i1:i2, j1:j2) + CMPLX(blk(:, :), 0.0_dp, KIND=dp)
         IF (is_sym .AND. irow /= icol) THEN
            pmat(j1:j2, i1:i2) = pmat(j1:j2, i1:i2) + CMPLX(TRANSPOSE(blk(:, :)), 0.0_dp, KIND=dp)
         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)

      ! Imag part — iterate stored blocks of cpmat (antisymmetric)
      IF (.NOT. real_only) THEN
         rmt = dbcsr_get_matrix_type(cpmat)
         is_sym = (rmt /= dbcsr_type_no_symmetry)
         CALL dbcsr_iterator_start(iter, cpmat)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, irow, icol, blk)
            i1 = atom_first(irow); i2 = i1 + SIZE(blk, 1) - 1
            j1 = atom_first(icol); j2 = j1 + SIZE(blk, 2) - 1
            pmat(i1:i2, j1:j2) = pmat(i1:i2, j1:j2) + CMPLX(0.0_dp, blk(:, :), KIND=dp)
            IF (is_sym .AND. irow /= icol) THEN
               pmat(j1:j2, i1:i2) = pmat(j1:j2, i1:i2) - CMPLX(0.0_dp, TRANSPOSE(blk(:, :)), KIND=dp)
            END IF
         END DO
         CALL dbcsr_iterator_stop(iter)
      END IF

      ! Write to text file (same format as previous dense dumps)
      kx = NINT(kpt(1)*10000.0_dp)
      ky = NINT(kpt(2)*10000.0_dp)
      kz = NINT(kpt(3)*10000.0_dp)
      WRITE (fname, '(A,A,SP,3I7.5,SS,A)') 'pmat_dump_', TRIM(tag), kx, ky, kz, '.txt'
      lun = 9876
      OPEN (unit=lun, file=fname, status='replace', action='write')
      WRITE (lun, '(A,A,A)') '# Sparsified P (', TRIM(tag), ') — dense complex (zeros = dropped by DBCSR sparsity)'
      WRITE (lun, '(A,3F12.6)') '# k = ', kpt
      WRITE (lun, '(A,I6)') '# n_basis_tot = ', n_basis_tot
      DO ii = 1, n_basis_tot
         DO jj = 1, n_basis_tot
            WRITE (lun, '(2I6,2E20.12)') ii, jj, &
               REAL(pmat(ii, jj), KIND=dp), AIMAG(pmat(ii, jj))
         END DO
      END DO
      CLOSE (lun)

      DEALLOCATE (pmat)
   END SUBROUTINE dump_sparsified_pmat

! **************************************************************************************************
!> \brief Pull a (block-symmetric) DBCSR pair (rpmat real-symmetric,
!>        cpmat imag-antisymmetric) into a dense complex AO matrix.
!>        For symmetric DBCSR storage only the upper-triangle blocks are
!>        present; the lower triangle is mirrored here:
!>             Re(P)(j,i) = + Re(P)(i,j)
!>             Im(P)(j,i) = - Im(P)(i,j)
!>        Diagonal blocks of cpmat must vanish for true antisymmetry.
! **************************************************************************************************
   SUBROUTINE dbcsr_to_dense_complex(rpmat, cpmat, real_only, atom_first, atom_size, pmat)
      TYPE(dbcsr_type), POINTER                          :: rpmat, cpmat
      LOGICAL, INTENT(IN)                                :: real_only
      INTEGER, DIMENSION(:), INTENT(IN)                  :: atom_first, atom_size
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: pmat

      CHARACTER(LEN=1)                                   :: rmt_r, rmt_c
      INTEGER                                            :: i1, i2, icol, irow, j1, j2
      LOGICAL                                            :: rsym, csym
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: blk
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(mp_comm_type)                                 :: group

      pmat = (0.0_dp, 0.0_dp)

      rmt_r = dbcsr_get_matrix_type(rpmat)
      rsym = (rmt_r /= dbcsr_type_no_symmetry)

      ! Real part
      CALL dbcsr_iterator_start(iter, rpmat)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, irow, icol, blk)
         i1 = atom_first(irow); i2 = i1 + SIZE(blk, 1) - 1
         j1 = atom_first(icol); j2 = j1 + SIZE(blk, 2) - 1
         pmat(i1:i2, j1:j2) = pmat(i1:i2, j1:j2) + CMPLX(blk(:, :), 0.0_dp, KIND=dp)
         IF (rsym .AND. irow /= icol) THEN
            pmat(j1:j2, i1:i2) = pmat(j1:j2, i1:i2) + CMPLX(TRANSPOSE(blk(:, :)), 0.0_dp, KIND=dp)
         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)

      IF (.NOT. real_only) THEN
         rmt_c = dbcsr_get_matrix_type(cpmat)
         csym = (rmt_c /= dbcsr_type_no_symmetry)
         CALL dbcsr_iterator_start(iter, cpmat)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, irow, icol, blk)
            i1 = atom_first(irow); i2 = i1 + SIZE(blk, 1) - 1
            j1 = atom_first(icol); j2 = j1 + SIZE(blk, 2) - 1
            pmat(i1:i2, j1:j2) = pmat(i1:i2, j1:j2) + CMPLX(0.0_dp, blk(:, :), KIND=dp)
            IF (csym .AND. irow /= icol) THEN
               pmat(j1:j2, i1:i2) = pmat(j1:j2, i1:i2) - CMPLX(0.0_dp, TRANSPOSE(blk(:, :)), KIND=dp)
            END IF
         END DO
         CALL dbcsr_iterator_stop(iter)
      END IF

      ! Replicate the dense matrix across the DBCSR MPI group so every rank
      ! has the full P(k_IBZ). Each rank contributed only its locally-owned
      ! blocks above; the sum collects them on every rank.
      CALL dbcsr_get_info(rpmat, group=group)
      IF (group%num_pe > 1) CALL group%sum(pmat)
   END SUBROUTINE dbcsr_to_dense_complex

! **************************************************************************************************
!> \brief Reverse of dbcsr_to_dense_complex: write the upper-triangle blocks
!>        of (Re(pmat), Im(pmat)) into srpmat (symmetric) and scpmat
!>        (antisymmetric). Existing block sparsity of (srpmat, scpmat) is
!>        preserved; missing blocks in the dense matrix are silently dropped.
! **************************************************************************************************
   SUBROUTINE dense_complex_to_dbcsr(pmat, atom_first, atom_size, real_only, srpmat, scpmat)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: pmat
      INTEGER, DIMENSION(:), INTENT(IN)                  :: atom_first, atom_size
      LOGICAL, INTENT(IN)                                :: real_only
      TYPE(dbcsr_type), POINTER                          :: srpmat, scpmat

      INTEGER                                            :: i1, i2, icol, irow, j1, j2
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: blk
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL dbcsr_set(srpmat, 0.0_dp)
      IF (.NOT. real_only) CALL dbcsr_set(scpmat, 0.0_dp)

      CALL dbcsr_iterator_start(iter, srpmat)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, irow, icol, blk)
         i1 = atom_first(irow); i2 = i1 + SIZE(blk, 1) - 1
         j1 = atom_first(icol); j2 = j1 + SIZE(blk, 2) - 1
         blk(:, :) = REAL(pmat(i1:i2, j1:j2), KIND=dp)
      END DO
      CALL dbcsr_iterator_stop(iter)

      IF (.NOT. real_only) THEN
         CALL dbcsr_iterator_start(iter, scpmat)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, irow, icol, blk)
            i1 = atom_first(irow); i2 = i1 + SIZE(blk, 1) - 1
            j1 = atom_first(icol); j2 = j1 + SIZE(blk, 2) - 1
            blk(:, :) = AIMAG(pmat(i1:i2, j1:j2))
         END DO
         CALL dbcsr_iterator_stop(iter)
      END IF
   END SUBROUTINE dense_complex_to_dbcsr

END MODULE kpoint_density_spglib
