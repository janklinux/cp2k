!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                  !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
MODULE qs_basis_rotation_methods
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE input_constants,                 ONLY: do_method_dftb
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: kpoint_sym_type,&
                                              kpoint_type
   USE orbital_pointers,                ONLY: nso
   USE orbital_transformation_matrices, ONLY: calculate_rotmat,&
                                              orbrotmat_type,&
                                              release_rotmat
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! Global parameters (only in this module)

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_basis_rotation_methods'

   ! Public subroutines

   PUBLIC :: qs_basis_rotation

CONTAINS

! **************************************************************************************************
!> \brief Cached check of CP2K_AIMS_DUMP_MOS env var. Mirrors dump_mos_active()
!>        in kpoint_methods.F so the per-rotation DBG WRITEs below stay quiet
!>        unless explicitly opted in.
! **************************************************************************************************
   LOGICAL FUNCTION dbg_spglib_active() RESULT(active)
      CHARACTER(LEN=8) :: env
      INTEGER          :: status
      INTEGER, SAVE    :: cached = -1
      IF (cached < 0) THEN
         CALL get_environment_variable('CP2K_AIMS_DUMP_MOS', env, status=status)
         IF (status == 0 .AND. (TRIM(env) == '1' .OR. TRIM(env) == 'true' .OR. TRIM(env) == 'TRUE')) THEN
            cached = 1
         ELSE
            cached = 0
         END IF
      END IF
      active = (cached == 1)
   END FUNCTION dbg_spglib_active

! **************************************************************************************************
!> \brief   Construct basis set rotation matrices
!> \param qs_env ...
!> \param kpoints ...
! **************************************************************************************************
   SUBROUTINE qs_basis_rotation(qs_env, kpoints)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(kpoint_type), POINTER                         :: kpoints

      INTEGER                                            :: ik, ikind, ir, ira, irot, jr, lval, &
                                                            nkind, nrot, i, j
      REAL(KIND=dp), DIMENSION(3, 3)                     :: rotmat, rcart, rtest
      REAL(KIND=dp), DIMENSION(5, 5)                     :: rtest5
      REAL(KIND=dp)                                      :: max_dev, trD2
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis
      TYPE(kpoint_sym_type), POINTER                     :: kpsym
      TYPE(orbrotmat_type), DIMENSION(:), POINTER        :: orbrot
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CPASSERT(ASSOCIATED(qs_env))
      CPASSERT(ASSOCIATED(kpoints))
      IF (ASSOCIATED(kpoints%kind_rotmat)) THEN
         CALL get_qs_env(qs_env, cell=cell)
         CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)
         CALL get_qs_kind_set(qs_kind_set, maxlgto=lval)
         nrot = SIZE(kpoints%kind_rotmat, 1)
         nkind = SIZE(kpoints%kind_rotmat, 2)
         IF (dbg_spglib_active()) THEN
            WRITE(0,'(A,2I5)') 'DBG qs_basis_rotation: nrot,nkind=', nrot, nkind
            WRITE(0,'(A,I5)') 'DBG lval=', lval
         END IF
         ! remove possible old rotation matrices
         DO irot = 1, nrot
            DO ikind = 1, nkind
               IF (ASSOCIATED(kpoints%kind_rotmat(irot, ikind)%rmat)) THEN
                  DEALLOCATE (kpoints%kind_rotmat(irot, ikind)%rmat)
               END IF
            END DO
         END DO
         ! check all rotations needed
         NULLIFY (orbrot)
         CALL get_qs_env(qs_env, dft_control=dft_control)
         DO ik = 1, kpoints%nkp
            kpsym => kpoints%kp_sym(ik)%kpoint_sym
            IF (kpsym%apply_symmetry) THEN
               DO irot = 1, SIZE(kpsym%rotp)
                  ir = kpsym%rotp(irot)
                  ira = 0
                  DO jr = 1, SIZE(kpoints%ibrot)
                     IF (ir == kpoints%ibrot(jr)) ira = jr
                  END DO
                  IF (ira > 0) THEN
                     IF (.NOT. ASSOCIATED(kpoints%kind_rotmat(ira, 1)%rmat)) THEN
                        ! kpsym%rot stores crys_sym%rt = REAL(spglib_rotations) where spglib_rotations
                        ! is C row-major, so in Fortran (column-major) it is R_spglib^T.
                        ! kpop convention: R_spglib[irot] maps k_star -> k_IBZ (backward).
                        ! The D-matrix needs R_cart_forward (k_IBZ -> k_star):
                        !   rcart_back = hmat * R_spglib * h_inv = hmat * TRANSPOSE(kpsym%rot) * h_inv
                        !   R_cart_forward = rcart_back^{-1} = TRANSPOSE(rcart_back)  (orthogonal matrix)
                        IF (dbg_spglib_active()) THEN
                           WRITE(0,'(A,2I4,A,9F8.4)') &
                              'DBG rot_raw ira=', ira, ik, ' rot=', ((kpsym%rot(i,j,irot),j=1,3),i=1,3)
                        END IF
                        rcart(1:3, 1:3) = MATMUL(cell%hmat, &
                                                  MATMUL(TRANSPOSE(kpsym%rot(:, :, irot)), cell%h_inv))
                        IF (dbg_spglib_active()) THEN
                           rtest = MATMUL(TRANSPOSE(rcart), rcart)
                           DO i = 1, 3; rtest(i, i) = rtest(i, i) - 1.0_dp; END DO
                           max_dev = MAXVAL(ABS(rtest))
                           WRITE(0,'(A,2I4,A,9F8.4,A,E10.3)') &
                              'DBG ira=', ira, ik, ' R_cart=', ((rcart(i,j),j=1,3),i=1,3), &
                              ' |R^TR-I|=', max_dev
                        END IF
                        ! Pass R_cart_forward = TRANSPOSE(rcart_back) to calculate_rotmat
                        rotmat(1:3, 1:3) = TRANSPOSE(rcart(1:3, 1:3))
                        CALL calculate_rotmat(orbrot, rotmat, lval)
                        IF (dbg_spglib_active()) THEN
                           IF (lval >= 1) THEN
                              WRITE(0,'(A,I3,A,9F7.3)') 'DBG D_orb(l=1) ira=', ira, &
                                 ' mat=', ((orbrot(1)%mat(i,j),j=1,3),i=1,3)
                              rtest(1:3,1:3) = MATMUL(TRANSPOSE(orbrot(1)%mat(1:3,1:3)), orbrot(1)%mat(1:3,1:3))
                              DO i=1,3; rtest(i,i)=rtest(i,i)-1.0_dp; END DO
                              WRITE(0,'(A,I3,A,E10.3)') 'DBG D_orb(l=1) ira=', ira, &
                                 ' |D^TD-I|=', MAXVAL(ABS(rtest(1:3,1:3)))
                              WRITE(0,'(A,I3,A,F8.4,A)') 'DBG D_orb(l=1) ira=', ira, &
                                 ' Tr=', orbrot(1)%mat(1,1)+orbrot(1)%mat(2,2)+orbrot(1)%mat(3,3), &
                                 ' (expect 0 for 120-deg, +3 for identity)'
                           END IF
                           IF (lval >= 2) THEN
                              WRITE(0,'(A,I3,A)') 'DBG D_orb(l=2) ira=', ira, ' full mat (rows 1-5):'
                              DO i = 1, 5
                                 WRITE(0,'(A,I1,A,5F8.4)') '  row', i, '=', (orbrot(2)%mat(i,j),j=1,5)
                              END DO
                              rtest5(1:5,1:5) = MATMUL(TRANSPOSE(orbrot(2)%mat(1:5,1:5)), orbrot(2)%mat(1:5,1:5))
                              DO i = 1, 5; rtest5(i,i) = rtest5(i,i) - 1.0_dp; END DO
                              WRITE(0,'(A,I3,A,E10.3)') 'DBG D_orb(l=2) ira=', ira, &
                                 ' |D^TD-I|=', MAXVAL(ABS(rtest5(1:5,1:5)))
                              trD2 = 0.0_dp
                              DO i = 1, 5; trD2 = trD2 + orbrot(2)%mat(i,i); END DO
                              WRITE(0,'(A,I3,A,F8.4,A)') 'DBG D_orb(l=2) ira=', ira, &
                                 ' Tr=', trD2, ' (expect -1 for 120-deg, +3 for identity)'
                           END IF
                        END IF
                        IF (dft_control%qs_control%method_id == do_method_dftb) THEN
                           CPABORT("ROTMAT")
                        ELSE
                           DO ikind = 1, nkind
                              CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis)
                              NULLIFY (kpoints%kind_rotmat(ira, ikind)%rmat)
                              CALL set_rotmat_basis(kpoints%kind_rotmat(ira, ikind)%rmat, orbrot, orb_basis)
                              IF (dbg_spglib_active()) THEN
                                 WRITE(0,'(A,I4,A,I4,A,I4,A,4F8.4)') &
                                    'DBG rmat ira=', ira, ' ikind=', ikind, &
                                    ' size=', SIZE(kpoints%kind_rotmat(ira,ikind)%rmat, 1), &
                                    ' diag(1:4)=', (kpoints%kind_rotmat(ira,ikind)%rmat(i,i), &
                                                    i=1,MIN(4,SIZE(kpoints%kind_rotmat(ira,ikind)%rmat,1)))
                              END IF
                           END DO
                        END IF
                     END IF
                  END IF
               END DO
            END IF
         END DO
         CALL release_rotmat(orbrot)
      END IF

   END SUBROUTINE qs_basis_rotation

! **************************************************************************************************
!> \brief ...
!> \param rmat ...
!> \param orbrot ...
!> \param basis ...
! **************************************************************************************************
   SUBROUTINE set_rotmat_basis(rmat, orbrot, basis)
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rmat
      TYPE(orbrotmat_type), DIMENSION(:), POINTER        :: orbrot
      TYPE(gto_basis_set_type), POINTER                  :: basis

      INTEGER                                            :: fs1, fs2, iset, ishell, l, nset, nsgf
      INTEGER, DIMENSION(:), POINTER                     :: nshell
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf, lshell

      CALL get_gto_basis_set(gto_basis_set=basis, nsgf=nsgf)
      ALLOCATE (rmat(nsgf, nsgf))
      rmat = 0.0_dp

      CALL get_gto_basis_set(gto_basis_set=basis, nset=nset, nshell=nshell, l=lshell, &
                             first_sgf=first_sgf)
      DO iset = 1, nset
         DO ishell = 1, nshell(iset)
            l = lshell(ishell, iset)
            fs1 = first_sgf(ishell, iset)
            fs2 = fs1 + nso(l) - 1
            rmat(fs1:fs2, fs1:fs2) = orbrot(l)%mat(1:nso(l), 1:nso(l))
         END DO
      END DO

   END SUBROUTINE set_rotmat_basis

END MODULE qs_basis_rotation_methods
