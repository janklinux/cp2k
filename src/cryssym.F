!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief K-points and crystal symmetry routines
!> \author jgh
! **************************************************************************************************
MODULE cryssym

   USE bibliography,                    ONLY: Togo2018,&
                                              cite_reference
   USE kinds,                           ONLY: dp
   USE spglib_f08,                      ONLY: spg_get_international,&
                                              spg_get_major_version,&
                                              spg_get_micro_version,&
                                              spg_get_minor_version,&
                                              spg_get_multiplicity,&
                                              spg_get_pointgroup,&
                                              spg_get_schoenflies,&
                                              spg_get_symmetry
   USE string_utilities,                ONLY: strip_control_codes
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: csym_type, release_csym_type, print_crys_symmetry, print_kp_symmetry
   PUBLIC :: crys_sym_gen, kpoint_gen

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cryssym'

! **************************************************************************************************
!> \brief CSM type
!> \par   Content:
!>
! **************************************************************************************************
   TYPE csym_type
      LOGICAL                                     :: symlib = .FALSE.
      LOGICAL                                     :: fullgrid = .FALSE.
      INTEGER                                     :: plevel = 0
      INTEGER                                     :: punit = -1
      INTEGER                                     :: istriz = -1
      REAL(KIND=dp)                               :: delta = 1.0e-8_dp
      REAL(KIND=dp), DIMENSION(3, 3)              :: hmat = 0.0_dp
      ! KPOINTS
      REAL(KIND=dp), DIMENSION(3)                 :: wvk0 = 0.0_dp
      INTEGER, DIMENSION(3)                       :: mesh = 0
      INTEGER                                     :: nkpoint = 0
      INTEGER                                     :: nat = 0
      INTEGER, DIMENSION(:), ALLOCATABLE          :: atype
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: scoord
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: xkpoint
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE    :: wkpoint
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: kpmesh
      INTEGER, DIMENSION(:, :), ALLOCATABLE       :: kplink
      INTEGER, DIMENSION(:), ALLOCATABLE          :: kpop
      !SPGLIB
      CHARACTER(len=11)                           :: international_symbol = ""
      CHARACTER(len=6)                            :: pointgroup_symbol = ""
      CHARACTER(len=10)                           :: schoenflies = ""
      INTEGER                                     :: n_operations = 0
      INTEGER, DIMENSION(:, :, :), ALLOCATABLE    :: rotations
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: translations
      !Symmetry group (populated from spglib data)
      REAL(KIND=dp), DIMENSION(3, 3, 48)          :: rt = 0.0_dp
      REAL(KIND=dp), DIMENSION(3, 48)             :: vt = 0.0_dp
      INTEGER, ALLOCATABLE, DIMENSION(:, :)       :: f0
      INTEGER                                     :: nrtot = 0
      INTEGER, DIMENSION(48)                      :: ibrot = 1
   END TYPE csym_type

CONTAINS

! **************************************************************************************************
!> \brief Release the CSYM type
!> \param csym  The CSYM type
! **************************************************************************************************
   SUBROUTINE release_csym_type(csym)
      TYPE(csym_type)                                    :: csym

      IF (ALLOCATED(csym%rotations)) THEN
         DEALLOCATE (csym%rotations)
      END IF
      IF (ALLOCATED(csym%translations)) THEN
         DEALLOCATE (csym%translations)
      END IF
      IF (ALLOCATED(csym%atype)) THEN
         DEALLOCATE (csym%atype)
      END IF
      IF (ALLOCATED(csym%scoord)) THEN
         DEALLOCATE (csym%scoord)
      END IF
      IF (ALLOCATED(csym%xkpoint)) THEN
         DEALLOCATE (csym%xkpoint)
      END IF
      IF (ALLOCATED(csym%wkpoint)) THEN
         DEALLOCATE (csym%wkpoint)
      END IF
      IF (ALLOCATED(csym%kpmesh)) THEN
         DEALLOCATE (csym%kpmesh)
      END IF
      IF (ALLOCATED(csym%kplink)) THEN
         DEALLOCATE (csym%kplink)
      END IF
      IF (ALLOCATED(csym%kpop)) THEN
         DEALLOCATE (csym%kpop)
      END IF
      IF (ALLOCATED(csym%f0)) THEN
         DEALLOCATE (csym%f0)
      END IF

   END SUBROUTINE release_csym_type

! **************************************************************************************************
!> \brief ...
!> \param csym ...
!> \param scoor ...
!> \param types ...
!> \param hmat ...
!> \param delta ...
!> \param iounit ...
! **************************************************************************************************
   SUBROUTINE crys_sym_gen(csym, scoor, types, hmat, delta, iounit)
      TYPE(csym_type)                                    :: csym
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: scoor
      INTEGER, DIMENSION(:), INTENT(IN)                  :: types
      REAL(KIND=dp), INTENT(IN)                          :: hmat(3, 3)
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: delta
      INTEGER, INTENT(IN), OPTIONAL                      :: iounit

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'crys_sym_gen'

      INTEGER                                            :: handle, ierr, major, micro, minor, nat, &
                                                            nop, tra_mat(3, 3)
      LOGICAL                                            :: spglib

      CALL timeset(routineN, handle)

      !..total number of atoms
      nat = SIZE(scoor, 2)
      csym%nat = nat

      ! output unit
      IF (PRESENT(iounit)) THEN
         csym%punit = iounit
      ELSE
         csym%punit = -1
      END IF

      ! accuracy for symmetry
      IF (PRESENT(delta)) THEN
         csym%delta = delta
      ELSE
         csym%delta = 1.e-6_dp
      END IF

      !..set cell values
      csym%hmat = hmat

      ! atom types
      ALLOCATE (csym%atype(nat))
      csym%atype(1:nat) = types(1:nat)

      ! scaled coordinates
      ALLOCATE (csym%scoord(3, nat))
      csym%scoord(1:3, 1:nat) = scoor(1:3, 1:nat)

      csym%n_operations = 0

      !..try spglib
      major = spg_get_major_version()
      minor = spg_get_minor_version()
      micro = spg_get_micro_version()
      IF (major == 0) THEN
         CALL cp_warn(__LOCATION__, "Symmetry library SPGLIB not available")
         spglib = .FALSE.
      ELSE
         spglib = .TRUE.
         CALL cite_reference(Togo2018)
         ierr = spg_get_international(csym%international_symbol, TRANSPOSE(hmat), scoor, types, nat, delta)
         IF (ierr == 0) THEN
            CALL cp_warn(__LOCATION__, "Symmetry Library SPGLIB failed")
            spglib = .FALSE.
         ELSE
            nop = spg_get_multiplicity(TRANSPOSE(hmat), scoor, types, nat, delta)
            ALLOCATE (csym%rotations(3, 3, nop), csym%translations(3, nop))
            csym%n_operations = nop
            ierr = spg_get_symmetry(csym%rotations, csym%translations, nop, &
                                    TRANSPOSE(hmat), scoor, types, nat, delta)
            ! Schoenflies Symbol
            csym%schoenflies = ' '
            ierr = spg_get_schoenflies(csym%schoenflies, TRANSPOSE(hmat), scoor, types, nat, delta)
            ! Point Group
            csym%pointgroup_symbol = ' '
            tra_mat = 0
            ierr = spg_get_pointgroup(csym%pointgroup_symbol, tra_mat, &
                                      csym%rotations, csym%n_operations)

            CALL strip_control_codes(csym%international_symbol)
            CALL strip_control_codes(csym%schoenflies)
            CALL strip_control_codes(csym%pointgroup_symbol)
         END IF
      END IF
      csym%symlib = spglib

      CALL timestop(handle)

   END SUBROUTINE crys_sym_gen

! **************************************************************************************************
!> \brief ...
!> \param csym ...
!> \param nk ...
!> \param symm ...
!> \param shift ...
!> \param full_grid ...
! **************************************************************************************************
   SUBROUTINE kpoint_gen(csym, nk, symm, shift, full_grid)
      TYPE(csym_type)                                    :: csym
      INTEGER, INTENT(IN)                                :: nk(3)
      LOGICAL, INTENT(IN), OPTIONAL                      :: symm
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: shift(3)
      LOGICAL, INTENT(IN), OPTIONAL                      :: full_grid

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'kpoint_gen'

      INTEGER                                            :: handle, i, ik, j, nkp, nkpts
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kpop, xptr
      LOGICAL                                            :: fullmesh
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: wkp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: xkp

      CALL timeset(routineN, handle)

      IF (PRESENT(shift)) THEN
         csym%wvk0 = shift
      ELSE
         csym%wvk0 = 0.0_dp
      END IF

      csym%istriz = -1
      IF (PRESENT(symm)) THEN
         IF (symm) csym%istriz = 1
      END IF

      IF (PRESENT(full_grid)) THEN
         fullmesh = full_grid
      ELSE
         fullmesh = .FALSE.
      END IF
      csym%fullgrid = fullmesh

      csym%nkpoint = 0
      csym%mesh(1:3) = nk(1:3)

      nkpts = nk(1)*nk(2)*nk(3)
      ALLOCATE (xkp(3, nkpts), wkp(nkpts), kpop(nkpts))
      ! kp: link
      ALLOCATE (csym%kplink(2, nkpts))
      csym%kplink = 0

      ! go through all the options
      IF (csym%symlib) THEN
         ! symmetry library is available
         IF (fullmesh) THEN
            ! full mesh requested
            CALL full_grid_gen(nk, xkp, wkp, shift)
            IF (csym%istriz == 1) THEN
               ! use inversion symmetry
               CALL inversion_symm(xkp, wkp, csym%kplink(1, :))
            ELSE
               ! full kpoint mesh is used
            END IF
         ELSE IF (csym%istriz /= 1) THEN
            ! use inversion symmetry
            CALL full_grid_gen(nk, xkp, wkp, shift)
            CALL inversion_symm(xkp, wkp, csym%kplink(1, :))
         ELSE
            ! use symmetry library to reduce k-points
            IF (SUM(ABS(csym%wvk0)) /= 0.0_dp) THEN
               CALL cp_abort(__LOCATION__, "MacDonald shifted k-point meshes are only "// &
                             "possible without symmetrization.")
            END IF

            CALL full_grid_gen(nk, xkp, wkp, shift)
            CALL kp_symmetry(csym, xkp, wkp, kpop)

         END IF
      ELSE
         ! no symmetry library is available
         CALL full_grid_gen(nk, xkp, wkp, shift)
         IF (csym%istriz /= 1 .AND. fullmesh) THEN
            ! full kpoint mesh is used
            DO i = 1, nkpts
               csym%kplink(1, i) = i
            END DO
         ELSE
            ! use inversion symmetry
            CALL inversion_symm(xkp, wkp, csym%kplink(1, :))
         END IF
      END IF
      ! count kpoints
      nkp = 0
      DO i = 1, nkpts
         IF (wkp(i) > 0.0_dp) nkp = nkp + 1
      END DO

      ! store reduced kpoint set
      csym%nkpoint = nkp
      ALLOCATE (csym%xkpoint(3, nkp), csym%wkpoint(nkp))
      ALLOCATE (xptr(nkp))
      j = 0
      DO ik = 1, nkpts
         IF (wkp(ik) > 0.0_dp) THEN
            j = j + 1
            csym%wkpoint(j) = wkp(ik)
            csym%xkpoint(1:3, j) = xkp(1:3, ik)
            xptr(j) = ik
         END IF
      END DO
      CPASSERT(j == nkp)

      ! kp: mesh
      ALLOCATE (csym%kpmesh(3, nkpts))
      csym%kpmesh(1:3, 1:nkpts) = xkp(1:3, 1:nkpts)

      ! kp: link
      DO ik = 1, nkpts
         i = csym%kplink(1, ik)
         DO j = 1, nkp
            IF (i == xptr(j)) THEN
               csym%kplink(2, ik) = j
               EXIT
            END IF
         END DO
      END DO
      DEALLOCATE (xptr)

      ! kp: operations
      ALLOCATE (csym%kpop(nkpts))
      IF (csym%symlib .AND. csym%istriz == 1 .AND. .NOT. fullmesh) THEN
         ! atomic symmetry operations possible
         csym%kpop(1:nkpts) = kpop(1:nkpts)
         DO ik = 1, nkpts
            CPASSERT(csym%kpop(ik) /= 0)
         END DO
      ELSE
         ! only time reversal symmetry
         DO ik = 1, nkpts
            IF (wkp(ik) > 0.0_dp) THEN
               csym%kpop(ik) = 1
            ELSE
               csym%kpop(ik) = 2
            END IF
         END DO
      END IF

      DEALLOCATE (xkp, wkp, kpop)

      CALL timestop(handle)

   END SUBROUTINE kpoint_gen

! **************************************************************************************************
!> \brief IBZ reduction by manually applying spglib rotation matrices to the full k-grid
!> Also populates csym%rt/vt/f0/nrtot/ibrot.
!> \param csym ...
!> \param xkp   full k-point grid from full_grid_gen (fractional reciprocal coordinates)
!> \param wkp   output multiplicities: >0 for IBZ representatives, 0 otherwise
!> \param kpop  output rotation index mapping each k-point to its representative
! **************************************************************************************************
   SUBROUTINE kp_symmetry(csym, xkp, wkp, kpop)
      TYPE(csym_type)                                    :: csym
      REAL(KIND=dp), DIMENSION(:, :)                     :: xkp
      REAL(KIND=dp), DIMENSION(:)                        :: wkp
      INTEGER, DIMENSION(:)                              :: kpop

      INTEGER                                            :: i, ia, irot, j, nat, nkpts, cp2k_rep
      REAL(KIND=dp)                                      :: dk(3), rk(3), det
      REAL(KIND=dp)                                      :: rt_invt(3, 3, 48)
      REAL(KIND=dp)                                      :: A(3, 3)

      nat = csym%nat
      nkpts = csym%mesh(1)*csym%mesh(2)*csym%mesh(3)

      ! --- populate symmetry group arrays from spglib rotation data ---
      csym%nrtot = csym%n_operations
      CPASSERT(csym%n_operations <= 48)
      csym%rt = 0.0_dp
      csym%vt = 0.0_dp
      csym%ibrot = 0
      rt_invt = 0.0_dp
      DO irot = 1, csym%n_operations
         csym%rt(:, :, irot) = REAL(csym%rotations(:, :, irot), KIND=dp)
         csym%vt(:, irot) = csym%translations(:, irot)
         csym%ibrot(irot) = irot
         ! Precompute R^{-T} for k-point transformation.
         ! spglib rotations act on direct-lattice fractional coords; k-points in
         ! reciprocal fractional coords transform as R^{-T} (= cofactor / det).
         A = csym%rt(:, :, irot)
         det = A(1, 1)*(A(2, 2)*A(3, 3) - A(2, 3)*A(3, 2)) &
             - A(1, 2)*(A(2, 1)*A(3, 3) - A(2, 3)*A(3, 1)) &
             + A(1, 3)*(A(2, 1)*A(3, 2) - A(2, 2)*A(3, 1))
         rt_invt(1, 1, irot) = (A(2, 2)*A(3, 3) - A(2, 3)*A(3, 2))/det
         rt_invt(1, 2, irot) = -(A(2, 1)*A(3, 3) - A(2, 3)*A(3, 1))/det
         rt_invt(1, 3, irot) = (A(2, 1)*A(3, 2) - A(2, 2)*A(3, 1))/det
         rt_invt(2, 1, irot) = -(A(1, 2)*A(3, 3) - A(1, 3)*A(3, 2))/det
         rt_invt(2, 2, irot) = (A(1, 1)*A(3, 3) - A(1, 3)*A(3, 1))/det
         rt_invt(2, 3, irot) = -(A(1, 1)*A(3, 2) - A(1, 2)*A(3, 1))/det
         rt_invt(3, 1, irot) = (A(1, 2)*A(2, 3) - A(1, 3)*A(2, 2))/det
         rt_invt(3, 2, irot) = -(A(1, 1)*A(2, 3) - A(1, 3)*A(2, 1))/det
         rt_invt(3, 3, irot) = (A(1, 1)*A(2, 2) - A(1, 2)*A(2, 1))/det
      END DO

      ! atom permutation table: f0(ia, irot) = image of atom ia under operation irot
      IF (ALLOCATED(csym%f0)) DEALLOCATE (csym%f0)
      ALLOCATE (csym%f0(nat, csym%n_operations))
      csym%f0 = 0
      DO irot = 1, csym%n_operations
         DO ia = 1, nat
            rk = MATMUL(csym%rt(:, :, irot), csym%scoord(:, ia)) + csym%translations(:, irot)
            DO j = 1, nat
               dk = rk - csym%scoord(:, j)
               dk = dk - ANINT(dk)
               IF (ALL(ABS(dk) < csym%delta)) THEN
                  csym%f0(ia, irot) = j
                  EXIT
               END IF
            END DO
         END DO
      END DO

      ! --- IBZ reduction: manual rotation loop ---
      ! Each k-point starts as its own representative.
      ! Scanning in index order guarantees the lowest-index point becomes the representative.
      DO i = 1, nkpts
         csym%kplink(1, i) = i
      END DO
      DO i = 1, nkpts
         IF (csym%kplink(1, i) /= i) CYCLE   ! already mapped to a lower-index representative
         DO irot = 1, csym%n_operations
            ! standard rotation: R^{-T} * k_i (correct for reciprocal fractional coords)
            rk = MATMUL(rt_invt(:, :, irot), xkp(:, i))
            rk = rk - ANINT(rk)              ! fold to [-0.5, 0.5)
            j = kpoint_cp2k_index(rk, csym%mesh)
            IF (j > i .AND. csym%kplink(1, j) == j) csym%kplink(1, j) = i
            ! time-reversal: -R^{-T} * k_i
            rk = -MATMUL(rt_invt(:, :, irot), xkp(:, i))
            rk = rk - ANINT(rk)
            j = kpoint_cp2k_index(rk, csym%mesh)
            IF (j > i .AND. csym%kplink(1, j) == j) csym%kplink(1, j) = i
         END DO
      END DO

      ! --- fill wkp (multiplicities) from kplink ---
      wkp = 0.0_dp
      DO i = 1, nkpts
         wkp(csym%kplink(1, i)) = wkp(csym%kplink(1, i)) + 1.0_dp
      END DO
      DO i = 1, nkpts
         IF (csym%kplink(1, i) /= i) wkp(i) = 0.0_dp
      END DO

      ! --- find symmetry operation index for each k-point (kpop) ---
      ! kpop(i): index into csym%rt of the operation mapping representative -> i;
      !          positive = rotation only, negative = rotation + time-reversal
      DO i = 1, nkpts
         cp2k_rep = csym%kplink(1, i)
         IF (cp2k_rep == i) THEN
            kpop(i) = 1
         ELSE
            kpop(i) = 1   ! fallback: identity
            DO irot = 1, csym%n_operations
               rk = MATMUL(rt_invt(:, :, irot), xkp(:, cp2k_rep))
               dk = rk - xkp(:, i)
               dk = dk - ANINT(dk)
               IF (ALL(ABS(dk) < 1.0e-10_dp)) THEN
                  kpop(i) = irot; EXIT
               END IF
               rk = -MATMUL(rt_invt(:, :, irot), xkp(:, cp2k_rep))
               dk = rk - xkp(:, i)
               dk = dk - ANINT(dk)
               IF (ALL(ABS(dk) < 1.0e-10_dp)) THEN
                  kpop(i) = -irot; EXIT
               END IF
            END DO
         END IF
      END DO

   END SUBROUTINE kp_symmetry

! **************************************************************************************************
!> \brief Compute the CP2K xkp linear index for a fractional k-point.
!>        Inverse of full_grid_gen: k(j) = (2*i_j - n_j - 1) / (2*n_j)
!>        Returns 0 if rk does not correspond to a grid point.
!> \param rk    fractional k-point coordinates, folded to [-0.5, 0.5)
!> \param mesh  k-point mesh dimensions
! **************************************************************************************************
   FUNCTION kpoint_cp2k_index(rk, mesh) RESULT(idx)
      REAL(KIND=dp), INTENT(IN)                          :: rk(3)
      INTEGER, INTENT(IN)                                :: mesh(3)
      INTEGER                                            :: idx, ix, iy, iz

      ix = (NINT(rk(1)*2.0_dp*REAL(mesh(1), dp)) + mesh(1) + 1)/2
      iy = (NINT(rk(2)*2.0_dp*REAL(mesh(2), dp)) + mesh(2) + 1)/2
      iz = (NINT(rk(3)*2.0_dp*REAL(mesh(3), dp)) + mesh(3) + 1)/2
      IF (ix < 1 .OR. ix > mesh(1) .OR. &
          iy < 1 .OR. iy > mesh(2) .OR. &
          iz < 1 .OR. iz > mesh(3)) THEN
         idx = 0
      ELSE
         idx = (ix - 1)*mesh(2)*mesh(3) + (iy - 1)*mesh(3) + iz
      END IF

   END FUNCTION kpoint_cp2k_index
! **************************************************************************************************
!> \brief ...
!> \param nk ...
!> \param xkp ...
!> \param wkp ...
!> \param shift ...
! **************************************************************************************************
   SUBROUTINE full_grid_gen(nk, xkp, wkp, shift)
      INTEGER, INTENT(IN)                                :: nk(3)
      REAL(KIND=dp), DIMENSION(:, :)                     :: xkp
      REAL(KIND=dp), DIMENSION(:)                        :: wkp
      REAL(KIND=dp), INTENT(IN)                          :: shift(3)

      INTEGER                                            :: i, ix, iy, iz
      REAL(KIND=dp)                                      :: kpt_latt(3)

      wkp = 0.0_dp
      i = 0
      DO ix = 1, nk(1)
         DO iy = 1, nk(2)
            DO iz = 1, nk(3)
               i = i + 1
               kpt_latt(1) = REAL(2*ix - nk(1) - 1, KIND=dp)/(2._dp*REAL(nk(1), KIND=dp))
               kpt_latt(2) = REAL(2*iy - nk(2) - 1, KIND=dp)/(2._dp*REAL(nk(2), KIND=dp))
               kpt_latt(3) = REAL(2*iz - nk(3) - 1, KIND=dp)/(2._dp*REAL(nk(3), KIND=dp))
               xkp(1:3, i) = kpt_latt(1:3)
               wkp(i) = 1.0_dp
            END DO
         END DO
      END DO
      DO i = 1, nk(1)*nk(2)*nk(3)
         xkp(1:3, i) = xkp(1:3, i) + shift(1:3)
      END DO

   END SUBROUTINE full_grid_gen

! **************************************************************************************************
!> \brief ...
!> \param xkp ...
!> \param wkp ...
!> \param link ...
! **************************************************************************************************
   SUBROUTINE inversion_symm(xkp, wkp, link)
      REAL(KIND=dp), DIMENSION(:, :)                     :: xkp
      REAL(KIND=dp), DIMENSION(:)                        :: wkp
      INTEGER, DIMENSION(:)                              :: link

      INTEGER                                            :: i, j, nkpts

      nkpts = SIZE(wkp, 1)

      link(:) = 0
      DO i = 1, nkpts
         IF (link(i) == 0) link(i) = i
         DO j = i + 1, nkpts
            IF (wkp(j) == 0) CYCLE
            IF (ALL(xkp(:, i) == -xkp(:, j))) THEN
               wkp(i) = wkp(i) + wkp(j)
               wkp(j) = 0.0_dp
               link(j) = i
               EXIT
            END IF
         END DO
      END DO

   END SUBROUTINE inversion_symm


! **************************************************************************************************
!> \brief ...
!> \param csym ...
! **************************************************************************************************
   SUBROUTINE print_crys_symmetry(csym)
      TYPE(csym_type)                                    :: csym

      INTEGER                                            :: i, iunit, j, plevel

      iunit = csym%punit
      IF (iunit >= 0) THEN
         plevel = csym%plevel
         WRITE (iunit, "(/,T2,A)") "Crystal Symmetry Information"
         IF (csym%symlib) THEN
            WRITE (iunit, '(A,T71,A10)') "       International Symbol: ", ADJUSTR(TRIM(csym%international_symbol))
            WRITE (iunit, '(A,T71,A10)') "       Point Group Symbol: ", ADJUSTR(TRIM(csym%pointgroup_symbol))
            WRITE (iunit, '(A,T71,A10)') "       Schoenflies Symbol: ", ADJUSTR(TRIM(csym%schoenflies))
            !
            WRITE (iunit, '(A,T71,I10)') "       Number of Symmetry Operations: ", csym%n_operations
            IF (plevel > 0) THEN
               DO i = 1, csym%n_operations
                  WRITE (iunit, '(A,i4,T51,3I10,/,T51,3I10,/,T51,3I10)') &
                     "           Rotation #: ", i, (csym%rotations(j, :, i), j=1, 3)
                  WRITE (iunit, '(T36,3F15.7)') csym%translations(:, i)
               END DO
            END IF
         ELSE
            WRITE (iunit, "(T2,A)") "SPGLIB for Crystal Symmetry Information determination is not availale"
         END IF
      END IF

   END SUBROUTINE print_crys_symmetry

! **************************************************************************************************
!> \brief ...
!> \param csym ...
! **************************************************************************************************
   SUBROUTINE print_kp_symmetry(csym)
      TYPE(csym_type), INTENT(IN)                        :: csym

      INTEGER                                            :: i, iunit, nat, plevel

      iunit = csym%punit
      IF (iunit >= 0) THEN
         plevel = csym%plevel
         WRITE (iunit, "(/,T2,A)") "K-point Symmetry Information"
         WRITE (iunit, '(A,T67,I14)') "       Number of Special K-points: ", csym%nkpoint
         WRITE (iunit, '(T19,A,T74,A)') " Wavevector Basis ", " Weight"
         DO i = 1, csym%nkpoint
            WRITE (iunit, '(T2,i10,3F10.5,T71,I10)') i, csym%xkpoint(1:3, i), NINT(csym%wkpoint(i))
         END DO
         WRITE (iunit, '(/,A,T63,3I6)') "       K-point Mesh: ", csym%mesh(1), csym%mesh(2), csym%mesh(3)
         WRITE (iunit, '(T19,A,T54,A)') " Wavevector Basis ", " Special Points    Rotation"
         DO i = 1, csym%mesh(1)*csym%mesh(2)*csym%mesh(3)
            WRITE (iunit, '(T2,i10,3F10.5,T45,3I12)') i, csym%kpmesh(1:3, i), &
               csym%kplink(1:2, i), csym%kpop(i)
         END DO
         IF (csym%nrtot > 0) THEN
            WRITE (iunit, '(/,A)') "       Atom Transformation Table"
            nat = SIZE(csym%f0, 1)
            DO i = 1, csym%nrtot
               WRITE (iunit, '(T10,A,I3,(T21,12I5))') " Rot=", csym%ibrot(i), csym%f0(1:nat, i)
            END DO
         END IF
      END IF

   END SUBROUTINE print_kp_symmetry

END MODULE cryssym
