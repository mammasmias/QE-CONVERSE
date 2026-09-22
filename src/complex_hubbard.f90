! Copyright (C) 2026 QE-CONVERSE contributors
! Distributed under the GNU General Public License; see LICENSE.
!
MODULE complex_hubbard
  ! Use QE's complex U+V representation also for the on-site Dudarev
  ! functional. This keeps the QE module ABI and its complex Broyden mixer.
  USE kinds, ONLY : DP
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: check_complex_hubbard, prepare_complex_hubbard, report_complex_hubbard
CONTAINS
  SUBROUTINE check_complex_hubbard()
    USE ldaU, ONLY : lda_plus_u, lda_plus_u_kind, Hubbard_projectors, &
                    is_hubbard_back, orbital_resolved, hub_pot_fix, &
                    Hubbard_J0, Hubbard_beta
    USE ions_base, ONLY : ntyp => nsp
    USE control_flags, ONLY : gamma_only
    USE noncollin_module, ONLY : noncolin
    IMPLICIT NONE
    IF (.NOT. lda_plus_u) RETURN
    IF (gamma_only) CALL errore('complex_hubbard', &
      'Use an explicit complex k-point grid, not K_POINTS gamma', 1)
    IF (noncolin) CALL errore('complex_hubbard', &
      'Complex Hubbard support here requires collinear spins', 1)
    IF (Hubbard_projectors == 'pseudo') CALL errore('complex_hubbard', &
      'Pseudo Hubbard projectors are not supported', 1)
    IF (ANY(is_hubbard_back(1:ntyp))) CALL errore('complex_hubbard', &
      'Only one Hubbard manifold per atom is supported', 1)
    IF (lda_plus_u_kind /= 0 .AND. lda_plus_u_kind /= 2) &
      CALL errore('complex_hubbard', 'Only Dudarev U and U+V are supported', 1)
    IF (orbital_resolved .OR. hub_pot_fix) CALL errore('complex_hubbard', &
      'Orbital-resolved or fixed Hubbard potentials are not supported', 1)
    IF (ANY(Hubbard_J0(1:ntyp) /= 0._DP) .OR. &
        ANY(Hubbard_beta(1:ntyp) /= 0._DP)) CALL errore('complex_hubbard', &
      'This implementation supports U and V without J0 or beta', 1)
  END SUBROUTINE check_complex_hubbard

  SUBROUTINE prepare_complex_hubbard()
    ! Call AFTER potinit has read/initialized the original real U density,
    ! and BEFORE the first perturbed diagonalization or SCF mixer allocation.
    USE ldaU, ONLY : lda_plus_u, lda_plus_u_kind, Hubbard_U, Hubbard_V, &
                    init_hubbard, is_hubbard, ldim_u, ll, nwfcU, offsetU, &
                    nsg, nsgnew, v_nsg, eth, phase_fac
    USE ions_base, ONLY : nat, ityp, ntyp => nsp
    USE uspp_param, ONLY : upf
    USE lsda_mod, ONLY : nspin
    USE scf, ONLY : rho
    USE io_global, ONLY : stdout
    IMPLICIT NONE
    INTEGER :: na, nt, is, ld, viz, old_nwfc
    REAL(DP) :: initial_energy
    INTEGER, ALLOCATABLE :: old_offset(:)
    INTEGER, EXTERNAL :: find_viz

    CALL check_complex_hubbard()
    IF (.NOT. lda_plus_u) RETURN
    IF (lda_plus_u_kind == 2) RETURN

    IF (.NOT. ALLOCATED(rho%ns)) CALL errore('complex_hubbard', &
      'The initial on-site occupations have not been initialized', 1)
    old_nwfc = nwfcU
    initial_energy = eth
    ALLOCATE(old_offset(nat))
    old_offset = offsetU

    ! V(I,I,1) is U_I in QE's extended functional. All genuine intersite
    ! interactions remain exactly zero; no artificial V is introduced.
    Hubbard_V = 0._DP
    DO na = 1, nat
      Hubbard_V(na,na,1) = Hubbard_U(ityp(na))
    ENDDO
    lda_plus_u_kind = 2
    ! Rebuild neighborhoods, orbital labels, dimensions and complex arrays.
    ! Do not deallocate wfcU: existing projector buffers retain their layout.
    IF (ALLOCATED(ll)) DEALLOCATE(ll)
    CALL init_hubbard(upf(1:ntyp)%psd, nspin, .FALSE.)
    IF (nwfcU /= old_nwfc .OR. ANY(offsetU /= old_offset)) &
      CALL errore('complex_hubbard', 'U to U+V conversion changed projector ordering', 1)
    DEALLOCATE(old_offset)

    nsg = CMPLX(0._DP, 0._DP, KIND=DP)
    DO na = 1, nat
      nt = ityp(na)
      IF (.NOT. is_hubbard(nt)) CYCLE
      ld = ldim_u(nt)
      viz = find_viz(na,na)
      IF (viz <= 0) CALL errore('complex_hubbard', 'Missing self-site neighbor', 1)
      DO is = 1, nspin
        ! In array order both store conjg(p_row)*p_column. The initial
        ! new_ns matrix is real symmetric; subsequent nsg updates are complex.
        nsg(1:ld,1:ld,viz,na,is) = CMPLX(rho%ns(1:ld,1:ld,is,na), 0._DP, KIND=DP)
      ENDDO
    ENDDO
    nsgnew = nsg
    phase_fac = CMPLX(1._DP, 0._DP, KIND=DP)
    CALL v_hubbard_extended(nsg, v_nsg, eth)
    IF (ABS(eth-initial_energy) > 1.e-10_DP * MAX(1._DP, ABS(initial_energy))) &
      CALL errore('complex_hubbard', 'U to U+V conversion changed the initial Hubbard energy', 1)
    WRITE(stdout,'(/5X,A)') 'Complex Hubbard: on-site U initialized in the U+V representation (intersite V=0).'
  END SUBROUTINE prepare_complex_hubbard

  SUBROUTINE report_complex_hubbard()
    USE ldaU, ONLY : lda_plus_u, lda_plus_u_kind, nsg, ldim_u, is_hubbard
    USE ions_base, ONLY : nat, ityp
    USE lsda_mod, ONLY : nspin
    USE io_global, ONLY : stdout
    IMPLICIT NONE
    INTEGER :: na, is, ld, viz
    REAL(DP) :: hermitian_error
    INTEGER, EXTERNAL :: find_viz
    IF (.NOT. lda_plus_u) RETURN
    IF (lda_plus_u_kind /= 2) RETURN
    hermitian_error = 0._DP
    DO na = 1, nat
      IF (.NOT. is_hubbard(ityp(na))) CYCLE
      ld = ldim_u(ityp(na))
      viz = find_viz(na,na)
      DO is = 1, nspin
        hermitian_error = MAX(hermitian_error, MAXVAL(ABS(nsg(1:ld,1:ld,viz,na,is) - &
          CONJG(TRANSPOSE(nsg(1:ld,1:ld,viz,na,is))))))
      ENDDO
    ENDDO
    WRITE(stdout,'(5X,A,ES16.8)') 'Complex Hubbard: max |Im n| = ', MAXVAL(ABS(AIMAG(nsg)))
    WRITE(stdout,'(5X,A,ES16.8)') 'Complex Hubbard: onsite Hermiticity error = ', hermitian_error
    IF (hermitian_error > 1.e-10_DP) &
      CALL errore('complex_hubbard', 'Non-Hermitian on-site occupations after SCF', 1)
  END SUBROUTINE report_complex_hubbard
END MODULE complex_hubbard
