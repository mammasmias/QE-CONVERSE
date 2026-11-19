! Copyright (C) 2001-2026 Quantum ESPRESSO / QE-CONVERSE contributors
! Distributed under the GNU General Public License; see LICENSE.
!
! Modern form of the legacy velocity + Green-function algorithm.
! Solve (H_k - eps_n + alpha P) du_n = -Q (dH_k/dk) u_n.
! P spans the same occ bands as the covariant finite-difference method.
! All energies, including dH/dk and alpha, are in Ry (no Hartree conversion).
MODULE converse_kdotp
  USE kinds, ONLY : dp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dudk_kdotp
  INTEGER :: occupied
  REAL(dp) :: projector_shift
CONTAINS
  SUBROUTINE dudk_kdotp(ik, occ, dudk)
    USE wvfct, ONLY : npwx, npw, nbnd, g2kin, current_k
    USE wavefunctions, ONLY : evc
    USE klist, ONLY : ngk, xk
    USE lsda_mod, ONLY : lsda, current_spin, isk
    USE gipaw_module, ONLY : q_gipaw, conv_threshold, evq
    USE buffers, ONLY : save_buffer
    USE io_files, ONLY : iunwfc, nwordwfc
    USE cell_base, ONLY : tpiba
    USE becmod, ONLY : becp, allocate_bec_type, deallocate_bec_type
    USE uspp, ONLY : nkb, okvan
    USE noncollin_module, ONLY : npol
    USE control_flags, ONLY : gamma_only
    USE mp_bands, ONLY : intra_bgrp_comm
    USE mp, ONLY : mp_sum
    USE io_global, ONLY : stdout
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ik, occ
    COMPLEX(dp), INTENT(OUT) :: dudk(npwx,nbnd,3)
    COMPLEX(dp), ALLOCATABLE :: rhs(:,:), hpsi(:,:)
    REAL(dp), ALLOCATABLE :: precond(:,:), energies(:)
    REAL(dp) :: k0(3), shifted(3), zero(3), eprec, threshold, residual
    INTEGER :: ipol, sig, ibnd, iterations
    LOGICAL :: converged
    EXTERNAL :: cg_psi, cgsolve_all

    dudk = (0.0_dp,0.0_dp)
    IF (gamma_only .OR. npol /= 1 .OR. okvan) &
      CALL errore('dudk_kdotp','requires complex collinear norm-conserving wavefunctions',1)
    IF (q_gipaw <= 0.0_dp .OR. conv_threshold <= 0.0_dp) &
      CALL errore('dudk_kdotp','q_gipaw and conv_threshold must be positive',1)
    IF (occ == 0) RETURN
    occupied = occ
    current_k = ik
    npw = ngk(ik)
    IF (lsda) current_spin = isk(ik)
    k0 = xk(:,ik)
    zero = 0.0_dp
    ! Refresh the reference states with the same GIPAW Hamiltonian as covariant.
    CALL compute_u_kq(ik, zero)
    evc = evq
    CALL save_buffer(evc, nwordwfc, iunwfc, ik)
    CALL allocate_bec_type(nkb, nbnd, becp)
    ALLOCATE(rhs(npwx,occ), hpsi(npwx,occ), precond(npwx,occ), energies(occ))
    CALL prepare_k(ik,k0)
    CALL h_psi_gipaw(npwx,npw,occ,evc,hpsi)
    DO ibnd = 1, occ
      energies(ibnd) = REAL(SUM(CONJG(evc(1:npw,ibnd))*hpsi(1:npw,ibnd)),dp)
    ENDDO
    CALL mp_sum(energies,intra_bgrp_comm)
    projector_shift = MAX(2.0_dp*(MAXVAL(energies)-MINVAL(energies)),1.0_dp)
    precond = 0.0_dp
    DO ibnd = 1, occ
      eprec = 1.35_dp*SUM(g2kin(1:npw)*ABS(evc(1:npw,ibnd))**2)
      CALL mp_sum(eprec,intra_bgrp_comm)
      eprec = MAX(eprec,1.0_dp)
      precond(1:npw,ibnd) = 1.0_dp/MAX(1.0_dp,g2kin(1:npw)/eprec)
    ENDDO
    threshold = MIN(SQRT(conv_threshold),1.0e-8_dp)
    DO ipol = 1, 3
      ! Central difference of the FULL GIPAW Hamiltonian, not of eigenvectors.
      ! This includes k-dependent KB, SO, NMR and Hubbard projectors.
      rhs = (0.0_dp,0.0_dp)
      DO sig = -1, 1, 2
        shifted = k0
        shifted(ipol) = shifted(ipol) + sig*q_gipaw/tpiba
        CALL prepare_k(ik,shifted)
        hpsi = (0.0_dp,0.0_dp)
        CALL h_psi_gipaw(npwx,npw,occ,evc,hpsi)
        rhs(1:npw,:) = rhs(1:npw,:) + sig*hpsi(1:npw,:)/(2.0_dp*q_gipaw)
      ENDDO
      CALL prepare_k(ik,k0)
      CALL project_out(npw,occ,rhs)
      rhs = -rhs
      CALL cgsolve_all(sternheimer_operator,cg_psi,energies,rhs,dudk(:,:,ipol), &
           precond,npwx,npw,threshold,ik,iterations,converged,residual,occ,1)
      WRITE(stdout,'(5X,"kdotp direction ",I1,": iterations=",I5," residual=",ES12.4)') &
           ipol,iterations,residual
      IF (.NOT. converged) CALL errore('dudk_kdotp','Sternheimer solver did not converge',ipol)
      CALL project_out(npw,occ,dudk(:,:,ipol))
    ENDDO
    ! Return projectors, kinetic energy, Hubbard phases and xk at the original k.
    CALL prepare_k(ik,k0)
    DEALLOCATE(rhs,hpsi,precond,energies)
    CALL deallocate_bec_type(becp)
  END SUBROUTINE dudk_kdotp

  SUBROUTINE prepare_k(ik,k)
    USE wvfct, ONLY : npw, npwx, nbnd, g2kin
    USE klist, ONLY : xk, igk_k
    USE gvect, ONLY : g
    USE cell_base, ONLY : tpiba2
    USE uspp, ONLY : nkb, vkb
    USE uspp_init, ONLY : init_us_2
    USE gipaw_module, ONLY : lambda_so
    USE paw_gipaw, ONLY : paw_vkb
    USE ldaU, ONLY : lda_plus_u, lda_plus_u_kind, Hubbard_projectors
    USE basis, ONLY : wfcatom, swfcatom, natomwfc
    USE becmod, ONLY : becp, allocate_bec_type, deallocate_bec_type
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ik
    REAL(dp), INTENT(IN) :: k(3)
    xk(:,ik) = k
    g2kin(1:npw) = ((k(1)+g(1,igk_k(1:npw,ik)))**2 + &
                   (k(2)+g(2,igk_k(1:npw,ik)))**2 + &
                   (k(3)+g(3,igk_k(1:npw,ik)))**2)*tpiba2
    IF (nkb > 0) CALL init_us_2(npw,igk_k(1,ik),xk(1,ik),vkb)
    IF (ANY(lambda_so /= 0.0_dp)) CALL init_gipaw_2(npw,igk_k(1,ik),xk(1,ik),paw_vkb)
    IF (lda_plus_u .AND. lda_plus_u_kind == 2) CALL phase_factor(ik)
    IF (lda_plus_u .AND. Hubbard_projectors /= 'pseudo') THEN
      ! orthoUwfc_k owns the global becp temporarily (natomwfc columns).
      ! Release our response workspace, then restore it with nbnd columns.
      CALL deallocate_bec_type(becp)
      ALLOCATE(wfcatom(npwx,natomwfc),swfcatom(npwx,natomwfc))
      CALL orthoUwfc_k(ik,.FALSE.)
      DEALLOCATE(wfcatom,swfcatom)
      CALL allocate_bec_type(nkb,nbnd,becp)
    ENDIF
  END SUBROUTINE prepare_k

  SUBROUTINE project_out(n,m,psi)
    USE wvfct, ONLY : npwx
    USE wavefunctions, ONLY : evc
    USE mp_bands, ONLY : intra_bgrp_comm
    USE mp, ONLY : mp_sum
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n,m
    COMPLEX(dp), INTENT(INOUT) :: psi(npwx,m)
    COMPLEX(dp) :: overlap(occupied,m)
    CALL zgemm('C','N',occupied,m,n,(1.0_dp,0.0_dp),evc,npwx,psi,npwx, &
               (0.0_dp,0.0_dp),overlap,occupied)
    CALL mp_sum(overlap,intra_bgrp_comm)
    CALL zgemm('N','N',n,m,occupied,(-1.0_dp,0.0_dp),evc,npwx,overlap,occupied, &
               (1.0_dp,0.0_dp),psi,npwx)
    psi(n+1:npwx,:) = (0.0_dp,0.0_dp)
  END SUBROUTINE project_out

  SUBROUTINE sternheimer_operator(n,psi,hpsi,energies,ik,m)
    USE wvfct, ONLY : npwx
    USE wavefunctions, ONLY : evc
    USE mp_bands, ONLY : intra_bgrp_comm
    USE mp, ONLY : mp_sum
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n,ik,m
    REAL(dp), INTENT(IN) :: energies(m)
    COMPLEX(dp), INTENT(IN) :: psi(npwx,m)
    COMPLEX(dp), INTENT(OUT) :: hpsi(npwx,m)
    COMPLEX(dp) :: overlap(occupied,m)
    INTEGER :: ibnd
    hpsi = (0.0_dp,0.0_dp)
    IF (m == 0) RETURN
    CALL h_psi_gipaw(npwx,n,m,psi,hpsi)
    DO ibnd = 1,m
      hpsi(1:n,ibnd) = hpsi(1:n,ibnd)-energies(ibnd)*psi(1:n,ibnd)
    ENDDO
    CALL zgemm('C','N',occupied,m,n,(1.0_dp,0.0_dp),evc,npwx,psi,npwx, &
               (0.0_dp,0.0_dp),overlap,occupied)
    CALL mp_sum(overlap,intra_bgrp_comm)
    CALL zgemm('N','N',n,m,occupied,CMPLX(projector_shift,0.0_dp,dp),evc,npwx, &
               overlap,occupied,(1.0_dp,0.0_dp),hpsi,npwx)
    hpsi(n+1:npwx,:) = (0.0_dp,0.0_dp)
  END SUBROUTINE sternheimer_operator
END MODULE converse_kdotp
