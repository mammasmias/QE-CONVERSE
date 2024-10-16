!
! Copyright (C) 2001-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


!-----------------------------------------------------------------------
! Compute the covariant derivative
!-----------------------------------------------------------------------
SUBROUTINE compute_dudk_new(dudk_method)
  USE kinds,         ONLY : dp
  USE wvfct,         ONLY : npw, nbnd, npwx
  USE klist,         ONLY : nks, lgauss 
  USE cell_base,     ONLY : tpiba
  USE io_global,     ONLY : stdout
  USE mp_global,     ONLY : my_pool_id, me_pool, root_pool
  USE gipaw_module,  ONLY : q_gipaw, alpha_pv, evq
  USE io_files,      ONLY : nwordwfc, diropn
  implicit none
  character(80), intent(in) :: dudk_method
  complex(dp), allocatable :: dudk(:,:,:)
  real(dp):: emin, emax
  logical :: exst
  integer :: ik, ipol, occ
  integer, parameter :: iundudk1 = 75, iundudk2 = 76, iundudk3 = 77
  
  call start_clock ('compute_dudk')

  call diropn(iundudk1, 'dudk1', 2*nwordwfc, exst)
  call diropn(iundudk2, 'dudk2', 2*nwordwfc, exst)
  call diropn(iundudk3, 'dudk3', 2*nwordwfc, exst)
  
  ! allocate covariant derivatives
    allocate (dudk(npwx,nbnd,3) )
    allocate(evq(npwx,nbnd))

  if (ik == 1) write(stdout,*)
  write(stdout,'(5X,''Computing du/dk '',$)')
  if (trim(dudk_method) == 'read') then
    write(stdout,'(''(read from disk)'')')

  elseif (trim(dudk_method) == 'covariant') then
    write(stdout,'(''(covariant derivative)'')') 
    do ik = 1, nks
    call find_nbnd_occ(ik, occ, emin, emax)
    write(stdout,'(5X,''k-point:'',I4,4X,''occ='',I3)') ik, occ
    call dudk_covariant(ik, occ, dudk)
    do ipol = 1, 3
        call davcio( dudk(1,1,ipol), 2*nwordwfc, iundudk1 + ipol - 1, ik, +1 )
    enddo
    enddo


  elseif (trim(dudk_method) == 'singlepoint') then
    write(stdout,'(''(single point)'')') 
    do ik = 1, nks
    call find_nbnd_occ(ik, occ, emin, emax)
    write(stdout,'(5X,''k-point:'',I4,4X,''occ='',I3)') ik, occ
    call dudk_covariant_single_point(ik, occ, dudk)
!    call dudk_covariant_single_point_old(ik, occ, dudk)
    do ipol = 1, 3
        call davcio( dudk(1,1,ipol), 2*nwordwfc, iundudk1 + ipol - 1, ik, +1 )
    enddo
    enddo

    
!  elseif (trim(dudk_method) == 'kdotp') then
!    write(stdout,'(''(k \dot p perturbation)'')') 
!    call find_nbnd_occ(ik, occ, emin, emax)
   ! determine alpha_pv and store it via gipaw_module
!    alpha_pv = emax - emin
!    if (.not. lgauss) alpha_pv = 2.d0 * (emax - emin)
!    alpha_pv = max(alpha_pv, 1.d-2)
!    if ( ( me_pool == root_pool ) .and. ( my_pool_id == 0) ) &
!        write(*,'(5X,''k-point_ik:'',I4,4X,''pool:'',I4,4X,''occ='',I3, 4X, &
!                     ''alpha_pv='',F10.4)') ik, my_pool_id+1, occ, alpha_pv
!    call dudk_kdotp(ik, occ, dudk)

  else
    write(stdout,*)
    call errore('compute_dudk', 'unknown du/dk method: '//trim(dudk_method), 1)
  endif

   write(stdout,'(5X,''done with du/dk'',3X,''q_gipaw = '',F8.6)') q_gipaw
   deallocate(dudk)

  call stop_clock ('compute_dudk')
END SUBROUTINE compute_dudk_new





!-----------------------------------------------------------------------
! covariant derivative
!-----------------------------------------------------------------------
  SUBROUTINE dudk_covariant(ik, occ, dudk)
  USE kinds,                ONLY : dp  
  USE cell_base,            ONLY : tpiba
  USE wvfct,                ONLY : npw, nbnd, npwx
  USE wavefunctions,        ONLY : evc
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE buffers,              ONLY : get_buffer, save_buffer
  USE mp_global,            ONLY : inter_pool_comm, intra_pool_comm
  USE mp,                   ONLY : mp_sum
  USE mp_bands,             ONLY : intra_bgrp_comm, inter_bgrp_comm
  USE gipaw_module,         ONLY : q_gipaw, evq
  IMPLICIT NONE
  !
  complex(dp), external :: zdotc
  complex(dp), allocatable :: overlap(:,:) 
  complex(dp) :: dudk(npwx,nbnd,3)
  real(dp) :: q_gipaw3(3), delta_k
  integer :: ik, i, occ, sig 
  integer :: ibnd, jbnd, ipol
  integer :: ibnd_start, ibnd_end

  if (occ == 0) return

  ! allocate/initialize  overlap matrix and dudk
  allocate(overlap(occ,occ))
  overlap (:,:) = (0.0_dp, 0.0_dp)   
  dudk(:,:,:) = (0.0_dp, 0.0_dp)

  ! read the wavefunction
  q_gipaw3(:) = 0.d0
  call compute_u_kq(ik, q_gipaw3)
  ! save the refreshed (now SO including) result on disk ( as evc(:,:) )
  call save_buffer(evc, nwordwfc, iunwfc, ik)
  delta_k = q_gipaw/tpiba

  ! loop over crystal directions
  do ipol = 1, 3
    dudk(1:npwx,1:occ,ipol) = (0.d0,0.d0)
    !
    ! loop over +/-1
    do sig = -1, 1, 2
      !
      ! set the k+q point
      q_gipaw3(:) = 0.d0
      q_gipaw3(ipol) = delta_k * sig
      call compute_u_kq(ik, q_gipaw3) ! gives evq(:,:) !!!

      ! compute overlaps <k+q|q> = <evc|evq>
      do ibnd = 1, occ
        do jbnd = 1, occ
          overlap(ibnd,jbnd) = zdotc(npw, evc(1,ibnd), 1, evq(1,jbnd), 1)
        enddo
      enddo

#ifdef __MPI
      call mp_sum( overlap, intra_bgrp_comm )   
#endif
      call invert_matrix(occ, overlap) 

      do ibnd = 1, occ
        do jbnd = 1, occ
          dudk(1:npw,ibnd,ipol) = dudk(1:npw,ibnd,ipol) + &
                        sig * 0.5d0/(delta_k*tpiba) * &
                        overlap(jbnd,ibnd) * evq(1:npw,jbnd)
        enddo
      enddo

    enddo ! sig
  enddo ! ipol
  deallocate(overlap)

END SUBROUTINE dudk_covariant





!-----------------------------------------------------------------------
! (k dot p) perturbation
!-----------------------------------------------------------------------
!SUBROUTINE dudk_kdotp(ik, occ, dudk)
!  USE kinds,                ONLY : dp  
!  USE cell_base,            ONLY : tpiba
!  USE wvfct,                ONLY : nbnd, npwx, current_k, igk 
!  USE wavefunctions_module, ONLY : evc
!  USE buffers,              ONLY : get_buffer, save_buffer
!  USE gvect,                ONLY : ngm, gcutm, g
!  USE gipaw_module,         ONLY : q_gipaw, evq, nbnd_occ
!  USE paw_gipaw,            ONLY : lambda_so
!  IMPLICIT NONE
!  complex(dp) :: dudk(npwx,nbnd,3), v_evc(npwx,nbnd)
!  real(dp) :: q_gipaw3(3)
!  integer :: ik, i, occ, jpol
!  integer :: ibnd, jbnd, ipol

!  if (occ == 0) return
  !
!  dudk(:,:,:) = (0.d0,0.d0)
!  q_gipaw3(:) = 0.d0
  !
  ! loop over crystal directions
!  !
!  do ipol = 1, 3
    ! apply the velocity operator
!    v_evc(:,:) = (0.d0,0.d0)
!    call apply_vel(evc, v_evc, ik, ipol, q_gipaw3)
!   call apply_vel_H(evc, v_evc, ik, ipol, q_gipaw3)
    ! solve (k dot p) system
!    call greenfunction(ik, v_evc, dudk(:,:,ipol), q_gipaw3)
!  enddo
   
!END SUBROUTINE dudk_kdotp





  
!-----------------------------------------------------------------------
! covariant derivative (single point)
!-----------------------------------------------------------------------
SUBROUTINE dudk_covariant_single_point(ik, occ, dudk)
  USE kinds,     ONLY : dp  
  USE cell_base, ONLY : bg, tpiba
  USE gvect,     ONLY : ngm, g
  USE klist,     ONLY : igk_k
  USE wvfct,     ONLY : npw, nbnd, npwx 
  USE wavefunctions,    ONLY : evc
!  USE gvecs  ,          ONLY : nls, nlsm
  USE fft_base,         ONLY : dffts
  USE fft_interfaces,   ONLY : fwfft, invfft
  USE mp_global,        ONLY : intra_pool_comm 
  USE mp_bands,         ONLY : intra_bgrp_comm, inter_bgrp_comm
  USE mp,               ONLY : mp_sum
  USE gipaw_module,     ONLY : evq
  USE io_files,         ONLY : nwordatwfc, iunsat, iunwfc, nwordwfc
  USE buffers,   ONLY : get_buffer, save_buffer
  implicit none
  complex(dp), external :: zdotc
  integer :: ik, i, occ, sig
  complex(dp), allocatable :: psir(:,:), aux(:), aux2(:,:)
  complex(dp) :: dudk(npwx,nbnd,3)
  integer :: ibnd, jbnd, ipol
  integer :: nls, nlm
  complex(dp), allocatable :: overlap(:,:)
  real(dp) :: bmod, q_gipaw3(3)
  
   
!!! nls = dffts%nl
!!! nlm = dffts%nlm
  if (occ == 0) return

!  allocate real space wfcs
  allocate ( psir(dffts%nnr,occ), aux(dffts%nnr), aux2(dffts%nnr,occ))
  call get_buffer(evc, nwordwfc, iunwfc, ik)
  ! transform the wfcs to real space
  do ibnd = 1, occ
    aux(:) = (0.d0, 0.d0)
    aux(dffts%nl(igk_k(1:npw,ik))) = evc(1:npw,ibnd)
    CALL invfft ('Wave', aux, dffts)
    psir(:,ibnd) = aux(:)
  enddo

  ! allocate/initialize  overlap matrix and dudk
  allocate( overlap(occ,occ) )
  overlap(:,:) = (0.0_dp, 0.0_dp)   
  dudk(:,:,:) = (0.0_dp, 0.0_dp)

  ! loop over crystal directions
  do ipol = 1, 3
    dudk(1:npwx,1:nbnd,ipol) = (0.d0,0.d0)

    ! loop over +/-1
    do sig = -1, 1, 2

      ! compute the covariant derivative
      do ibnd = 1, occ
        call multiply_iqr(sig, dffts, bg(:,ipol), psir(:,ibnd), aux2(:,ibnd))
        call fwfft ('Wave', aux2(:,ibnd), dffts)
        evq(:,ibnd) = aux2(dffts%nl(igk_k(1:npw,ik)),ibnd)
      enddo
  
      ! compute overlaps <evc|evq>
      do ibnd = 1, occ
        do jbnd = 1, occ
          overlap(ibnd,jbnd) = zdotc(npw, evc(1,ibnd), 1, evq(1,jbnd), 1)
        enddo
      enddo

#ifdef __MPI
      call mp_sum( overlap, intra_bgrp_comm )   
#endif
      call invert_matrix(occ, overlap) 
      ! 
      bmod = sqrt( sum(bg(1:3,ipol)**2.d0) ) * tpiba  
      !
      do ibnd = 1, occ
        do jbnd = 1, occ
          dudk(1:npw,ibnd,ipol) = dudk(1:npw,ibnd,ipol) + &
                        sig * 0.5d0/bmod * &
                        overlap(jbnd,ibnd) * evq(1:npw,jbnd)
        enddo
      enddo

    enddo ! sig
  enddo ! ipol
  deallocate(psir, aux, aux2, overlap)

END SUBROUTINE dudk_covariant_single_point



  
!-----------------------------------------------------------------------
! covariant derivative (single point)
!-----------------------------------------------------------------------
SUBROUTINE dudk_covariant_single_point_old(ik, occ, dudk)
  USE kinds,     ONLY : dp  
  USE cell_base, ONLY : bg, tpiba
  USE gvect,     ONLY : ngm, g
  USE klist,     ONLY : igk_k
  USE wvfct,     ONLY : npw, nbnd, npwx 
  USE wavefunctions, ONLY : evc
!  USE gvecs  ,   ONLY : nls, nlsm
  USE fft_base,         ONLY : dffts, dfftp
  USE fft_interfaces,   ONLY : fwfft, invfft
  USE io_files,             ONLY : nwordatwfc, iunsat, iunwfc, nwordwfc
!  USE ldaU,                 ONLY : lda_plus_u, swfcatom
  USE mp_global,            ONLY : intra_pool_comm  !, inter_pool_comm, nproc, npool
  USE mp_bands,             ONLY : intra_bgrp_comm, inter_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE buffers,   ONLY : get_buffer, save_buffer
  USE gipaw_module,     ONLY : evq
  implicit none
  integer :: ik, i, occ, sig, nrxxs 
  complex(dp), allocatable :: psir(:,:), aux(:), aux2(:,:)
  complex(dp) :: dudk(npwx,nbnd,3)
  integer :: ibnd, jbnd, ipol
  complex(dp), allocatable :: overlap(:,:)
  real(dp) :: bmod, q_gipaw3(3)

  if (occ == 0) return

  nrxxs = dffts%nnr

  ! allocate real space wfcs
  allocate ( psir(dffts%nnr,occ), aux(dffts%nnr), aux2(dffts%nnr,occ) )

  ! read the wavefunction (not necessary here, wfc for q=0 already stored in evc(:,:) 
   call get_buffer(evc, nwordwfc, iunwfc, ik)
  ! IF ( lda_plus_u ) CALL davcio( swfcatom, nwordatwfc, iunsat, ik, -1 )
  
  ! transform to real space
  do ibnd = 1, occ
    aux(1:nrxxs) = (0.d0, 0.d0)
    aux(dffts%nl(igk_k(1:npw,ik))) = evc(1:npw,ibnd)
    CALL invfft ('Wave', aux, dffts)
    psir(1:nrxxs,ibnd) = aux(1:nrxxs)
  enddo

  ! allocate overlap matrix
  allocate( overlap(occ,occ) )

  ! loop over crystal directions
  do ipol = 1, 3
    dudk(1:npwx,1:nbnd,ipol) = (0.d0,0.d0)

    ! loop over +/-1
    do sig = -1, 1, 2

      ! compute overlaps
      do ibnd = 1, occ
        do jbnd = 1, occ
!          call apply_eigr(sig, bg(:,ipol), psir(:,jbnd), aux)
          call multiply_iqr(sig, dfftp, bg(:,ipol), psir(:,jbnd), aux)
          overlap(ibnd,jbnd) = dot_product(psir(:,ibnd), aux(:))/ &
                               dble(dffts%nr1x*dffts%nr2x*dffts%nr3x)
        enddo
      enddo
#ifdef __MPI
      call mp_sum( overlap, intra_bgrp_comm )   
#endif
      
      call invert_matrix(occ, overlap)
      
      ! compute the covariant derivative
      do ibnd = 1, occ
!         call apply_eigr(sig, bg(:,ipol), psir(:,ibnd), aux2(:,ibnd))
        call multiply_iqr(sig, dffts, bg(:,ipol), psir(:,ibnd), aux2(:,ibnd))
        CALL fwfft ('Wave', aux2(:,ibnd), dffts)
      enddo
      ! 
      bmod = sum(bg(1:3,ipol)**2.d0 ) * tpiba      
      !
      do ibnd = 1, occ
        do jbnd = 1, occ
          dudk(1:npw,ibnd,ipol) = dudk(1:npw,ibnd,ipol) + &
                        sig * 0.5d0/bmod * &
                        overlap(jbnd,ibnd) * aux2(dffts%nl(igk_k(1:npw,ik)),jbnd)
        enddo
      enddo
      !
    enddo ! sig
  enddo ! ipol
  deallocate(psir, aux, aux2, overlap)

END SUBROUTINE dudk_covariant_single_point_old



!-----------------------------------------------------------------------
SUBROUTINE invert_matrix(n, a)
!-----------------------------------------------------------------------
  USE kinds, ONLY : dp
  implicit none
  integer :: n, ipiv(n), info, lwork
  complex(dp) :: a(n,n)  
  complex(dp), allocatable, dimension(:) :: work

  if (n == 0) return
  if (n == 1) then
    a(1,1) = 1.d0/a(1,1)
    return
  endif

  call ZGETRF(n, n, a, n, ipiv, info)
  if (info /= 0) call errore('invert_matrix', 'ZGETRF info=', info)

  allocate(work(1))
  lwork = -1
  call ZGETRI(n, a, n, ipiv, work, lwork, info)
  if (info /= 0) call errore('invert_matrix', 'ZGETRI(1) info=', info)

  lwork = nint(real(work(1)))
  deallocate(work)
  allocate(work(lwork))
  call ZGETRI(n, a, n, ipiv, work, lwork, info)
  if (info /= 0) call errore('invert_matrix', 'ZGETRI(2) info=', info)
END SUBROUTINE invert_matrix




!-----------------------------------------------------------------------
! apply e^(i G r) to a wavefunction in real space
!-----------------------------------------------------------------------
SUBROUTINE apply_eigr(sig, gv, psi, res)
  USE kinds,     ONLY : dp
  USE constants, ONLY : tpi
  USE fft_base,  ONLY : dffts, dfftp
  USE mp_global, ONLY : me_pool
  implicit none
  real(dp) :: gv(3), phase
  complex(dp) :: psi(dffts%nnr), res(dffts%nnr)
  integer :: sig, i1, i2, i3, itmp, izub, izlb, ir
  integer :: nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s

  nr1s = dffts%nr1;   nr2s = dffts%nr2;   nr3s = dffts%nr3
  nrx1s = dffts%nr1x; nrx2s = dffts%nr2x; nrx3s = dffts%nr3x
 
#ifdef __MPI
     izub=0
     do itmp=1,me_pool + 1
        izlb=izub+1
        izub=izub+dffts%nr3p(itmp)
     enddo
#else
     izlb=1
     izub=nr3s
#endif

  do i1 = 1, nr1s
    do i2 = 1, nr2s
      do i3 = izlb, izub
        phase = sig* tpi * ( gv(1) * dble(i1)/dble(nr1s) + &
                             gv(2) * dble(i2)/dble(nr2s) + &
                             gv(3) * dble(i3)/dble(nr3s) )
        itmp = i3 - izlb + 1
        ir = i1 + (i2-1)*nrx1s + (itmp-1)*nrx1s*nrx2s
        res(ir) = cmplx(cos(phase),sin(phase)) * psi(ir) 
      enddo
    enddo
  enddo

END SUBROUTINE apply_eigr
!-----------------------------------------------------------------------
SUBROUTINE multiply_iqr(sig, dfft, xq, psi, res)
!!
  !! Multiply exp(i*q*r) to func.
  !! Real space indexing is adapted from Modules/compute_dipole.f90
  !!
!----------------------------------------------------------------------------
    USE kinds,               ONLY : DP
    USE constants,           ONLY : tpi
    USE cell_base,           ONLY : at
    USE fft_types,           ONLY : fft_type_descriptor, fft_index_to_3d
    !
    IMPLICIT NONE
    !
    TYPE (fft_type_descriptor), INTENT(IN) :: dfft
    ! fft_type_descriptor. dffts or dfftp
    REAL(DP), INTENT(IN) :: xq(3)
    ! q vector in cartesian coordinate
    COMPLEX(DP) :: psi(dfft%nnr), res(dfft%nnr)
    ! input, output: real-space function
    !
    LOGICAL :: offrange
    INTEGER :: ir, i, j, k, ir_end, sig
    REAL(DP) :: arg, xq_cry(3)
    COMPLEX(DP) :: phase
    !
    xq_cry = xq
    CALL cryst_to_cart(1, xq_cry, at, -1)
    !
#if defined (__MPI)
    ir_end = MIN(dfft%nnr, dfft%nr1x*dfft%my_nr2p*dfft%my_nr3p)
#else
    ir_end = dfft%nnr
#endif
    !
    DO ir = 1, ir_end
      !
      CALL fft_index_to_3d(ir, dfft, i, j, k, offrange)
      IF ( offrange ) CYCLE
      !
      ! (i,j,k) is the zero-based coordinate of the real-space grid
      arg = sig *  tpi * (  xq_cry(1) * REAL(i, DP) / REAL(dfft%nr1, DP) &
                   + xq_cry(2) * REAL(j, DP) / REAL(dfft%nr2, DP) &
                   + xq_cry(3) * REAL(k, DP) / REAL(dfft%nr3, DP)  )
      phase = CMPLX( COS(arg), SIN(arg), kind=DP )
      !
      res(ir) = psi(ir) * phase
      !
    END DO ! ir
!----------------------------------------------------------------------------
  END SUBROUTINE multiply_iqr
  !----------------------------------------------------------------------------
  !



!-----------------------------------------------------------------------
SUBROUTINE find_nbnd_occ(ik, nbnd_occ, emin, emax)
!-----------------------------------------------------------------------
  USE kinds,     ONLY : dp
  USE wvfct,     ONLY : wg, nbnd
  USE klist,     ONLY : wk, lgauss, degauss, ngauss, two_fermi_energies
  USE constants, ONLY : pi
  USE pwcom
  implicit none
  integer, intent(in) :: ik
  integer, intent(out):: nbnd_occ
  real(dp), intent(out) :: emin, emax
  real(dp) :: small, xmax, fac, e_target
  integer :: ibnd

  IF ( two_fermi_energies ) THEN
     ef = (ef_up + ef_dw) / 2
  ENDIF

  if (lgauss) then ! metallic case
    small = 6.9626525973374d-5
    xmax = sqrt(-log(sqrt(pi)*small))
    if (ngauss == -99) then
      fac = 1.d0 / sqrt(small)
      xmax = 2.d0 * log(0.5d0*(fac+sqrt(fac*fac-4.d0)))
    endif
    e_target = ef + xmax * degauss
  endif

  nbnd_occ = 0
  emin = 1d6
  emax = -1d6
  do ibnd = 1, nbnd
    if (lgauss) then ! metallic case
      emin = min(emin, et(ibnd,ik))
      emax = e_target
      if (et(ibnd,ik) < e_target) nbnd_occ = ibnd
    elseif (wg(ibnd,ik) > 1d-4*wk(ik)) then ! insulator
      emin = min(emin, et(ibnd,ik))
      emax = max(emax, et(ibnd,ik))
      nbnd_occ = nbnd_occ + 1
    endif
  enddo

END SUBROUTINE find_nbnd_occ


!-----------------------------------------------------------------------
!SUBROUTINE apply_vel_H_old(psi, vel_psi, ik, ipol)
  !-----------------------------------------------------------------------
  !
  ! ... Apply the velocity operator
  ! ...   v = p + dV^{NL}_{k,k}/dk = dH_k/dk
  ! ...
  ! ... Here we use Hartree atomic units, so that:
  ! ...   V^{NL} => V^{NL} * ryd_to_hartree
  !-----------------------------------------------------------------------
!  USE kinds,                ONLY : DP
!  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
!  USE klist,                ONLY : xk
!  USE wvfct,                ONLY : nbnd, npwx, npw, igk, g2kin, ecutwfc, current_k
!  USE lsda_mod,             ONLY : current_spin, isk, lsda
!  USE becmod,               ONLY : bec_type, becp, calbec, &
!                                   allocate_bec_type, deallocate_bec_type
!  USE cell_base,            ONLY : tpiba, tpiba2
!  USE gipaw_module,         ONLY : q_gipaw
!  USE paw_gipaw,            ONLY : paw_vkb, lambda_so
!  USE uspp,                 ONLY : vkb, nkb
!  USE gvect,                ONLY : g, ngm

  !-- paramters ----------------------------------------------------------
!  IMPLICIT NONE
!  INTEGER, INTENT(IN) :: ipol       ! cartesian direction (1..3)
!  INTEGER, INTENT(IN) :: ik         ! k-point
!  COMPLEX(DP), INTENT(IN) :: psi(npwx,nbnd)
!  COMPLEX(DP), INTENT(OUT) :: vel_psi(npwx,nbnd)

  !-- local variables ----------------------------------------------------
!  real(dp), parameter :: ryd_to_hartree = 0.5d0
!  complex(dp), allocatable :: aux(:,:), vkb_save(:,:)
!  real(dp) :: dk, oldk(3)
!  integer :: isign

!  call start_clock('apply_vel')
!  vel_psi = (0.d0,0.d0)

  ! set dk (= delta_k ?)
!  dk = q_gipaw/tpiba  ! original value (taken from _so_gipaw)
 ! dk = q_gipaw/2.d0 
  
  ! allocate temporary arrays, save old NL-potential
!  allocate(aux(npwx,nbnd), vkb_save(npwx,nkb))
!  vkb_save = vkb
!  oldk(:) = xk(:,ik)

  !====================================================================
  ! compute (1/2|dk|) ( H_{k+dk} |psi> - H_{k-dk}|psi> )
  !====================================================================
!  call allocate_bec_type(nkb, nbnd, becp)
!  do isign = -1,1,2
!    xk(ipol,ik) = oldk(ipol) + isign * dk     ! k \pm dk

    ! set up H_{k \pm dk}
!    call init_us_2(npw, igk, xk(1,ik), vkb) !no_phase has no inflence!
!    if (any(lambda_so /= 0.d0)) &           !no_phase has no inflence!
!            call init_gipaw_2(npw, igk, xk(1,ik), paw_vkb)  
!    call g2_kin(ik)
!    call h_psi (npwx, npw, nbnd, psi, aux)
    !
!    vel_psi = vel_psi + dble(isign) * ryd_to_hartree * aux/(2.d0*dk*tpiba)
!  enddo

!  call deallocate_bec_type(becp)

  ! restore NL-potential at k
!  xk(:,ik) = oldk(:)
!  vkb = vkb_save
  
  ! free memory
!  deallocate(aux, vkb_save)

!END SUBROUTINE apply_vel_H_old

