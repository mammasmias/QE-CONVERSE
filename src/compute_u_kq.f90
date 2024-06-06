!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE compute_u_kq(ik, q)
  !----------------------------------------------------------------------------
  !
  ! ... diagonalize at k+q
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : RytoeV, tpi
  USE io_global,            ONLY : stdout
  USE io_files,             ONLY : nwordwfcU, iunhub, iunwfc, nwordwfc
  USE mp,                   ONLY : mp_sum
  USE mp_pools,             ONLY : inter_pool_comm, intra_pool_comm, me_pool
#ifdef __BANDS
  USE mp_bands,             ONLY : me_bgrp, inter_bgrp_comm
#endif
  USE klist,                ONLY : nkstot, nks, xk, ngk, igk_k
  USE uspp,                 ONLY : vkb, nkb
  USE wvfct,                ONLY : et, nbnd, npwx, g2kin, &
                                   current_k, nbndx, btype
  USE gvecw,                ONLY : gcutw
  USE control_flags,        ONLY : ethr, io_level, lscf, istep, max_cg_iter
  USE control_flags,        ONLY : cntrl_isolve => isolve, iverbosity
  USE ldaU,                 ONLY : lda_plus_u, wfcU
  USE lsda_mod,             ONLY : current_spin, lsda, isk
  USE noncollin_module,     ONLY : noncolin, npol
  USE wavefunctions,        ONLY : evc  
  USE gvect,                ONLY : g, ngm, ngl
  USE cell_base,            ONLY : at, bg, omega, tpiba, tpiba2
  USE bp,                   ONLY : lelfield
  USE becmod,               ONLY : becp
  USE random_numbers,       ONLY : randy
  USE buffers
  USE paw_gipaw,            ONLY : paw_vkb
  USE uspp_init,            ONLY : init_us_2
  USE gipaw_module
  IMPLICIT NONE
  INTEGER :: ik, iter       ! k-point, current iterations
  REAL(DP) :: q(3)          ! q-vector
  REAL(DP) :: avg_iter
  INTEGER :: ig, i, jk, npw
  REAL(DP) :: xkold(3)
  REAL(DP), allocatable :: et_old(:)
  REAL(DP) :: rr, arg
  logical :: q_is_zero

  CALL start_clock( 'c_bands' )

  !
  !
  !
  IF (.NOT.allocated(btype)) THEN
     !
     ! Initialize the diagonalization (done only once...)
     !
     if (isolve == 1) then
        nbndx = nbnd ! CG
     elseif (isolve == 0) then
        nbndx = 6*nbnd ! Davidson
     else
        call errore('compute_u_kq', &
                'Don''t even try to use this isolve!', abs(isolve))
     endif
     ! UWG: speed up the diagonalization: 
     !      if (allocated(btype)) deallocate(btype)
     !      allocate(btype(nbndx,nkstot))
     !      btype(1:nbndx,:) = 1
     ! btype(:,:) = 1, only if the corresponding state has to be
     !                 converged to full accuracy, zero otherwise!            
     allocate(btype(nbnd,nks))
     btype(:,:) = 0
     do jk = 1, nks
        btype(1:nbnd_occ(jk), jk) = 1
     enddo

!     btype(:,:) = 1

  ENDIF

  ! check if |q| is zero 
  q_is_zero = .false.
  if (sqrt(q(1)*q(1)+q(2)*q(2)+q(3)*q(3)) < 1d-8) q_is_zero = .true.

  !iter = 2   
  ! UWG: at least for q=0, we should use iter = 1 ...
  ! if ( q_is_zero ) iter = 1
  iter = 1
  istep = 0
  max_cg_iter = 200
  ethr = 1d-13
  lscf = .false.

  ! save eigenvalues
  allocate( et_old(nbnd) )
  et_old(:) = et(:,ik)

  !! debug
  if (iverbosity > 10) &
    WRITE(stdout, '(5X,"compute_u_kq: q = (",F10.4,",",F10.4,",",F10.4,")")') q

  avg_iter = 0.D0

  current_k = ik
  IF ( lsda ) current_spin = isk(ik)
  npw = ngk(ik)
  
  ! same sorting of G-vector at k+q
!  call gk_sort(xk(1,ik),ngm,g,gcutw,npw,igk_k(1,ik),g2kin)
  
  ! set the k-point
  xkold(:) = xk(:,ik)
  xk(:,ik) = xk(:,ik) + q(:)
  g2kin(1:npw) = ( ( xk(1,ik) + g(1,igk_k(1:npw,ik)) )**2 + &
                   ( xk(2,ik) + g(2,igk_k(1:npw,ik)) )**2 + &
                   ( xk(3,ik) + g(3,igk_k(1:npw,ik)) )**2 ) * tpiba2

  ! various initializations
  IF ( nkb > 0 ) CALL init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb )

!  !<UWG>                                                                                             
  if (any(lambda_so /= 0.d0)) &
      call init_gipaw_2(npw, igk_k(1,ik), xk(1,ik), paw_vkb)
!  !</UWG>                                                       

  ! read in wavefunctions from the previous iteration
  CALL get_buffer( evc, nwordwfc, iunwfc, ik)

  ! Needed for LDA+U
!  IF ( lda_plus_u ) CALL get_buffer( wfcU, nwordwfcU, iunhub, ik )

  ! randomize a little bit in case of CG diagonalization
!  if ( isolve == 1 ) then
!#ifdef __BANDS
!    rr = randy(ik+nks*me_bgrp) ! starting from a well defined k-dependent seed
!#else
!    rr = randy(ik+nks*me_pool) ! starting from a well defined k-dependent seed
!#endif
!    do i = 1, nbnd
!      do ig = 1, npw
!        rr = r_rand*(2.d0*randy() - 1.d0)
!        arg = tpi * randy()
!        evc(ig,i) = evc(ig,i)*CMPLX(1.d0+rr*cos(arg),rr*sin(arg),kind=DP)
!      enddo
!    enddo
!  endif


  ! diagonalization of bands for k-point ik
  call diag_bands_gipaw ( iter, ik, avg_iter )

  ! UWG: divison by nkstot not necessary (wrong!) here:
  ! avg_iter = avg_iter / nkstot 
  ! instead we have to collect and average the iterations with respect to pools   
  ! call mp_sum( avg_iter, inter_pool_comm )
  ! avg_iter = avg_iter / npool
  ! alternatively, we can leave out both...

  !! debug
  if (iverbosity > 20) &
    write(stdout,'(5X,"ethr = ",1PE9.2,",  avg # of iterations =",0PF5.1)') &
         ethr, avg_iter

  ! check if diagonalization was ok
  if (iverbosity > 20) then
    write(stdout,'(5X,''eigenvalues at k:'')')
    write(stdout,'(8F9.4)') et_old(1:nbnd)*RytoeV
    write(stdout,'(5X,''eigenvalues at k+q:'')')
    write(stdout,'(8F9.4)') et(1:nbnd,ik)*RytoeV
  endif

  do i = 1, nbnd
    if (abs(et(i,ik) - et_old(i))*RytoeV > 0.2d0) then
      write(stdout,'(5X,''ATTENTION: ik='',I4,''  ibnd='',I3,$)') ik, i
      write(stdout,'(2X,''eigenvalues differ too much!'')')
      write(stdout,'(5X,2(F10.4,2X))') et_old(i)*RytoeV, et(i,ik)*RytoeV
    endif
  enddo

  ! restore the k-point and eigenvalues
  xk(:,ik) = xkold(:)
  et(:,ik) = et_old(:)
  deallocate(et_old)


  ! set wfct{k+q} and restore wavefunctions ( q = 0 )
  evq = evc
  call start_clock('buffer_IO')
  CALL get_buffer(evc, nwordwfc, iunwfc, ik)
  call stop_clock('buffer_IO')

  ! Needed for LDA+U
  IF ( lda_plus_u ) CALL get_buffer( wfcU, nwordwfcU, iunhub, ik )


  CALL stop_clock( 'c_bands' )  
  RETURN

END SUBROUTINE compute_u_kq
