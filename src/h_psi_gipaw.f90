!
! Copyright (C) 2002-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE h_psi_gipaw( lda, n, m, psi, hpsi )
  !----------------------------------------------------------------------------
  !! This routine computes the product of the Hamiltonian matrix with m 
  !! wavefunctions contained in psi.
  !
  !! \(\textit{Wrapper routine}\): performs bgrp parallelization on 
  !! non-distributed bands. If suitable and required, calls old H\psi 
  !! routine h_psi_ .
  !
  USE kinds,              ONLY: DP
  USE noncollin_module,   ONLY: npol
  USE xc_lib,             ONLY: exx_is_active
  USE mp_bands,           ONLY: use_bgrp_in_hpsi, inter_bgrp_comm, nbgrp
  USE mp,                 ONLY: mp_allgather, mp_size, &
                                mp_type_create_column_section, mp_type_free
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda
  !! leading dimension of arrays psi, spsi, hpsi
  INTEGER, INTENT(IN) :: n
  !! true dimension of psi, spsi, hpsi
  INTEGER, INTENT(IN) :: m
  !! number of states psi
  COMPLEX(DP), INTENT(IN) :: psi(lda*npol,m)
  !! the wavefunction
  COMPLEX(DP), INTENT(OUT) :: hpsi(lda*npol,m)
  !! Hamiltonian dot psi
  !
  ! ... local variables
  !
  INTEGER :: m_start, m_end
  INTEGER :: column_type
  INTEGER, ALLOCATABLE :: recv_counts(:), displs(:)
  !
  !
  CALL start_clock( 'h_psi_bgrp' ); !write (*,*) 'start h_psi_bgrp'; FLUSH(6)
  !
     CALL h_psi_gipaw_( lda, n, m, psi, hpsi )
  !
  CALL stop_clock( 'h_psi_bgrp' )
  !
  !
  RETURN
  !
END SUBROUTINE h_psi_gipaw
!
!----------------------------------------------------------------------------
SUBROUTINE h_psi_gipaw_ ( lda, n, m, psi, hpsi )
  !----------------------------------------------------------------------------
  !
  ! ... This routine computes the product of the Hamiltonian
  ! ... matrix with m wavefunctions contained in psi
  !
  ! ... input:
  ! ...    lda   leading dimension of arrays psi, spsi, hpsi
  ! ...    n     true dimension of psi, spsi, hpsi
  ! ...    m     number of states psi
  ! ...    psi
  !
  ! ... output:
  ! ...    hpsi  H*psi
  !
  USE kinds,    ONLY : DP
  USE uspp,     ONLY : vkb, nkb
  USE klist,    ONLY : igk_k
  USE fft_base, ONLY : dffts, dfftp, dfftb
!  USE ldaU,     ONLY : lda_plus_u
  USE lsda_mod, ONLY : current_spin
  USE scf,      ONLY : vrs
  USE orbital_magnetization,      ONLY : dvrs
  USE gvect,    ONLY : gstart, ngm
  USE gvecs,    ONLY : ngms
  USE gipaw_module,  ONLY : lambda_so
  USE nmr_mod
#ifdef EXX
  USE exx,      ONLY : vexx
  USE funct,    ONLY : exx_is_active
  USE realus,   ONLY: real_space
  USE wvfct,    ONLY : g2kin, nbndx, nbnd
#endif

  !
  IMPLICIT NONE
  !
  ! ... input/output arguments
  !
  INTEGER     :: lda, n, m
  COMPLEX(DP) :: psi(lda,m) 
  COMPLEX(DP) :: hpsi(lda,m)
  integer     :: nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s,nrxxs

  ! various initializations:
  nrxxs = dffts%nnr
  nr1s=dffts%nr1
  nr2s=dffts%nr2
  nr3s=dffts%nr3
  nrx1s=dffts%nr1x
  nrx2s=dffts%nr2x
  nrx3s=dffts%nr3x
  !
  CALL start_clock( 'h_psi' )
  !
  if ( any(m_0 /= 0.d0))  then
          call h_psi_k_nmr( )
  else  
  CALL h_psi_k( )
  endif
  !
  CALL stop_clock( 'h_psi' )
  !
  RETURN
  !
  CONTAINS
!-----------------------------------------------------------------------
     SUBROUTINE h_psi_k( )
!-----------------------------------------------------------------------
       !
       ! ... k-points version
       !
       USE wavefunctions,        ONLY : psic
       USE gipaw_module,         ONLY : lambda_so
       USE wvfct,                ONLY : current_k
       USE klist,                ONLY : nelup, neldw
       USE fft_interfaces,       ONLY : invfft, fwfft
       USE becmod,                  ONLY: bec_type, becp, calbec
       USE orbital_magnetization,  ONLY : dvrs
       USE wvfct,    ONLY : g2kin, nbndx, nbnd
       USE nmr_mod

       !
       IMPLICIT NONE
       !
       INTEGER :: ibnd, j, ik
       complex(dp), allocatable, dimension(:,:) :: p_psic
       !
       ik = current_k
       CALL start_clock( 'init' )
       !
       ! ... Here we apply the kinetic energy (k+G)^2 psi
       !
       DO ibnd = 1, m
          !
          hpsi(1:n,ibnd) = g2kin(1:n) * psi(1:n,ibnd)
  !        print*, 'hpsi', hpsi(n,ibnd), ibnd
          !
       END DO
 !      stop
       !
  ! new part from h_psi PW
       CALL vloc_psi_k_gipaw( lda, n, m, psi, vrs(1,current_spin), hpsi )

  ! ... Here the product with the non local potential V_NL psi
  ! ... (not in the real-space case: it is done together with V_loc)
  !
  IF ( nkb > 0) THEN
     !
     CALL start_clock( 'h_psi:calbec' )
     CALL calbec( n, vkb, psi, becp, m )
     CALL stop_clock( 'h_psi:calbec' )
     CALL add_vuspsi( lda, n, m, hpsi )
     !
  ENDIF
  ! ... Add the paramagnetic term to the Hamiltonian
  if (any(lambda_so /= 0.d0)) then
          call add_so_Fnl (lda, n, m, 1.d0, psi, hpsi)
  endif
  if (any(m_0 /= 0.d0)) then
          call add_nmr_Fnl (lda, n, m, psi, hpsi)
  endif


       
     END SUBROUTINE h_psi_k     
     !
!-----------------------------------------------------------------------
SUBROUTINE vloc_psi_k_gipaw( lda, n, m, psi, v, hpsi )
  !-----------------------------------------------------------------------
  !! Calculation of Vloc*psi using dual-space technique - k-points:
  !
  !! * fft to real space;
  !! * product with the potential v on the smooth grid;
  !! * back to reciprocal space;
  !! * addition to the hpsi.
  !
  USE parallel_include
  USE kinds,                  ONLY : DP
  USE wvfct,                  ONLY : current_k
  USE klist,                  ONLY : igk_k
  USE mp_bands,               ONLY : me_bgrp
  USE fft_base,               ONLY : dffts
  USE fft_wave
  USE fft_helper_subroutines, ONLY : fftx_ntgrp, tg_get_nnr, tg_get_group_nr3
  USE wavefunctions,          ONLY : psic
  USE gipaw_module,           ONLY : lambda_so
  USE nmr_mod
  USE orbital_magnetization  
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda
  !! leading dimension of arrays psi, hpsi
  INTEGER, INTENT(IN) :: n
  !! true dimension of psi, hpsi
  INTEGER, INTENT(IN) :: m
  !! number of states psi
  COMPLEX(DP), INTENT(IN) :: psi(lda,m)
  !! the wavefunction
  COMPLEX(DP), INTENT(INOUT) :: hpsi(lda,m)
  !! Hamiltonian dot psi
  REAL(DP), INTENT(IN) :: v(dffts%nnr)
  !! the total pot. in real space (smooth grid) for current spin

  INTEGER :: ibnd, j, incr, ctv
  INTEGER :: i, iin, right_nnr, right_nr3, right_inc
  COMPLEX(DP), ALLOCATABLE :: vpsi(:,:)
  ! ... chunking parameters
  INTEGER, PARAMETER :: blocksize = 256
  INTEGER :: numblock
  ! ... Task Groups
  LOGICAL :: use_tg
  REAL(DP), ALLOCATABLE :: tg_v(:)
  COMPLEX(DP), ALLOCATABLE :: tg_psic(:), tg_vpsi(:,:), p_psic(:,:)
  INTEGER :: idx, brange, dffts_nnr
  !
  CALL start_clock( 'vloc_psi' )
  use_tg = dffts%has_task_groups
  !
  IF( use_tg ) THEN
     CALL start_clock( 'vloc_psi:tg_gather' )
     dffts_nnr = dffts%nnr_tg
     incr = fftx_ntgrp(dffts)
     ALLOCATE( tg_v(dffts_nnr) )
     ALLOCATE( tg_psic(dffts_nnr), tg_vpsi(lda,incr) )
     CALL tg_gather( dffts, v, tg_v )
     CALL stop_clock( 'vloc_psi:tg_gather' )
  ELSE
     dffts_nnr = dffts%nnr
     ALLOCATE( vpsi(lda,1) )
  ENDIF
  !
  IF ( use_tg ) THEN
     !
     CALL tg_get_nnr( dffts, right_nnr )
     !
     ! ... compute the number of chuncks
     numblock = (n+blocksize-1)/blocksize
     !
     DO ibnd = 1, m, fftx_ntgrp(dffts)
             !
        CALL tgwave_g2r( psi(:,ibnd:m), tg_psic, dffts, n, igk_k(:,current_k) )
        !
!        write (6,*) 'wfc R '
!        write (6,99) (tg_psic(i), i=1,400)
        !
        CALL tg_get_group_nr3( dffts, right_nr3 )
        !
        DO j = 1, dffts%nr1x*dffts%nr2x*right_nr3
           tg_psic(j) = tg_psic(j) * tg_v(j)
        ENDDO
        !
!        write (6,*) 'v psi R '
!        write (6,99) (tg_psic(i), i=1,400)
        !
        brange = m-ibnd+1
        !
        CALL tgwave_r2g( tg_psic, tg_vpsi(:,1:brange), dffts, n, igk_k(:,current_k) )
        !
        DO idx = 0, MIN(fftx_ntgrp(dffts)-1, m-ibnd)
           DO j = 1, numblock
              DO iin = (j-1)*blocksize+1, MIN(j*blocksize,n)
                 hpsi(iin,ibnd+idx) = hpsi(iin,ibnd+idx) + tg_vpsi(iin,idx+1)
              ENDDO
           ENDDO
        ENDDO
        !
     ENDDO
     !
  ELSE
     !
       if (current_k <= 0) call errore('h_psi_k', 'current_k??', 1)
       if (current_spin == 0) call errore('h_psi_k', 'current_spin??', 1)
       if (any(lambda_so /= 0.d0) .or. any(m_0 /= 0.d0)) &
         allocate( p_psic(1:dffts%nnr,3) )
     !
     DO ibnd = 1, m
        !
!        psic(1:dffts%nnr) = ( 0.D0, 0.D0 )
!        psic(dffts%nl(igk_k(1:n,current_k))) = psi(1:n,ibnd)
        CALL wave_g2r( psi(1:n,ibnd:ibnd), psic, dffts, igk=igk_k(:,current_k) )

        !
        !        write (6,*) 'wfc R '
!        write (6,99) (psic(i), i=1,400)
        !
        if (any(m_0 /= 0.d0)) then
!                print*, 'chiama nmr_valance'
            call add_nmr_valence(current_k, n, psi(1:n,ibnd), p_psic)
        else
          DO j = 1, dffts_nnr
             psic(j) = psic(j) * v(j)
          ENDDO
        endif
        !
!        write (6,*) 'v psi R '
!        write(6,*) (psic(i), i=1,400i)
        !
        ! Add the SO term to the Hamiltonian
        if (any(lambda_so /= 0.d0)) then
            call add_so_valence(current_k, n, 1.d0, psi(1:n,ibnd), p_psic)
        endif
        !
        CALL wave_r2g( psic(1:dffts_nnr), vpsi(1:n,:), dffts, igk=igk_k(:,current_k) )
        !
        ! Add the SO term to the Hamiltonian
     !   if (any(lambda_so /= 0.d0)) then
     !       call add_so_valence(current_k, n, 1.d0, psi(1:n,ibnd), p_psic)
     !   endif
        !
        !
        DO i = 1, n
           hpsi(i,ibnd) = hpsi(i,ibnd) + vpsi(i,1)
        ENDDO
        !
!        write (6,*) 'v psi G ', ibnd
!        write (6,99) (psic(i), i=1,400)
     ENDDO
!        do ibnd =1, m
!          print*, 'hpsi', hpsi(n,ibnd), ibnd
!        enddo
!        stop
  ENDIF
  !
  if (any(lambda_so /= 0.d0) .or. any(m_0 /= 0.d0))  then
          deallocate( p_psic )
!  if (any(m_0 /= 0.d0)) then
!          print*, 'chiama nmr_Fnl'
!          call add_nmr_Fnl (lda, n, m, psi, hpsi)
!         do ibnd =1, m
!          print*, 'hpsi', hpsi(n,ibnd), ibnd
!        enddo
!  endif
!  stop
  endif
  !
  IF ( use_tg ) THEN
     DEALLOCATE( tg_psic, tg_vpsi )
     DEALLOCATE( tg_v )
  ELSE
     DEALLOCATE( vpsi )
  ENDIF
  !
  CALL stop_clock( 'vloc_psi' )
  !
99 format ( 20 ('(',2f12.9,')') )
  !
  RETURN
  !
END SUBROUTINE vloc_psi_k_gipaw

!---------------------------------------------------------------
 ! add the SO valence term
!---------------------------------------------------------------
  SUBROUTINE add_so_valence(ik, n, alpha, psi, p_psic)
  USE wavefunctions,          ONLY : psic
  USE wvfct,                ONLY : current_k
  USE gvect,                ONLY :  ngm, g
  USE cell_base,            ONLY : tpiba
  USE lsda_mod,             ONLY : current_spin
  USE klist,                ONLY : xk
  USE klist,                  ONLY : igk_k
  USE fft_base,             ONLY : dffts, dfftp, dfftb
  USE fft_interfaces,       ONLY : invfft, fwfft
  USE orbital_magnetization, ONLY : a2gp8
  IMPLICIT NONE
  integer :: ik, n, ipol, ig, i, nrxxs_l
  complex(dp) :: p_psic(dffts%nnr,3)
  complex(dp) :: psi(n)
  real(dp) :: gk, sigma, alpha
  ! index for the cross product
  integer :: ind(2,3), ii, jj
  ind(:,1) = (/ 2, 3 /)
  ind(:,2) = (/ 3, 1 /)
  ind(:,3) = (/ 1, 2 /)
  call start_clock( 'add_so_valence' )

  !!RETURN
  ! compute (p+k)|psi> in real space
  p_psic(1:dffts%nnr,1:3) = (0.d0,0.d0)
  ik = current_k
  nrxxs_l = dffts%nnr
  do ipol = 1, 3
    do ig = 1, n
      gk = xk(ipol,ik) + g(ipol,igk_k(ig,ik))
      p_psic(dffts%nl(igk_k(ig,ik)),ipol) = -gk * tpiba * psi(ig)
    enddo
     CALL invfft ('Wave' , p_psic(:,ipol), dffts)
  enddo
   
  ! compute lambda_so ( sigma \cdot dvrs \times p_psi)
  sigma = 1.d0; if (current_spin == 2) sigma = -1.d0
  do ipol = 1, 3
    if (lambda_so(ipol) == 0.d0) cycle
    ii = ind(1,ipol)
    jj = ind(2,ipol)
    do i = 1, nrxxs_l
        psic(i) = psic(i) + alpha * lambda_so(ipol) * a2gp8 * sigma * ( &
        dvrs(i,current_spin,ii)*p_psic(i,jj) - &
        dvrs(i,current_spin,jj)*p_psic(i,ii) )
    enddo
  enddo

  CALL stop_clock( 'add_so_valence' )
  END SUBROUTINE add_so_valence

!---------------------------------------------------------------
! add the F_R^{NL} term of the SO reconstruction
!---------------------------------------------------------------
  SUBROUTINE add_so_Fnl(lda, n, m, alpha, psi, hpsi)
  USE kinds,      ONLY : DP
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp
  USE lsda_mod,   ONLY : current_spin
  USE paw_gipaw,       ONLY : paw_vkb, paw_nkb, paw_recon, paw_lmaxkb, paw_becp
  USE wvfct,           ONLY : current_k
  USE klist,           ONLY : xk, igk_k
  USE fft_base, ONLY : dffts, dfftp, dfftb
  USE wvfct,    ONLY : g2kin, nbndx, nbnd
  USE becmod,          ONLY: bec_type, becp, calbec
  USE orbital_magnetization, ONLY : a2gp8
  USE gipaw_module,  ONLY : radial_integral_paramagnetic_so
  USE gipaw_module,  ONLY : lx, ly, lz
  USE gipaw_module,  ONLY : lambda_so
  IMPLICIT NONE
  INTEGER :: lda, n, m
  INTEGER :: ibnd, ijkb0, nt, na, ikb, jkb, ih, jh
  INTEGER :: l1, m1, lm1, l2, m2, lm2, nbs1, nbs2
  integer :: nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s,nrxxs
  COMPLEX(DP), ALLOCATABLE :: ps (:,:), paw_becp4 (:,:)
  COMPLEX(DP) :: psi(lda,n), hpsi(lda,n)
  real(dp) :: sigma, alpha

  CALL start_clock( 'add_so_Fnl' )

  nrxxs = dffts%nnr
  nr1s=dffts%nr1
  nr2s=dffts%nr2
  nr3s=dffts%nr3
  nrx1s=dffts%nr1x
  nrx2s=dffts%nr2x
  nrx3s=dffts%nr3x

  if (m > nbndx) call errore('add_so_Fnl', 'm > nbndx ???', m)
!  !!if (m > nbnd) call errore('add_so_Fnl', 'm > nbnd ???', m)
  ALLOCATE (ps(paw_nkb,m))
  allocate (paw_becp4(paw_nkb,m))
  ps(:,:) = (0.D0, 0.D0)
  paw_becp4 = (0d0, 0d0)

  call init_gipaw_2(n, igk_k(1,current_k), xk(1,current_k), paw_vkb)
  CALL calbec( n, paw_vkb, psi, paw_becp4, m )
  sigma = 1.d0; if (current_spin == 2) sigma = -1.d0

 ijkb0 = 0
  do nt = 1, ntyp
    do na = 1, nat
      if (ityp(na) .eq. nt) then
        do ih = 1, paw_recon(nt)%paw_nh
          ikb = ijkb0 + ih
          nbs1 = paw_recon(nt)%paw_indv(ih)
          l1 = paw_recon(nt)%paw_nhtol(ih)
          m1 = paw_recon(nt)%paw_nhtom(ih)
          lm1 = m1 + l1**2
         
          do jh = 1, paw_recon(nt)%paw_nh
            jkb = ijkb0 + jh
            nbs2 = paw_recon(nt)%paw_indv(jh)
            l2 = paw_recon(nt)%paw_nhtol(jh)
            m2 = paw_recon(nt)%paw_nhtom(jh)
            lm2 = m2 + l2**2

            if (l1 /= l2) cycle
            if (l1 == 0) cycle

            do ibnd = 1, m
              ps(ikb,ibnd) = ps(ikb,ibnd) + lambda_so(1) * a2gp8 * &
                sigma * lx(lm1,lm2) * &
                radial_integral_paramagnetic_so(nbs1,nbs2,nt) * &
                paw_becp4(jkb,ibnd)
            
             ps(ikb,ibnd) = ps(ikb,ibnd) + lambda_so(2) * a2gp8 * &
                sigma * ly(lm1,lm2) * &
                radial_integral_paramagnetic_so(nbs1,nbs2,nt) * &
                paw_becp4(jkb,ibnd)

              ps(ikb,ibnd) = ps(ikb,ibnd) + lambda_so(3) * a2gp8 * &
                sigma * lz(lm1,lm2) * &
                radial_integral_paramagnetic_so(nbs1,nbs2,nt) * &
                paw_becp4(jkb,ibnd)
            enddo   ! ibnd
          enddo   ! jh
        enddo   ! ih
        ijkb0 = ijkb0 + paw_recon(nt)%paw_nh
      endif    ! ityp(na) .eq. nt
    enddo   ! na
  enddo   ! nt
  ps = ps * alpha

  CALL ZGEMM( 'N', 'N', n, m, paw_nkb, ( 0.D0, 1.0D0 ) , paw_vkb, &
              lda, ps, paw_nkb, ( 1.D0, 0.D0 ) , hpsi, lda )

  deallocate (ps)
  deallocate (paw_becp4)

  CALL stop_clock( 'add_so_Fnl' )
  END SUBROUTINE add_so_Fnl

  !-----------------------------------------------------------------------
     SUBROUTINE h_psi_k_nmr()
!-----------------------------------------------------------------------
       !
       ! ... k-points version
       !
       USE wavefunctions,        ONLY : psic
       USE gvect,                ONLY : ngm, g
       USE gipaw_module,         ONLY : lambda_so
       USE wvfct,                ONLY : current_k, npw
       USE klist,                ONLY : nelup, neldw
       USE klist,                ONLY : xk, nks, igk_k, ngk
       USE fft_interfaces,       ONLY : invfft, fwfft
       USE becmod,                  ONLY: bec_type, becp, calbec
       USE orbital_magnetization,  ONLY : dvrs
       USE wvfct,    ONLY : g2kin, nbndx, nbnd
       USE gvecw,                ONLY : gcutw
       USE cell_base,            ONLY : tpiba2
       USE nmr_mod

       !
       IMPLICIT NONE
       !
       INTEGER :: ibnd, j, ik
       complex(dp), allocatable, dimension(:,:) :: p_psic
       !
       allocate( p_psic(1:nrxxs,3) )
       ik = current_k
       ! ... Here we apply the kinetic energy (k+G)^2 psi
       !
       DO ibnd = 1, m
          !
          hpsi(1:n,ibnd) = g2kin(1:n) * psi(1:n,ibnd)
          !
       END DO
       ! ... the local potential V_Loc psi. First the psi in real sp
       DO ibnd = 1, m
          !
          !
          psic(1:nrxxs) = ( 0.D0, 0.D0 )
          !
          psic(dffts%nl(igk_k(1:n,current_k))) = psi(1:n,ibnd)
          !
        CALL invfft ('Wave' , psic, dffts)
        if (any(m_0 /= 0.d0)) then
                call add_nmr_valence(current_k, n, psi(1:n,ibnd), p_psic)
        endif
          !
          ! ... back to reciprocal space
          !
          CALL fwfft ('Wave' , psic, dffts)
          !
          ! ... addition to the total product
          !
          hpsi(1:n,ibnd) = hpsi(1:n,ibnd) + psic(dffts%nl(igk_k(1:n,current_k)))
          !
        END DO
        IF ( nkb > 0) THEN
        !
          CALL start_clock( 'h_psi:calbec' )
          CALL calbec( n, vkb, psi, becp, m )
          CALL stop_clock( 'h_psi:calbec' )
          CALL add_vuspsi( lda, n, m, hpsi )
        !
        ENDIF
        ! ... Add the paramagnetic term to the Hamiltonian
        if (any(m_0 /= 0.d0)) then
          call add_nmr_Fnl (lda, n, m, psi, hpsi)
        endif
        deallocate( p_psic )

       END SUBROUTINE h_psi_k_nmr



END SUBROUTINE h_psi_gipaw_
