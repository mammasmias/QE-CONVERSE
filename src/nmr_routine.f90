!!
!! BEGIN_NMR: Changed by Timo Thonhauser for NMR shielding tensors
!!
!! Copyright (C) Timo Thonhauser 2007
!!
!! The nmr part of the code first gets initialized in "init_nmr()".
!! All necessary variables are set and the vector potential A_vec
!! is calculated in reciprocal space. The routine h_psi.f90 is then
!! modified to implement the vector potential into the Hamiltonian.
!! The subroutine init_us_2.f90 has been replaced by init_us_3.f90
!! in order to correct the non-local projectors of the pseudopotential.
!! Finally, in "calc_moment()", the induced current is computed and
!! from that the shielding is calculated. A correction with the first
!! moment of div(I_err) is also included.
!!
!!-----------------------------------------------------------------------
!
!  !-----------------------------------------------------------------------
  SUBROUTINE init_nmr()
!  !-----------------------------------------------------------------------
!  !
!  ! This subroutine initializes NMR variables ...
!  ! Rydberg atomic units are used: m=1/2, e^2=2, hbar=1, E_0=1
!  !
  USE  constants    , ONLY : pi, tpi, sqrt2
  USE  ions_base    , ONLY : tau, nat
  USE  io_global    , ONLY : stdout
  USE  cell_base    , ONLY : tpiba, tpiba2, omega, alat
  USE fft_base,         ONLY : dffts, dfftp
  USE fft_interfaces,   ONLY : fwfft, invfft
  USE  gvect        , ONLY : ngm, g, gg
  USE  gvecs        , ONLY : ngms
  USE kinds,     ONLY : dp
  USE  control_flags,  ONLY : gamma_only
  USE nmr_mod
  USE gipaw_module
!  !
  IMPLICIT NONE
!  !
  INTEGER                 :: ii, jj
  REAL(DP)                :: r_diff(3), r_length
  COMPLEX(DP)             :: prefac
  integer :: nrxx, nrxxs
  integer :: nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s
!  !
    nrxxs = dffts%nnr
    nrxx= dfftp%nnr 
    nr1s = dffts%nr1;   nr2s = dffts%nr2;   nr3s = dffts%nr3
    nrx1s = dffts%nr1x; nrx2s = dffts%nr2x; nrx3s = dffts%nr3x
!  !
!  ! Simple tests before running NMR ...
!  !
  write(stdout,'(4A)') cr, cr,'    Initializing NMR ...', cr
!  !
  if(ngm /= ngms .OR. nrxx/=nrxxs) CALL errore('init_nmr', &
  &   'Fine grid and smooth grid are NOT identical', 1)
  if( gamma_only ) CALL errore('init_nmr','gamma_only = TRUE', 2)
!  !
!  !
!  ! Allocating arrays for the vector potential and grad psi.
!  ! Everything should be cleaned up later in clean_pw.f90 ...
!  !
  write(stdout,'(A,i8)') '    Allocate A_vec(nrxxs,3)    with dffts%nnr =', nrxxs
  write(stdout,'(A,i8)') '    Allocate grad_psi(nrxxs,3) with dffts%nnr =', nrxxs
  write(stdout,'(A,i8)') '    Allocate A_0(3,nat)        with nat   =', nat
!
  allocate( A_vec(nrxxs,3), grad_psi(nrxxs,3), A_0(3,nat))
  A_vec(1:nrxxs,1:3)    = (0.0D0, 0.0D0)
  grad_psi(1:nrxxs,1:3) = (0.0D0, 0.0D0)
  A_0(1:3,1:nat)        = (0.0D0, 0.0D0)
!  !
!  !
!  ! Defining basic variables ...
!  !
  write(stdout,'(2A,i4)')  cr, '    Magnetic dipole at atom ', m_0_atom
  write(stdout,'(A,3f6.2,A)')  '    m_0 = (', m_0(1), m_0(2), m_0(3), ')  mu_b'
!  !
#if 0
  n_nlpp = 0
  do ii  = 1, nlpp_max
     if( r_nlpp(ii) > 0.0D0 ) then
        n_nlpp = n_nlpp + 1
        write(stdout,'(A,i4,A,e10.4,A)' )'    Non-local projector at atom ', &
          ii, ' with radius ', r_nlpp(ii), ' Bohr'
     end if
  end do
#endif
!  !
  if (m_0(1) /= 0.0D0) m_0_comp = 1
  if (m_0(2) /= 0.0D0) m_0_comp = 2
  if (m_0(3) /= 0.0D0) m_0_comp = 3
  m_0(1:3) = m_0(1:3) * sqrt2         ! Bohr magneton in Rydberg units
!  !
!  !
!  ! Initializing the vector potential in G space ...
!  !
  DO ii = 1, ngms
!     !
     IF (gg(ii) == 0.0D0) CYCLE
     prefac = -2.0D0 * IM * pi * fine_struct/(gg(ii)*tpiba*omega) ! Rydberg units,
     prefac = prefac * exp(-gg(ii)*tpiba2*sigma**2/4.0D0)         ! Gaussian
     prefac = prefac * exp(-IM * tpi *( &                         ! convention.
              g(1,ii)*tau(1,m_0_atom) + &                         ! Convolution
              g(2,ii)*tau(2,m_0_atom) + &                         ! with delta
              g(3,ii)*tau(3,m_0_atom)  ))                         ! function
     A_vec(dffts%nl(ii),1) = prefac * ( m_0(2)*g(3,ii)-m_0(3)*g(2,ii) )
     A_vec(dffts%nl(ii),2) = prefac * (-m_0(1)*g(3,ii)+m_0(3)*g(1,ii) )
     A_vec(dffts%nl(ii),3) = prefac * ( m_0(1)*g(2,ii)-m_0(2)*g(1,ii) )
     !
  END DO
!  !
!  !
!  ! Calculate and write out the vector potential at the center of 
!  ! all atoms for the phase twist of the projectors ...
!  !
  write(stdout, '(A)') cr
!  !
  DO jj = 1, nat
!     !
     r_diff(:) = alat * ( tau(:,jj) - tau(:,m_0_atom) )
     r_length  = sqrt( r_diff(1)**2 + r_diff(2)**2 + r_diff(3)**2 )
     A_0(1,jj) =  m_0(2)*r_diff(3) - m_0(3)*r_diff(2)
     A_0(2,jj) = -m_0(1)*r_diff(3) + m_0(3)*r_diff(1)
     A_0(3,jj) =  m_0(1)*r_diff(2) - m_0(2)*r_diff(1)
!     !
     if (r_length == 0.0D0) then
        A_0(:,jj) = (0.0D0, 0.0D0)
     else
        A_0(:,jj) = fine_struct/2.0D0 * A_0(:,jj) / r_length**3
     end if
!     !
     write(stdout, '(A,i4,A,3e14.6)') '    A_0(:,', jj, '): ', REAL(A_0(1:3,jj))
!     !
  END DO
!  !
  write(stdout, '(A)') cr
!  !
!  !
!  ! FFT A_vec into real space ...
!  !
  call invfft ('Rho',A_vec(:,1),dffts)
  call invfft ('Rho',A_vec(:,2),dffts)
  call invfft ('Rho',A_vec(:,3),dffts)
!  !
!  !
!  ! Make A_vec purely real ...
!  !
!  print*,'nrxxs', nrxxs
  A_vec(1:nrxxs,1:3) = CMPLX( DBLE(A_vec(1:nrxxs,1:3)), 0.0D0, kind=DP )
!  print*,'A', A_vec(1000,1:3)
!  stop
!  !
!  !
  END SUBROUTINE init_nmr
!  !-----------------------------------------------------------------------
!
!
!  !---------------------------------------------------------------
!  ! add the NMR valence term
!  !---------------------------------------------------------------
SUBROUTINE add_nmr_valence(ik, n, psi, p_psic)
USE wavefunctions,            ONLY : psic
USE wvfct,                ONLY : current_k
!  !USE gvect,                ONLY : nlm, ngm, g
  USE  gvect,                 ONLY : ngm, g, gg
  USE cell_base,            ONLY : tpiba
  USE lsda_mod,             ONLY : current_spin
  USE klist,                ONLY : xk, igk_k
  USE scf,                  ONLY : vrs
  USE fft_base,         ONLY : dffts, dfftp
  USE fft_interfaces,   ONLY : fwfft, invfft
  USE kinds,     ONLY : dp
  USE nmr_mod
  USE gipaw_module
!  !USE wvfct,                ONLY : igk
!  !USE gsmooth,              ONLY : nls, nr1s, nr2s, nr3s, &
!  !                                 nrx1s, nrx2s, nrx3s, nrxxs
  IMPLICIT NONE
  integer :: ik, n, ipol, ig, i
!  complex(dp) :: p_psic(nrxxs,3)
  complex(dp) :: p_psic(dffts%nnr,3)
  complex(dp) :: psi(n)
  real(dp) :: gk, Asquare
  integer :: nrxx, nrxxs
  integer :: nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s
    nrxxs = dffts%nnr
    nrxx= dfftp%nnr 
    nr1s = dffts%nr1;   nr2s = dffts%nr2;   nr3s = dffts%nr3
    nrx1s = dffts%nr1x; nrx2s = dffts%nr2x; nrx3s = dffts%nr3x
   !! psic(1:nrxxs) = vrs(1:nrxxs,current_spin)*psic(1:nrxxs)
!  !!RETURN
!  ! compute (p+k)|psi> in real space
 p_psic(1:dffts%nnr,1:3) = (0.d0,0.d0)
 ik = current_k
 do ipol = 1, 3
   do ig = 1, n
     gk = xk(ipol,ik) + g(ipol,igk_k(ig,ik))
     p_psic(dffts%nl(igk_k(ig,ik)),ipol) = gk * tpiba * psi(ig) * IM
   enddo
   CALL invfft ('Wave' , p_psic(:,ipol), dffts)
 enddo
!
!  ! BEGIN_NMR: Changed by Timo Thonhauser for NMR shielding tensors
!  !
!  ! Rydberg atomic units are used (m=1/2, e^2=2, hbar=1, E_0=1).
!  !
!print*, 'psic', psic(1), psic(100), psic(1000)
!stop
 do i = 1, nrxxs
   Asquare = A_vec(i,1)**2 + A_vec(i,2)**2 + A_vec(i,3)**2
   Asquare = Asquare * fine_struct**2/2.0d0                   ! A^2 term
   !!Asquare = 0.d0
   psic(i) = ( Asquare + vrs(i,current_spin) ) * psic(i)
!    ! negative sign according to Giovanni
   psic(i) = psic(i) - IM*sqrt(2.d0)*fine_struct* ( &
                    A_vec(i,1) * p_psic(i,1) +               & ! 2 p.(A psi)=  
                    A_vec(i,2) * p_psic(i,2) +               & ! 2 A.(p psi)
                    A_vec(i,3) * p_psic(i,3) )
 enddo
 END SUBROUTINE add_nmr_valence
!
!
!  !---------------------------------------------------------------
!  ! add the F_R^{NL} term of the NMR reconstruction
!  !---------------------------------------------------------------
 SUBROUTINE add_nmr_Fnl(lda, n, m, psi, hpsi)
 USE kinds,      ONLY : DP
 USE ions_base,  ONLY : nat, ntyp => nsp, ityp
 USE lsda_mod,   ONLY : current_spin
 USE paw_gipaw,       ONLY : paw_vkb, paw_nkb, paw_recon, paw_lmaxkb, paw_becp
!!  USE paw,        ONLY : paw_vkb, paw_nkb, paw_recon, paw_lmaxkb
 USE wvfct,           ONLY : current_k
 USE klist,           ONLY : xk, igk_k
 USE fft_base, ONLY : dffts, dfftp, dfftb
 USE wvfct,    ONLY : g2kin, nbndx, nbnd
 USE becmod,          ONLY: bec_type, becp, calbec
!!  USE wvfct,    ONLY : igk, g2kin, nbndx, nbnd
!!  USE gsmooth,  ONLY : nls, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, nrxxs
!!  USE paw,      ONLY : paw_becp
  USE gipaw_module
!!  USE gipaw_module, ONLY : radial_integral_paramagnetic, &
!!                              radial_integral_diamagnetic, &
!!                              lx, ly, lz
 USE nmr_mod
  IMPLICIT NONE
  INTEGER :: lda, n, m, i, j
  INTEGER :: ibnd, ijkb0, nt, na, ikb, jkb, ih, jh
  INTEGER :: l1, m1, lm1, l2, m2, lm2, nbs1, nbs2
  integer :: nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s,nrxxs
  COMPLEX(DP), ALLOCATABLE :: ps (:,:), paw_becp5 (:,:)
  COMPLEX(DP) :: psi(lda,n), hpsi(lda,n)
    
  nrxxs = dffts%nnr
  nr1s=dffts%nr1
  nr2s=dffts%nr2
  nr3s=dffts%nr3
  nrx1s=dffts%nr1x
  nrx2s=dffts%nr2x
  nrx3s=dffts%nr3x
  if (paw_nkb == 0) return
  if (m > nbndx) call errore('add_nmr_Fnl', 'm > nbndx ???', -m)
  if (m > nbnd) call errore('add_nmr_Fnl', 'm > nbnd ???', -m)
  ALLOCATE (ps(paw_nkb,m))
  allocate (paw_becp5(paw_nkb,m))
  ps(:,:) = (0.D0, 0.D0)
  paw_becp5 = (0d0, 0d0)
!
!!  call init_paw_2(n, igk, xk(1,current_k), paw_vkb)
!!  call ccalbec(paw_nkb, lda, n, m, paw_becp, paw_vkb, psi)
!debug

  call init_gipaw_2(n, igk_k(1,current_k), xk(1,current_k), paw_vkb)
  CALL calbec( n, paw_vkb, psi, paw_becp5, m )
!do jkb=1, paw_nkb
!  do ibnd = 1, m
!     print*, 'paw_becp5(jkb,ibnd)', paw_becp5(jkb,ibnd), jkb, ibnd
!  enddo
!enddo
!stop
!
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
!print*, 'paw_recon(nt)%paw_nh', paw_recon(nt)%paw_nh
!print*, 'nbs1, l1, m1, lm1', nbs1, l1, m1, lm1
!stop
          do jh = 1, paw_recon(nt)%paw_nh
            jkb = ijkb0 + jh
            nbs2 = paw_recon(nt)%paw_indv(jh)
            l2 = paw_recon(nt)%paw_nhtol(jh)
            m2 = paw_recon(nt)%paw_nhtom(jh)
            lm2 = m2 + l2**2
!print*, 'paw_recon(nt)%paw_nh', paw_recon(nt)%paw_nh
!print*, 'nbs1, l1, m1, lm1', nbs2, l2, m2, lm2
!stop
            if (l1 /= l2) cycle
            if (l1 == 0) cycle
            if (na /= m_0_atom) cycle
            do ibnd = 1, m
              ps(ikb,ibnd) = ps(ikb,ibnd) - fine_struct**2 / sqrt(2.d0) * &
                m_0(1) * lx(lm1,lm2) * &
                radial_integral_paramagnetic(nbs1,nbs2,nt) * &
                paw_becp5(jkb,ibnd)
!             print*, 'radial_integral_paramagnetic', radial_integral_paramagnetic(nbs1,nbs2,nt)
!             print*, 'paw_becp5(jkb,ibnd)', paw_becp5(jkb,ibnd)
!              print*, 'ps', ps(ikb,ibnd), ikb,ibnd

              ps(ikb,ibnd) = ps(ikb,ibnd) - fine_struct**2 / sqrt(2.d0) * &
                m_0(2) * ly(lm1,lm2) * &
                radial_integral_paramagnetic(nbs1,nbs2,nt) * &
                paw_becp5(jkb,ibnd)
!             print*, 'radial_integral_paramagnetic', radial_integral_paramagnetic(nbs1,nbs2,nt)

              ps(ikb,ibnd) = ps(ikb,ibnd) - fine_struct**2 / sqrt(2.d0) * &
                m_0(3) * lz(lm1,lm2) * &
                radial_integral_paramagnetic(nbs1,nbs2,nt) * &
                paw_becp5(jkb,ibnd)
!             print*, 'radial_integral_paramagnetic', radial_integral_paramagnetic(nbs1,nbs2,nt)
            enddo   ! ibnd
          enddo   ! jh
        enddo   ! ih
        ijkb0 = ijkb0 + paw_recon(nt)%paw_nh
      endif    ! ityp(na) .eq. nt
    enddo   ! na
  enddo   ! nt
! do ikb=1, paw_nkb
!   do ibnd=1, m
!    print*, 'ps',ps (ikb, paw_nkb)
!   enddo
!  enddo


  CALL ZGEMM( 'N', 'N', n, m, paw_nkb, ( 0.D0, 1.0D0 ) , paw_vkb, &
              lda, ps, paw_nkb, ( 1.D0, 0.D0 ) , hpsi, lda )
!     print*,'hpsi', hpsi(1,1)

  deallocate (ps)
  END SUBROUTINE add_nmr_Fnl 
!
!
!  !-----------------------------------------------------------------------
  SUBROUTINE calc_chemical_shift(orb_magn_tot, nmr_shift_core)
!  !-----------------------------------------------------------------------
  USE kinds,                ONLY : dp
  USE io_global,            ONLY : stdout
  USE ions_base,            ONLY : nat, atm, ityp
  USE parameters,           ONLY : ntypx
  USE nmr_mod
  implicit none
  real(dp) :: orb_magn_tot(3), chemical_shift(3), m_0_mod
  real(dp) :: nmr_shift_core(ntypx)

  m_0_mod = sqrt(sum(m_0(:)**2.d0))/sqrt(2.d0)
  chemical_shift = -orb_magn_tot/m_0_mod * 1d6
  write(stdout,*)
  write(stdout,'(5X,''NUCLEAR DIPOLE ON ATOM'',I4,'' ('',A3,''):'')') &
         m_0_atom, atm(ityp(m_0_atom))
  write(stdout,'(5X,''m_0                = '',3(F14.6))') m_0
  write(stdout,'(5X,''Chemical shift (ppm):'',3(F14.4))') chemical_shift/2.d0
  write(stdout,'(5X,''Core shift     (ppm):'',F14.4)') &
        nmr_shift_core(ityp(m_0_atom))
  END SUBROUTINE calc_chemical_shift


!#if 0
!!-----------------------------------------------------------------------
!!SUBROUTINE calc_moment()
!  !
!  ! This subroutine calculates the magnetic moment
!  ! m = 1/2 int dr r x j(r), where the current j(r) is
!  ! calculated by (* indicates complex conjugate):
!  !
!  !    j(r) = j_paramagnetic(r) + j_diamagnetic(r)
!  !
!  !    j_par(r) = i e hbar/(2 m) (psi (grad psi*) - psi* (grad psi))
!  !             = e hbar/m Im( psi* (grad psi) )
!  !
!  !    j_dia(r) = -e^2/mc A_vec psi* psi
!  !
!  ! The output is in multiples of mu_b = e hbar / 2 m.
!  ! A correction according to the first moment of div(I_err) is
!  ! also applied.
!  !
!  ! Rydberg atomic units are used: m=1/2, e^2=2, hbar=1, E_0=1
!  !
!  !
!!  USE  nmr_mod
!!  USE buffers      , ONLY : get_buffer
!!  USE  io_global   , ONLY : stdout
!!  USE  io_files    , ONLY : iunwfc, nwordwfc, iunigk
!!  USE  klist       , ONLY : nks, xk, ngk
!!  USE  cell_base   , ONLY : at, alat
!!  USE  ions_base   , ONLY : tau
!!  USE  constants   , ONLY : sqrt2, tpi
!!  USE  wvfct       , ONLY : npw, nbnd, igk
!!  USE  gsmooth     , ONLY : nrxxs,nls,nr1s,nr2s,nr3s,nrx1s,nrx2s,nrx3s,ngms
!!  USE  cell_base   , ONLY : tpiba
!!  USE  gvect       , ONLY : g
!!  USE  wavefunctions_module, ONLY : evc, psic
!!  USE  mp_global   , ONLY : me_pool
!!  USE  sticks      , ONLY : dfftp
!  !
!!  IMPLICIT NONE
!  !
!!  INTEGER                :: band, ii, jj, kk, ll, ind, ik
!!  REAL(DP)               :: rx, ry, rz, pos_x, pos_y, pos_z
!!  REAL(DP)               :: I_par(3), I_dia(3), sigma_band, sigma_tot
!!  REAL(DP)               :: I_corr(nlpp_max,3), M_corr(nlpp_max,3)
!!  REAL(DP),  ALLOCATABLE :: I_tot(:,:), I_par_tot(:,:), I_dia_tot(:,:)
!!  REAL(DP),  ALLOCATABLE :: m_tot(:,:), m_par(:,:), m_dia(:,:)
!!  REAL(DP),  ALLOCATABLE :: I_err(:,:)
!  !
!  !
!!  allocate( I_tot(nbnd,3), I_par_tot(nbnd,3), I_dia_tot(nbnd,3) )
!!  allocate( m_tot(nbnd,3), m_par(nbnd,3), m_dia(nbnd,3) )
!!  if (n_nlpp > 0) allocate( I_err(nrxxs,3) )
!  !
!!  I_tot(:,:)     = 0.0D0
!!  I_par_tot(:,:) = 0.0D0
!!  I_dia_tot(:,:) = 0.0D0
!!  m_tot(:,:)     = 0.0D0
!!  m_par(:,:)     = 0.0D0
!!  m_dia(:,:)     = 0.0D0
!!  I_corr(:,:)    = 0.0D0
!!  M_corr(:,:)    = 0.0D0
!!  if (n_nlpp > 0) I_err(:,:) = 0.0D0
!  !
!  !
!  ! Calculate the current and the moment for each k-point and band ...
!  !
!!  if ( nks > 1 ) REWIND( iunigk )
!  !
!!  DO ik = 1, nks
!     !
!     ! Read the wave functions ...
!     !
!!     IF ( nks > 1 ) READ( iunigk ) igk
!!     IF ( nks > 1 ) CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
!!     npw = ngk(ik)
!     !
!!     DO band = 1, nbnd
!        !
!!        psic(1:nrxxs)               = (0.0D0, 0.0D0)
!!        grad_psi(1:nrxxs,1:3)       = (0.0D0, 0.0D0)
!!        psic    (nls(igk(1:npw)))   = evc(1:npw, band)
!!        grad_psi(nls(igk(1:npw)),1) = IM*(xk(1,ik)+g(1,igk(1:npw)))*tpiba*evc(1:npw, band)
!!        grad_psi(nls(igk(1:npw)),2) = IM*(xk(2,ik)+g(2,igk(1:npw)))*tpiba*evc(1:npw, band)
!!        grad_psi(nls(igk(1:npw)),3) = IM*(xk(3,ik)+g(3,igk(1:npw)))*tpiba*evc(1:npw, band)
!        !
!!        call cft3s (grad_psi(:,1), nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2)
!!        call cft3s (grad_psi(:,2), nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2)
!!        call cft3s (grad_psi(:,3), nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2)
!!        call cft3s (psic         , nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2)
!        !
!!        DO kk=1, dfftp%npp(me_pool+1)
!!        DO jj=1, nrx2s
!!        DO ii=1, nrx1s
!           !
!!           ind        = ii + (jj-1)*nrx1s + (kk-1)*nrx2s*nrx1s
!           !
!           ! negative sign according to Giovanni
!!           I_par(1:3) = -2.0D0 * DIMAG( dconjg(psic(ind)) * grad_psi(ind,1:3) )
!!           I_dia(1:3) = -sqrt2*fine_struct*DBLE(A_vec(ind,1:3)*dconjg(psic(ind))*psic(ind))
!           !
!!           I_par(1:3) = I_par(1:3) * 2.0D0          ! assume double occupied bands
!!           I_dia(1:3) = I_dia(1:3) * 2.0D0
!           !
!!           I_par_tot(band,1:3) = I_par_tot(band,1:3) + I_par(1:3)
!!           I_dia_tot(band,1:3) = I_dia_tot(band,1:3) + I_dia(1:3)
!!           !
!!           pos_x = dble( ii-1                      )/dble(nrx1s) - tau(1,m_0_atom)
!!           pos_y = dble( jj-1                      )/dble(nrx2s) - tau(2,m_0_atom)
!!           pos_z = dble( kk+dfftp%ipp(me_pool+1)-1 )/dble(nrx3s) - tau(3,m_0_atom)
!           !
!!           if(pos_x >  0.5D0) pos_x = pos_x - 1.0D0
!!           if(pos_x < -0.5D0) pos_x = pos_x + 1.0D0
!!           if(pos_y >  0.5D0) pos_y = pos_y - 1.0D0
!!           if(pos_y < -0.5D0) pos_y = pos_y + 1.0D0
!!           if(pos_z >  0.5D0) pos_z = pos_z - 1.0D0
!!           if(pos_z < -0.5D0) pos_z = pos_z + 1.0D0
!           !
!!           rx    = alat * ( pos_x*at(1,1) + pos_y*at(1,2) + pos_z*at(1,3) )
!!           ry    = alat * ( pos_x*at(2,1) + pos_y*at(2,2) + pos_z*at(2,3) )
!!           rz    = alat * ( pos_x*at(3,1) + pos_y*at(3,2) + pos_z*at(3,3) )
!           !
!!           m_par(band,1) = m_par(band,1) + 0.5D0 * ( ry*I_par(3) - rz*I_par(2) )
!!           m_par(band,2) = m_par(band,2) + 0.5D0 * (-rx*I_par(3) + rz*I_par(1) )
!!           m_par(band,3) = m_par(band,3) + 0.5D0 * ( rx*I_par(2) - ry*I_par(1) )
!           !
!!           m_dia(band,1) = m_dia(band,1) + 0.5D0 * ( ry*I_dia(3) - rz*I_dia(2) )
!!           m_dia(band,2) = m_dia(band,2) + 0.5D0 * (-rx*I_dia(3) + rz*I_dia(1) )
!!           m_dia(band,3) = m_dia(band,3) + 0.5D0 * ( rx*I_dia(2) - ry*I_dia(1) )
!           !
!!           if(n_nlpp > 0) I_err(ind,1:3) = I_err(ind,1:3)+I_par(1:3)+I_dia(1:3)
!           !
!!        END DO
!!        END DO
!!        END DO
!        !
!!     END DO
!     !
!!  END DO
!  !
!!  I_par_tot(1:nbnd,1:3) = I_par_tot(1:nbnd,1:3) / dble(nr1s*nr2s*nr3s * nks)
!!  I_dia_tot(1:nbnd,1:3) = I_dia_tot(1:nbnd,1:3) / dble(nr1s*nr2s*nr3s * nks)
!!  I_tot(1:nbnd,1:3)     = I_par_tot(1:nbnd,1:3) + I_dia_tot(1:nbnd,1:3)
!!  !
!!  m_par(1:nbnd,1:3)     = m_par(1:nbnd,1:3)     / dble(nr1s*nr2s*nr3s * nks)
!!  m_dia(1:nbnd,1:3)     = m_dia(1:nbnd,1:3)     / dble(nr1s*nr2s*nr3s * nks)
!!  m_tot(1:nbnd,1:3)     = m_par(1:nbnd,1:3)     + m_dia(1:nbnd,1:3)
!  !
!!  CALL reduce( 3*nbnd, I_par_tot )
!!  CALL reduce( 3*nbnd, I_dia_tot )
!!  CALL reduce( 3*nbnd, I_tot )
!!  CALL reduce( 3*nbnd, m_par )
!!  CALL reduce( 3*nbnd, m_dia )
!!  CALL reduce( 3*nbnd, m_tot )
!  !
!  !
!  ! Calculate the first moment of div(I_err) ...
!  !
!!  if (n_nlpp > 0) then
!     !
!!     grad_psi(1:nrxxs,1:3)  = I_err(1:nrxxs,1:3)
!!     call cft3s (grad_psi(:,1), nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s,-1)
!!     call cft3s (grad_psi(:,2), nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s,-1)
!!     call cft3s (grad_psi(:,3), nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s,-1)
!!     psic(1:nrxxs) = (0.0D0, 0.0D0)
!!     psic(nls(1:ngms)) = IM*tpiba * ( &
!!                         g(1,1:ngms)*grad_psi(nls(1:ngms),1) + &
!!                         g(2,1:ngms)*grad_psi(nls(1:ngms),2) + &
!!                         g(3,1:ngms)*grad_psi(nls(1:ngms),3) )
!!     call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 1)
!     !
!!     DO kk=1, dfftp%npp(me_pool+1)
!!     DO jj=1, nrx2s
!!     DO ii=1, nrx1s
!        !
!!        ind = ii + (jj-1)*nrx1s + (kk-1)*nrx2s*nrx1s
!        !
!!        DO ll = 1, nlpp_max
!           !
!!           if( r_nlpp(ll) == 0.0D0 ) cycle
!           !
!!           pos_x = dble(ii-1                     )/dble(nrx1s) - tau(1,ll)
!           pos_y = dble(jj-1                     )/dble(nrx2s) - tau(2,ll)
!           pos_z = dble(kk+dfftp%ipp(me_pool+1)-1)/dble(nrx3s) - tau(3,ll)
!           !
!           if(pos_x >  0.5D0) pos_x = pos_x - 1.0D0
!           if(pos_x < -0.5D0) pos_x = pos_x + 1.0D0
!           if(pos_y >  0.5D0) pos_y = pos_y - 1.0D0
!           if(pos_y < -0.5D0) pos_y = pos_y + 1.0D0
!           if(pos_z >  0.5D0) pos_z = pos_z - 1.0D0
!           if(pos_z < -0.5D0) pos_z = pos_z + 1.0D0
!           !
!           rx    = alat * ( pos_x*at(1,1) + pos_y*at(1,2) + pos_z*at(1,3) )
!           ry    = alat * ( pos_x*at(2,1) + pos_y*at(2,2) + pos_z*at(2,3) )
!           rz    = alat * ( pos_x*at(3,1) + pos_y*at(3,2) + pos_z*at(3,3) )
!           !
!           if(sqrt(rx**2+ry**2+rz**2) > r_nlpp(ll)) cycle
!           !
!           I_corr(ll,1) = I_corr(ll,1) + rx * real(psic(ind))
!           I_corr(ll,2) = I_corr(ll,2) + ry * real(psic(ind))
!           I_corr(ll,3) = I_corr(ll,3) + rz * real(psic(ind))
!           !
!           EXIT
!           !
!        END DO
!        !
!     END DO
!     END DO
!     END DO
!     !
!     I_corr(:,:) = I_corr(:,:) / dble(nr1s*nr2s*nr3s * nks)
!     CALL reduce( 3*nlpp_max, I_corr )
!     !
!     write (stdout, '(2A)') cr, '  Corrections to the currents from the nlpp projectors ...'
!     !
!     DO ll = 1, nlpp_max
!        !
!        if( r_nlpp(ll) == 0.0D0 ) cycle
!        !
!        write (stdout,'(A,i4,A,3e14.6)') '  Correction atom ', &
!           ll, ': ', I_corr(ll,1), I_corr(ll,2), I_corr(ll,3)
!        !
!        pos_x = tau(1,ll) - tau(1,m_0_atom)
!        pos_y = tau(2,ll) - tau(2,m_0_atom)
!        pos_z = tau(3,ll) - tau(3,m_0_atom)
!        !
!        if(pos_x >  0.5D0) pos_x = pos_x - 1.0D0
!        if(pos_x < -0.5D0) pos_x = pos_x + 1.0D0
!        if(pos_y >  0.5D0) pos_y = pos_y - 1.0D0
!        if(pos_y < -0.5D0) pos_y = pos_y + 1.0D0
!        if(pos_z >  0.5D0) pos_z = pos_z - 1.0D0
!        if(pos_z < -0.5D0) pos_z = pos_z + 1.0D0
!        !
!        rx    = alat * ( pos_x*at(1,1) + pos_y*at(1,2) + pos_z*at(1,3) )
!        ry    = alat * ( pos_x*at(2,1) + pos_y*at(2,2) + pos_z*at(2,3) )
!        rz    = alat * ( pos_x*at(3,1) + pos_y*at(3,2) + pos_z*at(3,3) )
!        M_corr(ll,1) =  0.5D0*( ry*I_corr(ll,3)-rz*I_corr(ll,2) )
!        M_corr(ll,2) =  0.5D0*(-rx*I_corr(ll,3)+rz*I_corr(ll,1) )
!        M_corr(ll,3) =  0.5D0*( rx*I_corr(ll,2)-ry*I_corr(ll,1) )
!        !
!     END DO
!     !
!  END IF
!  !
!  !
!  !
!  ! Output of results ...
!  !
!  write (stdout, '(3A)') cr, cr, '  Analyzing the currents assuming double occupied bands ...'
!  DO band = 1, nbnd
!     write (stdout,'(A,i4,A,3e14.6)') &
!     '  Paramagnetic current band ', band, ': ', I_par_tot(band,1), I_par_tot(band,2), I_par_tot(band,3)
!  END DO
!  !
!  DO band = 1, nbnd
!     write (stdout,'(A,i4,A,3e14.6)') &
!     '  Diamagnetic  current band ', band, ': ', I_dia_tot(band,1), I_dia_tot(band,2), I_dia_tot(band,3)
!  END DO
!  !
!  DO band = 1, nbnd
!     write (stdout,'(A,i4,A,3e14.6)') &
!     '  Total        current band ', band, ': ', I_tot(band,1), I_tot(band,2), I_tot(band,3)
!  END DO
!  !
!  !
!  !
!  write (stdout, '(2A)') cr, '  Analyzing the moments ...'
!  DO band = 1, nbnd
!     write (stdout,'(A,i4,A,3e14.6)') &
!     '  Paramagnetic moment  band ', band, ': ', m_par(band,1), m_par(band,2), m_par(band,3)
!  END DO
!  !
!  DO band = 1, nbnd
!     write (stdout,'(A,i4,A,3e14.6)') &
!     '  Diamagnetic  moment  band ', band, ': ', m_dia(band,1), m_dia(band,2), m_dia(band,3)
!  END DO
!  !
!  DO band = 1, nbnd
!     write (stdout,'(A,i4,A,3e14.6)') &
!     '  Total        moment  band ', band, ': ', m_tot(band,1), m_tot(band,2), m_tot(band,3)
!  END DO
!  !
!  !
!  !
!  if (m_0(m_0_comp) == 0.0D0) m_0(m_0_comp) = sqrt2
!  !
!  write (stdout, '(2A,i4)') cr, '  Analyzing the NMR shielding assuming a dipole in direction ', m_0_comp
!  !
!  sigma_tot = 0.0D0
!  DO band = 1, nbnd
!     sigma_band = -1.0D6 * sqrt2 * m_par(band, m_0_comp)/m_0(m_0_comp)
!     sigma_tot  = sigma_tot + sigma_band
!     write (stdout,'(A,i4,A,f16.4,A,f16.4)') &
!     '  Paramagnetic shielding band ', band, ': ', sigma_band, '   Total: ', sigma_tot
!  END DO
!  !   
!  sigma_tot = 0.0D0
!  DO band = 1, nbnd
!     sigma_band = -1.0D6 * sqrt2 * m_dia(band, m_0_comp)/m_0(m_0_comp)
!     sigma_tot  = sigma_tot + sigma_band
!     write (stdout,'(A,i4,A,f16.4,A,f16.4)') &
!     '  Diamagnetic  shielding band ', band, ': ', sigma_band, '   Total: ', sigma_tot
!  END DO
!  !
!  sigma_tot = 0.0D0
!  DO band = 1, nbnd
!     sigma_band = -1.0D6 * sqrt2 * m_tot(band, m_0_comp)/m_0(m_0_comp)
!     sigma_tot  = sigma_tot + sigma_band
!     write (stdout,'(A,i4,A,f16.4,A,f16.4)') &
!     '  Total        shielding band ', band, ': ', sigma_band, '   Total: ', sigma_tot
!  END DO
!  !
!  !
!  !
!  IF (n_nlpp > 0) THEN
!     !
!     write (stdout, '(2A)') cr, '  Non-local projector corrected NMR shieldings ...'
!     !
!     DO ll = 1, nlpp_max
!        !
!        if( r_nlpp(ll) == 0.0D0 ) cycle
!        !
!        sigma_band = -1.0D6 * sqrt2 * M_corr(ll, m_0_comp)/m_0(m_0_comp)
!        sigma_tot  = sigma_tot + sigma_band
!        !
!        write (stdout,'(A,i4,A,f16.4,A,f16.4)') &
!          '  Total corr.  shielding atom ', ll, &
!          ': ', sigma_band, '   Total: ', sigma_tot
!        !
!     END DO
!     !
!  END IF
!  !
!  !
!  !
!  DEALLOCATE( I_tot, I_par_tot, I_dia_tot, m_tot, m_par, m_dia)
!  if (n_nlpp > 0) DEALLOCATE( I_err )
!  !
!END SUBROUTINE calc_moment
!#endif
!
!!-----------------------------------------------------------------------
!!
!!
!!
!! END_NMR
!!
