!-----------------------------------------------------------------------
SUBROUTINE newscf
!-----------------------------------------------------------------------
!! Recalculate self-consistent
  USE constants, ONLY : bohr_radius_angs, fpi, rytoev
  USE io_global, ONLY : stdout
  USE cell_base, ONLY : omega
  USE ions_base, ONLY : nat, tau, ityp, atm
  USE fft_base,  ONLY : dfftp, dffts
  USE scf,       ONLY : rho, rho_core
  USE xc_lib,    ONLY : xclib_set_threshold, dmxc
  USE basis, ONLY: starting_wfc, starting_pot, natomwfc
  USE cellmd,ONLY: lmovecell
  USE gvecs, ONLY: doublegrid
  USE gvect, ONLY: gstart
  USE wvfct, ONLY: btype
  USE klist, ONLY: nkstot, lgauss, ltetra, nks
  USE wvfct, ONLY: nbnd, nbndx
  USE noncollin_module, ONLY: report
  USE check_stop,    ONLY : check_stop_init
  USE fft_base,      ONLY : dfftp, dffts
  USE symm_base,     ONLY : nsym
  USE io_files,      ONLY : iunwfc, prefix, tmp_dir, postfix
  USE ldaU,          ONLY : lda_plus_u
  USE control_flags, ONLY : restart, io_level, lscf, iprint, &
                            david, max_cg_iter, nexxiter, &
                            isolve, tr2, ethr, mixing_beta, nmix, niter, &
                            iverbosity, do_makov_payne, ldftd3, noinv
  USE extrapolation, ONLY : extrapolate_charge
  USE kinds,         ONLY : dp
  USE input_parameters, ONLY : startingpot, startingwfc, restart_mode
  USE dfunct,        ONLY : newd
  USE scf,        ONLY: rho, rho_core, v, vltot, vrs, kedtau
  USE lsda_mod,   ONLY: nspin, current_spin
  USE orbital_magnetization, ONLY : dvrs
  USE wvfct,    ONLY : g2kin, nbndx, nbnd
  USE gipaw_module, ONLY : conv_threshold, assume_isolated
  USE mp_pools,        ONLY : intra_pool_comm, inter_pool_comm
  USE ener,                 ONLY : ef, ef_up, ef_dw, ef_cond
  USE martyna_tuckerman, ONLY : do_comp_mt
  USE nmr_mod,    ONLY : m_0
  !
  IMPLICIT NONE
  !
  INTEGER :: i, j, nrxxs, ik, occ
  CHARACTER(LEN=256) :: dirname
  INTEGER :: iter
  LOGICAL               :: exst
  REAL(DP), DIMENSION(3,3) :: chi(3,3)
  REAL(DP), DIMENSION(dfftp%nnr,1) ::  rhotot, sign_r
  REAL(DP) :: exxen
  REAL(dp) :: ehomo, elumo ! highest occupied and lowest unoccupied levels
  real(dp) :: emin, emax
  allocate (dvrs (dfftp%nnr, nspin, 3))

  
  
  !  dft='starting from scratch'
  restart  =.false.
  io_level = 0
  lscf=.true.
  lda_plus_u=.false.
  doublegrid=.false.
  lmovecell=.false.
  iprint=10000
!  read wfc and potential from preav scf. Generally is not suggested
!  starting_wfc='file'  ! read wfc from preav. scf
!  starting_pot='file'  ! read potential from preav scf
! Initialization of wfc
  starting_wfc='atomic'
  report=0
  CALL check_stop_init()
  CALL setup_para ( dfftp%nr3, 1, nbnd )
  CALL export_gstart_2_solvers(gstart)
  if ( .not. allocated (btype) ) then
     allocate( btype( nbnd, nkstot ) )
     btype(:,:) = 1
  end if
  !
  !
  SELECT CASE( trim( assume_isolated ) )
      !
    CASE( 'makov-payne', 'm-p', 'mp' )
      !
      do_makov_payne = .true.
      !
    CASE( 'martyna-tuckerman', 'm-t', 'mt' )
      !
      do_comp_mt     = .true.
      !
      CONTINUE
      !
    CASE DEFAULT
      !
      do_makov_payne = .false.
      do_comp_mt     = .false.
      !
  END SELECT
  !  we don't need symmetries
  !
  nsym=1
  noinv=.true.
  !
  ! these must be tuned for fast convergence
  !
  david = 6
  nbndx = max (nbndx, david*nbnd)
  isolve=0
  ethr = conv_threshold
  !
  nmix=8
  niter=100
  nexxiter=100
  !
  !call openfil !DFT+U(V) implementation
  call summary ( )
  call hinit0 ( )
  call potinit ( )
  
  CALL set_dvrs( dvrs, vrs, dfftp%nnr, nspin )

  call newd ( )
  call wfcinit_gipaw ( )
  CALL electrons_gipaw ( )
  ! 
  ! WARNING in NMR calculation : for non-metallic system (occupation='fixed') at least 1 unoccupied of Kohn-Sham states (nbnd) must be added
  !
  ! compute the Fermi Energy for not a metal 
  IF ( .NOT. lgauss .OR. ltetra ) THEN
          ! ... presumably not a metal
          if (any(m_0 /= 0.d0)) then
           do ik = 1, nks
             call find_nbnd_occ(ik, occ, emin, emax)
              if ( occ .eq. nbnd ) then
               print*, 'occupied bands =',occ, 'number of total bands =',nbnd
               print*, 'ERROR in NMR calculation! for non-metallic system &
                       at least 1 unoccupied of Kohn-Sham state must be added! nbnd =',occ+1
               stop
               endif
            enddo
           endif
          CALL get_homo_lumo (ehomo, elumo)
          ef = (ehomo + elumo)/2.d0
  ENDIF
  !
  CLOSE(unit=iunwfc, status='keep')
  !
  !
  RETURN
END SUBROUTINE newscf




