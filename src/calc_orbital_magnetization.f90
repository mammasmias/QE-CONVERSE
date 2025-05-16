!-----------------------------------------------------------------------
  SUBROUTINE calc_orbital_magnetization
  !-----------------------------------------------------------------------
  USE kinds,                 ONLY : dp
  USE io_files,              ONLY : nwordwfc, iunwfc
  USE wvfct,                 ONLY : npw, g2kin, nbnd, npwx, wg, et, &
                                   current_k
  USE klist,                 ONLY : xk, nks, lgauss, igk_k, ngk
  USE fft_base,              ONLY : dffts, dfftp
  USE cell_base,             ONLY : tpiba2, tpiba
  USE gvecw,                 ONLY : ecutwfc
  USE gvect,                 ONLY : ngm, g
  USE io_global,             ONLY : stdout
  USE io_files,              ONLY : iunsat, nwordatwfc
  USE uspp,                  ONLY : nkb, vkb
  USE becmod,                ONLY : becp, allocate_bec_type, calbec
  USE becmod,                ONLY : deallocate_bec_type
  USE wavefunctions,         ONLY : evc
  USE lsda_mod,              ONLY : current_spin, lsda, isk, nspin
  USE gipaw_module,          ONLY : lambda_so, dudk_method
  USE paw_gipaw,             ONLY : paw_vkb, paw_nkb, paw_becp
  USE ener,                  ONLY : ef
!  USE ldaU,                 ONLY : swfcatom, lda_plus_u
  USE buffers,               ONLY : get_buffer
  USE scf,                   ONLY : vrs
  USE mp_global,             ONLY : my_pool_id
  USE mp_bands,              ONLY : intra_bgrp_comm, inter_bgrp_comm
  USE mp_world,              ONLY : mpime, world_comm
  USE mp,                    ONLY : mp_barrier, mp_sum
  USE mp_pools,              ONLY : my_pool_id, me_pool, root_pool,  &
                                    inter_pool_comm, intra_pool_comm
  USE mp_images,             ONLY : my_image_id, inter_image_comm, nimage, &
                                    intra_image_comm
  USE gipaw_module,          ONLY : q_gipaw, nmr_shift_core
  USE uspp_init,             ONLY : init_us_2
  USE gvecw,                 ONLY : gcutw
  USE control_flags,         ONLY : iverbosity
  USE orbital_magnetization
  USE nmr_mod

  implicit none
  complex(dp), external :: ZDOTC
  integer, parameter :: iundudk1 = 75, iundudk2 = 76, iundudk3 = 77
  real(dp), parameter :: rydtohar = 0.5d0
  complex(dp), allocatable :: dudk_bra(:,:), dudk_ket(:,:), hpsi(:)
  complex(dp), allocatable :: vkb_save(:,:), aux(:,:)
  complex(dp) :: braket
  real(dp) :: kp_berry(3), kp_berry2(3), kp_M_LC(3), kp_M_IC(3), tmp1(3), tmp2(3)
  integer :: ik, ibnd, jbnd, kk, ii, jj, occ, nrxxs, nr1, nr2, nr3
  real(dp) :: tmp(3), emin, emax
  ! index for the cross product
  integer :: ind(2,3)
  
  !debug mpi 
  integer :: ierr, num_procs, rank, root, i, dd, comm_rank
  integer :: ibnd_start, ibnd_end
  real :: partial_value, gathered_kp_berry(6)

  ind(:,1) = (/ 2, 3 /)
  ind(:,2) = (/ 3, 1 /)
  ind(:,3) = (/ 1, 2 /)

  ! various initializations:
  nrxxs = dffts%nnr
  nr1=dfftp%nr1
  nr2=dfftp%nr2
  nr3=dfftp%nr3

  do ik = 1, nks
    npw = ngk (ik) !number of plane waves for each k point
  enddo
  
  delta_k = q_gipaw/tpiba
   
!  compute the covariant derivatives
  call compute_dudk_new (dudk_method)
  ! allocate memory
  allocate (dudk_bra(npwx,nbnd), dudk_ket(npwx,nbnd), hpsi(npwx))
  
  call start_clock ('orbital_magnetization')
  
  write(stdout,*)
  write(stdout,'(5X,''Computing the orbital magnetization (bohr mag/cell):'')')
  berry_curvature = 0.d0
  orb_magn_LC = 0.d0
  orb_magn_IC = 0.d0
  delta_M_bare = 0.d0
  delta_M_para = 0.d0
  delta_M_dia = 0.d0

  ! allocate the derivatives of the projectors
  call allocate_bec_type(nkb, nbnd, becp)
  allocate(dbecp(nkb,nbnd,3), paw_dbecp(paw_nkb,nbnd,3))
  allocate(vkb_save(npwx,nkb), aux(nkb,nbnd))
#define __USE_BARRIER
   
  CALL set_dvrs( dvrs, vrs, dfftp%nnr, nspin ) 
   ! loop over k-points
  do ik = 1, nks
    npw = ngk(ik)
    call find_nbnd_occ(ik, occ, emin, emax)
    if (lgauss) occ = nbnd
    !
    !band parallelization
    call divide(inter_bgrp_comm, occ, ibnd_start, ibnd_end)
    !
    ! setup the hamiltonian
    current_k = ik
    current_spin = 1
    if (lsda) current_spin = isk(ik)
     CALL gk_sort( xk(1,ik), ngm, g, gcutw, ngk(ik), igk_k(1,ik), g2kin )
     g2kin(1:ngk(ik)) = g2kin(1:ngk(ik)) * tpiba2
    call get_buffer(evc, nwordwfc, iunwfc, ik)
    if (nkb > 0) then
      call init_us_2(ngk(ik), igk_k(1,ik), xk(1,ik), vkb)
      CALL calbec( npw, vkb, evc, becp, nbnd )
    endif
!    if (lda_plus_u) call davcio(swfcatom, nwordatwfc, iunsat, ik, -1)
   
    ! compute the diamagnetic terms
    call init_gipaw_2(ngk(ik), igk_k(1,current_k), xk(1,current_k), paw_vkb)
    call calbec( npw, paw_vkb, evc, paw_becp, nbnd )
    if (any(m_0 /= 0.d0))       call calc_delta_M_dia_nmr
    if (any(lambda_so /= 0.d0)) call calc_delta_M_dia_so

    call compute_dbecp  ! for deltaM bare
    call compute_paw_dbecp ! for delta_M_para_so

    ! loop over the magnetization directions
    do kk =  1, 3
      ii = ind(1,kk)
      jj = ind(2,kk)
      ! read the bra and the ket
      call davcio(dudk_bra, 2*nwordwfc, iundudk1 + ii - 1, ik, -1)
      call davcio(dudk_ket, 2*nwordwfc, iundudk1 + jj - 1, ik, -1)
      
      ! compute the orbital magnetization
      kp_berry(kk) = 0.d0
      kp_M_IC(kk) = 0.d0
      kp_M_LC(kk) = 0.d0
      do ibnd = ibnd_start, ibnd_end
       ! IC term and Berry curvature
        braket = zdotc(ngk(ik), dudk_bra(1,ibnd), 1, dudk_ket(1,ibnd), 1)
        kp_berry(kk) = kp_berry(kk) + 2.d0*wg(ibnd,ik)*imag(braket)
        berry_curvature(kk) = berry_curvature(kk) + &
                              2.d0*wg(ibnd,ik)*imag(braket)
        kp_M_IC(kk) = kp_M_IC(kk) + 2.d0*wg(ibnd,ik)*et(ibnd,ik)*imag(braket)
        orb_magn_IC(kk) = orb_magn_IC(kk) + &
                         2.d0*wg(ibnd,ik)*et(ibnd,ik)*imag(braket)

        ! this is the LC term
        call h_psi_gipaw(npwx, ngk(ik), 1, dudk_ket(1:npwx,ibnd), hpsi)
        braket = zdotc(ngk(ik), dudk_bra(1,ibnd), 1, hpsi, 1)
        kp_M_LC(kk) = kp_M_LC(kk) + wg(ibnd,ik)*imag(braket)
        orb_magn_LC(kk) = orb_magn_LC(kk) + wg(ibnd,ik)*imag(braket)
        
        call h_psi_gipaw(npwx, ngk(ik), 1, dudk_bra(1:npwx,ibnd), hpsi)
        braket = zdotc(ngk(ik), dudk_ket(1,ibnd), 1, hpsi, 1)
        kp_M_LC(kk) = kp_M_LC(kk) - wg(ibnd,ik)*imag(braket)
        orb_magn_LC(kk) = orb_magn_LC(kk) - wg(ibnd,ik)*imag(braket)

      enddo
     ! compute the GIPAW corrections
      call calc_delta_M_bare
      if (any(lambda_so /= 0.d0)) call calc_delta_M_para_so
      if (any(m_0 /= 0.d0))       call calc_delta_M_para_nmr
    enddo ! kk
    ! Parallel reductions
#ifdef __MPI
  call mp_sum( kp_berry, intra_bgrp_comm )
  call mp_sum( kp_M_LC, intra_bgrp_comm )
  call mp_sum( kp_M_IC, intra_bgrp_comm )
  call mp_sum( kp_berry, inter_bgrp_comm )
  call mp_sum( kp_M_LC, inter_bgrp_comm )
  call mp_sum( kp_M_IC, inter_bgrp_comm )
#endif
  
    if (me_pool == root_pool) then
      write(*,'(''BC: k-point:'',I5,2X,''pool:'',I4,4X,9F12.6)') ik, my_pool_id+1, kp_berry
      write(*,'(''LC: k-point:'',I5,2X,''pool:'',I4,4X,9F12.6)') ik, my_pool_id+1, kp_M_LC*rydtohar
      write(*,'(''IC: k-point:'',I5,2X,''pool:'',I4,4X,9F12.6)') ik, my_pool_id+1, kp_M_IC*rydtohar
!      write(*,'(5X,''k-point:'',I4,4X,''pool:'',I4,4X,3(F10.4),/,5X,(9F12.6))') &
!            ik, my_pool_id+1, xk(:,ik), kp_berry, kp_M_LC*rydtohar, kp_M_IC*rydtohar
    endif
  enddo !ik
  
#ifdef __MPI
 ! no reduction for delta_M_bare and delta_M_para and delta_M_dia
  call mp_sum(orb_magn_LC, intra_bgrp_comm )
  call mp_sum(orb_magn_IC, intra_bgrp_comm )
  call mp_sum(berry_curvature, intra_bgrp_comm )
  call mp_sum(orb_magn_LC, inter_bgrp_comm )
  call mp_sum(orb_magn_IC, inter_bgrp_comm )
  call mp_sum(berry_curvature, inter_bgrp_comm )
#endif
  
  ! no reduction for delta_M_bare and delta_M_para and delta_M_dia 


  ! close files
  close(unit=iundudk1, status='keep')
  close(unit=iundudk2, status='keep')
  close(unit=iundudk3, status='keep')

  orb_magn_tot = orb_magn_LC + orb_magn_IC + &
                 delta_M_bare + delta_M_dia + delta_M_para

  write(stdout,*)
 ! print results
  write(stdout,'(5X,''SPIN-ORBIT:'')')
  write(stdout,'(5X,''lambda_so          = '',3(F14.6))') lambda_so
  write(stdout,*)

  write(stdout,'(5X,''Berry curvature    = '',3(F14.6))') berry_curvature
  write(stdout,'(5X,''Fermi energy       = '',F14.6,'' rydberg'')') ef
  write(stdout,*)

  write(stdout,'(5X,''NUCLEAR DIPOLE ON ATOM'',I4,'':'')') m_0_atom
  write(stdout,'(5X,''m_0                = '',3(F14.6))') m_0
  write(stdout,*)

  if (iverbosity > 0) then
    write(stdout,'(5X,''(without Berry curvature term)'')')
    write(stdout,'(5X,''M_LC               = '',3(F14.6))') orb_magn_LC
    write(stdout,'(5X,''M_IC               = '',3(F14.6))') orb_magn_IC
    write(stdout,'(5X,''Delta_M_bare       = '',3(F14.6))') delta_M_bare
    write(stdout,'(5X,''Delta_M_para       = '',3(F14.6))') delta_M_para
    write(stdout,'(5X,''Delta_M_dia        = '',3(F14.6))') delta_M_dia
    write(stdout,'(5X,''M_tot              = '',3(F14.6))') orb_magn_tot
    write(stdout,*)
    write(stdout,'(5X,''(with Berry curvature term)'')')
   endif
  orb_magn_LC = orb_magn_LC - ef*berry_curvature
  orb_magn_IC = orb_magn_IC - ef*berry_curvature
  orb_magn_tot = orb_magn_tot - 2.d0*ef*berry_curvature
  write(stdout,'(5X,''M_LC               = '',3(F14.6))') orb_magn_LC
  write(stdout,'(5X,''M_IC               = '',3(F14.6))') orb_magn_IC
  write(stdout,'(5X,''Delta_M            = '',3(F14.6))') &
        delta_M_bare + delta_M_para + delta_M_dia
  write(stdout,'(5X,''M_tot              = '',3(F14.6))') orb_magn_tot

  ! free memory
  CALL deallocate_bec_type ( becp )
  deallocate( dudk_bra, dudk_ket, hpsi )

  ! close files
  close(unit=iundudk1)
  close(unit=iundudk2)
  close(unit=iundudk3)
  
  call stop_clock ('orbital_magnetization')
  ! go on, reporting the g-tensor
  if (any(lambda_so /= 0.d0)) call calc_g_tensor(orb_magn_tot)
  ! go on, reporting the chemical shift
  if (any(m_0 /= 0.d0)) call calc_chemical_shift(orb_magn_tot, nmr_shift_core)




  CONTAINS
    !------------------------------------------------------------------
    ! derivative of the beta's
    !------------------------------------------------------------------
    SUBROUTINE compute_dbecp
    USE cell_base, ONLY : tpiba
    USE uspp_init, ONLY : init_us_2
    implicit none
    integer :: ipol, sig
    real(dp) :: kq(3)
    if (nkb == 0) return
    vkb_save = vkb
    do ipol = 1, 3
      dbecp(:,:,ipol) = (0.d0,0.d0)
      do sig = -1, 1, 2
        kq(:) = xk(:,ik)
        kq(ipol) = kq(ipol) + sig * delta_k
        call init_us_2_no_phase (npw, igk_k(1,ik), kq, vkb)
        call calbec( npw, vkb, evc, aux, nbnd )
        dbecp(:,:,ipol) = dbecp(:,:,ipol) + &
                          0.5d0*sig/(delta_k*tpiba) * aux(:,:)
      enddo
!      print*, 'dbecp', dbecp(:,:,ipol), ipol
    enddo
    vkb = vkb_save
    END SUBROUTINE compute_dbecp

    !------------------------------------------------------------------
    ! derivative of GIPAW projectors
    !------------------------------------------------------------------
    SUBROUTINE compute_paw_dbecp
    USE cell_base, ONLY : tpiba
!    USE g_tensor_module, ONLY : init_paw_2_no_phase
    implicit none
    integer :: ipol, sig
    real(dp) :: kq(3)
    if (paw_nkb == 0) return
    do ipol = 1, 3
      paw_dbecp(:,:,ipol) = (0.d0,0.d0)
      do sig = -1, 1, 2
        kq(:) = xk(:,ik)
        kq(ipol) = kq(ipol) + sig * delta_k
        !!!!call init_paw_2(npw, igk, kq, paw_vkb)
        call init_gipaw_2_no_phase(ngk(ik), igk_k(1,ik), kq, paw_vkb)
        call calbec( npw, paw_vkb, evc, paw_becp, nbnd )
        paw_dbecp(:,:,ipol) = paw_dbecp(:,:,ipol) + &
                              0.5d0*sig/(delta_k*tpiba) * paw_becp(:,:)
      enddo
    enddo
    END SUBROUTINE compute_paw_dbecp

    !------------------------------------------------------------------
    ! GIPAW correction (delta_M_bare)
    !------------------------------------------------------------------
    SUBROUTINE calc_delta_M_bare
    USE ions_base,  ONLY : nat, ntyp => nsp, ityp
    USE uspp_param, ONLY : nh
    USE uspp,       ONLY : deeq, indv, nhtol
    implicit none
    complex(dp) :: tmp, becp_product
    integer :: ibnd, ijkb0, nt, na, jh, jkb, ih, ikb
    integer :: nbs1, nbs2, l1, l2
    if (nkb == 0) return
    tmp = (0.d0,0.d0)
    ijkb0 = 0
    do nt = 1, ntyp
      do na = 1, nat
        if ( ityp(na) == nt ) then
          do jh = 1, nh(nt)
            jkb = ijkb0 + jh
            nbs1 = indv(jh,nt)
            l1 = nhtol(jh,nt)
            do ih = 1, nh(nt)
              ikb = ijkb0 + ih
              nbs2 = indv(ih,nt)
              l2 = nhtol(ih,nt)
              if (l1 /= l2) cycle
              do ibnd = 1, occ
                becp_product = conjg(dbecp(ikb,ibnd,ii)) * dbecp(jkb,ibnd,jj)
                tmp = tmp + wg(ibnd,ik) * deeq(ih,jh,na,current_spin) * &
                            becp_product
              enddo
            enddo
          enddo
          ijkb0 = ijkb0 + nh(nt)
        endif
      enddo
    enddo
    !PRINT*, mpime, kk, -2.d0*tmp
    ! check the sign and real or imag!!
    delta_M_bare(kk) = delta_M_bare(kk) - 2.d0*imag(tmp)
    END SUBROUTINE calc_delta_M_bare


    !------------------------------------------------------------------
    ! GIPAW correction (Delta_M_para), SO case
    !------------------------------------------------------------------
    SUBROUTINE calc_delta_M_para_so
    USE constants,  ONLY : e2
    USE ions_base,  ONLY : nat, ntyp => nsp, ityp
    USE uspp_param, ONLY : nh
    USE uspp,       ONLY : deeq
    USE paw_gipaw,  ONLY : paw_vkb, paw_nkb, paw_recon, paw_lmaxkb
    USE gipaw_module, ONLY : radial_integral_paramagnetic_so
    USE gipaw_module, ONLY : lx, ly, lz
    USE orbital_magnetization, ONLY : a2gp4, a2gp8
    implicit none
    complex(dp) :: tmp, becp_product
    integer :: ibnd, ijkb0, nt, na, jh, jkb, ih, ikb
    integer :: l1, m1, lm1, l2, m2, lm2, nbs1, nbs2
    real(dp) :: sigma
    if (paw_nkb == 0) return
    sigma = 1.d0; if (current_spin == 2) sigma = -1.d0
    tmp = (0.d0,0.d0)
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
              if ( l1 /= l2 ) cycle
              if (l1 == 0) cycle
              do ibnd = 1, occ
                becp_product = conjg(paw_dbecp(ikb,ibnd,ii)) * &
                                     paw_dbecp(jkb,ibnd,jj)
                tmp = tmp + lambda_so(1) * a2gp8 * &
                      sigma * lx(lm1,lm2) * (0.d0,1.d0) * &
                      radial_integral_paramagnetic_so(nbs1,nbs2,nt) * &
                      wg(ibnd,ik) * becp_product
                tmp = tmp + lambda_so(2) * a2gp8 * &
                      sigma * ly(lm1,lm2) * (0.d0,1.d0) * &
                      radial_integral_paramagnetic_so(nbs1,nbs2,nt) * &
                      wg(ibnd,ik) * becp_product
                tmp = tmp + lambda_so(3) * a2gp8 * &
                      sigma * lz(lm1,lm2) * (0.d0,1.d0) * &
                      radial_integral_paramagnetic_so(nbs1,nbs2,nt) * &
                      wg(ibnd,ik) * becp_product
              enddo   ! ibnd
            enddo   ! jh
          enddo   ! ih
          ijkb0 = ijkb0 + paw_recon(nt)%paw_nh
        endif   ! if (ityp...
      enddo   ! na
    enddo  ! nt
    !PRINT*, kk, tmp
    ! check the sign and real or imag!!
    delta_M_para(kk) = delta_M_para(kk) - 2.d0*imag(tmp)
    END SUBROUTINE calc_delta_M_para_so

    !------------------------------------------------------------------
    ! GIPAW correction (delta_M_dia), SO case
    !------------------------------------------------------------------
    SUBROUTINE calc_delta_M_dia_so
    USE ions_base,  ONLY : nat, ntyp => nsp, ityp
    USE upf_params, ONLY : lmaxx !lmaxx = 3?
    USE constants,  ONLY : pi
    USE uspp_param, ONLY : nh
    USE uspp,       ONLY : deeq, ap
    USE paw_gipaw,        ONLY : paw_vkb, paw_nkb, paw_recon, paw_lmaxkb
    USE gipaw_module, ONLY : lx, ly, lz, alpha
    USE orbital_magnetization, ONLY : a2gp4, a2gp8
    USE gipaw_module, ONLY : gprime, radial_integral_diamagnetic_so
    IMPLICIT NONE
    integer :: l1, m1, lm1, l2, m2, lm2, ih, ikb, nbs1, jh, jkb, nbs2
    integer :: nt, ibnd, na, lm, nrc, ijkb0, lmaxx2
    complex(dp) :: bec_product
    complex(dp), allocatable :: dia_corr(:)
    real(dp) :: diamagnetic_tensor(3,3), sigma
    lmaxx2 = 3
    if (paw_nkb == 0) return
    sigma = 1.d0; if (current_spin == 2) sigma = -1.d0

    allocate(dia_corr(lmaxx2**2))
    dia_corr = 0.0_dp
    diamagnetic_tensor = 0.d0

    ijkb0 = 0
    do nt = 1, ntyp
      do na = 1, nat
        if (ityp (na) == nt) then
          do ih = 1, paw_recon(nt)%paw_nh
            ikb = ijkb0 + ih
                       nbs1 = paw_recon(nt)%paw_indv(ih)
            l1 = paw_recon(nt)%paw_nhtol(ih)
            m1 = paw_recon(nt)%paw_nhtom(ih)
            lm1 = m1+l1**2
            do jh = 1, paw_recon(nt)%paw_nh
              jkb = ijkb0 + jh
              nbs2 = paw_recon(nt)%paw_indv(jh)
              l2 = paw_recon(nt)%paw_nhtol(jh)
              m2 = paw_recon(nt)%paw_nhtom(jh)
              lm2 = m2 + l2**2

              do ibnd = 1, nbnd
                bec_product = paw_becp(jkb,ibnd) * conjg(paw_becp(ikb,ibnd))
                !!!PRINT*, ikb,jkb,bec_product
                !<apsi> s/non-trace-zero component
                ! 2/3 to separate the non-trace vanishing component
                ! 1/(2c^2) from the equation (59) in PM-PRB
                if ( l1 == l2 .AND. m1 == m2 ) then
                  diamagnetic_tensor(1,1) = diamagnetic_tensor(1,1) &
                    + 2.0_dp / 3.0_dp * bec_product &
                    * radial_integral_diamagnetic_so(nbs1,nbs2,nt) &
                    * wg(ibnd,ik) * alpha**2.d0
                  diamagnetic_tensor(2,2) = diamagnetic_tensor(2,2) &
                    + 2.0_dp / 3.0_dp * bec_product &
                    * radial_integral_diamagnetic_so(nbs1,nbs2,nt) &
                    * wg(ibnd,ik) * alpha**2.d0
                  diamagnetic_tensor(3,3) = diamagnetic_tensor(3,3) &
                    + 2.0_dp / 3.0_dp * bec_product &
                    * radial_integral_diamagnetic_so(nbs1,nbs2,nt) &
                    * wg(ibnd,ik) * alpha**2.d0
                endif

                ! 2/3 to separate the non-trace vanishing component
                do lm = 5, 9
                  dia_corr(lm) = dia_corr(lm) + bec_product / 3.0_dp &
                    * radial_integral_diamagnetic_so(nbs1,nbs2,nt) &
                    * ap(lm,lm1,lm2) * wg(ibnd,ik) * alpha**2.d0
                enddo

              enddo  ! ibnd
            enddo  ! jh
          enddo  ! ih
          ijkb0 = ijkb0 + paw_recon(nt)%paw_nh
        endif
      enddo  ! na
    enddo  ! nt

    !  transform in cartesian coordinates
    dia_corr(5:9) = - sqrt(4.0_dp * pi/5.0_dp) * dia_corr(5:9)
    diamagnetic_tensor(1,1) = diamagnetic_tensor(1,1) &
      + sqrt(3.0_dp) * dia_corr(8) - dia_corr(5)
    diamagnetic_tensor(2,2) = diamagnetic_tensor(2,2) &
      - sqrt(3.0_dp) * dia_corr(8) - dia_corr(5)
    diamagnetic_tensor(3,3) = diamagnetic_tensor(3,3) &
      + dia_corr(5) * 2.0_dp
    diamagnetic_tensor(1,2) = diamagnetic_tensor(1,2) &
      +  dia_corr(9) * sqrt(3.0_dp)
    diamagnetic_tensor(2,1) = diamagnetic_tensor(1,2)
    diamagnetic_tensor(1,3) = diamagnetic_tensor(1,3) &
      - dia_corr(6) * sqrt(3.0_dp)
    diamagnetic_tensor(3,1) = diamagnetic_tensor(1,3)
    diamagnetic_tensor(2,3) = diamagnetic_tensor(2,3) &
      - dia_corr(7) * sqrt(3.0_dp)
    diamagnetic_tensor(3,2) = diamagnetic_tensor(2,3)
    deallocate (dia_corr)

    !!!WRITE(*,'(''dia:'',/,3(3(E14.4,2X),/))') diamagnetic_tensor
    !!! </uwg>  (gprime/4.d0) correction factor added, compare GIPAW
    delta_M_dia = delta_M_dia + (gprime/4.d0)*0.5d0*sigma &
                                 *matmul(diamagnetic_tensor, lambda_so)

    END SUBROUTINE calc_delta_M_dia_so

  !-----------------------------------------------------------------------
  SUBROUTINE calc_g_tensor(orb_magn_tot)
  !-----------------------------------------------------------------------
  USE kinds,                ONLY : dp
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE wvfct,                ONLY : npw, g2kin, nbnd, npwx, wg, et, &
                                   current_k
  USE klist,                ONLY : xk, nks, igk_k, tot_magnetization
  USE cell_base,            ONLY : tpiba2
  USE gvect,                ONLY : ngm, g
  USE io_global,            ONLY : stdout
  USE wavefunctions,        ONLY : evc
  USE lsda_mod,             ONLY : current_spin, lsda, isk
  USE control_flags,        ONLY : iverbosity
  USE buffers,              ONLY : get_buffer
  USE gvecw,                ONLY : gcutw
  USE paw_gipaw,            ONLY : paw_vkb, paw_nkb, paw_recon, paw_lmaxkb, &
                                   paw_becp
  implicit none
  real(dp), parameter :: rydtohar = 0.5d0
  real(dp) :: orb_magn_tot(3)
  real(dp) :: delta_rmc, delta_rmc_gipaw, lambda_mod, delta_g_orb_magn(3)
  real(dp) :: tmp, delta_g_total(3)
  integer :: ik, ibnd, ipol, s_weight, unp_ele

  !--------------------------------------------------
  ! Relativistic Mass Correction and diamagnetic term
  !--------------------------------------------------
  delta_rmc = 0.d0
  delta_rmc_gipaw = 0.d0
  unp_ele = nint(tot_magnetization)

  do ik = 1, nks
    call get_buffer(evc, nwordwfc, iunwfc, ik)

    current_k = ik
    if (lsda) current_spin = isk(ik)
    s_weight = 1; if (current_spin == 2) s_weight = -1

    CALL gk_sort( xk(1,ik), ngm, g, gcutw, ngk(ik), igk_k(1,ik), g2kin )
    g2kin(1:npw) = g2kin(1:npw) * tpiba2
    do ibnd = 1, nbnd
      delta_rmc = delta_rmc - s_weight &
            * wg(ibnd,ik) * sum(g2kin(1:npw) &
            * conjg(evc(1:npw,ibnd)) * evc(1:npw,ibnd))
    end do

    call init_gipaw_2(ngk(ik), igk_k(1,current_k), xk(1,current_k), paw_vkb)
    call calbec( ngk(ik), paw_vkb, evc, paw_becp, nbnd )


    call relativistic_mass_correction (ik, tmp)
    delta_rmc_gipaw = delta_rmc_gipaw + s_weight * tmp

  enddo
! no reduction for delta_rmc_gipaw
#if defined(__MPI)
  CALL mp_sum( delta_rmc, intra_bgrp_comm )
#endif


  delta_rmc = delta_rmc * rydtohar * alpha**2.d0 * g_e * 1d6 / unp_ele
  delta_rmc_gipaw = delta_rmc_gipaw * rydtohar * alpha**2.d0 * g_e * 1d6 / unp_ele

  ! report results
  lambda_mod = sqrt(sum(lambda_so(:)**2.d0))
  delta_g_orb_magn = orb_magn_tot/lambda_mod * 1d6/ unp_ele
  write(stdout,*)
  write(stdout,'(5X,''SPIN-ORBIT:'')')
  write(stdout,'(5X,''lambda_so          = '',3(F14.6))') lambda_so
  write(stdout,'(5X,''unpaired electrons = '',3(I10))') unp_ele
  write(stdout,'(5X,''Contributions to the g-tensor (ppm):'')')
  write(stdout,'(5X,''delta_g RMC        = '',F14.4)') delta_rmc
  write(stdout,'(5X,''delta_g RMC(GIPAW) = '',F14.4)') delta_rmc_gipaw
  write(stdout,'(5X,''delta_g SO         = '',3(F14.4))') delta_g_orb_magn

  ! put all terms together
  !
  ! delta_g_total = delta_g_orb_magn + delta_g_dia
  ! delta_g_dia still allready in delta_g_orb_magn included
  !
  delta_g_total = delta_g_orb_magn
  do ipol = 1, 3
    if (lambda_so(ipol) == 0.d0) cycle
    delta_g_total(ipol) = delta_g_total(ipol) + delta_rmc+delta_rmc_gipaw
  enddo
  write(stdout,*)
  write(stdout,'(5X,''delta_g total      = '',3(F14.4))') delta_g_total

  END SUBROUTINE calc_g_tensor

  !====================================================================
  ! GIPAW contribution to the RMC term
  !====================================================================
  SUBROUTINE relativistic_mass_correction (ik, rmc_gipaw)
 ! USE atom,              ONLY : r, rab
  USE ions_base,         ONLY : nat, ityp, ntyp => nsp
  USE paw_gipaw,         ONLY : paw_vkb, paw_nkb, paw_recon, paw_lmaxkb, &
                                paw_becp
  USE wvfct,             ONLY : wg, nbnd
  USE gipaw_module,      ONLY : radial_integral_rmc
  IMPLICIT NONE
  real(dp), intent(inout):: rmc_gipaw
  integer :: ik, l1, m1, lm1, l2, m2, lm2, ih, ikb, nbs1, jh, jkb, nbs2
  integer :: nt, ibnd, na, lm, nrc, ijkb0
  complex(dp) :: rmc_corr
  complex(dp) :: bec_product

  rmc_corr = 0.0_dp

  do ibnd = 1, nbnd
    ijkb0 = 0
    do nt = 1, ntyp
      do na = 1, nat
        if (ityp (na) == nt) then
          do ih = 1, paw_recon(nt)%paw_nh
            ikb = ijkb0 + ih
            nbs1=paw_recon(nt)%paw_indv(ih)
            l1=paw_recon(nt)%paw_nhtol(ih)
            m1=paw_recon(nt)%paw_nhtom(ih)
            lm1=m1+l1**2
            do jh = 1, paw_recon(nt)%paw_nh
              jkb = ijkb0 + jh
              nbs2=paw_recon(nt)%paw_indv(jh)
              l2=paw_recon(nt)%paw_nhtol(jh)
              m2=paw_recon(nt)%paw_nhtom(jh)
              lm2=m2+l2**2
              if ( l1 /= l2 ) cycle
              if ( m1 /= m2 ) cycle
              bec_product = paw_becp(jkb,ibnd) * conjg(paw_becp(ikb,ibnd))
              rmc_corr = rmc_corr + bec_product &
                       * radial_integral_rmc(nbs1,nbs2,nt) * wg(ibnd,ik)
             enddo
            enddo
          ijkb0 = ijkb0 + paw_recon(nt)%paw_nh
        endif
      enddo
    enddo
  enddo
  rmc_gipaw = -real(rmc_corr,dp)
  END SUBROUTINE relativistic_mass_correction
  !
  !------------------------------------------------------------------
  ! GIPAW correction (Delta_M_dia), NMR case
  !------------------------------------------------------------------
  !
    SUBROUTINE calc_delta_M_dia_nmr
    USE ions_base,  ONLY : nat, ntyp => nsp, ityp
    USE upf_params, ONLY : lmaxx 
    USE constants,  ONLY : pi
    USE uspp_param, ONLY : nh
    USE uspp,       ONLY : deeq, ap
    USE paw_gipaw,        ONLY : paw_vkb, paw_nkb, paw_recon, paw_lmaxkb
    USE gipaw_module, ONLY : lx, ly, lz, alpha
    USE nmr_mod,    ONLY : m_0, m_0_atom, fine_struct
    USE gipaw_module, ONLY : radial_integral_diamagnetic
    IMPLICIT NONE
    complex(dp) :: tmp, bec_product
    integer :: ibnd, ijkb0, nt, na, jh, jkb, ih, ikb
    integer :: l1, m1, lm1, l2, m2, lm2, nbs1, nbs2, lm
    complex(dp), allocatable :: dia_corr(:)
    real(dp) :: diamagnetic_tensor(3,3)

    allocate (dia_corr(lmaxx**2))
    diamagnetic_tensor = 0.d0
    dia_corr = 0.0_dp

    if (paw_nkb == 0) return
    tmp = (0.d0,0.d0)
    ijkb0 = 0
    do nt = 1, ntyp
      do na = 1, nat
        if ( ityp(na) == nt ) then
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
              lm2=m2+l2**2
              if (na /= m_0_atom) cycle
              do ibnd = 1, nbnd
                bec_product = paw_becp(jkb,ibnd) * CONJG( paw_becp(ikb,ibnd) )
                !<apsi> s/non-trace-zero component
                ! 2/3 to separate the non-trace vanishing component
                ! 1/(2c^2) from the equation (59) in PM-PRB
                IF ( l1 == l2 .AND. m1 == m2 ) THEN
                  diamagnetic_tensor(1,1) = diamagnetic_tensor(1,1) &
                              + 2.0_dp / 3.0_dp * bec_product &
                              * radial_integral_diamagnetic(nbs1,nbs2,nt) &
                              * wg(ibnd,ik) * alpha ** 2 / 2.0_dp
                  diamagnetic_tensor(2,2) = diamagnetic_tensor(2,2) &
                              + 2.0_dp / 3.0_dp * bec_product &
                              * radial_integral_diamagnetic(nbs1,nbs2,nt) &
                              * wg(ibnd,ik) * alpha ** 2 / 2.0_dp
                  diamagnetic_tensor(3,3) = diamagnetic_tensor(3,3) &
                              + 2.0_dp / 3.0_dp * bec_product &
                              * radial_integral_diamagnetic(nbs1,nbs2,nt) &
                              * wg(ibnd,ik) * alpha ** 2 / 2.0_dp
                END IF
                               ! 2/3 to separate the non-trace vanishing component
                do lm = 5, 9
                   dia_corr(lm) =  dia_corr(lm) + bec_product / 3.0_dp &
                              * radial_integral_diamagnetic(nbs1,nbs2,nt) &
                              * ap(lm,lm1,lm2) * wg(ibnd,ik) * alpha ** 2 &
                              / 2.0_dp
                enddo
              enddo  ! ibnd
            enddo  ! jh
          enddo  ! ih
          ijkb0 = ijkb0 + paw_recon(nt)%paw_nh
        endif
      enddo  ! na
    enddo  ! bt

    !
    !  transform in cartesian coordinates
    !

    dia_corr(5:9) = - sqrt(4.0_dp * pi/5.0_dp) * dia_corr(5:9)

    diamagnetic_tensor(1,1) = diamagnetic_tensor(1,1) &
         + sqrt(3.0_dp) * dia_corr(8) - dia_corr(5)
    diamagnetic_tensor(2,2) = diamagnetic_tensor(2,2) &
         - sqrt(3.0_dp) * dia_corr(8) - dia_corr(5)
    diamagnetic_tensor(3,3) = diamagnetic_tensor(3,3) &
         + dia_corr(5) * 2.0_dp
    diamagnetic_tensor(1,2) = diamagnetic_tensor(1,2) &
         +  dia_corr(9) * sqrt(3.0_dp)
    diamagnetic_tensor(2,1) = diamagnetic_tensor(1,2)
    diamagnetic_tensor(1,3) = diamagnetic_tensor(1,3) &
         - dia_corr(6) * sqrt(3.0_dp)
    diamagnetic_tensor(3,1) = diamagnetic_tensor(1,3)
    diamagnetic_tensor(2,3) = diamagnetic_tensor(2,3) &
         - dia_corr(7) * sqrt(3.0_dp)
    diamagnetic_tensor(3,2) = diamagnetic_tensor(2,3)
    deallocate(dia_corr)

    delta_M_dia = delta_M_dia - sqrt(2.d0)*matmul(diamagnetic_tensor, m_0)



    END SUBROUTINE calc_delta_M_dia_nmr
    !------------------------------------------------------------------
    ! GIPAW correction (Delta_M_para), NMR case
    !------------------------------------------------------------------
    SUBROUTINE calc_delta_M_para_nmr
    USE ions_base,  ONLY : nat, ntyp => nsp, ityp
    USE uspp_param, ONLY : nh
    USE uspp,       ONLY : deeq
    USE paw_gipaw,  ONLY : paw_vkb, paw_nkb, paw_recon, paw_lmaxkb
    USE gipaw_module, ONLY : radial_integral_paramagnetic
    USE gipaw_module, ONLY : lx, ly, lz
    USE nmr_mod,    ONLY : m_0, m_0_atom, fine_struct
    implicit none
    complex(dp) :: tmp, becp_product
    integer :: ibnd, ijkb0, nt, na, jh, jkb, ih, ikb
    integer :: l1, m1, lm1, l2, m2, lm2, nbs1, nbs2
    if (paw_nkb == 0) return
    tmp = (0.d0,0.d0)
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
              if ( l1 /= l2 ) cycle
              !!if (l1 == 0) cycle
              if (na /= m_0_atom) cycle
              do ibnd = 1, occ
                becp_product = conjg(paw_dbecp(ikb,ibnd,ii)) * &
                                     paw_dbecp(jkb,ibnd,jj)

                tmp = tmp - fine_struct**2 * &
                      m_0(1) * lx(lm1,lm2) * (0.d0,1.d0) * &
                      radial_integral_paramagnetic(nbs1,nbs2,nt) * &
                      wg(ibnd,ik) * becp_product
                tmp = tmp - fine_struct**2 * &
                      m_0(2) * ly(lm1,lm2) * (0.d0,1.d0) * &
                      radial_integral_paramagnetic(nbs1,nbs2,nt) * &
                      wg(ibnd,ik) * becp_product
                tmp = tmp - fine_struct**2 * &
                      m_0(3)* lz(lm1,lm2) * (0.d0,1.d0) * &
                      radial_integral_paramagnetic(nbs1,nbs2,nt) * &
                      wg(ibnd,ik) * becp_product
              enddo   ! ibnd
            enddo   ! jh
          enddo   ! ih
          ijkb0 = ijkb0 + paw_recon(nt)%paw_nh
        endif   ! if (ityp...
      enddo   ! na
    enddo  ! nt
    !!PRINT*, kk, tmp
    ! check the sign and real or imag!!
    delta_M_para(kk) = delta_M_para(kk) - 2.d0*imag(tmp)/sqrt(2.d0)
    END SUBROUTINE calc_delta_M_para_nmr


  END SUBROUTINE calc_orbital_magnetization



