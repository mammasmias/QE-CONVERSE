!
! Copyright (C) 2001-2009 GIPAW and Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
PROGRAM qe_converse
  !-----------------------------------------------------------------------
  !
  ! ... This is the main driver of the QE-CONVERSE  program. 
  ! ... It controls the initialization routines.
  ! ... Features: NMR chemical shifts
  !...            EPR g-tensor
  ! ... References:
  ! ... S. Fioccola, L. Giacomazzi, D. Ceresoli, N. Richard, A. Hemeryck, L. Martin-Samos Comp.Phys.Comm.,XXX, (2024)
  ! ...
  USE gipaw_module
  USE constants, ONLY : eps12
  USE control_flags, ONLY : gamma_only, io_level
  USE environment, ONLY : environment_start, environment_end
  USE io_files, ONLY : prefix, tmp_dir, nwordwfc, iunwfc
  USE io_global, ONLY : ionode, ionode_id
  USE kinds, ONLY : DP
  USE lsda_mod, ONLY : nspin, starting_magnetization !FZ: added starting_magnetization for spinors
  USE mp, ONLY : mp_bcast
  USE mp_world, ONLY : world_comm
  USE mp_global, ONLY : mp_startup
  USE mp_pools,        ONLY : intra_pool_comm, inter_pool_comm
  USE mp_bands,        ONLY : intra_bgrp_comm, inter_bgrp_comm
  USE mp_bands,        ONLY : nbgrp
  USE mp_pools,        ONLY : nproc_pool
  USE mp_images,        ONLY : my_image_id, nimage
  USE paw_variables, ONLY : okpaw
  USE scf, ONLY : rho_core, rhog_core
  USE uspp, ONLY : okvan
  USE xc_lib, ONLY : xclib_dft_is  !FZ: for qe6.7, for hybrid functional, metaGGA, add dft_is_gradient for spinors
  USE ldaU, ONLY : lda_plus_u
  USE scf,   ONLY : rho, rho_core, rhog_core, v, vrs!FZ: test spinors added vrs 
  USE fft_rho,  ONLY : rho_g2r  
  USE fft_base,  ONLY : dfftp, dffts    
  USE noncollin_module,   ONLY : noncolin, npol, m_loc, ux, angle1, angle2  !FZ: test for spinors
  USE ions_base,          ONLY : nat, tau, ntyp => nsp, ityp, zv  !FZ: test for spinors
  USE control_flags, ONLY : restart, io_level, lscf, iprint, &
                            david, max_cg_iter, &
                            tr2, ethr, mixing_beta, nmix, niter, &
                            iverbosity
  USE gvecs, ONLY: doublegrid
  USE gvect, ONLY: gstart
  USE wvfct, ONLY: btype
  USE wvfct, ONLY: nbnd, nbndx
  USE symm_base,     ONLY : nsym
  USE noncollin_module, ONLY: report
  USE basis, ONLY: starting_wfc
  USE buffers,          ONLY : open_buffer, close_buffer, save_buffer
  USE dfunct,             ONLY : newd
  USE orbital_magnetization, ONLY : dvrs
  USE cellmd,           ONLY : cell_factor
  USE nmr_mod

  IMPLICIT NONE

  character(len=11) :: codename = 'QE-CONVERSE'
  character ( len = 256 ) :: outdir
  character (len=256), external :: trimcheck
  integer :: ios, na, ipol
  logical :: exst

  NAMELIST / input_qeconverse / prefix, outdir, &
                        diagonalization, isolve, iverbosity, &
                        verbosity, q_gipaw, dudk_method, &
                        diago_thr_init, conv_threshold, &
                        tr2, mixing_beta, assume_isolated, &
                        lambda_so, m_0, m_0_atom

! begin with the initialization part                
#if defined(__MPI)
  call mp_startup (start_images=.true., images_only=.false. )
#else
  call mp_startup(start_images=.false.)
#endif
  call environment_start (codename)

if (.not. ionode .or. my_image_id > 0) goto 400

  call input_from_file()

! define input default values
  prefix = 'prefix'
  CALL get_environment_variable ( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM ( outdir ) == ' ' ) outdir = './'
  diagonalization = 'david'
  verbosity = 'high'
  isolve = -1
  iverbosity = -1 
  q_gipaw = 0.01d0
  dudk_method = 'covariant'
  diago_thr_init = 1d-7
  conv_threshold = 1e-8
  mixing_beta = 0.5
  tr2 = 1e-8
  lambda_so(:) = 0.d0
  m_0(:) = 0.d0

    read ( 5, input_qeconverse, iostat = ios )
  tmp_dir = outdir

  ! further checks
  if (isolve /= -1) &
     call infomsg('qe_converse', '*** isolve is obsolete, use diagonalization instead ***')
  if (iverbosity /= -1) &
     call infomsg('qe_converse', '*** iverbosity is obsolete, use verbosity instead ***')

  select case (diagonalization)
     case('david')
       isolve = 0
     case('cg')
       isolve = 1
     case default
       call errore('qe_converse', 'diagonalization can be ''david'' or ''cg''', 1)
  end select

  select case (verbosity)
     case('low')
       iverbosity = 1
     case('medium')
       iverbosity = 11
     case('high')
       iverbosity = 21
     case default
       call errore('qe_converse', 'verbosity can be ''low'', ''medium'' or ''high''', 1)
  end select
400 continue


#ifdef __MPI
  ! broadcast input variables  
  call gipaw_bcast_input
#endif

  io_level = 1
  cell_factor = 1.1d0
  call start_clock ('read_file')
  !read ground state wavefunctions  
  CALL read_file ( )
  call stop_clock ('read_file')
  call gipaw_setup ( )
  if (any(m_0 /= 0.d0)) then
     call init_nmr
  endif

  call newscf

  call calc_orbital_magnetization ( )
  
  call environment_end(codename)
  call print_clock_gipaw
  call stop_code( .true. )

  STOP


END PROGRAM qe_converse

#ifdef __MPI
!-----------------------------------------------------------------------
SUBROUTINE gipaw_bcast_input
!-----------------------------------------------------------------------
  !
  ! ... Broadcast input data to all processors
  !
  USE gipaw_module
  USE nmr_mod
  USE control_flags, ONLY : tr2, mixing_beta, iverbosity 
  USE mp_world,      ONLY : world_comm
  USE mp,            ONLY : mp_bcast
  USE io_files, ONLY : prefix, tmp_dir

  implicit none
  integer :: root = 0


  CALL mp_bcast ( tmp_dir, root, world_comm )
  CALL mp_bcast ( prefix, root, world_comm )
  CALL mp_bcast ( diagonalization, root, world_comm )
  CALL mp_bcast ( verbosity, root, world_comm )
  CALL mp_bcast (isolve, root, world_comm)
  CALL mp_bcast (iverbosity, root, world_comm)
  CALL mp_bcast ( q_gipaw, root, world_comm )
  CALL mp_bcast ( dudk_method, root, world_comm )
  CALL mp_bcast ( diago_thr_init, root, world_comm )
  CALL mp_bcast ( conv_threshold, root, world_comm )
  CALL mp_bcast ( mixing_beta, root, world_comm )
  CALL mp_bcast ( tr2, root, world_comm )
  CALL mp_bcast ( lambda_so, root, world_comm )
  CALL mp_bcast ( m_0, root, world_comm )
  CALL mp_bcast ( m_0_atom, root, world_comm )
  CALL mp_bcast ( assume_isolated, root, world_comm )


END SUBROUTINE gipaw_bcast_input
#endif
!-----------------------------------------------------------------------
SUBROUTINE print_clock_gipaw
  !-----------------------------------------------------------------------
  !
  ! ... Print clocks
  !
  USE io_global,  ONLY : stdout
  IMPLICIT NONE


  write(stdout,*) '    Initialization:'
  write(stdout,*)
  call print_clock ('read_file')
  call print_clock ('gipaw_setup')
  call print_clock ('wfcinit')
  call print_clock ('hinit0')
  call print_clock ('potinit')
  write(stdout,*)
  write(stdout,*) '    SCF calculation:'
  write(stdout,*)
  call print_clock ('electrons')
  call print_clock ('h_psi')
  CALL print_clock( 'add_so_valence' )
  call print_clock( 'add_so_Fnl' )
  write(stdout,*)
  write(stdout,*) '    g-tensor:'
  write(stdout,*)
  call print_clock ('compute_dudk')
  call print_clock ('orbital_magnetization')

END SUBROUTINE print_clock_gipaw

