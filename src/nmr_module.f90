!-----------------------------------------------------------------------
MODULE nmr_mod
  !
  ! Module for variables related to the calculation of the NMR shielding.
  !
  USE kinds            , ONLY : DP
  !
  IMPLICIT NONE
  !
  !
  CHARACTER(LEN=1), PARAMETER :: cr          = char(10)
  REAL(DP)        , PARAMETER :: fine_struct = 7.2973525376D-3
  COMPLEX(DP)     , PARAMETER :: IM          = (0.0D0, 1.0D0)
  !
  INTEGER                     :: m_0_atom = 1, m_0_comp = 1, n_nlpp = 0
  INTEGER, PARAMETER :: nlpp_max            = 1000
  REAL(DP)                    :: m_0(3), sigma = 0.0D0, toggle = 1.0D0
  REAL(DP)                    :: r_nlpp(nlpp_max)
  COMPLEX(DP),    ALLOCATABLE :: A_vec(:,:), grad_psi(:,:), A_0(:,:)
  CHARACTER(len=1)            :: m_0_dir
  LOGICAL                     :: delete_dudk_files

!-----------------------------------------------------------------------
  CONTAINS
!-----------------------------------------------------------------------
  subroutine set_m_0_dir()
    if (m_0(1) > 0.5) m_0_dir = 'x'
    if (m_0(2) > 0.5) m_0_dir = 'y'
    if (m_0(3) > 0.5) m_0_dir = 'z'
  end subroutine

  function dudk_name_x()    
    character(len=4) :: atom_index
    character(len=80) :: dudk_name_x
    call set_m_0_dir
    write(atom_index, '(I4.4)') m_0_atom
    dudk_name_x = 'atom.' // atom_index // m_0_dir // 'dudk1'
  end function

  function dudk_name_y()    
    character(len=4) :: atom_index
    character(len=80) :: dudk_name_y
    call set_m_0_dir
    write(atom_index, '(I4.4)') m_0_atom
    dudk_name_y = 'atom.' // atom_index // m_0_dir // 'dudk2'
  end function

  function dudk_name_z()    
    character(len=4) :: atom_index
    character(len=80) :: dudk_name_z
    call set_m_0_dir
    write(atom_index, '(I4.4)') m_0_atom
    dudk_name_z = 'atom.' // atom_index // m_0_dir // 'dudk3'
  end function

#if 0
!-----------------------------------------------------------------------
!  SUBROUTINE calc_chemical_shift(orb_magn_tot, nmr_shift_core)
  !-----------------------------------------------------------------------
!  USE kinds,                ONLY : dp
  USE io_global,            ONLY : stdout
  USE ions_base,            ONLY : nat, atm, ityp
  USE parameters,           ONLY : ntypx
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
#endif
!-----------------------------------------------------------------------
END MODULE nmr_mod
!-----------------------------------------------------------------------
